#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_chi_b_taylordiff_v2.csv")
const CEP_REFERENCE_PATH = joinpath(
    PROJECT_ROOT, "data", "reference", "pnjl", "issue130_phase_reference_v2",
    "accepted", "tables", "cep_boundary_accepted_phase_map_v1.csv",
)

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using Main.Constants_PNJL: ħc_MeV_fm

function parse_args(args::Vector{String})
    output = DEFAULT_OUTPUT
    cep_reference = CEP_REFERENCE_PATH
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--output"
            i == length(args) && error("missing value for --output")
            i += 1
            output = args[i]
        elseif arg == "--cep-reference"
            i == length(args) && error("missing value for --cep-reference")
            i += 1
            cep_reference = args[i]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/export_pnjl_chi_b_taylordiff_baseline.jl [--output <path>] [--cep-reference <path>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return (output=output, cep_reference=cep_reference)
end

function load_xi0_cep(path::String)
    lines = readlines(path)
    isempty(lines) && error("CEP reference is empty: $path")
    header = split(strip(lines[1]), ',')
    idx_xi = findfirst(==("xi"), header)
    idx_T = findfirst(==("T_CEP_MeV"), header)
    idx_T === nothing && (idx_T = findfirst(==("T_midpoint_MeV"), header))
    idx_muB = findfirst(==("muB_CEP_MeV"), header)
    idx_muq = findfirst(==("muq_CEP_proxy_MeV"), header)
    idx_muq === nothing && (idx_muq = findfirst(==("mu_CEP_proxy_MeV"), header))
    idx_muq === nothing && (idx_muq = findfirst(==("muq_CEP_MeV"), header))
    (idx_xi === nothing || idx_T === nothing || (idx_muB === nothing && idx_muq === nothing)) &&
        error("CEP reference is missing xi/T/mu columns: $path")
    for line in lines[2:end]
        cols = split(strip(line), ',')
        parse(Float64, cols[idx_xi]) == 0.0 || continue
        muB = idx_muB === nothing ? 3.0 * parse(Float64, cols[idx_muq]) : parse(Float64, cols[idx_muB])
        return (T_CEP_MeV=parse(Float64, cols[idx_T]), muB_CEP_MeV=muB)
    end
    error("xi=0 CEP reference row not found: $path")
end

function baseline_points(cep_reference::String=CEP_REFERENCE_PATH)
    cep = load_xi0_cep(cep_reference)
    muB_CEP_fm = cep.muB_CEP_MeV / ħc_MeV_fm
    return [
        (label="fixed_p8_t4", T_fm=0.55, muB_fm=0.30, xi=0.0, p_num=8, t_num=4, orders=1:4),
        (label="cep_true_p8_t4", T_fm=cep.T_CEP_MeV / ħc_MeV_fm, muB_fm=muB_CEP_fm, xi=0.0, p_num=8, t_num=4, orders=1:4),
        (label="cep_minus_0p5_p8_t4", T_fm=(cep.T_CEP_MeV - 0.5) / ħc_MeV_fm, muB_fm=muB_CEP_fm, xi=0.0, p_num=8, t_num=4, orders=1:4),
        (label="cep_plus_0p5_p8_t4", T_fm=(cep.T_CEP_MeV + 0.5) / ħc_MeV_fm, muB_fm=muB_CEP_fm, xi=0.0, p_num=8, t_num=4, orders=1:4),
        (label="cep_true_p24_t8", T_fm=cep.T_CEP_MeV / ħc_MeV_fm, muB_fm=muB_CEP_fm, xi=0.0, p_num=24, t_num=8, orders=1:3),
    ]
end

function main(args::Vector{String})
    parsed = parse_args(args)
    output = parsed.output
    mkpath(dirname(output))
    PNJL = Models.pnjl_module()

    open(output, "w") do io
        println(io, "label,T_fm,muB_fm,xi,p_num,t_num,order,chi_B_taylordiff")
        for pt in baseline_points(parsed.cep_reference)
            kwargs = (; xi=pt.xi, p_num=pt.p_num, t_num=pt.t_num)
            for order in pt.orders
                value = PNJL.chi_B(pt.T_fm, pt.muB_fm; order=order, derivative_backend=:taylordiff, kwargs...)
                @printf(io, "%s,%.16e,%.16e,%.6f,%d,%d,%d,%.16e\n",
                    pt.label,
                    pt.T_fm,
                    pt.muB_fm,
                    pt.xi,
                    pt.p_num,
                    pt.t_num,
                    order,
                    value,
                )
            end
        end
    end

    println("PNJL chi_B TaylorDiff baseline written to: " * output)
end

main(ARGS)
