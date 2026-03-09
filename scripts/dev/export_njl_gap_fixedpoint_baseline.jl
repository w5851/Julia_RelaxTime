#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_NJL_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "njl", "baseline_njl_gap_fixedpoints_v1.csv")
const DEFAULT_NJL2_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "njl", "baseline_njl2_gap_fixedpoints_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))

const GAP_POINTS = [
    (label="cold_mu0", T_MeV=50.0, mu_MeV=0.0),
    (label="mid_mu0", T_MeV=100.0, mu_MeV=0.0),
    (label="hot_mu0", T_MeV=150.0, mu_MeV=0.0),
    (label="mid_dense", T_MeV=100.0, mu_MeV=200.0),
    (label="hot_dense", T_MeV=150.0, mu_MeV=200.0),
]

function parse_args(args::Vector{String})
    njl_output = DEFAULT_NJL_OUTPUT
    njl2_output = DEFAULT_NJL2_OUTPUT
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--njl-output"
            i == length(args) && error("missing value for --njl-output")
            i += 1
            njl_output = args[i]
        elseif arg == "--njl2-output"
            i == length(args) && error("missing value for --njl2-output")
            i += 1
            njl2_output = args[i]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/export_njl_gap_fixedpoint_baseline.jl [--njl-output <path>] [--njl2-output <path>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return njl_output, njl2_output
end

function export_model(kind::Symbol, output::String)
    mkpath(dirname(output))
    model = Models.create_model(kind)

    open(output, "w") do io
        println(io, "label,T_MeV,mu_MeV,phi_u,phi_d,phi_s,Phi,PhiBar,omega")
        for pt in GAP_POINTS
            T_fm = pt.T_MeV / Constants_PNJL.ħc_MeV_fm
            mu_fm = pt.mu_MeV / Constants_PNJL.ħc_MeV_fm
            st = Models.solve_gap(model, T_fm, mu_fm; p_num=24, t_num=6, xi=0.0, residual_norm_max=1e-6)
            x = Models.state_vector(st)
            omega = Models.omega(model, st, T_fm, mu_fm; p_num=24, t_num=6, xi=0.0)
            @printf(io, "%s,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                pt.label, pt.T_MeV, pt.mu_MeV, x[1], x[2], x[3], x[4], x[5], omega)
        end
    end
end

function main(args::Vector{String})
    njl_output, njl2_output = parse_args(args)
    export_model(:NJL, njl_output)
    export_model(:NJL2, njl2_output)
    println("NJL gap baseline written to: " * njl_output)
    println("NJL2 gap baseline written to: " * njl2_output)
end

main(ARGS)