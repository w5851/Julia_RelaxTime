"""
Scan Mott transition temperature and meson masses vs xi at muB=0.

Usage:
  julia --project=. scripts/analysis/scan_mott_meson_vs_xi_mu0.jl [--smoke] [--output <path>]
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "analysis", "mott_reference_mapping.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .MottReferenceMapping: estimate_mott_temperature
using .ScanCSV: ScanCSV

const _MESON_WORKFLOW_MOD = Main.Models.meson_workflow_module()

struct Options
    output::String
    smoke::Bool
    xi_values::Vector{Float64}
    t_values_mev::Vector{Float64}
    mesons::Vector{Symbol}
end

function _frange(lo::Float64, hi::Float64, step::Float64)
    vals = Float64[]
    x = lo
    while x <= hi + 1e-9
        push!(vals, round(x; digits=10))
        x += step
    end
    return vals
end

function parse_args(args::Vector{String})
    output = joinpath("data", "outputs", "results", "relaxtime", "scan", "mott_meson_vs_xi_mu0.csv")
    smoke = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--output"
            i == length(args) && error("missing value for --output")
            i += 1
            output = args[i]
        elseif arg == "--smoke"
            smoke = true
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/analysis/scan_mott_meson_vs_xi_mu0.jl [--smoke] [--output <path>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    if smoke
        xi_values = [-0.1, 0.0, 0.1]
        t_values = [150.0, 170.0, 190.0]
        mesons = [:pi, :K]
    else
        xi_values = _frange(-0.4, 0.4, 0.05)
        t_values = _frange(140.0, 240.0, 10.0)
        mesons = [:pi, :K, :eta, :eta_prime, :sigma_pi, :sigma_K, :sigma, :sigma_prime]
    end

    return Options(output, smoke, xi_values, t_values, mesons)
end

function _row_values(cols::Vector{String}, row::Dict{String,Any})
    return [string(get(row, c, "")) for c in cols]
end

function _calc_point(T_MeV::Float64, xi::Float64, mesons::Vector{Symbol})
    T_fm = T_MeV / ħc_MeV_fm
    mu_fm = 0.0
    return _MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
        T_fm,
        mu_fm;
        xi=xi,
        mesons=Tuple(mesons),
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=25,),
        mass_kwargs=(iterations=25,),
    )
end

function main()
    opts = parse_args(ARGS)
    mkpath(dirname(opts.output))

    cols = [
        "xi", "meson", "T_Mott_MeV", "method", "approx", "status",
        "points_used", "T_min_MeV", "T_max_MeV",
    ]

    open(opts.output, "w") do io
        ScanCSV.write_metadata(io, Dict(
            "format" => "scan_csv_v1",
            "script" => "scripts/analysis/scan_mott_meson_vs_xi_mu0.jl",
            "muB_MeV" => "0.0",
            "smoke" => string(opts.smoke),
        ))
        ScanCSV.write_header(io, cols)

        for xi in opts.xi_values
            per_meson = Dict(m => (Float64[], Float64[]) for m in opts.mesons)
            status = "success"

            for T_MeV in opts.t_values_mev
                local res
                try
                    res = _calc_point(T_MeV, xi, opts.mesons)
                catch
                    status = "fallback"
                    continue
                end

                for m in opts.mesons
                    mr = res.meson_results[m]
                    gap = if hasproperty(mr, :gap)
                        Float64(mr.gap)
                    elseif hasproperty(mr, :gaps)
                        Float64(mr.gaps.min)
                    else
                        NaN
                    end
                    if isfinite(gap)
                        push!(per_meson[m][1], T_MeV)
                        push!(per_meson[m][2], gap)
                    end
                end
            end

            for m in opts.mesons
                Ts, gaps = per_meson[m]
                if length(Ts) < 2
                    row = Dict{String,Any}(
                        "xi" => xi,
                        "meson" => String(m),
                        "T_Mott_MeV" => NaN,
                        "method" => "insufficient_data",
                        "approx" => true,
                        "status" => "fail",
                        "points_used" => length(Ts),
                        "T_min_MeV" => NaN,
                        "T_max_MeV" => NaN,
                    )
                    println(io, join(_row_values(cols, row), ','))
                    continue
                end

                est = estimate_mott_temperature(Ts, gaps)
                row = Dict{String,Any}(
                    "xi" => xi,
                    "meson" => String(m),
                    "T_Mott_MeV" => est.T_mott_MeV,
                    "method" => String(est.method),
                    "approx" => Bool(est.approx),
                    "status" => status,
                    "points_used" => length(Ts),
                    "T_min_MeV" => minimum(Ts),
                    "T_max_MeV" => maximum(Ts),
                )
                println(io, join(_row_values(cols, row), ','))
            end
        end
    end

    println("Wrote scan CSV: " * opts.output)
end

main()
