#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_scan_fixedpoints_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
using .PNJL: run_tmu_scan, run_trho_scan

function parse_args(args::Vector{String})
    output = DEFAULT_OUTPUT
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--output"
            i == length(args) && error("missing value for --output")
            i += 1
            output = args[i]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/export_pnjl_scan_fixedpoint_baseline.jl [--output <path>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return output
end

function _parse_single_row(path::String)
    lines = readlines(path)
    length(lines) >= 2 || error("invalid csv output: $path")
    header = split(lines[1], ',')
    cols = split(lines[2], ',')
    return Dict(header[i] => cols[i] for i in eachindex(header))
end

function _f(row::AbstractDict, key::String)
    haskey(row, key) || error("missing column '$key'")
    return parse(Float64, row[key])
end

function main(args::Vector{String})
    output = parse_args(args)
    mkpath(dirname(output))

    tmp_dir = mktempdir()
    tmu_out = joinpath(tmp_dir, "tmu_single.csv")
    trho_out = joinpath(tmp_dir, "trho_single.csv")

    run_tmu_scan(
        T_values=[150.0],
        mu_values=[0.0],
        xi_values=[0.0],
        output_path=tmu_out,
        overwrite=true,
        resume=false,
        use_phase_aware=false,
        thermo_backend=:legacy,
        p_num=8,
        t_num=4,
        iterations=80,
    )

    run_trho_scan(
        T_values=[150.0],
        rho_values=[0.2],
        xi_values=[0.0],
        output_path=trho_out,
        overwrite=true,
        resume=false,
        reverse_rho=false,
        thermo_backend=:legacy,
        p_num=8,
        t_num=4,
        iterations=80,
    )

    tmu_row = _parse_single_row(tmu_out)
    trho_row = _parse_single_row(trho_out)

    open(output, "w") do io
        println(io, "kind,T,mu,rho,xi,pressure_fm4,rho_out,mu_avg_MeV,entropy_fm3,energy_fm4")
        @printf(io, "tmu,%.6f,%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e\n",
            150.0, 0.0, NaN, 0.0,
            _f(tmu_row, "pressure_fm4"),
            _f(tmu_row, "rho"),
            NaN,
            _f(tmu_row, "entropy_fm3"),
            _f(tmu_row, "energy_fm4"),
        )
        @printf(io, "trho,%.6f,%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e\n",
            150.0, NaN, 0.2, 0.0,
            _f(trho_row, "pressure_fm4"),
            NaN,
            _f(trho_row, "mu_avg_MeV"),
            _f(trho_row, "entropy_fm3"),
            _f(trho_row, "energy_fm4"),
        )
    end

    println("baseline exported to: " * output)
end

main(ARGS)
