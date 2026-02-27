#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Models: run_trho_scan
using DelimitedFiles
using Printf

const DEFAULT_BASELINE_DIR = joinpath(PROJECT_ROOT, "data", "outputs", "results", "debug", "baseline_compare")
const DEFAULT_RAW = joinpath(DEFAULT_BASELINE_DIR, "julia_baseline_raw.csv")
const DEFAULT_UNIFIED = joinpath(DEFAULT_BASELINE_DIR, "julia_baseline_unified.csv")
const DEFAULT_FORTRAN = joinpath(DEFAULT_BASELINE_DIR, "fortran_baseline_unified.csv")

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :profile => "a1",
        :fortran_csv => DEFAULT_FORTRAN,
        :raw_output => DEFAULT_RAW,
        :output => DEFAULT_UNIFIED,
        :xi => 0.0,
        :ud_ratio => 0.876,
        :s_target => 0.0,
        :model_kind => :RPNJL,
        :thermo_backend => :models,
        :solver_backend => :models,
        :seed_policy => :hybrid_continuity,
        :p_num => 24,
        :t_num => 8,
        :iterations => 140,
        :max_points => 200,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            val = args[i + 1]
            i += 1
            return val
        end
        if arg == "--profile"
            opts[:profile] = lowercase(require_value())
        elseif arg == "--fortran-csv"
            opts[:fortran_csv] = require_value()
        elseif arg == "--raw-output"
            opts[:raw_output] = require_value()
        elseif arg == "--output"
            opts[:output] = require_value()
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--ud-ratio"
            opts[:ud_ratio] = parse(Float64, require_value())
        elseif arg == "--s-target"
            opts[:s_target] = parse(Float64, require_value())
        elseif arg == "--model-kind"
            opts[:model_kind] = Symbol(uppercase(require_value()))
        elseif arg == "--seed-policy"
            opts[:seed_policy] = Symbol(lowercase(require_value()))
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--iterations"
            opts[:iterations] = parse(Int, require_value())
        elseif arg == "--max-points"
            opts[:max_points] = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return opts
end

function print_usage()
    println("Usage: julia scripts/debug/export_julia_baseline.jl [options]\n")
    println("Options:")
    println("  --profile a1|fortran         Point source profile (default a1)")
    println("  --fortran-csv <path>         Fortran unified CSV for profile=fortran")
    println("  --raw-output <path>          Raw TrhoScan output path")
    println("  --output <path>              Unified Julia output path")
    println("  --xi <value>                 ξ value (default 0.0)")
    println("  --ud-ratio <value>           ρu/ρd target (default 0.876)")
    println("  --s-target <value>           ρs target (default 0.0)")
    println("  --model-kind PNJL|RPNJL      Model kind (default RPNJL)")
    println("  --seed-policy <symbol>       Trho seed policy (default hybrid_continuity)")
    println("  --p-num / --t-num <int>      Integration grids")
    println("  --iterations <int>           Nonlinear iterations (default 140)")
    println("  --max-points <int>           Limit points for profile=fortran (default 200)")
end

function a1_points()
    return [
        (3.0, 0.010),
        (5.0, 0.020),
        (8.0, 0.050),
        (20.0, 0.100),
        (50.0, 0.300),
        (100.0, 0.500),
        (150.0, 0.800),
        (200.0, 1.200),
    ]
end

function read_csv_lines(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("csv file is empty: $path")
    header = split(lines[1], ',')
    rows = [split(line, ',') for line in lines[2:end] if !isempty(strip(line))]
    return header, rows
end

function col_index_map(header::Vector{<:AbstractString})
    return Dict(name => idx for (idx, name) in enumerate(header))
end

function parse_float(cell::AbstractString)
    s = strip(cell)
    if isempty(s)
        return NaN
    end
    t = lowercase(s)
    if t == "nan"
        return NaN
    end
    return parse(Float64, s)
end

function points_from_fortran(path::AbstractString; max_points::Int=200)
    header, rows = read_csv_lines(path)
    cmap = col_index_map(header)
    iT = get(cmap, "T_MeV", 0)
    irho = get(cmap, "rho_target_over_rho0", 0)
    iT == 0 && error("missing column T_MeV in $path")
    irho == 0 && error("missing column rho_target_over_rho0 in $path")

    pts = Tuple{Float64, Float64}[]
    seen = Set{Tuple{Int, Int}}()
    for row in rows
        T = parse_float(row[iT])
        rho = parse_float(row[irho])
        key = (round(Int, T * 1000), round(Int, rho * 1000))
        if !(key in seen)
            push!(pts, (T, rho))
            push!(seen, key)
        end
        length(pts) >= max_points && break
    end
    return pts
end

function select_points(opts)
    profile = opts[:profile]
    if profile == "a1"
        return a1_points()
    elseif profile == "fortran"
        return points_from_fortran(String(opts[:fortran_csv]); max_points=Int(opts[:max_points]))
    else
        error("unsupported profile: $profile")
    end
end

function convert_raw_to_unified(raw_csv::AbstractString, out_csv::AbstractString)
    header, rows = read_csv_lines(raw_csv)
    cmap = col_index_map(header)

    required = [
        "T_MeV", "rho", "mu_u_MeV", "mu_d_MeV", "mu_s_MeV", "mu_B_MeV", "mu_Q_MeV", "mu_S_MeV",
        "pressure_fm4", "rho_u_fm3", "rho_d_fm3", "rho_s_fm3", "phi_u", "phi_d", "phi_s", "Phi1", "Phi2",
        "residual_norm", "converged", "message",
    ]
    for name in required
        haskey(cmap, name) || error("missing required raw column: $name")
    end

    mkpath(dirname(out_csv))

    function sanitize_cell(value)
        txt = string(value)
        txt = replace(txt, '\r' => ' ')
        txt = replace(txt, '\n' => ' ')
        txt = replace(txt, ',' => ';')
        txt = replace(txt, '"' => "'")
        return txt
    end

    open(out_csv, "w") do io
        println(io, join([
            "row_id",
            "T_MeV",
            "rho_target_over_rho0",
            "mu_B_MeV",
            "mu_Q_MeV",
            "mu_S_MeV",
            "pressure_fm4",
            "phi_u",
            "phi_d",
            "phi_s",
            "Phi",
            "PhiBar",
            "mu_u_MeV",
            "mu_d_MeV",
            "mu_s_MeV",
            "rho_u_fm3",
            "rho_d_fm3",
            "rho_s_fm3",
            "converged",
            "residual_norm",
            "message",
        ], ","))

        row_id = 0
        for row in rows
            row_id += 1
            rec = [
                string(row_id),
                row[cmap["T_MeV"]],
                row[cmap["rho"]],
                row[cmap["mu_B_MeV"]],
                row[cmap["mu_Q_MeV"]],
                row[cmap["mu_S_MeV"]],
                row[cmap["pressure_fm4"]],
                row[cmap["phi_u"]],
                row[cmap["phi_d"]],
                row[cmap["phi_s"]],
                row[cmap["Phi1"]],
                row[cmap["Phi2"]],
                row[cmap["mu_u_MeV"]],
                row[cmap["mu_d_MeV"]],
                row[cmap["mu_s_MeV"]],
                row[cmap["rho_u_fm3"]],
                row[cmap["rho_d_fm3"]],
                row[cmap["rho_s_fm3"]],
                row[cmap["converged"]],
                row[cmap["residual_norm"]],
                sanitize_cell(row[cmap["message"]]),
            ]
            println(io, join(rec, ","))
        end
    end
end

function run_export(opts)
    points = select_points(opts)
    isempty(points) && error("no points selected")

    raw_output = String(opts[:raw_output])
    out_csv = String(opts[:output])
    mkpath(dirname(raw_output))

    for (idx, (T, rho)) in enumerate(points)
        println(@sprintf("[%d/%d] T=%.3f MeV rho=%.6f", idx, length(points), T, rho))
        stats = run_trho_scan(
            ;
            T_values=[Float64(T)],
            rho_values=[Float64(rho)],
            xi_values=[Float64(opts[:xi])],
            output_path=raw_output,
            overwrite=(idx == 1),
            resume=false,
            reverse_rho=false,
            constraint_mode=:fixed_asymmetric_rho,
            asym_ud_ratio_target=Float64(opts[:ud_ratio]),
            asym_s_target=Float64(opts[:s_target]),
            thermo_backend=Symbol(opts[:thermo_backend]),
            solver_backend=Symbol(opts[:solver_backend]),
            model_kind=Symbol(opts[:model_kind]),
            seed_policy=Symbol(opts[:seed_policy]),
            p_num=Int(opts[:p_num]),
            t_num=Int(opts[:t_num]),
            iterations=Int(opts[:iterations]),
        )
        println(@sprintf("  -> success=%d failure=%d", stats.success, stats.failure))
    end

    convert_raw_to_unified(raw_output, out_csv)
    println("Unified Julia baseline written to: $out_csv")
end

function main()
    opts = parse_args(copy(ARGS))
    run_export(opts)
end

main()
