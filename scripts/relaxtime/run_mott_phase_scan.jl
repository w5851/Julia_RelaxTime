#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))
include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflows", "MesonMassWorkflow.jl"))

using .ScanCSV: ScanCSV
using .Constants_PNJL: ħc_MeV_fm
using .MesonMassWorkflow: solve_gap_and_meson_point
using TOML
using Dates
using SHA

struct ScanOptions
    outdir::String
    output_csv::String
    config_path::String
    xi_list::Vector{Float64}
    T_min_MeV::Float64
    T_max_MeV::Float64
    T_step_MeV::Float64
    muB_MeV::Float64
    resume::Bool
    overwrite::Bool
    p_num::Int
    t_num::Int
    max_iter::Int
end

@inline function _json_escape(s::AbstractString)
    out = IOBuffer()
    for c in s
        if c == '"'
            print(out, "\\\"")
        elseif c == '\\'
            print(out, "\\\\")
        elseif c == '\n'
            print(out, "\\n")
        elseif c == '\r'
            print(out, "\\r")
        elseif c == '\t'
            print(out, "\\t")
        else
            print(out, c)
        end
    end
    return String(take!(out))
end

function _to_json(x)
    if x === nothing
        return "null"
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x isa Integer || x isa AbstractFloat
        return string(x)
    elseif x isa AbstractString
        return "\"$(_json_escape(x))\""
    elseif x isa AbstractVector
        return "[" * join((_to_json(v) for v in x), ",") * "]"
    elseif x isa AbstractDict
        pairs_sorted = sort(collect(pairs(x)); by=kv -> String(kv.first))
        parts = String[]
        for (k, v) in pairs_sorted
            push!(parts, "\"$(_json_escape(String(k)))\":" * _to_json(v))
        end
        return "{" * join(parts, ",") * "}"
    else
        return "\"$(_json_escape(string(x)))\""
    end
end

function _write_json(path::String, x)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, _to_json(x))
    end
end

@inline function _fmt(x)
    x isa Bool && return x ? "true" : "false"
    return string(x)
end

function _print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_mott_phase_scan.jl [options]")
    println("Options:")
    println("  --config <path>      TOML config (default config/workflows/relaxtime/profiles/mott_phase_muB0_xi_scan.toml)")
    println("  --outdir <path>      Output directory")
    println("  --output <name>      Output CSV filename under outdir (default mott_phase_scan.csv)")
    println("  --tmin <MeV>         Temperature min (default 120)")
    println("  --tmax <MeV>         Temperature max (default 260)")
    println("  --tstep <MeV>        Temperature step (default 2)")
    println("  --xi-list v1,v2,...  Xi list override")
    println("  --resume             Resume by existing key")
    println("  --overwrite          Overwrite output csv")
    println("  --p-num <int>        Gap solver momentum nodes")
    println("  --t-num <int>        Gap solver angle nodes")
    println("  --max-iter <int>     Solver iteration cap")
    println("  -h, --help           Show help")
end

function _parse_xi_list(raw::AbstractString)
    vals = Float64[]
    for seg in split(raw, ',')
        s = strip(seg)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && throw(ArgumentError("--xi-list cannot be empty"))
    return unique(sort(vals))
end

function _default_cfg_dict()
    return Dict{String,Any}(
        "scan" => Dict{String,Any}(
            "mott_phase" => Dict{String,Any}(
                "muB_MeV" => 0.0,
                "xi_list" => Any[-0.3, -0.15, 0.0, 0.15, 0.3],
                "T_min_MeV" => 120.0,
                "T_max_MeV" => 260.0,
                "T_step_MeV" => 2.0,
                "p_num" => 12,
                "t_num" => 6,
                "max_iter" => 40,
                "resume" => true,
                "overwrite" => false,
            ),
        ),
    )
end

function _build_options(args::Vector{String})
    default_cfg_path = joinpath(PROJECT_ROOT, "config", "workflows", "relaxtime", "profiles", "mott_phase_muB0_xi_scan.toml")
    outdir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "mott_phase", "default")
    output_csv = "mott_phase_scan.csv"
    config_path = default_cfg_path

    cli = Dict{String,Any}()

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $arg"))
            i += 1
            return args[i]
        end

        if arg == "--config"
            config_path = require_value()
        elseif arg == "--outdir"
            outdir = require_value()
        elseif arg == "--output"
            output_csv = require_value()
        elseif arg == "--tmin"
            cli["T_min_MeV"] = parse(Float64, require_value())
        elseif arg == "--tmax"
            cli["T_max_MeV"] = parse(Float64, require_value())
        elseif arg == "--tstep"
            cli["T_step_MeV"] = parse(Float64, require_value())
        elseif arg == "--xi-list"
            cli["xi_list"] = _parse_xi_list(require_value())
        elseif arg == "--resume"
            cli["resume"] = true
        elseif arg == "--overwrite"
            cli["overwrite"] = true
        elseif arg == "--p-num"
            cli["p_num"] = parse(Int, require_value())
        elseif arg == "--t-num"
            cli["t_num"] = parse(Int, require_value())
        elseif arg == "--max-iter"
            cli["max_iter"] = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            _print_usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end

    cfg = _default_cfg_dict()
    if isfile(config_path)
        tcfg = TOML.parsefile(config_path)
        if haskey(tcfg, "scan") && haskey(tcfg["scan"], "mott_phase")
            merge!(cfg["scan"]["mott_phase"], tcfg["scan"]["mott_phase"])
        end
    end
    merge!(cfg["scan"]["mott_phase"], cli)

    mp = cfg["scan"]["mott_phase"]
    xi_list = unique(sort(Float64.(mp["xi_list"])))
    T_min_MeV = Float64(mp["T_min_MeV"])
    T_max_MeV = Float64(mp["T_max_MeV"])
    T_step_MeV = Float64(mp["T_step_MeV"])
    muB_MeV = Float64(get(mp, "muB_MeV", 0.0))
    resume = Bool(get(mp, "resume", true))
    overwrite = Bool(get(mp, "overwrite", false))
    p_num = Int(get(mp, "p_num", 12))
    t_num = Int(get(mp, "t_num", 6))
    max_iter = Int(get(mp, "max_iter", 40))

    T_step_MeV > 0 || throw(ArgumentError("T_step_MeV must be positive"))
    T_max_MeV >= T_min_MeV || throw(ArgumentError("T_max_MeV must be >= T_min_MeV"))

    return ScanOptions(
        outdir,
        output_csv,
        config_path,
        xi_list,
        T_min_MeV,
        T_max_MeV,
        T_step_MeV,
        muB_MeV,
        resume,
        overwrite,
        p_num,
        t_num,
        max_iter,
    ), cfg
end

function _write_run_artifacts(opts::ScanOptions, cfg::Dict{String,Any}, out_csv::String)
    run_id = string(Dates.format(now(UTC), dateformat"yyyymmddTHHMMSS"), "_", bytes2hex(rand(UInt8, 4)))
    effective = Dict{String,Any}(
        "schema_version" => "v1",
        "profile_name" => "mott_phase_muB0_xi_scan",
        "scan" => Dict{String,Any}(
            "mott_phase" => Dict{String,Any}(
                "muB_MeV" => opts.muB_MeV,
                "xi_list" => opts.xi_list,
                "T_min_MeV" => opts.T_min_MeV,
                "T_max_MeV" => opts.T_max_MeV,
                "T_step_MeV" => opts.T_step_MeV,
                "p_num" => opts.p_num,
                "t_num" => opts.t_num,
                "max_iter" => opts.max_iter,
                "resume" => opts.resume,
                "overwrite" => opts.overwrite,
            ),
        ),
    )
    effective_json = _to_json(effective)
    config_hash = bytes2hex(sha256(effective_json))

    manifest = Dict{String,Any}(
        "run_id" => run_id,
        "timestamp_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ"),
        "script" => "scripts/relaxtime/run_mott_phase_scan.jl",
        "config_path" => opts.config_path,
        "config_hash" => config_hash,
        "output_csv" => out_csv,
        "schema_version" => "v1",
    )

    _write_json(joinpath(opts.outdir, "effective_config.json"), effective)
    _write_json(joinpath(opts.outdir, "run_manifest.json"), manifest)
    return run_id
end

function main()
    opts, _cfg = _build_options(ARGS)

    mkpath(opts.outdir)
    out_csv = joinpath(opts.outdir, opts.output_csv)

    if opts.overwrite && isfile(out_csv)
        rm(out_csv)
    end

    key_cols = ["T_MeV", "muB_MeV", "xi"]
    existing_keys = (isfile(out_csv) && opts.resume && !opts.overwrite) ?
        ScanCSV.read_existing_keys(out_csv, key_cols) : Set{Tuple{Vararg{Float64}}}()

    run_id = _write_run_artifacts(opts, _cfg, out_csv)

    cols = [
        "run_id",
        "T_MeV", "muB_MeV", "xi",
        "M_pi", "M_K", "Gamma_pi", "Gamma_K",
        "residual_pi", "residual_K",
        "root_quality_pi", "root_quality_K",
        "selected_method_pi", "selected_method_K",
        "m_u", "m_d", "m_s",
        "status", "error_code", "error_message", "timestamp_utc",
    ]

    is_new = !isfile(out_csv)
    open(out_csv, is_new ? "w" : "a") do io
        if is_new
            ScanCSV.write_metadata(io, Dict(
                "format" => "scan_csv_v1",
                "script" => "scripts/relaxtime/run_mott_phase_scan.jl",
                "x_label" => "Temperature",
                "x_unit" => "MeV",
                "y_unit.M_pi" => "fm^-1",
                "y_unit.M_K" => "fm^-1",
                "y_unit.Gamma_pi" => "fm^-1",
                "y_unit.Gamma_K" => "fm^-1",
                "y_unit.m_u" => "fm^-1",
                "y_unit.m_d" => "fm^-1",
                "y_unit.m_s" => "fm^-1",
            ))
            ScanCSV.write_header(io, cols)
        end

        for xi in opts.xi_list
            meson_seed_state = nothing
            mixed_seed_tracking_state = nothing

            T = opts.T_min_MeV
            while T <= opts.T_max_MeV + 1e-9
                key = (Float64(T), Float64(opts.muB_MeV), Float64(xi))
                if key in existing_keys
                    T += opts.T_step_MeV
                    continue
                end

                timestamp = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ")
                row = Dict{String,Any}(
                    "run_id" => run_id,
                    "T_MeV" => T,
                    "muB_MeV" => opts.muB_MeV,
                    "xi" => xi,
                    "M_pi" => NaN,
                    "M_K" => NaN,
                    "Gamma_pi" => NaN,
                    "Gamma_K" => NaN,
                    "residual_pi" => NaN,
                    "residual_K" => NaN,
                    "root_quality_pi" => "",
                    "root_quality_K" => "",
                    "selected_method_pi" => "",
                    "selected_method_K" => "",
                    "m_u" => NaN,
                    "m_d" => NaN,
                    "m_s" => NaN,
                    "status" => "ok",
                    "error_code" => "",
                    "error_message" => "",
                    "timestamp_utc" => timestamp,
                )

                try
                    T_fm = T / ħc_MeV_fm
                    mu_fm = (opts.muB_MeV / ħc_MeV_fm) / 3.0
                    res = solve_gap_and_meson_point(
                        T_fm,
                        mu_fm;
                        xi=xi,
                        mesons=(:pi, :K),
                        meson_seed_state=meson_seed_state,
                        mixed_seed_tracking_state=mixed_seed_tracking_state,
                        mixed_branch_align=:identity_track_label_output,
                        p_num=opts.p_num,
                        t_num=opts.t_num,
                        solver_kwargs=(; iterations=opts.max_iter),
                        mass_kwargs=(; iterations=opts.max_iter),
                    )

                    qp = res.quark_params
                    mpi = res.meson_results[:pi]
                    mk = res.meson_results[:K]

                    row["M_pi"] = mpi.mass
                    row["M_K"] = mk.mass
                    row["Gamma_pi"] = mpi.gamma
                    row["Gamma_K"] = mk.gamma
                    row["residual_pi"] = mpi.residual
                    row["residual_K"] = mk.residual
                    row["root_quality_pi"] = String(mpi.root_quality)
                    row["root_quality_K"] = String(mk.root_quality)
                    row["selected_method_pi"] = String(mpi.root_diagnostics.selected_method)
                    row["selected_method_K"] = String(mk.root_diagnostics.selected_method)
                    row["m_u"] = qp.m.u
                    row["m_d"] = qp.m.d
                    row["m_s"] = qp.m.s

                    meson_seed_state = res.meson_seed_state
                    mixed_seed_tracking_state = res.mixed_seed_tracking
                catch e
                    row["status"] = "error"
                    row["error_code"] = "E_SOLVE"
                    row["error_message"] = sprint(showerror, e)
                end

                println(io, join([_fmt(get(row, c, "")) for c in cols], ','))
                push!(existing_keys, key)
                T += opts.T_step_MeV
            end
        end
    end

    println("Wrote Mott phase scan CSV: ", out_csv)
end

main()
