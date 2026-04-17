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
    include_mixed::Bool
    force_global_fallback::Bool
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

@inline function _norm_slash(path::AbstractString)
    return replace(String(path), '\\' => '/')
end

function _repo_relpath(path::AbstractString)
    abs_path = normpath(abspath(String(path)))
    rel = try
        relpath(abs_path, PROJECT_ROOT)
    catch
        nothing
    end
    if rel !== nothing
        rel_s = String(rel)
        if rel_s == "." || !(startswith(rel_s, "..") || startswith(rel_s, string("..", Base.Filesystem.path_separator)))
            return _norm_slash(rel_s)
        end
        return _norm_slash(rel_s)
    end
    return _norm_slash(abs_path)
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
    println("  --include-mixed      Also compute/output mixed mesons (eta, eta_prime, sigma, sigma_prime)")
    println("  --force-global-fallback  Force global fallback path in meson workflow (experimental)")
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
                "include_mixed" => false,
                "force_global_fallback" => false,
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
        elseif arg == "--include-mixed"
            cli["include_mixed"] = true
        elseif arg == "--force-global-fallback"
            cli["force_global_fallback"] = true
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
    include_mixed = Bool(get(mp, "include_mixed", false))
    force_global_fallback = Bool(get(mp, "force_global_fallback", false))

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
        include_mixed,
        force_global_fallback,
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
                "include_mixed" => opts.include_mixed,
                "force_global_fallback" => opts.force_global_fallback,
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
        "config_path" => _repo_relpath(opts.config_path),
        "config_hash" => config_hash,
        "output_csv" => _repo_relpath(out_csv),
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
        "governance_candidate_count_pi", "governance_candidate_count_K",
        "governance_selection_reason_pi", "governance_selection_reason_K",
        "second_pass_triggered_pi", "second_pass_triggered_K",
        "second_pass_candidate_count_pi", "second_pass_candidate_count_K",
        "m_u", "m_d", "m_s",
        "status", "error_code", "error_message", "timestamp_utc",
    ]

    if opts.include_mixed
        append!(cols, [
            "M_eta", "M_eta_prime", "M_sigma", "M_sigma_prime",
            "Gamma_eta", "Gamma_eta_prime", "Gamma_sigma", "Gamma_sigma_prime",
            "residual_eta", "residual_eta_prime", "residual_sigma", "residual_sigma_prime",
            "root_quality_eta", "root_quality_eta_prime", "root_quality_sigma", "root_quality_sigma_prime",
            "selected_method_eta", "selected_method_eta_prime", "selected_method_sigma", "selected_method_sigma_prime",
            "governance_candidate_count_eta", "governance_candidate_count_eta_prime", "governance_candidate_count_sigma", "governance_candidate_count_sigma_prime",
            "governance_selection_reason_eta", "governance_selection_reason_eta_prime", "governance_selection_reason_sigma", "governance_selection_reason_sigma_prime",
            "second_pass_triggered_eta", "second_pass_triggered_eta_prime", "second_pass_triggered_sigma", "second_pass_triggered_sigma_prime",
            "second_pass_candidate_count_eta", "second_pass_candidate_count_eta_prime", "second_pass_candidate_count_sigma", "second_pass_candidate_count_sigma_prime",
        ])
    end

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
                "y_unit.M_eta" => "fm^-1",
                "y_unit.M_eta_prime" => "fm^-1",
                "y_unit.M_sigma" => "fm^-1",
                "y_unit.M_sigma_prime" => "fm^-1",
                "y_unit.Gamma_eta" => "fm^-1",
                "y_unit.Gamma_eta_prime" => "fm^-1",
                "y_unit.Gamma_sigma" => "fm^-1",
                "y_unit.Gamma_sigma_prime" => "fm^-1",
                "y_unit.m_u" => "fm^-1",
                "y_unit.m_d" => "fm^-1",
                "y_unit.m_s" => "fm^-1",
            ))
            ScanCSV.write_header(io, cols)
        end

        for xi in opts.xi_list
            equilibrium_seed_state = nothing
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
                    "governance_candidate_count_pi" => 0,
                    "governance_candidate_count_K" => 0,
                    "governance_selection_reason_pi" => "",
                    "governance_selection_reason_K" => "",
                    "second_pass_triggered_pi" => false,
                    "second_pass_triggered_K" => false,
                    "second_pass_candidate_count_pi" => 0,
                    "second_pass_candidate_count_K" => 0,
                    "m_u" => NaN,
                    "m_d" => NaN,
                    "m_s" => NaN,
                    "status" => "ok",
                    "error_code" => "",
                    "error_message" => "",
                        "timestamp_utc" => timestamp,
                    )

                    if opts.include_mixed
                        row["M_eta"] = NaN
                        row["M_eta_prime"] = NaN
                        row["M_sigma"] = NaN
                        row["M_sigma_prime"] = NaN
                        row["Gamma_eta"] = NaN
                        row["Gamma_eta_prime"] = NaN
                        row["Gamma_sigma"] = NaN
                        row["Gamma_sigma_prime"] = NaN
                        row["residual_eta"] = NaN
                        row["residual_eta_prime"] = NaN
                        row["residual_sigma"] = NaN
                        row["residual_sigma_prime"] = NaN
                        row["root_quality_eta"] = ""
                        row["root_quality_eta_prime"] = ""
                        row["root_quality_sigma"] = ""
                        row["root_quality_sigma_prime"] = ""
                        row["selected_method_eta"] = ""
                        row["selected_method_eta_prime"] = ""
                        row["selected_method_sigma"] = ""
                        row["selected_method_sigma_prime"] = ""
                        row["governance_candidate_count_eta"] = 0
                        row["governance_candidate_count_eta_prime"] = 0
                        row["governance_candidate_count_sigma"] = 0
                        row["governance_candidate_count_sigma_prime"] = 0
                        row["governance_selection_reason_eta"] = ""
                        row["governance_selection_reason_eta_prime"] = ""
                        row["governance_selection_reason_sigma"] = ""
                        row["governance_selection_reason_sigma_prime"] = ""
                        row["second_pass_triggered_eta"] = false
                        row["second_pass_triggered_eta_prime"] = false
                        row["second_pass_triggered_sigma"] = false
                        row["second_pass_triggered_sigma_prime"] = false
                        row["second_pass_candidate_count_eta"] = 0
                        row["second_pass_candidate_count_eta_prime"] = 0
                        row["second_pass_candidate_count_sigma"] = 0
                        row["second_pass_candidate_count_sigma_prime"] = 0
                    end

                try
                    T_fm = T / ħc_MeV_fm
                    mu_fm = (opts.muB_MeV / ħc_MeV_fm) / 3.0
                    mesons = opts.include_mixed ? (:pi, :K, :eta, :eta_prime, :sigma, :sigma_prime) : (:pi, :K)
                    res = solve_gap_and_meson_point(
                        T_fm,
                        mu_fm;
                        xi=xi,
                        seed_state=(equilibrium_seed_state === nothing ? Main.Models.HADRON_SEED_5 : equilibrium_seed_state),
                        mesons=mesons,
                        meson_seed_state=meson_seed_state,
                        mixed_seed_tracking_state=mixed_seed_tracking_state,
                        mixed_branch_align=:strict_sign_binding,
                        p_num=opts.p_num,
                        t_num=opts.t_num,
                        solver_kwargs=(; iterations=opts.max_iter),
                        mass_kwargs=(; iterations=opts.max_iter),
                        force_global_fallback=opts.force_global_fallback,
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
                    row["governance_candidate_count_pi"] = Int(getproperty(mpi.root_diagnostics, :governance_candidate_count))
                    row["governance_candidate_count_K"] = Int(getproperty(mk.root_diagnostics, :governance_candidate_count))
                    row["governance_selection_reason_pi"] = String(getproperty(mpi.root_diagnostics, :governance_selection_reason))
                    row["governance_selection_reason_K"] = String(getproperty(mk.root_diagnostics, :governance_selection_reason))
                    row["second_pass_triggered_pi"] = Bool(getproperty(mpi.root_diagnostics, :second_pass_triggered))
                    row["second_pass_triggered_K"] = Bool(getproperty(mk.root_diagnostics, :second_pass_triggered))
                    row["second_pass_candidate_count_pi"] = Int(getproperty(mpi.root_diagnostics, :second_pass_candidate_count))
                    row["second_pass_candidate_count_K"] = Int(getproperty(mk.root_diagnostics, :second_pass_candidate_count))
                    row["m_u"] = qp.m.u
                    row["m_d"] = qp.m.d
                    row["m_s"] = qp.m.s

                    if opts.include_mixed
                        eta = res.meson_results[:eta]
                        eta_prime = res.meson_results[:eta_prime]
                        sigma = res.meson_results[:sigma]
                        sigma_prime = res.meson_results[:sigma_prime]

                        row["M_eta"] = eta.mass
                        row["M_eta_prime"] = eta_prime.mass
                        row["M_sigma"] = sigma.mass
                        row["M_sigma_prime"] = sigma_prime.mass
                        row["Gamma_eta"] = eta.gamma
                        row["Gamma_eta_prime"] = eta_prime.gamma
                        row["Gamma_sigma"] = sigma.gamma
                        row["Gamma_sigma_prime"] = sigma_prime.gamma
                        row["residual_eta"] = eta.residual
                        row["residual_eta_prime"] = eta_prime.residual
                        row["residual_sigma"] = sigma.residual
                        row["residual_sigma_prime"] = sigma_prime.residual
                        row["root_quality_eta"] = String(eta.root_quality)
                        row["root_quality_eta_prime"] = String(eta_prime.root_quality)
                        row["root_quality_sigma"] = String(sigma.root_quality)
                        row["root_quality_sigma_prime"] = String(sigma_prime.root_quality)
                        row["selected_method_eta"] = String(eta.root_diagnostics.selected_method)
                        row["selected_method_eta_prime"] = String(eta_prime.root_diagnostics.selected_method)
                        row["selected_method_sigma"] = String(sigma.root_diagnostics.selected_method)
                        row["selected_method_sigma_prime"] = String(sigma_prime.root_diagnostics.selected_method)
                        row["governance_candidate_count_eta"] = Int(getproperty(eta.root_diagnostics, :governance_candidate_count))
                        row["governance_candidate_count_eta_prime"] = Int(getproperty(eta_prime.root_diagnostics, :governance_candidate_count))
                        row["governance_candidate_count_sigma"] = Int(getproperty(sigma.root_diagnostics, :governance_candidate_count))
                        row["governance_candidate_count_sigma_prime"] = Int(getproperty(sigma_prime.root_diagnostics, :governance_candidate_count))
                        row["governance_selection_reason_eta"] = String(getproperty(eta.root_diagnostics, :governance_selection_reason))
                        row["governance_selection_reason_eta_prime"] = String(getproperty(eta_prime.root_diagnostics, :governance_selection_reason))
                        row["governance_selection_reason_sigma"] = String(getproperty(sigma.root_diagnostics, :governance_selection_reason))
                        row["governance_selection_reason_sigma_prime"] = String(getproperty(sigma_prime.root_diagnostics, :governance_selection_reason))
                        row["second_pass_triggered_eta"] = Bool(getproperty(eta.root_diagnostics, :second_pass_triggered))
                        row["second_pass_triggered_eta_prime"] = Bool(getproperty(eta_prime.root_diagnostics, :second_pass_triggered))
                        row["second_pass_triggered_sigma"] = Bool(getproperty(sigma.root_diagnostics, :second_pass_triggered))
                        row["second_pass_triggered_sigma_prime"] = Bool(getproperty(sigma_prime.root_diagnostics, :second_pass_triggered))
                        row["second_pass_candidate_count_eta"] = Int(getproperty(eta.root_diagnostics, :second_pass_candidate_count))
                        row["second_pass_candidate_count_eta_prime"] = Int(getproperty(eta_prime.root_diagnostics, :second_pass_candidate_count))
                        row["second_pass_candidate_count_sigma"] = Int(getproperty(sigma.root_diagnostics, :second_pass_candidate_count))
                        row["second_pass_candidate_count_sigma_prime"] = Int(getproperty(sigma_prime.root_diagnostics, :second_pass_candidate_count))
                    end

                    equilibrium_seed_state = collect(res.equilibrium.x_state)
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
