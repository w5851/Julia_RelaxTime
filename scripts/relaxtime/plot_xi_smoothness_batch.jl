#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PLOT_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "plot_scan_csv.py")
const DEFAULT_FIG_ROOT = joinpath(PROJECT_ROOT, "data", "outputs", "figures", "relaxtime", "plan_c")
const DEFAULT_YS = "tau_u,tau_s,eta_over_s,zeta_over_s,sigma_over_T"

using JSON3

struct Options
    manifest::String
    fig_root::String
end

@inline function _is_subpath(child::String, parent::String)
    c = normpath(abspath(child))
    p = normpath(abspath(parent))
    rel = try
        relpath(c, p)
    catch
        return false
    end
    return !(rel == ".." || startswith(rel, "..") || startswith(rel, "../") || startswith(rel, "..\\"))
end

function _usage()
    println("Usage: julia --project=. scripts/relaxtime/plot_xi_smoothness_batch.jl --manifest <run_manifest.json> [--fig-root <dir>]")
    println("Options:")
    println("  --manifest <path>   batch manifest json (required)")
    println("  --fig-root <path>   figure root dir (default data/outputs/figures/relaxtime/plan_c)")
    println("  -h, --help          show help")
end

function _parse_args(args::Vector{String})
    manifest = nothing
    fig_root = DEFAULT_FIG_ROOT
    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $arg"))
            i += 1
            return args[i]
        end

        if arg == "--manifest"
            manifest = require_value()
        elseif arg == "--fig-root"
            fig_root = require_value()
        elseif arg in ("-h", "--help")
            _usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end

    manifest === nothing && throw(ArgumentError("--manifest is required"))
    return Options(normpath(abspath(String(manifest))), normpath(abspath(String(fig_root))))
end

function _resolve_result_csv(path_like::AbstractString, manifest_dir::String)
    candidate = String(path_like)
    roots = (normpath(abspath(manifest_dir)), normpath(abspath(PROJECT_ROOT)))

    if isabspath(candidate)
        resolved = normpath(candidate)
        any(root -> _is_subpath(resolved, root), roots) || throw(ArgumentError("result_csv is outside allowed roots (manifest_dir/project_root): $(resolved)"))
        return resolved
    end

    rel_from_manifest = normpath(joinpath(manifest_dir, candidate))
    if isfile(rel_from_manifest)
        any(root -> _is_subpath(rel_from_manifest, root), roots) || throw(ArgumentError("result_csv is outside allowed roots (manifest_dir/project_root): $(rel_from_manifest)"))
        return rel_from_manifest
    end

    rel_from_repo = normpath(joinpath(PROJECT_ROOT, candidate))
    any(root -> _is_subpath(rel_from_repo, root), roots) || throw(ArgumentError("result_csv is outside allowed roots (manifest_dir/project_root): $(rel_from_repo)"))
    return rel_from_repo
end

function _python_bin()
    py = Sys.which("python")
    py !== nothing && return String(py)
    py3 = Sys.which("python3")
    py3 !== nothing && return String(py3)
    throw(ArgumentError("python interpreter not found (need python or python3 in PATH)"))
end

function _sample_id_from_entry(sample::Any, idx::Int)
    raw = haskey(sample, "sample_id") ? sample["sample_id"] : string("sample_", idx)
    sid_raw = strip(String(raw))
    isempty(sid_raw) && return string("sample_", idx)
    sid = join([c for c in sid_raw if isletter(c) || isdigit(c) || c == '_' || c == '-'])
    isempty(sid) && return string("sample_", idx)
    return sid
end

function _run_plot_for_sample(py::String, result_csv::String, out_dir::String)
    cmd = Cmd([
        py,
        PLOT_SCRIPT,
        "--mode", "lines",
        "--csv", result_csv,
        "--x", "xi",
        "--ys", DEFAULT_YS,
        "--out-dir", out_dir,
        "--check",
    ])
    run(Cmd(cmd; dir=PROJECT_ROOT))
end

function main(args::Vector{String}=copy(ARGS))
    opts = _parse_args(args)
    isfile(PLOT_SCRIPT) || throw(ArgumentError("plot script not found: $PLOT_SCRIPT"))
    isfile(opts.manifest) || throw(ArgumentError("manifest not found: $(opts.manifest)"))

    manifest_dir = dirname(opts.manifest)
    payload = JSON3.read(read(opts.manifest, String))
    haskey(payload, "samples") || throw(ArgumentError("manifest missing samples"))

    samples = payload["samples"]
    mkpath(opts.fig_root)
    py = _python_bin()

    total = 0
    plotted = 0
    skipped = 0
    for (idx, sample) in enumerate(samples)
        total += 1
        sid = _sample_id_from_entry(sample, idx)
        if !haskey(sample, "result_csv")
            @warn "skip sample without result_csv" sample_id=sid
            skipped += 1
            continue
        end
        result_csv = _resolve_result_csv(String(sample["result_csv"]), manifest_dir)
        if !isfile(result_csv)
            @warn "skip sample with missing result_csv" sample_id=sid result_csv=result_csv
            skipped += 1
            continue
        end

        sample_fig_dir = joinpath(opts.fig_root, sid)
        mkpath(sample_fig_dir)
        _run_plot_for_sample(py, result_csv, sample_fig_dir)
        plotted += 1
    end

    println("xi smoothness plotting done: total=$(total), plotted=$(plotted), skipped=$(skipped)")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
