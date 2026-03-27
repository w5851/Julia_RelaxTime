#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "config", "WorkflowConfig.jl"))

using .WorkflowConfig: normalize_merge_validate
using TOML
using Dates
using SHA

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

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_relaxtime_orchestrator.jl <transport|cross-section> [options]")
    println("Options:")
    println("  --config <path>             Config path (default config/workflows/relaxtime/default.toml)")
    println("  --outdir <path>             Output directory (default data/outputs/results/relaxtime/orchestrated/default)")
    println("  --resume                    Override resume=true")
    println("  --overwrite                 Override overwrite=true")
    println("  --fail-on-fallback          Reserve switch (recorded in manifest)")
    println("  --processes p1,p2,...       cross-section only; overrides process list")
    println("  -h, --help                  Show help")
end

function _deep_set!(cfg::Dict{String,Any}, keys::Vector{String}, value)
    cur = cfg
    for k in keys[1:end-1]
        if !haskey(cur, k) || !(cur[k] isa Dict{String,Any})
            cur[k] = Dict{String,Any}()
        end
        cur = cur[k]
    end
    cur[keys[end]] = value
    return cfg
end

function _canonicalize_processes(raw::AbstractString)
    vals = String[]
    for v in split(raw, ',')
        s = strip(v)
        isempty(s) && continue
        push!(vals, s)
    end
    isempty(vals) && throw(ArgumentError("--processes cannot be empty"))
    return vals
end

function parse_args(args::Vector{String})
    isempty(args) && (print_usage(); exit(1))
    cmd = args[1]
    cmd in ("transport", "cross-section") || error("unknown subcommand: $(cmd)")

    opts = Dict{String,Any}(
        "config_path" => joinpath(PROJECT_ROOT, "config", "workflows", "relaxtime", "default.toml"),
        "outdir" => joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "orchestrated", "default"),
        "resume" => nothing,
        "overwrite" => nothing,
        "fail_on_fallback" => false,
        "processes" => nothing,
    )

    i = 2
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            i += 1
            return args[i]
        end

        if arg == "--config"
            opts["config_path"] = require_value()
        elseif arg == "--outdir"
            opts["outdir"] = require_value()
        elseif arg == "--resume"
            opts["resume"] = true
        elseif arg == "--overwrite"
            opts["overwrite"] = true
        elseif arg == "--fail-on-fallback"
            opts["fail_on_fallback"] = true
        elseif arg == "--processes"
            opts["processes"] = _canonicalize_processes(require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return cmd, opts
end

function _build_cli_cfg(cmd::String, opts::Dict{String,Any})
    cli = Dict{String,Any}()
    if opts["resume"] !== nothing
        _deep_set!(cli, ["scan", "transport", "resume"], Bool(opts["resume"]))
    end
    if opts["overwrite"] !== nothing
        _deep_set!(cli, ["scan", "transport", "overwrite"], Bool(opts["overwrite"]))
    end
    if cmd == "cross-section" && opts["processes"] !== nothing
        _deep_set!(cli, ["scan", "cross_section", "processes"], opts["processes"])
    end
    return cli
end

function run_orchestrator(cmd::String, opts::Dict{String,Any})
    default_cfg = TOML.parsefile(joinpath(PROJECT_ROOT, "config", "workflows", "relaxtime", "default.toml"))
    toml_cfg = TOML.parsefile(String(opts["config_path"]))
    aliases = TOML.parsefile(joinpath(PROJECT_ROOT, "config", "workflows", "relaxtime", "schema", "aliases_v1.toml"))
    cli_cfg = _build_cli_cfg(cmd, opts)

    merged = normalize_merge_validate(default_cfg, toml_cfg, cli_cfg, aliases)
    effective = merged.effective

    outdir = String(opts["outdir"])
    mkpath(outdir)

    effective_json = joinpath(outdir, "effective_config.json")
    _write_json(effective_json, effective)

    cfg_hash = bytes2hex(sha256(_to_json(effective)))
    run_id = string("relaxtime-orch-", Dates.format(now(UTC), "yyyymmddTHHMMSS"), "-", first(cfg_hash, 8))
    manifest = Dict{String,Any}(
        "schema_version" => get(effective, "schema_version", "v1"),
        "config_path" => String(opts["config_path"]),
        "config_hash" => cfg_hash,
        "run_id" => run_id,
        "subcommand" => cmd,
        "timestamp_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS") * "Z",
        "fail_on_fallback" => Bool(opts["fail_on_fallback"]),
        "fallback_events" => Any[],
        "trace" => String.(merged.trace),
    )
    _write_json(joinpath(outdir, "run_manifest.json"), manifest)

    println("[orchestrator] command=$(cmd) outdir=$(outdir) run_id=$(run_id)")
end

function main()
    cmd, opts = parse_args(copy(ARGS))
    run_orchestrator(cmd, opts)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
