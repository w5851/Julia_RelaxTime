#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
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

@inline function _cmd_string_to_symbol(cmd::String)
    if cmd == "transport"
        return :transport
    end
    if cmd == "cross-section"
        return :cross_section
    end
    throw(ArgumentError("unknown subcommand: $(cmd)"))
end

function run_orchestrator(cmd::String, opts::Dict{String,Any})
    cmd_symbol = _cmd_string_to_symbol(cmd)
    result = Main.Models.run_relaxtime_orchestrator_pipeline(
        cmd_symbol;
        config_path=String(opts["config_path"]),
        outdir=String(opts["outdir"]),
        processes=opts["processes"],
        resume=opts["resume"],
        overwrite=opts["overwrite"],
        fail_on_fallback=Bool(opts["fail_on_fallback"]),
    )
    println("[orchestrator] command=$(cmd) outdir=$(result.outdir) run_id=$(result.run_id)")
end

function main()
    cmd, opts = parse_args(copy(ARGS))
    run_orchestrator(cmd, opts)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
