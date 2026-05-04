#!/usr/bin/env julia

using Dates
using JSON3
using Libdl

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_TAG = Dates.format(now(), "yyyymmdd_HHMMSS")
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "data", "outputs", "perf", "compile_breakdown")
const SYSIMAGE_PATH = joinpath(PROJECT_ROOT, "build", "JuliaRelaxTime." * Libdl.dlext)

Base.@kwdef mutable struct CliConfig
    tag::String = DEFAULT_TAG
    build_sysimage_if_missing::Bool = false
    modes::Set{Symbol} = Set([:cold_cli, :hot_in_process, :sysimage_cli])
end

function parse_bool(raw::AbstractString)
    lowered = lowercase(strip(raw))
    lowered in ("1", "true", "yes", "on") && return true
    lowered in ("0", "false", "no", "off") && return false
    throw(ArgumentError("invalid boolean value: $(raw)"))
end

function parse_args(args::Vector{String})
    cfg = CliConfig()
    for arg in args
        if arg in ("-h", "--help")
            print_help()
            exit(0)
        elseif startswith(arg, "--tag=")
            cfg.tag = split(arg, "="; limit=2)[2]
        elseif startswith(arg, "--build-sysimage-if-missing=")
            cfg.build_sysimage_if_missing = parse_bool(split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--modes=")
            raw = split(arg, "="; limit=2)[2]
            cfg.modes = Set(Symbol.(filter(!isempty, strip.(split(raw, ",")))))
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end
    for mode in cfg.modes
        mode in (:cold_cli, :hot_in_process, :sysimage_cli) || throw(ArgumentError("unsupported mode: $(mode)"))
    end
    return cfg
end

function print_help()
    println("""
Compile-time breakdown runner for stable scan entrypoints

Usage:
  julia --project=. scripts/perf/pnjl/compile_time_breakdown.jl [options]

Options:
  --tag=yyyymmdd_HHMMSS                 output directory suffix
  --build-sysimage-if-missing=true     build repository sysimage before comparison
  --modes=cold_cli,hot_in_process      choose subset from cold_cli, hot_in_process, sysimage_cli
  --help                               show help
""")
end

function ensure_sysimage!(cfg::CliConfig)
    if isfile(SYSIMAGE_PATH)
        return SYSIMAGE_PATH
    end
    cfg.build_sysimage_if_missing || return nothing
    build_script = joinpath(PROJECT_ROOT, "scripts", "dev", "build_sysimage.jl")
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(build_script)`
    run(Cmd(cmd; dir=PROJECT_ROOT))
    return isfile(SYSIMAGE_PATH) ? SYSIMAGE_PATH : nothing
end

function cli_command(kind::Symbol, output_path::String; with_sysimage::Union{Nothing, String}=nothing)
    cmd = String[joinpath(Sys.BINDIR, Base.julia_exename())]
    if with_sysimage !== nothing
        push!(cmd, "--sysimage=$(with_sysimage)")
    end
    append!(cmd, [
        "--project=.",
        "scripts/models/run_unified_scan.jl",
        "scan",
        String(kind),
        "--model_kind=PNJL",
        "--T_values=150",
        "--xi_values=0.0",
        "--output_path=$(output_path)",
        "--overwrite=true",
    ])
    if kind === :tmu
        push!(cmd, "--mu_values=0,100")
    elseif kind === :trho
        push!(cmd, "--rho_values=0.1,0.2")
    else
        throw(ArgumentError("unsupported kind: $(kind)"))
    end
    return cmd
end

function run_cli_measurement(kind::Symbol, stage_dir::String; with_sysimage::Union{Nothing, String}=nothing)
    output_path = joinpath(stage_dir, "$(String(kind)).csv")
    log_path = joinpath(stage_dir, "command.log")
    mkpath(stage_dir)
    command_vec = cli_command(kind, output_path; with_sysimage=with_sysimage)
    elapsed_ms = open(log_path, "w") do io
        println(io, "cwd=$(PROJECT_ROOT)")
        println(io, "command=$(join(command_vec, ' '))")
        println(io, "started_at=$(Dates.format(now(), Dates.ISODateTimeFormat))")
        cmd = Cmd(Cmd(command_vec); dir=PROJECT_ROOT)
        t0 = time_ns()
        run(pipeline(cmd; stdout=io, stderr=io))
        return (time_ns() - t0) / 1.0e6
    end
    point_count = 2
    return (
        mode = with_sysimage === nothing ? "cold_cli" : "sysimage_cli",
        kind = String(kind),
        point_count = point_count,
        total_ms = elapsed_ms,
        mean_ms_per_point = elapsed_ms / point_count,
        output_path = output_path,
        output_exists = isfile(output_path),
        log_path = log_path,
        command = command_vec,
    )
end

function include_models()
    if !isdefined(@__MODULE__, :Constants_PNJL)
        Base.include(@__MODULE__, joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
    end
    if !isdefined(@__MODULE__, :Models)
        Base.include(@__MODULE__, joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
    end
    return getfield(@__MODULE__, :Models)
end

function hot_scan_kwargs(kind::Symbol, output_path::String)
    common = (
        model_kind = :PNJL,
        T_values = [150.0],
        xi_values = [0.0],
        output_path = output_path,
        overwrite = true,
    )
    if kind === :tmu
        return merge(common, (mu_values=[0.0, 100.0],))
    elseif kind === :trho
        return merge(common, (rho_values=[0.1, 0.2],))
    end
    throw(ArgumentError("unsupported kind: $(kind)"))
end

function run_hot_measurement(kind::Symbol, stage_dir::String)
    Models = include_models()
    mkpath(stage_dir)
    output_path = joinpath(stage_dir, "$(String(kind)).csv")
    kwargs = hot_scan_kwargs(kind, output_path)

    # Warm the exact in-process path once to absorb JIT and method specialization.
    Base.invokelatest(Models.run_scan_pipeline, kind; kwargs...)

    measured_output = joinpath(stage_dir, "$(String(kind))_measured.csv")
    measured_kwargs = hot_scan_kwargs(kind, measured_output)
    t0 = time_ns()
    stats = Base.invokelatest(Models.run_scan_pipeline, kind; measured_kwargs...)
    elapsed_ms = (time_ns() - t0) / 1.0e6
    point_count = 2
    return (
        mode = "hot_in_process",
        kind = String(kind),
        point_count = point_count,
        total_ms = elapsed_ms,
        mean_ms_per_point = elapsed_ms / point_count,
        output_path = measured_output,
        output_exists = isfile(measured_output),
        manifest_path = hasproperty(stats, :manifest_path) ? String(getproperty(stats, :manifest_path)) : "",
        command = ["Models.run_scan_pipeline($(kind); ...) after one warmup in the same Julia process"],
    )
end

function write_reports(run_dir::String, sysimage_used, results)
    machine = (
        os = string(Sys.KERNEL, " / ", Sys.MACHINE),
        cpu_threads = Sys.CPU_THREADS,
        julia_version = string(VERSION),
        sysimage_path = sysimage_used === nothing ? nothing : String(sysimage_used),
    )
    payload = (
        generated_at = Dates.format(now(), Dates.ISODateTimeFormat),
        machine = machine,
        interpretation = (
            cold_cli = "includes Julia process startup, project activation, JIT, script execution, and output writing",
            hot_in_process = "same process, exact scan pipeline path, one unmeasured warmup then one measured run",
            sysimage_cli = "fresh Julia process with repository sysimage; still includes process startup and script execution",
        ),
        results = results,
    )
    json_path = joinpath(run_dir, "compile_breakdown.json")
    md_path = joinpath(run_dir, "compile_breakdown.md")
    open(json_path, "w") do io
        JSON3.pretty(io, payload)
    end
    open(md_path, "w") do io
        println(io, "# Compile-Time Breakdown")
        println(io)
        println(io, "- Generated at: `$(payload.generated_at)`")
        println(io, "- OS / arch: `$(machine.os)`")
        println(io, "- CPU threads: `$(machine.cpu_threads)`")
        println(io, "- Julia: `$(machine.julia_version)`")
        println(io, "- Sysimage: `$(something(machine.sysimage_path, "not used"))`")
        println(io)
        println(io, "## Summary")
        println(io)
        println(io, "| Kind | Mode | Points | Total (ms) | Mean (ms/pt) |")
        println(io, "| --- | --- | ---: | ---: | ---: |")
        for row in results
            println(io, "| `$(row.kind)` | `$(row.mode)` | $(row.point_count) | $(round(row.total_ms, digits=1)) | $(round(row.mean_ms_per_point, digits=1)) |")
        end
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "- `cold_cli`: $(payload.interpretation.cold_cli)")
        println(io, "- `hot_in_process`: $(payload.interpretation.hot_in_process)")
        println(io, "- `sysimage_cli`: $(payload.interpretation.sysimage_cli)")
        println(io)
        println(io, "## Details")
        for row in results
            println(io)
            println(io, "### `$(row.kind)` / `$(row.mode)`")
            println(io)
            println(io, "- total: `$(round(row.total_ms, digits=1)) ms`")
            println(io, "- mean: `$(round(row.mean_ms_per_point, digits=1)) ms/pt`")
            println(io, "- output exists: `$(row.output_exists)`")
            if hasproperty(row, :log_path)
                println(io, "- log: `$(row.log_path)`")
            end
            if hasproperty(row, :manifest_path) && !isempty(row.manifest_path)
                println(io, "- manifest: `$(row.manifest_path)`")
            end
            println(io, "- command:")
            println(io, "  ```text")
            println(io, "  $(join(row.command, ' '))")
            println(io, "  ```")
        end
    end
    return json_path, md_path
end

function main(args=ARGS)
    cfg = parse_args(args)
    run_dir = joinpath(OUTPUT_ROOT, cfg.tag)
    mkpath(run_dir)

    sysimage = ensure_sysimage!(cfg)
    results = NamedTuple[]
    for kind in (:tmu, :trho)
        if :cold_cli in cfg.modes
            push!(results, run_cli_measurement(kind, joinpath(run_dir, String(kind), "cold_cli")))
            write_reports(run_dir, sysimage, results)
        end
        if :hot_in_process in cfg.modes
            push!(results, run_hot_measurement(kind, joinpath(run_dir, String(kind), "hot_in_process")))
            write_reports(run_dir, sysimage, results)
        end
        if :sysimage_cli in cfg.modes && sysimage !== nothing
            push!(results, run_cli_measurement(kind, joinpath(run_dir, String(kind), "sysimage_cli"); with_sysimage=sysimage))
            write_reports(run_dir, sysimage, results)
        end
    end

    json_path, md_path = write_reports(run_dir, sysimage, results)
    println("Saved compile breakdown reports:")
    println("  JSON -> $(json_path)")
    println("  Markdown -> $(md_path)")
end

main()
