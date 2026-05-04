#!/usr/bin/env julia

using Dates
using JSON3
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "data", "outputs", "perf", "stable_entrypoints")

Base.@kwdef struct BaselineSpec
    key::String
    label::String
    category::String
    description::String
    workload_shape::String
    point_count::Int
    command_args::Vector{String}
    output_filename::String
end

Base.@kwdef mutable struct CliConfig
    samples::Int = 1
    warmups::Int = 0
    only_keys::Union{Nothing, Set{String}} = nothing
    tag::String = Dates.format(now(), "yyyymmdd_HHMMSS")
    keep_artifacts::Bool = true
end

const BASELINE_SPECS = [
    BaselineSpec(
        key = "scan_tmu_minimal",
        label = "Unified scan tmu minimal grid",
        category = "stable-cli",
        description = "End-to-end wall-clock baseline for the Models unified T-mu scan entrypoint.",
        workload_shape = "len(T)=1, len(mu)=2, len(xi)=1 => 1*2*1",
        point_count = 2,
        command_args = [
            "--project=.",
            "scripts/models/run_unified_scan.jl",
            "scan",
            "tmu",
            "--model_kind=PNJL",
            "--T_values=150",
            "--mu_values=0,100",
            "--xi_values=0.0",
            "--output_path=__OUTPUT_PATH__",
            "--overwrite=true",
        ],
        output_filename = "tmu.csv",
    ),
    BaselineSpec(
        key = "scan_trho_minimal",
        label = "Unified scan trho minimal grid",
        category = "stable-cli",
        description = "End-to-end wall-clock baseline for the Models unified T-rho scan entrypoint.",
        workload_shape = "len(T)=1, len(rho)=2, len(xi)=1 => 1*2*1",
        point_count = 2,
        command_args = [
            "--project=.",
            "scripts/models/run_unified_scan.jl",
            "scan",
            "trho",
            "--model_kind=PNJL",
            "--T_values=150",
            "--rho_values=0.1,0.2",
            "--xi_values=0.0",
            "--output_path=__OUTPUT_PATH__",
            "--overwrite=true",
        ],
        output_filename = "trho.csv",
    ),
]

function parse_args(args::Vector{String})
    cfg = CliConfig()
    for arg in args
        if arg in ("-h", "--help")
            print_help()
            exit(0)
        elseif startswith(arg, "--samples=")
            cfg.samples = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--warmups=")
            cfg.warmups = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--only=")
            raw = split(arg, "="; limit=2)[2]
            cfg.only_keys = Set(filter(!isempty, strip.(split(raw, ","))))
        elseif startswith(arg, "--tag=")
            cfg.tag = split(arg, "="; limit=2)[2]
        elseif startswith(arg, "--keep-artifacts=")
            cfg.keep_artifacts = parse(Bool, lowercase(split(arg, "="; limit=2)[2]))
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end
    cfg.samples > 0 || throw(ArgumentError("--samples must be > 0"))
    cfg.warmups >= 0 || throw(ArgumentError("--warmups must be >= 0"))
    return cfg
end

function print_help()
    println("""
Stable CLI entrypoint baseline runner

Usage:
  julia --project=. scripts/perf/pnjl/stable_entrypoint_baseline.jl [options]

Options:
  --samples=1                measured runs per entrypoint
  --warmups=0                unrecorded warmup runs per entrypoint
  --only=scan_tmu_minimal    comma-separated spec keys
  --tag=yyyymmdd_HHMMSS      output directory suffix
  --keep-artifacts=true      keep per-run output CSV/logs
  --help                     show help
""")
end

function selected_specs(cfg::CliConfig)
    specs = BASELINE_SPECS
    if cfg.only_keys !== nothing
        specs = filter(spec -> spec.key in cfg.only_keys, specs)
    end
    isempty(specs) && throw(ArgumentError("no baseline specs selected"))
    return specs
end

function stats(values::Vector{Float64})
    sorted = sort(values)
    n = length(sorted)
    mid = cld(n, 2)
    median = isodd(n) ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2
    return (
        min_ms = minimum(sorted),
        median_ms = median,
        mean_ms = sum(sorted) / n,
        max_ms = maximum(sorted),
    )
end

function build_machine_info()
    return (
        os = string(Sys.KERNEL, " / ", Sys.MACHINE),
        cpu_threads = Sys.CPU_THREADS,
        julia_version = string(VERSION),
        project_root = PROJECT_ROOT,
        measurement_scope = "end-to-end wall-clock time including Julia process startup, project activation, script execution, and output writing",
    )
end

function resolve_command(spec::BaselineSpec, artifact_dir::String)
    output_path = joinpath(artifact_dir, spec.output_filename)
    args = map(spec.command_args) do token
        replace(token, "__OUTPUT_PATH__" => output_path)
    end
    julia_exe = joinpath(Sys.BINDIR, Base.julia_exename())
    return vcat([julia_exe], args), output_path
end

function run_once(spec::BaselineSpec, run_dir::String; sample_index::Int, stage::String)
    artifact_dir = joinpath(run_dir, spec.key, stage, "sample_$(lpad(sample_index, 2, '0'))")
    mkpath(artifact_dir)
    command_vec, output_path = resolve_command(spec, artifact_dir)
    log_path = joinpath(artifact_dir, "command.log")

    elapsed_ms = open(log_path, "w") do io
        println(io, "cwd=$(PROJECT_ROOT)")
        println(io, "command=$(join(command_vec, ' '))")
        println(io, "started_at=$(Dates.format(now(), Dates.ISODateTimeFormat))")
        cmd = Cmd(Cmd(command_vec); dir=PROJECT_ROOT)
        t0 = time_ns()
        run(pipeline(cmd; stdout=io, stderr=io))
        return (time_ns() - t0) / 1.0e6
    end

    return (
        elapsed_ms = elapsed_ms,
        output_path = output_path,
        output_exists = isfile(output_path),
        log_path = log_path,
    )
end

function maybe_cleanup!(cfg::CliConfig, spec_dir::String)
    cfg.keep_artifacts && return
    isdir(spec_dir) || return
    rm(spec_dir; recursive=true, force=true)
end

function run_spec(cfg::CliConfig, spec::BaselineSpec, run_dir::String)
    println()
    println("==> $(spec.label) [$(spec.key)]")
    println("    warmups=$(cfg.warmups) samples=$(cfg.samples)")
    println("    workload=$(spec.workload_shape), points=$(spec.point_count)")

    warmup_records = NamedTuple[]
    for idx in 1:cfg.warmups
        print("    warmup $idx/$((cfg.warmups)): ")
        record = run_once(spec, run_dir; sample_index=idx, stage="warmup")
        println(@sprintf("%.1f ms", record.elapsed_ms))
        push!(warmup_records, record)
    end

    sample_records = NamedTuple[]
    for idx in 1:cfg.samples
        print("    sample $idx/$((cfg.samples)): ")
        record = run_once(spec, run_dir; sample_index=idx, stage="measured")
        println(@sprintf("%.1f ms", record.elapsed_ms))
        push!(sample_records, record)
    end

    measured = [record.elapsed_ms for record in sample_records]
    summary = stats(measured)
    println(@sprintf("    median=%.1f ms mean=%.1f ms", summary.median_ms, summary.mean_ms))

    spec_dir = joinpath(run_dir, spec.key)
    return (
        key = spec.key,
        label = spec.label,
        category = spec.category,
        description = spec.description,
        workload_shape = spec.workload_shape,
        point_count = spec.point_count,
        command = resolve_command(spec, joinpath(spec_dir, "measured", "sample_01"))[1],
        output_filename = spec.output_filename,
        warmups = cfg.warmups,
        samples = cfg.samples,
        summary = summary,
        per_point = (
            min_ms = summary.min_ms / spec.point_count,
            median_ms = summary.median_ms / spec.point_count,
            mean_ms = summary.mean_ms / spec.point_count,
            max_ms = summary.max_ms / spec.point_count,
        ),
        measured_runs = sample_records,
        warmup_runs = warmup_records,
        artifact_dir = spec_dir,
    )
end

function write_json_report(path::String, payload)
    open(path, "w") do io
        JSON3.pretty(io, payload)
    end
end

function write_markdown_report(path::String, cfg::CliConfig, results, machine_info, run_dir::String)
    open(path, "w") do io
        println(io, "# Stable Entrypoint Performance Baseline")
        println(io)
        println(io, "- Generated at: `$(Dates.format(now(), Dates.ISODateTimeFormat))`")
        println(io, "- Tag: `$(cfg.tag)`")
        println(io, "- Output root: `$(run_dir)`")
        println(io, "- Samples per entrypoint: `$(cfg.samples)`")
        println(io, "- Warmups per entrypoint: `$(cfg.warmups)`")
        println(io)
        println(io, "## Machine Assumptions")
        println(io)
        println(io, "- OS / arch: `$(machine_info.os)`")
        println(io, "- CPU threads: `$(machine_info.cpu_threads)`")
        println(io, "- Julia: `$(machine_info.julia_version)`")
        println(io, "- Measurement scope: $(machine_info.measurement_scope)")
        println(io)
        println(io, "## Summary")
        println(io)
        println(io, "| Key | Points | Samples | Min (ms) | Median (ms) | Mean (ms) | Mean (ms/pt) | Max (ms) |")
        println(io, "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for result in results
            s = result.summary
            println(io, "| `$(result.key)` | $(result.point_count) | $(result.samples) | $(round(s.min_ms, digits=1)) | $(round(s.median_ms, digits=1)) | $(round(s.mean_ms, digits=1)) | $(round(result.per_point.mean_ms, digits=1)) | $(round(s.max_ms, digits=1)) |")
        end
        println(io)
        println(io, "## Entrypoint Details")
        for result in results
            println(io)
            println(io, "### `$(result.key)`")
            println(io)
            println(io, "- Label: $(result.label)")
            println(io, "- Category: `$(result.category)`")
            println(io, "- Description: $(result.description)")
            println(io, "- Workload shape: `$(result.workload_shape)`")
            println(io, "- Point count: `$(result.point_count)`")
            println(io, "- Command:")
            println(io, "  ```text")
            println(io, "  $(join(result.command, ' '))")
            println(io, "  ```")
            println(io, "- Artifact directory: `$(result.artifact_dir)`")
            println(io, "- Expected output filename: `$(result.output_filename)`")
            println(io, "- Aggregate timing:")
            println(io, "  - min: `$(round(result.summary.min_ms, digits=1)) ms`")
            println(io, "  - median: `$(round(result.summary.median_ms, digits=1)) ms`")
            println(io, "  - mean: `$(round(result.summary.mean_ms, digits=1)) ms`")
            println(io, "- Per-point timing:")
            println(io, "  - min: `$(round(result.per_point.min_ms, digits=1)) ms/pt`")
            println(io, "  - median: `$(round(result.per_point.median_ms, digits=1)) ms/pt`")
            println(io, "  - mean: `$(round(result.per_point.mean_ms, digits=1)) ms/pt`")
            println(io, "- Measured runs:")
            for (idx, record) in enumerate(result.measured_runs)
                println(io, "  - sample $(idx): `$(round(record.elapsed_ms, digits=1)) ms`, output_exists=`$(record.output_exists)`, log=`$(record.log_path)`")
            end
        end
        println(io)
        println(io, "## Reuse Template")
        println(io)
        println(io, "未来追加新的稳定入口时，沿用以下字段：")
        println(io)
        println(io, "1. `key` / `label` / `description`")
        println(io, "2. 真实 CLI 命令与输出工件路径")
        println(io, "3. 机器假设（OS、CPU threads、Julia 版本）")
        println(io, "4. `samples` / `warmups`")
        println(io, "5. `min/median/mean/max` wall-clock 时间")
        println(io, "6. 每次样本的日志与输出文件存在性")
    end
end

function main(args=ARGS)
    cfg = parse_args(args)
    specs = selected_specs(cfg)
    run_dir = joinpath(OUTPUT_ROOT, cfg.tag)
    mkpath(run_dir)

    machine_info = build_machine_info()
    results = map(spec -> run_spec(cfg, spec, run_dir), specs)

    if !cfg.keep_artifacts
        for result in results
            maybe_cleanup!(cfg, result.artifact_dir)
        end
    end

    json_path = joinpath(run_dir, "baseline_summary.json")
    md_path = joinpath(run_dir, "baseline_summary.md")
    payload = (
        generated_at = Dates.format(now(), Dates.ISODateTimeFormat),
        tag = cfg.tag,
        machine = machine_info,
        config = (
            samples = cfg.samples,
            warmups = cfg.warmups,
            keep_artifacts = cfg.keep_artifacts,
        ),
        results = results,
    )
    write_json_report(json_path, payload)
    write_markdown_report(md_path, cfg, results, machine_info, run_dir)

    println()
    println("Saved baseline reports:")
    println("  JSON -> $(json_path)")
    println("  Markdown -> $(md_path)")
end

main()
