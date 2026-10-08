#!/usr/bin/env julia

"""Repeatable timing probe for fixed ``(T, μB)`` charged-GBU points.

The benchmark separates process/JIT effects from steady-state work:

* the first background call is recorded as a cold-in-process sample;
* no-seed and seed calls are then repeated in the same Julia process;
* density has a separate first-call and repeated sample group.

Julia process startup is outside this script and is reported explicitly as
excluded. The density workload uses the current infinite-thermal kernel at a
small resolution and is a screening probe by default. When
``BENCH_CGBU_WORKLOAD=production`` is selected, the same entrypoint calls the
full ``Workflow.channel_density`` path and records the
local/Mott/topology/eta/q-order/tail gates; that workload is intentionally
separate because it is orders of magnitude slower.
"""

const CGBU_BENCH_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

isdefined(Main, :Models) || Base.include(Main, joinpath(CGBU_BENCH_ROOT, "src", "models", "Models.jl"))
isdefined(Main, :ChargedGBUResearchWorkflow) || Base.include(Main,
    joinpath(CGBU_BENCH_ROOT, "src", "models", "workflow_apps", "ChargedGBUResearchWorkflow.jl"))

module ChargedGBUContourBenchmark

using Dates
using JSON3
using Printf
using SHA
using TOML

include("charged_gbu_benchmark_utils.jl")
const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const Models = Main.Models
const Workflow = Main.ChargedGBUResearchWorkflow
const BenchmarkUtils = ChargedGBUBenchmarkUtils
const CHANNELS = Workflow.R.CHANNELS

_env_int(name, default) = parse(Int, get(ENV, name, string(default)))
_env_float(name, default) = parse(Float64, get(ENV, name, string(default)))
_env_string(name, default) = get(ENV, name, default)
_env_bool(name, default) = lowercase(get(ENV, name, default ? "true" : "false")) in
    ("1", "true", "yes", "on")
elapsed_seconds(start_ns::UInt64) = (Base.time_ns() - start_ns) / 1.0e9

function screening_settings()
    settings = (mesh=_env_int("BENCH_CGBU_MESH", 64),
        cut_nodes=_env_int("BENCH_CGBU_CUT_NODES", 32),
        tail_nodes=_env_int("BENCH_CGBU_TAIL_NODES", 32),
        omega_nodes=_env_int("BENCH_CGBU_OMEGA_NODES", 64),
        q_nodes=_env_int("BENCH_CGBU_Q_NODES", 8),
        qmax=_env_float("BENCH_CGBU_QMAX", 8.0),
        repeats=_env_int("BENCH_CGBU_REPEATS", 5),
        density_repeats=_env_int("BENCH_CGBU_DENSITY_REPEATS", 3),
        include_density=_env_bool("BENCH_CGBU_INCLUDE_DENSITY", true),
        workload=_env_string("BENCH_CGBU_WORKLOAD", "screening"),
        channels=Symbol.(split(_env_string("BENCH_CGBU_CHANNELS", "pi_plus,pi_minus,K_plus,K_minus"), ',')))
    settings.workload in ("screening", "production") || throw(ArgumentError("workload must be screening or production"))
    !isempty(settings.channels) && all(in(CHANNELS), settings.channels) &&
        length(unique(settings.channels)) == length(settings.channels) || throw(ArgumentError("invalid channels"))
    return BenchmarkUtils.validate_settings(settings)
end

parse_points(spec::AbstractString) = BenchmarkUtils.parse_points(spec)

default_points() = parse_points(_env_string("BENCH_CGBU_POINTS",
    "165.923:23.525,164.238:112.304,151.593:315.980," *
    "118.524:553.066,79.957:719.076,40.000:800.000"))

"""The BQS background contract is exactly `(T_MeV, muB_MeV)`."""
point(T_MeV, muB_MeV) = (T_MeV=Float64(T_MeV), muB_MeV=Float64(muB_MeV))

function background(model, T_MeV, muB_MeV, config; seed=nothing)
    return Workflow.background(model, point(T_MeV, muB_MeV), config; seed=seed)
end

_timed_call(f) = BenchmarkUtils.measure_call(f)

function density_probe(bg, channel, settings)
    started = Base.time_ns()
    q_values, q_weights = Workflow.R.gauleg(0.0, settings.qmax, settings.q_nodes)
    rows = NamedTuple[]
    for (q, weight) in zip(q_values, q_weights)
        kernel = Workflow.I.kernel(bg, channel, q; cut_nodes=settings.cut_nodes,
            split_inv_fm=max(36.0, q + 28.0 * bg.T + 2.0))
        profile = Workflow.P.profile(kernel; mesh=settings.mesh, tail_nodes=settings.tail_nodes)
        shell = Workflow.Y.shell(profile; nodes=settings.omega_nodes)
        push!(rows, (q_inv_fm=q, q_weight_inv_fm=weight,
            shell_contribution_inv_fm2=Float64(shell.density),
            root_count=Int(shell.root_count), output_finite=isfinite(shell.density)))
    end
    # shell.density already contains q²/(2π²); apply the outer q quadrature here.
    value = sum(row.q_weight_inv_fm * row.shell_contribution_inv_fm2 for row in rows)
    return (elapsed_s=elapsed_seconds(started), density_inv_fm3=value,
        output_finite=isfinite(value) && all(row.output_finite for row in rows),
        negative_contribution=isfinite(value) && value < 0, rows=rows,
        status="screening_evaluated", production_gates_passed=false,
        gate_scope="finite_shell_output_only;full_production_gates_omitted")
end

function channel_probe(bg, settings, config)
    results = Dict{String,Any}()
    for ch in settings.channels
        channel_started = Base.time_ns()
        results[String(ch)] = try
            if settings.workload == "screening"
                density_probe(bg, ch, settings)
            else
                start = Base.time_ns()
                result = Workflow.channel_density(bg, ch, config, _ -> nothing)
                (elapsed_s=elapsed_seconds(start), density_inv_fm3=result.density,
                    output_finite=isfinite(result.density), status=result.status,
                    negative_contribution=isfinite(result.density) && result.density < 0,
                    production_gates_passed=result.passed,
                    final_order=get(result, :final_order, nothing),
                    root_count=isempty(result.rows) ? nothing : maximum(r.root_count for r in result.rows))
            end
        catch err
            err isa InterruptException && rethrow()
            (elapsed_s=elapsed_seconds(channel_started), density_inv_fm3=NaN, output_finite=false,
                negative_contribution=false, production_gates_passed=false,
                status="evaluation_failed", reason=sprint(showerror, err))
        end
    end
    return results
end

function _tracked_source_paths()
    output = read(`git -C $PROJECT_ROOT ls-files -z -- src config scripts/analysis/relaxtime benchmark/relaxtime Project.toml Manifest.toml`, String)
    paths = String[]
    for line in split(output, '\0')
        isempty(strip(line)) && continue
        path = joinpath(PROJECT_ROOT, strip(line))
        isfile(path) && push!(paths, path)
    end
    for path in (joinpath(PROJECT_ROOT, "benchmark", "relaxtime",
            "bench_charged_gbu_contour_point.jl"),
        joinpath(PROJECT_ROOT, "benchmark", "relaxtime",
            "charged_gbu_benchmark_utils.jl"))
        isfile(path) && push!(paths, path)
    end
    return unique(paths)
end

function source_snapshot(output)
    output === nothing && return
    for path in _tracked_source_paths()
        target = joinpath(dirname(abspath(output)), "source_snapshot", relpath(path, PROJECT_ROOT))
        mkpath(dirname(target)); cp(path, target; force=true)
    end
end

function source_hashes()
    return Dict(replace(relpath(path, PROJECT_ROOT), '\\' => '/') =>
        bytes2hex(SHA.sha256(read(path))) for path in _tracked_source_paths())
end

function _git_metadata()
    return (git_head=readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`),
        github_sha=get(ENV, "GITHUB_SHA", nothing),
        github_run_id=get(ENV, "GITHUB_RUN_ID", nothing),
        github_run_attempt=get(ENV, "GITHUB_RUN_ATTEMPT", nothing),
        runner_os=get(ENV, "RUNNER_OS", nothing),
        runner_arch=get(ENV, "RUNNER_ARCH", nothing))
end

jsonsafe(x::AbstractFloat) = isfinite(x) ? x : nothing
jsonsafe(x::NamedTuple) = Dict(String(k) => jsonsafe(v) for (k, v) in pairs(x))
jsonsafe(x::AbstractDict) = Dict(String(k) => jsonsafe(v) for (k, v) in pairs(x))
jsonsafe(x::AbstractArray) = [jsonsafe(v) for v in x]
jsonsafe(x) = x

function write_json(path, payload)
    path === nothing && return nothing
    mkpath(dirname(abspath(path)))
    temp = path * ".partial"
    open(temp, "w") do io
        JSON3.pretty(io, jsonsafe(payload)); write(io, '\n')
    end
    mv(temp, path; force=true)
end

function _background_summary(sample)
    sample.success || return BenchmarkUtils.timing_record(sample)
    value = sample.value
    return merge(BenchmarkUtils.timing_record(sample),
        (residual=Float64(value.bg.residual), seed_length=length(value.seed),
            bg=value.bg, seed=value.seed, omega=value.omega))
end

function _density_summary(sample)
    sample.success || return BenchmarkUtils.timing_record(sample)
    value = sample.value
    return merge(BenchmarkUtils.timing_record(sample),
        (channels=value,))
end

function _record_background_sample(model, config, T_MeV, muB_MeV; seed=nothing)
    return _timed_call() do
        background(model, T_MeV, muB_MeV, config; seed=seed)
    end
end

function channel_timings(samples, channel, workload)
    records = map(samples) do sample
        sample.success && haskey(sample.value, String(channel)) || return merge(
            BenchmarkUtils.timing_record(sample), (success=false, elapsed_s=sample.elapsed_s))
        value = sample.value[String(channel)]
        ok = value.output_finite && (workload == "screening" || value.production_gates_passed)
        return (success=ok, elapsed_s=value.elapsed_s, timed_s=nothing, bytes=nothing,
            gc_s=nothing, compile_s=nothing, recompile_s=nothing,
            error=ok ? nothing : value.status)
    end
    return BenchmarkUtils.summarize_timed(records)
end

function _point_benchmark(model, config, settings, T_MeV, muB_MeV;
    initial_record=nothing, continuation_seed=nothing, process_first=false,
    include_density=false)
    cold = initial_record === nothing ?
        _record_background_sample(model, config, T_MeV, muB_MeV) : initial_record
    cold_record = cold.success ? cold.value : nothing

    no_seed_samples = [_record_background_sample(model, config, T_MeV, muB_MeV)
        for _ in 1:settings.repeats]
    seed_base = continuation_seed === nothing ?
        (cold_record === nothing ? nothing : cold_record.seed) : continuation_seed
    seed_first = seed_base === nothing ? nothing :
        _record_background_sample(model, config, T_MeV, muB_MeV; seed=seed_base)
    seed_samples = seed_first === nothing || !seed_first.success ?
        NamedTuple[] : [_record_background_sample(model, config, T_MeV, muB_MeV;
            seed=seed_base) for _ in 1:settings.repeats]

    density = nothing
    if include_density && cold_record !== nothing
        density_call = () -> channel_probe(cold_record.bg, settings, config)
        density_first = _timed_call(density_call)
        density_samples = [_timed_call(density_call) for _ in 1:settings.density_repeats]
        accepted_samples = filter(density_samples) do sample
            sample.success && all(v.production_gates_passed for v in values(sample.value))
        end
        density = (first=_density_summary(density_first),
            background_source="first_no_seed_solution_at_this_point",
            repeats=[_density_summary(sample) for sample in density_samples],
            summary=BenchmarkUtils.summarize_timed(density_samples),
            accepted_timing_summary=BenchmarkUtils.summarize_timed(accepted_samples),
            channel_timing_summaries=Dict(String(ch) => channel_timings(density_samples, ch,
                settings.workload) for ch in settings.channels),
            gate_scope=settings.workload == "screening" ? "screening_only;full_production_gates_omitted" : "full_smoke_production_gates")
    end

    next_seed = if seed_first !== nothing && seed_first.success
        seed_first.value.seed
    elseif cold_record !== nothing
        cold_record.seed
    else
        nothing
    end
    report = (T_MeV=T_MeV, muB_MeV=muB_MeV,
        first_no_seed=merge(_background_summary(cold), (process_first=process_first,)),
        seed_source=continuation_seed === nothing ? "same_point_first_solution" : "prior_point_seed",
        no_seed=BenchmarkUtils.summarize_timed(no_seed_samples),
        no_seed_samples=[_background_summary(sample) for sample in no_seed_samples],
        seed_first=seed_first === nothing ? nothing : _background_summary(seed_first),
        seed=BenchmarkUtils.summarize_timed(seed_samples),
        seed_samples=[_background_summary(sample) for sample in seed_samples],
        density=density)
    return (report=report, next_seed=next_seed)
end

function run_single(model, config, settings, T_MeV, muB_MeV)
    first = _record_background_sample(model, config, T_MeV, muB_MeV)
    return _point_benchmark(model, config, settings, T_MeV, muB_MeV;
        initial_record=first, process_first=true, include_density=settings.include_density).report
end

function run_warm(model, config, settings; checkpoint=identity)
    points = default_points()
    records = NamedTuple[]
    continuation_seed = nothing
    for (index, (T_MeV, muB_MeV)) in enumerate(points)
        point_result = _point_benchmark(model, config, settings, T_MeV, muB_MeV;
            continuation_seed=continuation_seed, process_first=index == 1,
            include_density=settings.include_density)
        push!(records, point_result.report)
        continuation_seed = point_result.next_seed
        checkpoint(records)
        println("point $(index)/$(length(points)) complete: T=$(T_MeV), muB=$(muB_MeV)"); flush(stdout)
    end
    return (points=[(T_MeV=T, muB_MeV=μ) for (T, μ) in points], points_results=records,
        note="Only point 1 is process-first. Each point is prewarmed independently; repeated seed inputs are fixed to the preceding point solution.")
end

function host_metadata(config_path, source_before)
    return (julia_version=string(VERSION), cpu_name=string(Sys.CPU_NAME), cpu_threads=Sys.CPU_THREADS,
        julia_threads=Threads.nthreads(), kernel=string(Sys.KERNEL),
        process_startup_included=false, include_and_package_load_included=false,
        compiler_counter_scope="inside_at_timed;outer_inference_may_be_excluded",
        git=_git_metadata(), source_hashes_before=source_before,
        source_hashes_after=source_hashes(), config_path=replace(relpath(config_path, PROJECT_ROOT), '\\' => '/'))
end

function main()
    config_path = Workflow.DEFAULT_CONFIG
    config = TOML.parsefile(config_path)
    Workflow.validate_config(config)
    model = Models.create_model(:PNJL)
    settings = screening_settings()
    mode = lowercase(_env_string("BENCH_CGBU_MODE", "warm"))
    mode in ("single", "warm") || throw(ArgumentError("BENCH_CGBU_MODE must be single or warm"))
    output = let value = _env_string("BENCH_CGBU_OUTPUT", "")
        isempty(strip(value)) ? nothing : value
    end
    source_before = source_hashes()
    source_snapshot(output)
    write_json(output, (schema="charged_gbu_contour_benchmark_v4", status="running",
        settings=settings, resolved_production_config=config,
        host=host_metadata(config_path, source_before), source_hashes_before=source_before))
    points = default_points()
    println("Charged GBU infinite-thermal contour benchmark")
    @printf("mode=%s, kernel=charged_gbu_infinite_v1, points=%d, repeats=%d, density_repeats=%d\n",
        mode, length(points), settings.repeats, settings.density_repeats)
    @printf("screening: mesh=%d cut=%d tail=%d omega=%d q=%d qmax=%.3f\n", settings.mesh,
        settings.cut_nodes, settings.tail_nodes, settings.omega_nodes, settings.q_nodes, settings.qmax)
    result = if mode == "single"
        run_single(model, config, settings, points[1]...)
    else
        run_warm(model, config, settings; checkpoint=records -> write_json(output,
            (schema="charged_gbu_contour_benchmark_v4", status="running", settings=settings,
                resolved_production_config=config, host=host_metadata(config_path, source_before),
                source_hashes_before=source_before, points_results=records)))
    end
    source_after = source_hashes()
    payload = (schema="charged_gbu_contour_benchmark_v4", generated_at_utc=string(Dates.now(Dates.UTC)),
        mode=mode, status=source_before == source_after ? "complete" : "source_changed",
        source_unchanged=source_before == source_after,
        method=(route="charged_gbu_infinite_v1", thermal_target="infinity",
            background="FixedMuBConservedCharges", observable="fixed_quark_only_gbu_partial_yield",
            screening_only=settings.workload == "screening"), settings=settings,
        resolved_production_config=config,
        host=host_metadata(config_path, source_before), result=result)
    write_json(output, payload)
    output === nothing || println("benchmark JSON: $(abspath(output))")
    source_before == source_after || error("tracked source changed during benchmark")
    return payload
end

end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && ChargedGBUContourBenchmark.main()
