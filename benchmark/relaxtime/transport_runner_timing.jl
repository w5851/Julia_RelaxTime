#!/usr/bin/env julia

"""Measure the existing CLI without changing its arguments or numerical path."""
module TransportRunnerTiming

const SCRIPT_STARTED_NS = time_ns()

using Dates
using LinearAlgebra
using SHA
using TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const ENTRYPOINT = joinpath(ROOT, "scripts", "relaxtime", "run_phase_guided_transport_scan.jl")

function phase_metrics(measured)
    result = Dict{String,Any}(
        "elapsed_seconds" => measured.time,
        "gc_seconds" => measured.gctime,
        "allocated_bytes" => measured.bytes,
        "compile_timing_available" => hasproperty(measured, :compile_time),
    )
    if hasproperty(measured, :compile_time)
        result["compile_seconds"] = measured.compile_time
        result["recompile_seconds"] = measured.recompile_time
    end
    return result
end

function timed_entrypoint(output_path::AbstractString; entrypoint::AbstractString=ENTRYPOINT,
    namespace::Module=Main, script_started_ns::UInt64=SCRIPT_STARTED_NS)
    isempty(output_path) && throw(ArgumentError("output_path must not be empty"))
    output = abspath(output_path)
    ispath(output) && throw(ArgumentError("output_path already exists: $output"))
    isfile(entrypoint) || throw(ArgumentError("entrypoint does not exist: $entrypoint"))
    mkpath(dirname(output))
    started_ns = time_ns()
    report = Dict{String,Any}(
        "schema_version" => 1,
        "started_utc" => string(Dates.now(Dates.UTC)),
        "entrypoint" => abspath(entrypoint),
        "entrypoint_sha256" => bytes2hex(sha256(read(entrypoint))),
        "timing_script_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "argv" => copy(ARGS),
        "julia_version" => string(VERSION),
        "cpu_name" => Sys.CPU_NAME,
        "cpu_threads" => Sys.CPU_THREADS,
        "total_memory_bytes" => Int(Sys.total_memory()),
        "free_memory_bytes_at_start" => Int(Sys.free_memory()),
        "julia_threads" => Threads.nthreads(),
        "blas_threads_at_start" => BLAS.get_num_threads(),
        "calculation_sha" => get(ENV, "CALCULATION_SHA", ""),
        "workflow_sha" => get(ENV, "GITHUB_SHA", ""),
        "github_run_id" => get(ENV, "GITHUB_RUN_ID", ""),
        "wrapper_initialization_seconds" => Float64(started_ns - script_started_ns) / 1e9,
        "entrypoint_returned" => false,
        "compile_time_included_in_elapsed" => true,
        "allocation_semantics" => "Cumulative bytes allocated, not peak resident memory",
        "process_scope" => "GNU time supplies total process wall/CPU time and peak RSS; these phases exclude Julia startup and final report serialization",
        "acceptance_scope" => "Entry point return is not numerical acceptance; check the actual transport artifact separately",
    )
    stage = "runtime"
    try
        Threads.nthreads() == 1 && BLAS.get_num_threads() == 1 ||
            throw(ArgumentError("runner timing requires one Julia thread and one BLAS thread"))
        stage = "load"
        loaded = @timed Base.include(namespace, entrypoint)
        report["load"] = phase_metrics(loaded)
        report["blas_threads_after_load"] = BLAS.get_num_threads()
        BLAS.get_num_threads() == 1 || throw(ArgumentError("entrypoint loading changed the BLAS thread count"))
        stage = "scan"
        # The CLI owns parsing, provenance and computation. invokelatest is
        # required because the entrypoint was defined by include above.
        scanned = @timed Base.invokelatest(getfield(namespace, :main))
        report["scan"] = phase_metrics(scanned)
        report["entrypoint_returned"] = true
        return scanned.value
    catch err
        report["failed_stage"] = stage
        report["error"] = sprint(showerror, err)
        rethrow()
    finally
        report["measured_wrapper_seconds"] = Float64(time_ns() - started_ns) / 1e9
        report["blas_threads_at_end"] = BLAS.get_num_threads()
        report["finished_utc"] = string(Dates.now(Dates.UTC))
        open(output, "w") do io
            TOML.print(io, report; sorted=true)
        end
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    TransportRunnerTiming.timed_entrypoint(get(ENV, "TRANSPORT_RUNNER_TIMING_PATH", ""))
end
