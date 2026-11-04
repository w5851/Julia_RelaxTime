using Test
using TOML
using SHA
using LinearAlgebra: BLAS

if !isdefined(Main, :TransportRunnerTiming)
    include(joinpath(@__DIR__, "..", "..", "..", "benchmark", "relaxtime", "transport_runner_timing.jl"))
end

@testset "Runner timing preserves CLI execution and failure semantics" begin
    original_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        mktempdir() do folder
            entrypoint = joinpath(folder, "entrypoint.jl")
            output = joinpath(folder, "timing.toml")
            write(entrypoint, "const calls = Ref(0)\nmain() = (calls[] += 1; copy(ARGS))\n")
            namespace = Module(gensym(:TimedEntrypoint))
            arguments = copy(ARGS)
            returned = TransportRunnerTiming.timed_entrypoint(output; entrypoint, namespace)
            @test returned == ARGS == arguments
            @test getfield(namespace, :calls)[] == 1
            report = TOML.parsefile(output)
            @test report["entrypoint_returned"]
            @test report["argv"] == arguments
            @test report["entrypoint_sha256"] == bytes2hex(sha256(read(entrypoint)))
            @test report["julia_threads"] == report["blas_threads_after_load"] == report["blas_threads_at_end"] == 1
            @test report["total_memory_bytes"] > 0
            @test report["compile_time_included_in_elapsed"]
            for stage in ("load", "scan")
                @test isfinite(report[stage]["elapsed_seconds"]) && report[stage]["elapsed_seconds"] >= 0
                @test report[stage]["allocated_bytes"] >= 0
                @test report[stage]["compile_timing_available"] == (VERSION >= v"1.11")
                if report[stage]["compile_timing_available"]
                    @test isfinite(report[stage]["compile_seconds"]) && report[stage]["compile_seconds"] >= 0
                    @test isfinite(report[stage]["recompile_seconds"]) && report[stage]["recompile_seconds"] >= 0
                end
            end
            original = read(output)
            @test_throws ArgumentError TransportRunnerTiming.timed_entrypoint(output; entrypoint, namespace)
            @test read(output) == original
            @test getfield(namespace, :calls)[] == 1

            failing = joinpath(folder, "scan_failure.jl")
            failed_output = joinpath(folder, "scan_failure.toml")
            write(failing, "main() = error(\"intentional scan failure\")\n")
            @test_throws ErrorException TransportRunnerTiming.timed_entrypoint(failed_output;
                entrypoint=failing, namespace=Module(gensym(:ScanFailure)))
            failed = TOML.parsefile(failed_output)
            @test !failed["entrypoint_returned"] && failed["failed_stage"] == "scan"
            @test haskey(failed, "load") && !haskey(failed, "scan")
            @test occursin("intentional scan failure", failed["error"])
            @test ARGS == arguments

            failing_load = joinpath(folder, "load_failure.jl")
            failed_load_output = joinpath(folder, "load_failure.toml")
            write(failing_load, "error(\"intentional load failure\")\n")
            @test_throws LoadError TransportRunnerTiming.timed_entrypoint(failed_load_output;
                entrypoint=failing_load, namespace=Module(gensym(:LoadFailure)))
            failed_load = TOML.parsefile(failed_load_output)
            @test !failed_load["entrypoint_returned"] && failed_load["failed_stage"] == "load"
            @test !haskey(failed_load, "load") && !haskey(failed_load, "scan")
        end
    finally
        BLAS.set_num_threads(original_blas_threads)
    end
end
