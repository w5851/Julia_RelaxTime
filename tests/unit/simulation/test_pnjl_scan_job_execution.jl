using Test
using HTTP

# Exercise the real job lifecycle with a deterministic numerical stub.
module ScanJobExecutionTest
using HTTP
using JSON3
module Models
const RUN = Ref{Function}(() -> nothing)
run_tmu_scan(; kwargs...) = RUN[](; kwargs...)
end
include(joinpath(@__DIR__, "..", "..", "..", "src", "simulation", "fullserver", "shared.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "simulation", "fullserver", "pnjl_scan_jobs.jl"))
end

@testset "Scan worker terminal states and output paths" begin
    jobs = ScanJobExecutionTest
    params = Dict{Symbol,Any}(:mode => "point", :T_mev => 150.0, :mu_mev => 0.0, :xi => 0.0)

    @testset "cancelled attempts are not retried" begin
        id = jobs._new_job("tmu", Dict{String,Any}(); total_points=2)
        jobs._PNJL_SCAN_JOBS[id]["policy"] = Dict("max_retries" => 3, "timeout_seconds" => 0)
        calls = Ref(0)
        jobs.Models.RUN[] = function (; stop_requested, progress_cb, kwargs...)
            calls[] += 1
            @test !stop_requested()
            progress_cb(nothing, nothing)
            cancelled = jobs.handle_pnjl_scan_job_cancel(id)
            @test cancelled.status == 200
            @test jobs.JSON3.read(String(cancelled.body)).worker_active
            @test jobs._queue_snapshot().running == 1
            @test stop_requested()
            progress_cb(nothing, nothing)
            error("failure after cancellation")
        end
        fetch(jobs._launch_scan_job(id, "tmu", params))
        @test calls[] == 1
        @test jobs._PNJL_SCAN_JOBS[id]["status"] == "cancelled"
        @test jobs._PNJL_SCAN_JOBS[id]["progress"]["completed"] == 1
        @test !jobs._PNJL_SCAN_JOBS[id]["worker_active"]
        @test jobs._queue_snapshot().running == 0
        @test !any(event -> event["event_code"] == "retry", jobs._PNJL_SCAN_JOBS[id]["events"])
    end

    @testset "a deadline stops work at the next boundary" begin
        id = jobs._new_job("tmu", Dict{String,Any}())
        jobs._PNJL_SCAN_JOBS[id]["policy"] = Dict("max_retries" => 3, "timeout_seconds" => 120)
        jobs.Models.RUN[] = function (; stop_requested, kwargs...)
            @test !stop_requested()
            jobs._with_job(id) do job
                job["started_ns"] = time_ns() - UInt64(120_000_000_000)
            end
            @test stop_requested()
            @test jobs._queue_snapshot().running == 1
            error("failure after deadline")
        end
        fetch(jobs._launch_scan_job(id, "tmu", params))
        job = jobs._PNJL_SCAN_JOBS[id]
        @test job["reason_code"] == "TIMEOUT"
        @test count(event -> event["event_code"] == "ended", job["events"]) == 1
        @test !any(event -> event["event_code"] == "retry", job["events"])
        @test jobs._queue_snapshot().running == 0
    end

    @testset "a completed attempt still enforces its deadline" begin
        id = jobs._new_job("tmu", Dict{String,Any}())
        jobs._PNJL_SCAN_JOBS[id]["policy"] = Dict("max_retries" => 0, "timeout_seconds" => 120)
        jobs.Models.RUN[] = function (; output_path, kwargs...)
            jobs._with_job(id) do job
                job["started_ns"] = time_ns() - UInt64(120_000_000_000)
            end
            return (output=output_path, total=1, success=1, failure=0, skipped=0)
        end
        fetch(jobs._launch_scan_job(id, "tmu", params))
        job = jobs._PNJL_SCAN_JOBS[id]
        @test job["status"] == "failed"
        @test job["reason_code"] == "TIMEOUT"
        @test job["result"] === nothing
        @test last(job["events"])["reason_code"] == "TIMEOUT"
    end

    @testset "output requests rejected before registration" begin
        before = length(jobs._PNJL_SCAN_JOBS)
        bad = merge(params, Dict{Symbol,Any}(:output_path => joinpath(jobs._PROJECT_ROOT, "outside.csv")))
        @test_throws ArgumentError jobs._reserve_scan_job("tmu", bad, nothing, "fixture")
        @test length(jobs._PNJL_SCAN_JOBS) == before
        @test_throws ArgumentError jobs._ensure_safe_output_path(jobs._SCAN_OUTPUT_ROOT)
        good = joinpath(jobs._SCAN_OUTPUT_ROOT, "new_directory", "result.csv")
        @test jobs._ensure_safe_output_path(good) == good
        @test !jobs._is_within_output_root(jobs._SCAN_OUTPUT_ROOT * "_other", jobs._SCAN_OUTPUT_ROOT)
        if !Sys.iswindows()
            @test !jobs._is_within_output_root(uppercase(jobs._SCAN_OUTPUT_ROOT), jobs._SCAN_OUTPUT_ROOT)
            mktempdir(jobs._SCAN_OUTPUT_ROOT) do inside
                mktempdir() do outside
                    symlink(outside, joinpath(inside, "link"))
                    @test_throws ArgumentError jobs._ensure_safe_output_path(joinpath(inside, "link", "result.csv"))
                end
            end
        end
    end
end
