using Test
using HTTP
using JSON3

const PROJECT_ROOT_PSMD = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_PSMD, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_PSMD, "src", "simulation", "FullServerApp.jl"))
end

const PSMD = Main.FullServerApp

function _status_body(job_id::String)
    resp = PSMD.handle_pnjl_scan_job_status(job_id)
    return resp, JSON3.read(String(resp.body))
end

function _reset_metrics_fixture!()
    lock(PSMD._PNJL_SCAN_JOBS_LOCK) do
        empty!(PSMD._PNJL_SCAN_JOBS)
        PSMD._PNJL_SCAN_JOB_SEQ[] = 0
        if isdefined(PSMD, :_PNJL_SCAN_RUNTIME_METRICS)
            empty!(PSMD._PNJL_SCAN_RUNTIME_METRICS)
            PSMD._PNJL_SCAN_RUNTIME_METRICS["job_succeeded"] = 0
            PSMD._PNJL_SCAN_RUNTIME_METRICS["job_failed"] = 0
            PSMD._PNJL_SCAN_RUNTIME_METRICS["job_cancelled"] = 0
            PSMD._PNJL_SCAN_RUNTIME_METRICS["duration_le_10s"] = 0
            PSMD._PNJL_SCAN_RUNTIME_METRICS["duration_10s_60s"] = 0
            PSMD._PNJL_SCAN_RUNTIME_METRICS["duration_gt_60s"] = 0
        end
    end
    return nothing
end

@testset "PNJL scan metrics and diagnostics" begin
    _reset_metrics_fixture!()

    @testset "status payload exposes metrics snapshot" begin
        job_id = PSMD._new_job("tmu", Dict{String, Any}("kind" => "tmu"); total_points=1)
        lock(PSMD._PNJL_SCAN_JOBS_LOCK) do
            job = PSMD._PNJL_SCAN_JOBS[job_id]
            job["status"] = "failed"
            job["started_at"] = "2026-03-29T00:00:00"
            job["ended_at"] = "2026-03-29T00:00:05"
            job["reason_code"] = "EXECUTION_ERROR"
            PSMD._update_scan_runtime_metrics_on_terminal!(job)
        end

        status_resp, status_body = _status_body(job_id)
        @test status_resp.status == 200
        @test haskey(status_body, :metrics)
        @test haskey(status_body.metrics, :terminal)
        @test haskey(status_body.metrics, :duration_buckets)
        @test haskey(status_body.metrics, :queue)

        @test haskey(status_body.metrics.terminal, :succeeded)
        @test haskey(status_body.metrics.terminal, :failed)
        @test haskey(status_body.metrics.terminal, :cancelled)
        @test status_body.metrics.terminal.failed >= 1

        @test haskey(status_body.metrics.duration_buckets, :le_10s)
        @test haskey(status_body.metrics.duration_buckets, :between_10s_60s)
        @test haskey(status_body.metrics.duration_buckets, :gt_60s)
        @test status_body.metrics.duration_buckets.le_10s >= 1

        @test haskey(status_body.diagnostics, :reason_code)
        @test status_body.diagnostics.reason_code == "EXECUTION_ERROR"
    end

    _reset_metrics_fixture!()
end
