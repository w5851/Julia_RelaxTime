using Test
using HTTP
using JSON3

const PROJECT_ROOT_PSCT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_PSCT, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_PSCT, "src", "simulation", "FullServerApp.jl"))
end

const PSCT = Main.FullServerApp

function _json_body(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

function _post_job(payload::Dict)
    req = HTTP.Request("POST", "/api/modules/pnjl-scan/jobs", ["Content-Type" => "application/json"], JSON3.write(payload))
    return PSCT.handle_pnjl_scan_job_create(req)
end

function _reset_cancel_timeout_fixture!()
    lock(PSCT._PNJL_SCAN_JOBS_LOCK) do
        empty!(PSCT._PNJL_SCAN_JOBS)
        PSCT._PNJL_SCAN_JOB_SEQ[] = 0
        if isdefined(PSCT, :_PNJL_SCAN_IDEMPOTENCY_CACHE)
            empty!(PSCT._PNJL_SCAN_IDEMPOTENCY_CACHE)
        end
    end
    return nothing
end

@testset "PNJL scan cancel and timeout governance" begin
    _reset_cancel_timeout_fixture!()

    @testset "cancel queued/running allowed" begin
        payload = Dict(
            "kind" => "tmu",
            "params" => Dict(
                "T_values" => [150.0],
                "mu_values" => [0.0],
                "xi" => 0.0,
                "max_retries" => 0,
            ),
        )
        created = _post_job(payload)
        created_body = _json_body(created)
        @test created.status == 202

        cancel_resp = PSCT.handle_pnjl_scan_job_cancel(String(created_body.job_id))
        cancel_body = _json_body(cancel_resp)
        @test cancel_resp.status == 200
        @test cancel_body.job_status == "cancelled"
    end

    @testset "cancel terminal state rejected" begin
        lock(PSCT._PNJL_SCAN_JOBS_LOCK) do
            PSCT._PNJL_SCAN_JOBS["done_for_cancel"] = Dict{String, Any}(
                "job_id" => "done_for_cancel",
                "seq" => 100,
                "kind" => "tmu",
                "status" => "succeeded",
                "created_at" => "",
                "started_at" => "",
                "ended_at" => "",
                "progress" => Dict{String, Any}("total" => 1, "completed" => 1, "percent" => 100.0),
                "result" => Dict{String, Any}(),
                "error" => nothing,
                "internal_error" => nothing,
                "request" => Dict{String, Any}(),
            )
        end

        cancel_resp = PSCT.handle_pnjl_scan_job_cancel("done_for_cancel")
        cancel_body = _json_body(cancel_resp)
        @test cancel_resp.status == 409
        @test cancel_body.error_code == "JOB_NOT_CANCELLABLE"
    end

    @testset "timeout produces deterministic terminal code/state" begin
        job_id = PSCT._new_job("tmu", Dict{String, Any}("kind" => "tmu"); total_points=1)
        lock(PSCT._PNJL_SCAN_JOBS_LOCK) do
            job = PSCT._PNJL_SCAN_JOBS[job_id]
            job["started_at"] = "2026-03-29T00:00:00"
            job["status"] = "running"
            job["policy"] = Dict("timeout_seconds" => 1)
        end

        changed = PSCT._maybe_mark_job_timeout!(job_id; now_ts="2026-03-29T00:00:01")
        @test changed == true

        status_resp = PSCT.handle_pnjl_scan_job_status(job_id)
        status_body = _json_body(status_resp)
        @test status_resp.status == 200
        @test status_body.job_status == "failed"
        @test status_body.error == "scan execution timeout"
        @test status_body.diagnostics.reason_code == "TIMEOUT"
    end

    _reset_cancel_timeout_fixture!()
end
