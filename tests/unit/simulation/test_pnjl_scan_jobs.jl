using Test
using HTTP
using JSON3
using Dates: DateTime

const PROJECT_ROOT_PSJ = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_PSJ, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_PSJ, "src", "simulation", "FullServerApp.jl"))
end

const PSJ = Main.FullServerApp

function _post_jobs_request(payload::Dict)
    return HTTP.Request(
        "POST",
        "/api/modules/pnjl-scan/jobs",
        ["Content-Type" => "application/json"],
        JSON3.write(payload),
    )
end

function _body_dict(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

function _reset_jobs_state!()
    lock(PSJ._PNJL_SCAN_JOBS_LOCK) do
        empty!(PSJ._PNJL_SCAN_JOBS)
        PSJ._PNJL_SCAN_JOB_SEQ[] = 0
        if isdefined(PSJ, :_PNJL_SCAN_PRUNE_METRICS)
            empty!(PSJ._PNJL_SCAN_PRUNE_METRICS)
            PSJ._PNJL_SCAN_PRUNE_METRICS["total"] = 0
            PSJ._PNJL_SCAN_PRUNE_METRICS["ttl"] = 0
            PSJ._PNJL_SCAN_PRUNE_METRICS["keep_max"] = 0
        end
    end
    return nothing
end

@testset "PNJL scan jobs path guard" begin
    @testset "legal output path accepted and normalized" begin
        params = Dict{Symbol, Any}(
            :output_path => "data/outputs/results/pnjl/scan/tmu/unit_guard_ok.csv",
        )
        out = PSJ._safe_output_path("tmu", params, "job-guard-ok")
        @test startswith(out, joinpath(PROJECT_ROOT_PSJ, "data", "outputs"))
        @test endswith(out, "unit_guard_ok.csv")
    end

    @testset "escaping output path rejected" begin
        params = Dict{Symbol, Any}(
            :output_path => "..\\..\\Windows\\Temp\\pwn.csv",
        )
        @test_throws ArgumentError PSJ._safe_output_path("tmu", params, "job-guard-bad")
    end
end

@testset "PNJL scan jobs contract" begin
    _reset_jobs_state!()

    @testset "scan mode contract: point/scan" begin
        params_point_tmu = Dict{Symbol, Any}(
            :mode => "point",
            :T_mev => 150.0,
            :mu_mev => 0.0,
            :xi => 0.0,
        )
        @test PSJ._estimate_total_points("tmu", params_point_tmu) == 1

        params_point_trho = Dict{Symbol, Any}(
            :mode => "point",
            :T_mev => 150.0,
            :rho_value => 0.001,
            :xi => 0.0,
        )
        @test PSJ._estimate_total_points("trho", params_point_trho) == 1

        params_scan_tmu = Dict{Symbol, Any}(
            :mode => "scan",
            :T_values => [150.0, 160.0],
            :mu_values => [0.0, 100.0],
            :xi => 0.0,
        )
        @test PSJ._estimate_total_points("tmu", params_scan_tmu) == 4
    end

    @testset "create/status accepts point mode and persists policy mode" begin
        _reset_jobs_state!()
        output_rel = joinpath("data", "outputs", "results", "pnjl", "scan", "tmu", "unit_contract_point_$(time_ns()).csv")
        output_abs = joinpath(PROJECT_ROOT_PSJ, output_rel)

        req = _post_jobs_request(Dict(
            "kind" => "tmu",
            "params" => Dict(
                "mode" => "point",
                "T_mev" => 150.0,
                "mu_mev" => 0.0,
                "xi" => 0.0,
                "output_path" => output_rel,
                "overwrite" => true,
                "resume" => false,
                "max_retries" => 0,
                "p_num" => 12,
                "t_num" => 6,
            ),
        ))
        create_resp = PSJ.handle_pnjl_scan_job_create(req)
        create_body = _body_dict(create_resp)
        @test create_resp.status == 202
        @test create_body.status == "accepted"
        @test haskey(create_body, :job_id)

        job_id = String(create_body.job_id)
        status_resp = PSJ.handle_pnjl_scan_job_status(job_id)
        status_body = _body_dict(status_resp)
        @test status_resp.status == 200
        @test status_body.status == "ok"
        @test status_body.policy.mode == "point"

        final_status = get(status_body, :job_status, nothing)
        for _ in 1:120
            final_status in ("succeeded", "failed") && break
            sleep(0.25)
            latest_resp = PSJ.handle_pnjl_scan_job_status(job_id)
            latest_body = _body_dict(latest_resp)
            final_status = get(latest_body, :job_status, nothing)
        end
        @test final_status in ("succeeded", "failed")

        for _ in 1:20
            if !isfile(output_abs)
                break
            end
            try
                rm(output_abs; force=true)
                break
            catch
                sleep(0.1)
            end
        end
        @test !isfile(output_abs)
    end

    @testset "invalid mode rejected" begin
        req = _post_jobs_request(Dict(
            "kind" => "tmu",
            "params" => Dict(
                "mode" => "badmode",
                "T_values" => [150.0],
                "mu_values" => [0.0],
                "xi" => 0.0,
            ),
        ))
        resp = PSJ.handle_pnjl_scan_job_create(req)
        body = _body_dict(resp)
        @test resp.status == 400
        @test body.status == "error"
        @test body.error_code == "INVALID_REQUEST"
        @test occursin("mode", String(body.error))
    end

    @testset "max_retries gate" begin
        req = _post_jobs_request(Dict(
            "kind" => "tmu",
            "params" => Dict(
                "max_retries" => 4,
            ),
        ))
        resp = PSJ.handle_pnjl_scan_job_create(req)
        body = _body_dict(resp)
        @test resp.status == 400
        @test body.status == "error"
        @test body.error_code == "INVALID_REQUEST"
        @test occursin("max_retries", String(body.error))
        @test haskey(body, :message_id)
        @test haskey(body, :diagnostics)
        @test body.diagnostics.error_type in ("ErrorException", "ArgumentError")
    end

    @testset "xi strategy exclusivity gate" begin
        req = _post_jobs_request(Dict(
            "kind" => "tmu",
            "params" => Dict(
                "xi" => 0.1,
                "xi_values" => [0.0, 0.1],
            ),
        ))
        resp = PSJ.handle_pnjl_scan_job_create(req)
        body = _body_dict(resp)
        @test resp.status == 400
        @test body.status == "error"
        @test body.error_code == "INVALID_REQUEST"
        @test occursin("Use only one xi strategy", String(body.error))
        @test haskey(body, :message_id)
    end

    @testset "job lifecycle primitives" begin
        job_id = PSJ._new_job("tmu", Dict{String, Any}("kind" => "tmu"); total_points=2)
        pos = PSJ._queue_position(job_id)
        @test pos isa Int
        status_payload = PSJ.handle_pnjl_scan_job_status(job_id)
        status_body = _body_dict(status_payload)
        @test status_payload.status == 200
        @test status_body.job_status == "queued"
        @test haskey(status_body, :diagnostics)
        @test status_body.diagnostics.job_id == job_id
        @test status_body.diagnostics.kind == "tmu"
        @test status_body.diagnostics.job_status == "queued"

        result_payload = PSJ.handle_pnjl_scan_job_result(job_id)
        result_body = _body_dict(result_payload)
        @test result_payload.status == 409
        @test result_body.status == "error"
        @test result_body.error_code == "JOB_NOT_SUCCEEDED"
        @test haskey(result_body, :message_id)
        @test result_body.job_status == "queued"
        @test haskey(result_body, :diagnostics)
        @test result_body.diagnostics.job_id == job_id
    end

    @testset "queue limit gate" begin
        _reset_jobs_state!()
        lock(PSJ._PNJL_SCAN_JOBS_LOCK) do
            for i in 1:PSJ._PNJL_SCAN_MAX_PENDING
                id = "queued_$(i)"
                PSJ._PNJL_SCAN_JOBS[id] = Dict{String, Any}(
                    "job_id" => id,
                    "seq" => i,
                    "kind" => "tmu",
                    "status" => "queued",
                    "created_at" => "",
                    "started_at" => nothing,
                    "ended_at" => nothing,
                    "progress" => Dict{String, Any}("total" => nothing, "completed" => 0, "percent" => nothing),
                    "result" => nothing,
                    "error" => nothing,
                    "request" => Dict{String, Any}(),
                )
            end
            PSJ._PNJL_SCAN_JOB_SEQ[] = PSJ._PNJL_SCAN_MAX_PENDING
        end

        req = _post_jobs_request(Dict("kind" => "tmu", "params" => Dict()))
        resp = PSJ.handle_pnjl_scan_job_create(req)
        body = _body_dict(resp)
        @test resp.status == 429
        @test body.status == "error"
        @test body.error_code == "QUEUE_FULL"
        @test body.error == "queue is full"
        @test haskey(body, :message_id)
    end

    @testset "finished jobs pruning keeps active jobs" begin
        _reset_jobs_state!()
        lock(PSJ._PNJL_SCAN_JOBS_LOCK) do
            PSJ._PNJL_SCAN_JOBS["queued"] = Dict{String, Any}(
                "job_id" => "queued",
                "seq" => 1,
                "kind" => "tmu",
                "status" => "queued",
                "created_at" => "",
                "started_at" => nothing,
                "ended_at" => nothing,
                "progress" => Dict{String, Any}("total" => nothing, "completed" => 0, "percent" => nothing),
                "result" => nothing,
                "error" => nothing,
                "request" => Dict{String, Any}(),
            )
            PSJ._PNJL_SCAN_JOBS["running"] = Dict{String, Any}(
                "job_id" => "running",
                "seq" => 2,
                "kind" => "tmu",
                "status" => "running",
                "created_at" => "",
                "started_at" => nothing,
                "ended_at" => nothing,
                "progress" => Dict{String, Any}("total" => nothing, "completed" => 0, "percent" => nothing),
                "result" => nothing,
                "error" => nothing,
                "request" => Dict{String, Any}(),
            )

            for i in 1:5
                id = "done_$(i)"
                PSJ._PNJL_SCAN_JOBS[id] = Dict{String, Any}(
                    "job_id" => id,
                    "seq" => 2 + i,
                    "kind" => "tmu",
                    "status" => "succeeded",
                    "created_at" => "",
                    "started_at" => nothing,
                    "ended_at" => "",
                    "progress" => Dict{String, Any}("total" => 1, "completed" => 1, "percent" => 100.0),
                    "result" => Dict{String, Any}(),
                    "error" => nothing,
                    "request" => Dict{String, Any}(),
                )
            end

            PSJ._prune_finished_jobs_locked!(2)
        end

        lock(PSJ._PNJL_SCAN_JOBS_LOCK) do
            @test haskey(PSJ._PNJL_SCAN_JOBS, "queued")
            @test haskey(PSJ._PNJL_SCAN_JOBS, "running")
            @test !haskey(PSJ._PNJL_SCAN_JOBS, "done_1")
            @test !haskey(PSJ._PNJL_SCAN_JOBS, "done_2")
            @test !haskey(PSJ._PNJL_SCAN_JOBS, "done_3")
            @test haskey(PSJ._PNJL_SCAN_JOBS, "done_4")
            @test haskey(PSJ._PNJL_SCAN_JOBS, "done_5")
        end
    end

    @testset "finished jobs pruning honors ttl" begin
        _reset_jobs_state!()
        lock(PSJ._PNJL_SCAN_JOBS_LOCK) do
            PSJ._PNJL_SCAN_JOBS["done_old"] = Dict{String, Any}(
                "job_id" => "done_old",
                "seq" => 1,
                "kind" => "tmu",
                "status" => "succeeded",
                "created_at" => "",
                "started_at" => nothing,
                "ended_at" => "2026-03-01T00:00:00",
                "progress" => Dict{String, Any}("total" => 1, "completed" => 1, "percent" => 100.0),
                "result" => Dict{String, Any}(),
                "error" => nothing,
                "request" => Dict{String, Any}(),
            )
            PSJ._PNJL_SCAN_JOBS["done_new"] = Dict{String, Any}(
                "job_id" => "done_new",
                "seq" => 2,
                "kind" => "tmu",
                "status" => "failed",
                "created_at" => "",
                "started_at" => nothing,
                "ended_at" => "2026-03-28T00:00:15",
                "progress" => Dict{String, Any}("total" => 1, "completed" => 1, "percent" => 100.0),
                "result" => Dict{String, Any}(),
                "error" => "x",
                "request" => Dict{String, Any}(),
            )

            removed = PSJ._prune_finished_jobs_locked!(10; now=DateTime(2026, 3, 28, 0, 0, 20), ttl_seconds=10)
            @test removed == 1
        end

        lock(PSJ._PNJL_SCAN_JOBS_LOCK) do
            @test !haskey(PSJ._PNJL_SCAN_JOBS, "done_old")
            @test haskey(PSJ._PNJL_SCAN_JOBS, "done_new")
        end
    end

    @testset "status endpoint keeps internal errors private" begin
        _reset_jobs_state!()
        lock(PSJ._PNJL_SCAN_JOBS_LOCK) do
            PSJ._PNJL_SCAN_JOBS["failed_hidden"] = Dict{String, Any}(
                "job_id" => "failed_hidden",
                "seq" => 1,
                "kind" => "tmu",
                "status" => "failed",
                "created_at" => "",
                "started_at" => nothing,
                "ended_at" => "2026-03-28T00:00:15",
                "progress" => Dict{String, Any}("total" => 1, "completed" => 1, "percent" => 100.0),
                "result" => nothing,
                "error" => "scan execution failed",
                "internal_error" => "Stacktrace: secret",
                "request" => Dict{String, Any}(),
            )
        end

        resp = PSJ.handle_pnjl_scan_job_status("failed_hidden")
        body = _body_dict(resp)
        @test resp.status == 200
        @test body.status == "ok"
        @test body.job_status == "failed"
        @test body.error == "scan execution failed"
        @test !haskey(body, :internal_error)
    end

    @testset "not-found responses include message_id" begin
        status_resp = PSJ.handle_pnjl_scan_job_status("missing_job")
        status_body = _body_dict(status_resp)
        @test status_resp.status == 404
        @test status_body.status == "error"
        @test status_body.error_code == "JOB_NOT_FOUND"
        @test haskey(status_body, :message_id)

        result_resp = PSJ.handle_pnjl_scan_job_result("missing_job")
        result_body = _body_dict(result_resp)
        @test result_resp.status == 404
        @test result_body.status == "error"
        @test result_body.error_code == "JOB_NOT_FOUND"
        @test haskey(result_body, :message_id)
    end

    @testset "scan policy can be configured by env" begin
        withenv(
            "PNJL_SCAN_FINISHED_KEEP_MAX" => "11",
            "PNJL_SCAN_FINISHED_TTL_SECONDS" => "123",
        ) do
            policy = PSJ._scan_jobs_policy()
            @test policy.keep_max == 11
            @test policy.ttl_seconds == 123
        end
    end

    @testset "status exposes minimal governance metrics" begin
        _reset_jobs_state!()
        job_id = PSJ._new_job("tmu", Dict{String, Any}("kind" => "tmu"); total_points=1)
        resp = PSJ.handle_pnjl_scan_job_status(job_id)
        body = _body_dict(resp)
        @test resp.status == 200
        @test haskey(body, :governance)
        @test haskey(body.governance, :finished_keep_max)
        @test haskey(body.governance, :finished_ttl_seconds)
        @test haskey(body.governance, :pruned_total)
        @test haskey(body.governance, :pruned_ttl)
        @test haskey(body.governance, :pruned_keep_max)
    end

    @testset "prune metrics count ttl and keep-max removals" begin
        _reset_jobs_state!()
        lock(PSJ._PNJL_SCAN_JOBS_LOCK) do
            PSJ._PNJL_SCAN_JOBS["done_old"] = Dict{String, Any}(
                "job_id" => "done_old",
                "seq" => 1,
                "kind" => "tmu",
                "status" => "succeeded",
                "created_at" => "",
                "started_at" => nothing,
                "ended_at" => "2026-03-01T00:00:00",
                "progress" => Dict{String, Any}("total" => 1, "completed" => 1, "percent" => 100.0),
                "result" => Dict{String, Any}(),
                "error" => nothing,
                "internal_error" => nothing,
                "request" => Dict{String, Any}(),
            )
            PSJ._PNJL_SCAN_JOBS["done_new_1"] = Dict{String, Any}(
                "job_id" => "done_new_1",
                "seq" => 2,
                "kind" => "tmu",
                "status" => "failed",
                "created_at" => "",
                "started_at" => nothing,
                "ended_at" => "2026-03-28T00:00:25",
                "progress" => Dict{String, Any}("total" => 1, "completed" => 1, "percent" => 100.0),
                "result" => Dict{String, Any}(),
                "error" => "x",
                "internal_error" => "x",
                "request" => Dict{String, Any}(),
            )
            PSJ._PNJL_SCAN_JOBS["done_new_2"] = Dict{String, Any}(
                "job_id" => "done_new_2",
                "seq" => 3,
                "kind" => "tmu",
                "status" => "succeeded",
                "created_at" => "",
                "started_at" => nothing,
                "ended_at" => "2026-03-28T00:00:26",
                "progress" => Dict{String, Any}("total" => 1, "completed" => 1, "percent" => 100.0),
                "result" => Dict{String, Any}(),
                "error" => nothing,
                "internal_error" => nothing,
                "request" => Dict{String, Any}(),
            )

            _ = PSJ._prune_finished_jobs_locked!(1; now=DateTime(2026, 3, 28, 0, 0, 30), ttl_seconds=10)

            @test PSJ._PNJL_SCAN_PRUNE_METRICS["ttl"] == 1
            @test PSJ._PNJL_SCAN_PRUNE_METRICS["keep_max"] == 1
            @test PSJ._PNJL_SCAN_PRUNE_METRICS["total"] == 2
        end
    end

    _reset_jobs_state!()
end
