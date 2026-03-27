using Test
using HTTP
using JSON3

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
        @test occursin("max_retries", String(body.error))
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
        @test occursin("Use only one xi strategy", String(body.error))
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
        @test body.error == "queue is full"
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

    _reset_jobs_state!()
end
