using Test
using HTTP
using JSON3

const PROJECT_ROOT_SJLC = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_SJLC, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_SJLC, "src", "simulation", "FullServerApp.jl"))
end

const SJLC = Main.FullServerApp

function _body(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

function _reset_logging_fixture!()
    lock(SJLC._PNJL_SCAN_JOBS_LOCK) do
        empty!(SJLC._PNJL_SCAN_JOBS)
        SJLC._PNJL_SCAN_JOB_SEQ[] = 0
    end
    return nothing
end

@testset "scan job logging contract" begin
    _reset_logging_fixture!()

    @testset "event helper shape" begin
        job = Dict{String, Any}(
            "job_id" => "log-job-1",
            "kind" => "tmu",
            "status" => "queued",
        )
        evt = SJLC._new_job_event(job, "created")
        @test evt["job_id"] == "log-job-1"
        @test evt["kind"] == "tmu"
        @test evt["state"] == "queued"
        @test evt["event_code"] == "created"
        @test haskey(evt, "timestamp")
    end

    @testset "status endpoint includes lifecycle events" begin
        job_id = SJLC._new_job("tmu", Dict{String, Any}("kind" => "tmu"); total_points=1)
        status_resp = SJLC.handle_pnjl_scan_job_status(job_id)
        status_body = _body(status_resp)

        @test status_resp.status == 200
        @test haskey(status_body, :events)
        @test length(status_body.events) >= 1

        first_event = status_body.events[1]
        @test haskey(first_event, :job_id)
        @test haskey(first_event, :kind)
        @test haskey(first_event, :state)
        @test haskey(first_event, :timestamp)
        @test haskey(first_event, :event_code)
    end

    _reset_logging_fixture!()
end
