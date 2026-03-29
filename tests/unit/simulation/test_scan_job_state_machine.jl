using Test

const PROJECT_ROOT_SJSM = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_SJSM, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_SJSM, "src", "simulation", "FullServerApp.jl"))
end

const SJSM = Main.FullServerApp

@testset "scan job state transition table" begin
    @testset "legal transitions" begin
        @test SJSM._is_valid_job_status_transition("queued", "running")
        @test SJSM._is_valid_job_status_transition("queued", "cancelled")
        @test SJSM._is_valid_job_status_transition("running", "running")
        @test SJSM._is_valid_job_status_transition("running", "succeeded")
        @test SJSM._is_valid_job_status_transition("running", "failed")
        @test SJSM._is_valid_job_status_transition("running", "cancelled")
    end

    @testset "illegal transitions" begin
        @test !SJSM._is_valid_job_status_transition("queued", "succeeded")
        @test !SJSM._is_valid_job_status_transition("queued", "failed")
        @test !SJSM._is_valid_job_status_transition("succeeded", "running")
        @test !SJSM._is_valid_job_status_transition("failed", "running")
        @test !SJSM._is_valid_job_status_transition("cancelled", "running")
    end
end

@testset "scan job state mutation guard" begin
    job = Dict{String, Any}(
        "job_id" => "job-state-1",
        "status" => "queued",
        "started_at" => nothing,
        "ended_at" => nothing,
        "error" => nothing,
        "internal_error" => nothing,
        "result" => nothing,
    )

    SJSM._set_job_status!(job, "running")
    @test job["status"] == "running"
    @test job["started_at"] !== nothing
    @test job["ended_at"] === nothing

    SJSM._set_job_status!(job, "cancelled")
    @test job["status"] == "cancelled"
    @test job["ended_at"] !== nothing

    @test_throws ArgumentError SJSM._set_job_status!(job, "running")
    @test_throws ArgumentError SJSM._set_job_status!(job, "unknown")
end
