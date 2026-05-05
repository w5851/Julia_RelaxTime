using Test

const PROJECT_ROOT_FBCC = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :ServerLauncher)
    include(joinpath(PROJECT_ROOT_FBCC, "src", "simulation", "ServerLauncher.jl"))
end

const FBCC = Main.ServerLauncher

@testset "frontend-backend split config contract" begin
    @testset "backend runtime policy profiles" begin
        localhost = FBCC.server_runtime_policy("localhost")
        staging = FBCC.server_runtime_policy("staging")
        remote = FBCC.server_runtime_policy("remote")

        @test localhost.profile == "localhost"
        @test localhost.cors_allow_origins == "*"
        @test localhost.scan_max_running >= 1
        @test localhost.scan_max_pending >= localhost.scan_max_running
        @test localhost.scan_job_timeout_seconds == 0
        @test localhost.warmup_profile == "none"
        @test localhost.warmup_strict == false

        @test staging.profile == "staging"
        @test staging.cors_allow_origins != "*"
        @test staging.scan_job_timeout_seconds > 0
        @test staging.warmup_profile == "point"
        @test staging.warmup_strict == false

        @test remote.profile == "remote"
        @test remote.cors_allow_origins != "*"
        @test remote.scan_job_timeout_seconds > 0
        @test remote.warmup_profile == "service_core"
        @test remote.warmup_strict == false
    end

    @testset "backend policy can export env map" begin
        env_map = FBCC.runtime_policy_env("staging")
        @test haskey(env_map, "JRT_DEPLOY_PROFILE")
        @test haskey(env_map, "JRT_CORS_ALLOW_ORIGINS")
        @test haskey(env_map, "PNJL_SCAN_MAX_RUNNING")
        @test haskey(env_map, "PNJL_SCAN_MAX_PENDING")
        @test haskey(env_map, "PNJL_SCAN_JOB_TIMEOUT_SECONDS")
        @test haskey(env_map, "JRT_SERVER_WARMUP_PROFILE")
        @test haskey(env_map, "JRT_SERVER_WARMUP_STRICT")
        @test env_map["JRT_SERVER_WARMUP_PROFILE"] == "point"
    end

    @testset "backend warmup profiles expose stable step sets" begin
        @test FBCC.resolve_server_warmup_profile("none") === :none
        @test FBCC.resolve_server_warmup_profile("point") === :point
        @test FBCC.resolve_server_warmup_profile("service_core") === :service_core
        @test FBCC.list_server_warmup_steps(:none) == Symbol[]
        @test FBCC.list_server_warmup_steps(:point) == [:pnjl_point, :transport_point]
        @test FBCC.list_server_warmup_steps(:service_core) == [:pnjl_point, :transport_point, :phase_pipeline]
    end
end
