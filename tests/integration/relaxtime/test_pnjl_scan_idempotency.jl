using Test
using HTTP
using JSON3

const PROJECT_ROOT_PSI = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_PSI, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_PSI, "src", "simulation", "FullServerApp.jl"))
end

const PSI = Main.FullServerApp

function _post_idempotency_request(payload::Dict; idempotency_key::Union{Nothing, String}=nothing)
    headers = Pair{String, String}["Content-Type" => "application/json"]
    if idempotency_key !== nothing
        push!(headers, "Idempotency-Key" => idempotency_key)
    end
    return HTTP.Request("POST", "/api/modules/pnjl-scan/jobs", headers, JSON3.write(payload))
end

function _body(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

function _reset_idempotency_fixture!()
    lock(PSI._PNJL_SCAN_JOBS_LOCK) do
        empty!(PSI._PNJL_SCAN_JOBS)
        PSI._PNJL_SCAN_JOB_SEQ[] = 0
        if isdefined(PSI, :_PNJL_SCAN_IDEMPOTENCY_CACHE)
            empty!(PSI._PNJL_SCAN_IDEMPOTENCY_CACHE)
        end
    end
    return nothing
end

@testset "PNJL scan idempotency" begin
    _reset_idempotency_fixture!()

    @testset "same key + same payload replays same job" begin
        key = "idem-same-001"
        payload = Dict(
            "kind" => "tmu",
            "params" => Dict(
                "T_values" => [150.0],
                "mu_values" => [0.0],
                "xi" => 0.0,
                "max_retries" => 0,
            ),
        )

        first = PSI.handle_pnjl_scan_job_create(_post_idempotency_request(payload; idempotency_key=key))
        first_body = _body(first)
        second = PSI.handle_pnjl_scan_job_create(_post_idempotency_request(payload; idempotency_key=key))
        second_body = _body(second)

        @test first.status == 202
        @test second.status == 202
        @test first_body.status == "accepted"
        @test second_body.status == "accepted"
        @test second_body.job_id == first_body.job_id
        idem = get(second_body, :idempotency, Dict{Symbol, Any}())
        @test get(idem, :replayed, false) == true
    end

    @testset "same key + different payload conflicts" begin
        key = "idem-conflict-001"
        payload_a = Dict(
            "kind" => "tmu",
            "params" => Dict(
                "T_values" => [150.0],
                "mu_values" => [0.0],
                "xi" => 0.0,
                "max_retries" => 0,
            ),
        )
        payload_b = Dict(
            "kind" => "tmu",
            "params" => Dict(
                "T_values" => [160.0],
                "mu_values" => [0.0],
                "xi" => 0.0,
                "max_retries" => 0,
            ),
        )

        first = PSI.handle_pnjl_scan_job_create(_post_idempotency_request(payload_a; idempotency_key=key))
        conflict = PSI.handle_pnjl_scan_job_create(_post_idempotency_request(payload_b; idempotency_key=key))
        conflict_body = _body(conflict)

        @test first.status == 202
        @test conflict.status == 409
        @test conflict_body.status == "error"
        @test conflict_body.error_code == "IDEMPOTENCY_KEY_CONFLICT"
        @test haskey(conflict_body, :message_id)
    end

    _reset_idempotency_fixture!()
end
