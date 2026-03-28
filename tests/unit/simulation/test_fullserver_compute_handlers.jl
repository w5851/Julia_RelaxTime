using Test
using HTTP
using JSON3

const PROJECT_ROOT_FCH = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_FCH, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_FCH, "src", "simulation", "FullServerApp.jl"))
end

const FCH = Main.FullServerApp

function _make_compute_req(payload::Dict)
    return HTTP.Request(
        "POST",
        "/compute",
        ["Content-Type" => "application/json"],
        JSON3.write(payload),
    )
end

function _parse_body(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

@testset "FullServer compute handler" begin
    @testset "runtime error response is sanitized" begin
        req = _make_compute_req(Dict(
            "p1x" => 0.5,
            "p1y" => 0.0,
            "p1z" => 1.8,
            "p2x" => -0.5,
            "p2y" => 0.0,
            "p2z" => -1.8,
            "m1" => 1.52,
            "m2" => 1.52,
            "m3" => 1.52,
        ))

        resp = FCH.handle_compute(req)
        body = _parse_body(resp)

        @test resp.status == 500
        @test body.success == false
        @test body.data === nothing
        @test body.error_code == "COMPUTATION_ERROR"
        @test body.error == "Computation failed"
        @test haskey(body, :message_id)
        @test length(String(body.message_id)) > 0
    end
end
