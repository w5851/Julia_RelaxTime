using Test
using HTTP
using JSON3
using Logging: NullLogger, with_logger

const PROJECT_ROOT_FCH = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_FCH, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_FCH, "src", "simulation", "FullServerApp.jl"))
end

const FCH = Main.FullServerApp
if !isdefined(Main, :HTTPServer)
    include(joinpath(PROJECT_ROOT_FCH, "src", "simulation", "HTTPServer.jl"))
end

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
    @testset "both servers reject malformed and non-finite inputs" begin
        valid = Dict{String,Any}(
            "p1x" => 0.5, "p1y" => 0.0, "p1z" => 1.8,
            "p2x" => -0.5, "p2y" => 0.0, "p2z" => -1.8,
            "m1" => 1.52, "m2" => 1.52, "m3" => 1.52, "m4" => 1.52,
        )
        invalid_bodies = ["{", "null", "[]", ""]
        for value in (nothing, true, "1.52", -1.0, 0.0)
            push!(invalid_bodies, JSON3.write(merge(valid, Dict("m1" => value))))
        end
        push!(invalid_bodies, replace(JSON3.write(valid), "1.52" => "1e999"; count=1))
        for handler in (FCH.handle_compute, Main.HTTPServer.handle_compute)
            for raw in invalid_bodies
                response = with_logger(NullLogger()) do
                    handler(HTTP.Request("POST", "/compute", ["Content-Type" => "application/json"], raw))
                end
                body = _parse_body(response)
                @test response.status == 400
                @test body.success == false
                @test body.data === nothing
                @test body.error_code == "INVALID_INPUT"
            end
            response = handler(_make_compute_req(valid))
            @test response.status == 200
            @test _parse_body(response).data.validation.is_valid
        end
    end

    @testset "missing fields are input errors" begin
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

        @test resp.status == 400
        @test body.success == false
        @test body.data === nothing
        @test body.error_code == "INVALID_INPUT"
        @test body.error == "Invalid request parameters"
        @test haskey(body, :message_id)
        @test length(String(body.message_id)) > 0
    end
end
