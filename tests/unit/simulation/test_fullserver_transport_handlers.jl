using Test
using HTTP
using JSON3

const PROJECT_ROOT_FTH = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_FTH, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_FTH, "src", "simulation", "FullServerApp.jl"))
end

const FTH = Main.FullServerApp

function _transport_req(payload::Dict)
    return HTTP.Request(
        "POST",
        "/api/modules/transport-point/run",
        ["Content-Type" => "application/json"],
        JSON3.write(payload),
    )
end

function _transport_decode(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

@testset "FullServer transport-point handlers" begin
    @testset "malformed JSON and invalid parameter types return 400" begin
        for raw in ("{", "null", "[]", "1", "\"text\"", "{\"params\":null}", "{\"params\":[]}")
            resp = FTH.handle_transport_point(HTTP.Request("POST", "/api/modules/transport-point/run", [], raw))
            @test resp.status == 400
            @test _transport_decode(resp).error_code == "INVALID_REQUEST"
        end
        for (key, value) in (("T_mev", true), ("tau", []), ("tau", true), ("p_num", 1.5),
                             ("compute_bulk", "invalid"), ("transport", []))
            params = Dict{String,Any}("T_mev" => 150.0, "mu_mev" => 0.0, "tau" => 1.0)
            params[key] = value
            resp = FTH.handle_transport_point(_transport_req(Dict("params" => params)))
            @test resp.status == 400
            @test _transport_decode(resp).error_code == "INVALID_REQUEST"
        end
    end

    @testset "invalid request includes message_id" begin
        req = _transport_req(Dict(
            "params" => Dict(
                "T_mev" => 150.0,
                "mu_mev" => 0.0,
            ),
        ))
        resp = FTH.handle_transport_point(req)
        body = _transport_decode(resp)

        @test resp.status == 400
        @test body.status == "error"
        @test body.error_code == "INVALID_REQUEST"
        @test haskey(body, :message_id)
        @test length(String(body.message_id)) > 0
    end
end
