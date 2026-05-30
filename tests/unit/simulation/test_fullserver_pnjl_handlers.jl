using Test
using HTTP
using JSON3

const PROJECT_ROOT_FPH = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_FPH, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_FPH, "src", "simulation", "FullServerApp.jl"))
end

const FPH = Main.FullServerApp

function _pnjl_req(payload::Dict)
    return HTTP.Request(
        "POST",
        "/api/modules/pnjl-gap/run",
        ["Content-Type" => "application/json"],
        JSON3.write(payload),
    )
end

function _decode(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

@testset "FullServer PNJL handlers" begin
    @testset "happy path returns sync point payload" begin
        req = _pnjl_req(Dict(
            "params" => Dict(
                "T_mev" => 150.0,
                "mu_mev" => 0.0,
                "xi" => 0.0,
                "p_num" => 8,
                "t_num" => 4,
            ),
        ))
        resp = FPH.handle_pnjl_single_point(req)
        body = _decode(resp)

        @test resp.status == 200
        @test body.status == "ok"
        @test body.result.converged === true
        @test isfinite(body.result.pressure)
        @test body.result.residual_norm === nothing
        @test body.result.entropy === nothing
        @test body.result.energy === nothing
        @test length(body.result.masses) == 3
        @test all(isfinite, body.result.masses)
    end

    @testset "invalid request includes message_id" begin
        req = _pnjl_req(Dict(
            "params" => Dict(
                "T_mev" => 150.0,
            ),
        ))
        resp = FPH.handle_pnjl_single_point(req)
        body = _decode(resp)

        @test resp.status == 400
        @test body.status == "error"
        @test body.error_code == "INVALID_REQUEST"
        @test haskey(body, :message_id)
        @test length(String(body.message_id)) > 0
    end
end
