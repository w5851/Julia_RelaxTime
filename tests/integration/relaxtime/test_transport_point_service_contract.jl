using Test
using HTTP
using JSON3

const PROJECT_ROOT_TPSC = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :FullServerApp)
    if !isdefined(Main, :Models)
        include(joinpath(PROJECT_ROOT_TPSC, "src", "models", "Models.jl"))
    end
    include(joinpath(PROJECT_ROOT_TPSC, "src", "simulation", "FullServerApp.jl"))
end

const TPSC = Main.FullServerApp

function _transport_point_req(payload::Dict)
    return HTTP.Request(
        "POST",
        "/api/modules/transport-point/run",
        ["Content-Type" => "application/json"],
        JSON3.write(payload),
    )
end

function _transport_point_body(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

@testset "transport-point service contract" begin
    @testset "happy path returns sync point payload" begin
        req = _transport_point_req(Dict(
            "params" => Dict(
                "T_mev" => 150.0,
                "mu_mev" => 0.0,
                "xi" => 0.0,
                "tau" => 1.0,
                "compute_bulk" => false,
                "p_num" => 8,
                "t_num" => 4,
                "transport" => Dict(
                    "p_nodes" => 8,
                    "p_max" => 3.5,
                ),
            ),
        ))

        resp = TPSC.handle_transport_point(req)
        body = _transport_point_body(resp)

        @test resp.status == 200
        @test body.status == "ok"
        @test body.result.inputs.T_mev == 150.0
        @test body.result.inputs.mu_mev == 0.0
        @test body.result.inputs.tau.u == 1.0
        @test body.result.equilibrium.converged isa Bool
        @test length(body.result.equilibrium.masses) == 3
        @test isfinite(body.result.thermo_background.pressure)
        @test isfinite(body.result.transport.eta)
        @test isfinite(body.result.transport.sigma)
        @test body.result.reproducibility.physics_profile isa AbstractString
    end

    @testset "invalid tau payload returns INVALID_REQUEST" begin
        req = _transport_point_req(Dict(
            "params" => Dict(
                "T_mev" => 150.0,
                "mu_mev" => 0.0,
                "tau" => Dict(
                    "u" => 1.0,
                    "d" => 1.0,
                    "s" => 1.0,
                ),
            ),
        ))

        resp = TPSC.handle_transport_point(req)
        body = _transport_point_body(resp)

        @test resp.status == 400
        @test body.status == "error"
        @test body.error_code == "INVALID_REQUEST"
        @test haskey(body, :message_id)
    end
end
