using Test

include("../../src/pnjl/workflows/TransportWorkflow.jl")
using .TransportWorkflow

@testset "TransportWorkflow: gap -> transport (single point)" begin
    T = 0.15
    mu = 0.0
    xi = 0.0

    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)

    res = solve_gap_and_transport(
        T,
        mu;
        xi=xi,
        tau=tau,
        compute_tau=false,
        compute_bulk=false,   # 单测避免导数带来的多次求解
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=30,),
        transport_kwargs=(p_nodes=12, p_max=4.0,),
    )

    @test haskey(res, :equilibrium)
    @test res.equilibrium.converged isa Bool
    @test isfinite(res.thermo_params.T)
    @test isfinite(res.thermo_params.Φ)
    @test isfinite(res.thermo_params.Φbar)

    @test length(res.masses) == 3
    @test all(isfinite, res.masses)

    @test isfinite(res.transport.eta)
    @test isfinite(res.transport.sigma)
    @test res.transport.eta >= 0
    @test res.transport.sigma >= 0
end
