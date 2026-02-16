using Test

include("../../../src/pnjl/workflows/TransportWorkflow.jl")
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
        transport_config=TransportIntegrationConfig(p_nodes=12, p_max=4.0),
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

@testset "TransportWorkflow: densities legacy vs models" begin
    T = 0.15
    mu = 0.0
    xi = 0.0

    # Provide tau so workflow doesn't spend time computing it.
    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)

    # Solve equilibrium once; reuse it for both backends.
    base = TransportWorkflow.PNJL.solve(
        TransportWorkflow.PNJL.FixedMu(),
        T,
        mu;
        xi=xi,
        thermo_backend=:legacy,
        p_num=8,
        t_num=4,
        seed_strategy=TransportWorkflow.PNJL.DefaultSeed(phase_hint=:auto),
        iterations=30,
    )

    res_legacy = solve_gap_and_transport(
        T,
        mu;
        xi=xi,
        equilibrium=base,
        thermo_backend=:legacy,
        tau=tau,
        compute_tau=false,
        compute_bulk=false,
        p_num=8,
        t_num=4,
        transport_config=TransportIntegrationConfig(p_nodes=12, p_max=4.0),
    )

    res_models = solve_gap_and_transport(
        T,
        mu;
        xi=xi,
        equilibrium=base,
        thermo_backend=:models,
        tau=tau,
        compute_tau=false,
        compute_bulk=false,
        p_num=8,
        t_num=4,
        transport_config=TransportIntegrationConfig(p_nodes=12, p_max=4.0),
    )

    for k in (:u, :d, :s, :ubar, :dbar, :sbar)
        @test isfinite(getproperty(res_legacy.densities, k))
        @test isfinite(getproperty(res_models.densities, k))
        @test isapprox(
            getproperty(res_models.densities, k),
            getproperty(res_legacy.densities, k);
            rtol=1e-8,
            atol=1e-10,
        )
    end
end

@testset "TransportWorkflow: compute_bulk models backend smoke" begin
    # 选取 μ!=0 避免 n_B=0 导致的 0/0 敏感
    T = 0.15
    mu = 0.3
    xi = 0.0

    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)

    res = solve_gap_and_transport(
        T,
        mu;
        xi=xi,
        thermo_backend=:models,
        tau=tau,
        compute_tau=false,
        compute_bulk=true,
        p_num=6,
        t_num=3,
        solver_kwargs=(iterations=30,),
        transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
    )

    @test res.bulk_coeffs !== nothing
    @test isfinite(res.bulk_coeffs.v_n_sq)
    @test isfinite(res.bulk_coeffs.dμB_dT_sigma)
    @test all(isfinite, res.bulk_coeffs.masses)
end

