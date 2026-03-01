using Test

if !isdefined(Main, :Models)
    if !isdefined(Main, :Models)
    include("../../../src/models/Models.jl")
end
end
Main.Models.transport_workflow_module()
using .TransportWorkflow

@testset "TransportWorkflow smoke: solver_backend switch (legacy vs models)" begin
    T = 0.9
    mu = 0.0
    xi = 0.0

    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)

    common_kwargs = (
        xi=xi,
        thermo_backend=:models,
        tau=tau,
        compute_tau=false,
        compute_bulk=false,
        p_num=8,
        t_num=4,
        transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
    )

    res_legacy = TransportWorkflow.solve_gap_and_transport(T, mu; common_kwargs..., solver_backend=:legacy, solver_kwargs=(iterations=30,))

    models_solver = Main.Models.NLsolveGapSolver(method=:trust_region, jacobian=:finite, xtol=1e-10, ftol=1e-10)
    res_models = TransportWorkflow.solve_gap_and_transport(T, mu;
        common_kwargs...,
        solver_backend=:models,
        models_solver=models_solver,
        models_residual_norm_max=1e-4,
        seed_state=TransportWorkflow.HADRON_SEED_5,
    )

    @test res_legacy.equilibrium.converged isa Bool
    @test res_models.equilibrium.converged isa Bool

    @test all(isfinite, res_legacy.masses)
    @test all(isfinite, res_models.masses)

    @test isfinite(res_legacy.transport.eta)
    @test isfinite(res_models.transport.eta)
    @test isfinite(res_legacy.transport.sigma)
    @test isfinite(res_models.transport.sigma)

    # Comparability: omega at the two equilibria should be close (smoke-level tolerance).
    m = Main.Models.create_model(:PNJL)
    ωL = Main.Models.omega(m, res_legacy.equilibrium.x_state, T, Main.Models.normalize_mu_vec(mu); xi=xi, p_num=8, t_num=4)
    ωM = Main.Models.omega(m, res_models.equilibrium.x_state, T, Main.Models.normalize_mu_vec(mu); xi=xi, p_num=8, t_num=4)

    @test isfinite(ωL)
    @test isfinite(ωM)
    scale = max(1.0, abs(ωL))
    @test abs(ωM - ωL) / scale <= 1e-2
end
