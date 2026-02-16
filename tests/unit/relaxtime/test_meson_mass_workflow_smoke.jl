using Test

include("../../../src/pnjl/workflows/MesonMassWorkflow.jl")
using .MesonMassWorkflow

@testset "MesonMassWorkflow smoke: solve_gap_and_meson_point (single meson)" begin
    T = 0.15
    mu = 0.0
    xi = 0.0

    res = solve_gap_and_meson_point(
        T,
        mu;
        xi=xi,
        mesons=(:pi,),
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=30,),
        mass_kwargs=(iterations=30,),
    )

    @test haskey(res, :equilibrium)
    @test res.equilibrium.converged isa Bool

    @test length(res.equilibrium.masses) == 3
    @test all(isfinite, res.equilibrium.masses)

    @test haskey(res, :meson_results)
    @test haskey(res.meson_results, :pi)

    rpi = res.meson_results[:pi]
    @test rpi.converged isa Bool
    @test isfinite(rpi.threshold)
    @test rpi.threshold > 0

    # Prefer meaningful numerical output, but keep this smoke resilient.
    @test isfinite(rpi.mass) || !rpi.converged
    @test isfinite(rpi.gap) || !isfinite(rpi.mass)

    # Also validate backend keywords cheaply (no meson solve).
    res_models = solve_gap_and_meson_point(
        T,
        mu;
        xi=xi,
        thermo_backend=:models,
        solver_backend=:legacy,
        mesons=(),
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=30,),
    )

    @test res_models.equilibrium.converged isa Bool
    @test length(res_models.equilibrium.masses) == 3
    @test all(isfinite, res_models.equilibrium.masses)
    @test isempty(res_models.meson_results)
end
