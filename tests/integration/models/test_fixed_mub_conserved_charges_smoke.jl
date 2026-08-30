using Test

const _FIXED_MUB_INTEGRATION_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(_FIXED_MUB_INTEGRATION_ROOT, "src", "models", "Models.jl"))
end

@testset "FixedMuBConservedCharges unified joint solve smoke" begin
    hbarc_MeV_fm = 197.3269804
    T_fm = 170.0 / hbarc_MeV_fm
    muB_fm = 240.0 / hbarc_MeV_fm
    model = Models.create_model(:PNJL)
    mode = Models.FixedMuBConservedCharges(muB_fm, 0.4, 0.0)

    result = Models.solve(
        model,
        mode,
        T_fm;
        p_num=8,
        t_num=4,
        residual_norm_max=1e-5,
        iterations=200,
    )

    @test result isa Models.SolverResult
    @test result.converged
    @test result.residual_norm <= 1e-5
    @test length(result.solution) == 8
    @test length(result.x_state) == 5
    @test length(result.mu_vec) == 3

    conserved_mu = Models.conserved_mu_from_flavor(result.mu_vec...)
    @test conserved_mu.mu_B ≈ muB_fm atol=1e-8

    rho_vec = Models.model_rho(
        model,
        result.x_state,
        result.mu_vec,
        T_fm;
        p_num=8,
        t_num=4,
        xi=0.0,
    )
    densities = Models.conserved_densities_from_flavor(rho_vec)
    @test abs(densities.rho_S) <= 1e-8
    @test densities.rho_B > 1e-8
    @test densities.rho_Q / densities.rho_B ≈ 0.4 atol=1e-7

    # Public particle/antiparticle densities must reproduce the same net flavor
    # density; using `.quark` alone would be a different physical quantity.
    particle_density = Models.number_densities(
        model,
        result.x_state,
        T_fm,
        result.mu_vec;
        p_num=8,
        t_num=4,
        xi=0.0,
    )
    net_density = particle_density.quark - particle_density.antiquark
    @test net_density ≈ rho_vec rtol=1e-10 atol=1e-12
end
