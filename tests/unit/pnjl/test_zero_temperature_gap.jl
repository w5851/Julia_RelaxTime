using Test
using StaticArrays
using ForwardDiff
using LinearAlgebra: Symmetric, eigvals

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

@testset "PNJL zero-temperature self-consistent state" begin
    model = Models.create_model(:PNJL)
    solver = Models.NLsolveGapSolver(ftol=1e-12, xtol=1e-12)
    controls = (; solver, residual_norm_max=1e-10)

    @testset "Vacuum and Silver Blaze interval" begin
        vacuum = Models.solve_gap(model, 0.0, 0.0; controls...)
        dilute = Models.solve_gap(model, 0.0, 1.25; initial_guess=vacuum, xi=0.3, controls...)
        @test dilute.phi ≈ vacuum.phi rtol=1e-11
        @test vacuum.Phi == vacuum.PhiBar == 0.0
        @test Models.omega(model, dilute, 0.0, 1.25; xi=0.3) ≈
            Models.omega(model, vacuum, 0.0, 0.0) rtol=1e-12
        density = Models.number_densities(model, dilute, 0.0, 1.25; xi=0.3)
        @test all(iszero, density.quark)
        @test all(iszero, density.antiquark)
        @test maximum(abs, Models.gap_residual(model, vacuum, 0.0, 0.0)) < 1e-10
    end

    @testset "Dense state, conjugation, and postprocessing" begin
        dense = Models.solve_gap(model, 0.0, 2.2; controls...)
        residual = Models.gap_residual(model, dense, 0.0, 2.2)
        @test maximum(abs, residual) < 1e-10
        @test residual[4] == residual[5] == 0.0
        @test dense.Phi == dense.PhiBar == 0.0
        point = Models.solve_pnjl_point(T_mev=0.0,
            mu_mev=2.2 * Models.Constants_PNJL.ħc_MeV_fm)
        @test point.converged
        @test point.masses ≈ Models.calculate_mass_vec(model, dense.phi) rtol=1e-8
        @test point.pressure ≈ Models.model_pressure(model, dense, 2.2, 0.0) rtol=1e-10
        conjugate = Models.solve_gap(model, 0.0, -2.2; initial_guess=dense.phi, controls...)
        @test conjugate.phi ≈ dense.phi rtol=1e-10
        rho = Models.number_densities(model, dense, 0.0, 2.2)
        rho_negative = Models.number_densities(model, conjugate, 0.0, -2.2)
        @test rho.quark ≈ rho_negative.antiquark rtol=1e-10
        @test sum(rho.quark) > 0
        @test sum(rho.antiquark) == 0
        pressure_mu = ForwardDiff.derivative(mu -> -Models.omega(model, dense, 0.0, mu), 2.2)
        @test pressure_mu ≈ sum(rho.quark - rho.antiquark) rtol=1e-10
        @test Models.model_rho(model, dense, 2.2, 0.0) ≈ rho.quark - rho.antiquark rtol=1e-10
        arbitrary_polyakov = Models.MeanFieldState(dense.phi; Phi=0.2, PhiBar=0.7)
        @test Models.omega(model, arbitrary_polyakov, 0.0, 2.2) == Models.omega(model, dense, 0.0, 2.2)
        restarted = Models.solve_gap(model, 0.0, 2.2; initial_guess=arbitrary_polyakov, controls...)
        @test restarted.phi ≈ dense.phi rtol=1e-11
        @test restarted.Phi == restarted.PhiBar == 0.0
        asym_mu = SVector(2.2, 2.15, 2.05)
        asym = Models.solve_gap(model, 0.0, asym_mu; initial_guess=dense, xi=-0.3, controls...)
        @test maximum(abs, Models.gap_residual(model, asym, 0.0, asym_mu; xi=-0.3)) < 1e-10
        @test all(isfinite, Models.calculate_mass_vec(model, asym.phi))
    end

    @testset "Inputs and numerical policy" begin
        seed = [-1.84, -1.84, -2.23]
        for mu in (NaN, Inf)
            @test_throws ArgumentError Models.solve_gap(model, 0.0, mu)
        end
        for temperature in (-1.0, NaN, Inf)
            @test_throws ArgumentError Models.solve_gap(model, temperature, 0.0)
        end
        for xi in (-1.0, NaN, Inf)
            @test_throws ArgumentError Models.solve_gap(model, 0.0, 0.0; xi)
        end
        @test_throws ArgumentError Models.solve_gap(model, 0.0, [1.0, 1.0])
        @test_throws ArgumentError Models.solve_gap(model, 0.0, 0.0; initial_guess=[1.0, 2.0])
        @test_throws ArgumentError Models.solve_gap(model, 0.0, 0.0; initial_guess=[1.0, 2.0, 3.0, 4.0])
        @test_throws ArgumentError Models.solve_gap(model, 0.0, 0.0; initial_guess=[NaN, -1.0, -2.0])
        @test_throws ArgumentError Models.solve_gap(model, 0.0, 0.0; iterations=0)
        @test_throws ArgumentError Models.solve_gap(model, 0.0, 0.0; residual_norm_max=Inf)
        @test_throws ArgumentError Models.solve_gap(model, 0.0, 0.0; thermo_quadrature_policy=:tensor_gauss)
        @test_throws ArgumentError Models.solve_gap(model, 0.0, 0.0; thermo_quadrature_rtol=-1.0)
        finite_solver = Models.NLsolveGapSolver(jacobian=:finite, ftol=1e-11)
        finite = Models.solve_gap(model, 0.0, 0.0; initial_guess=seed, solver=finite_solver)
        @test maximum(abs, Models.gap_residual(model, finite, 0.0, 0.0)) < 1e-9
        @test_throws ErrorException Models.solve_gap(model, 0.0, 0.0;
            initial_guess=[-0.1, -0.1, -0.1], solver, iterations=1)
    end

    @testset "First-order branches and minimum selection" begin
        mu = 1.85
        low = Models.solve_gap(model, 0.0, mu;
            initial_guess=Models.SeedStrategies.HADRON_SEED_5, controls...)
        high = Models.solve_gap(model, 0.0, mu;
            initial_guess=Models.SeedStrategies.HIGH_DENSITY_SEED_5, controls...)
        saddle = Models.solve_gap(model, 0.0, mu; initial_guess=[-0.5, -0.5, -1.2], controls...)
        chosen = Models.solve_gap(model, 0.0, mu; controls...)
        @test maximum(abs, low.phi - high.phi) > 1.0
        @test Models.omega(model, high, 0.0, mu) < Models.omega(model, low, 0.0, mu)
        @test chosen.phi ≈ high.phi rtol=1e-9
        @test maximum(abs, Models.gap_residual(model, saddle, 0.0, mu)) < 1e-10
        potential = phi -> Models.omega(model, phi, 0.0, mu)
        @test minimum(eigvals(Symmetric(ForwardDiff.hessian(potential, saddle.phi)))) < 0
        @test minimum(eigvals(Symmetric(ForwardDiff.hessian(potential, chosen.phi)))) > 0
    end
end
