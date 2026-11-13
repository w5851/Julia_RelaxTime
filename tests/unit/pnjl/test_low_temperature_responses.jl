using Test
using StaticArrays
using ForwardDiff

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

@testset "Low-temperature Fermi surface responses" begin
    pi_kernel = Models.PNJLIntegrals
    td = Models.PNJLChiBTaylorDiff
    qopts = (; rtol=1e-9, atol=1e-12, maxevals=200000)

    @testset "Massless deconfined gas and derivatives" begin
        masses = SVector(0.0, 0.0, 0.0)
        mu0 = 1.8
        delta = td.taylor_variable(0.0, 4)
        for (T, xi) in ((0.0, 0.0), (0.00005, 0.3), (0.0005, -0.3), (0.1, 0.0))
            mu = mu0 + delta
            actual = -pi_kernel.calculate_log_sum_rs_reduced_adaptive(
                masses, 1.0, 1.0, SVector(mu, mu, mu), T, xi; qopts...)
            exact = 9 / (12pi^2 * sqrt(1 + xi)) *
                (mu^4 + 2pi^2 * mu^2 * T^2 + 7pi^4 * T^4 / 15)
            for n in 0:4
                @test td.nth_derivative_from_series(actual, n) ≈
                    td.nth_derivative_from_series(exact, n) rtol=2e-8 atol=2e-10
            end
        end
    end

    @testset "Massive sea zero-temperature limit" begin
        masses = SVector(1.5, 1.5, 2.5)
        mu0 = 1.8
        mu = mu0 + td.taylor_variable(0.0, 4)
        q = sqrt(mu0^2 - masses[1]^2)
        coeff = 6 / pi^2
        reference = (coeff * q^3 / 3, coeff * mu0 * q,
            coeff * (2mu0^2 - masses[1]^2) / q,
            coeff * mu0 * (2mu0^2 - 3masses[1]^2) / q^3)
        zero_sea = -pi_kernel.calculate_log_sum_rs_reduced_adaptive(
            masses, 0.1, 0.1, SVector(mu, mu, mu), 0.0, 0.0; qopts...)
        low_sea = -pi_kernel.calculate_log_sum_rs_reduced_adaptive(
            masses, 0.1, 0.1, SVector(mu, mu, mu), 0.01 / 197.327, 0.0; qopts...)
        for n in 1:4
            @test td.nth_derivative_from_series(zero_sea, n) ≈ reference[n] rtol=1e-11
            @test td.nth_derivative_from_series(low_sea, n) ≈ reference[n] rtol=1e-6 atol=1e-9
        end
        for mu0 in (1.25, 1.5)
            mu = mu0 + td.taylor_variable(0.0, 2)
            pressure = -pi_kernel.calculate_log_sum_rs_reduced_adaptive(
                masses, 0.1, 0.1, SVector(mu, mu, mu), 0.0, 0.0; qopts...)
            # At onset P ~ (mu-M)_+^(5/2): orders 0:2 exist and vanish.
            # Orders >= 3 at exactly mu=M are not defined by this test.
            for n in 0:2
                @test td.nth_derivative_from_series(pressure, n) == 0.0
            end
        end
    end

    @testset "AD refinement sees a zero primal integrand" begin
        derivative = ForwardDiff.derivative(0.0) do amplitude
            pi_kernel.integrate_rs_reduced_radial(q -> amplitude * exp(-q),
                0.3, 0.0; rtol=1e-10, atol=1e-12).value
        end
        @test derivative ≈ 1.0 rtol=1e-10
        zero_series = td.taylor_variable(0.0, 4)
        polynomial = zero_series * zero_series * zero_series * zero_series
        result = pi_kernel.integrate_rs_reduced_radial(q -> polynomial * exp(-q),
            0.3, 0.0; rtol=1e-10, atol=1e-12).value
        @test td.nth_derivative_from_series(result, 4) ≈ 24.0 rtol=1e-10
    end

    @testset "Charge conjugation and the degenerate-coordinate boundary" begin
        masses = SVector(1.5, 1.5, 2.5)
        delta = td.taylor_variable(0.0, 4)
        for ratio in (63.9, 64.1, 6000.0)
            T = 0.3 / ratio
            mu = 1.8 + delta
            positive = -pi_kernel.calculate_log_sum_rs_reduced_adaptive(
                masses, 0.12, 0.28, SVector(mu, mu, mu), T, 0.3; qopts...)
            negative = -pi_kernel.calculate_log_sum_rs_reduced_adaptive(
                masses, 0.28, 0.12, SVector(-mu, -mu, -mu), T, 0.3; qopts...)
            for n in 0:4
                @test td.nth_derivative_from_series(positive, n) ≈
                    td.nth_derivative_from_series(negative, n) rtol=1e-10 atol=1e-10
            end
            if ratio < 65
                # The original q-integral checks the identity on both sides of
                # its dispatch boundary, including unequal Polyakov fields.
                original = sum(masses) * zero(mu)
                for mass in masses
                    original += pi_kernel.integrate_rs_reduced_radial(
                        q -> 2T * q^2 * pi_kernel.calculate_log_term(
                            sqrt(q^2 + mass^2), mu, T, 0.12, 0.28),
                        mass, mu; thermal_scale=T, qopts...).value
                end
                original /= 2pi^2 * sqrt(1.3)
                for n in 0:4
                    @test td.nth_derivative_from_series(positive, n) ≈
                        td.nth_derivative_from_series(original, n) rtol=2e-6 atol=2e-9
                end
            end
        end
        @test_throws ArgumentError pi_kernel.integrate_rs_reduced_radial(
            exp, 1.5, 1.8; thermal_scale=0.001, maxevals=48)
    end
end

@testset "Adaptive derivatives preserve the base context quadrature" begin
    model = Models.create_model(:PNJL)
    td = Models.PNJLChiBTaylorDiff
    T = 60 / 197.327
    mu = SVector(2.2, 2.2, 2.2)
    quadrature = (; thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8, thermo_quadrature_atol=1e-11,
        thermo_quadrature_maxevals=200000)
    context = td.build_taylor_base_context(T, mu; model, p_num=12, t_num=4,
        thermo_quadrature_kwargs=quadrature)
    kwargs = (; model, p_num=12, t_num=4, base_context=context)
    pressure = td.pressure_series_direction(T, mu, SVector(1.0, 1.0, 1.0); order=2, kwargs...)
    density = Models.number_densities(model, context.x0, T, mu; quadrature...)
    @test td.nth_derivative_from_series(pressure, 0) ≈
        Models.model_pressure(model, context.x0, mu, T; quadrature...) rtol=1e-11
    @test td.nth_derivative_from_series(pressure, 1) ≈ sum(density.quark - density.antiquark) rtol=2e-7
    @test Models.chi2_B(T, 3mu[1]; kwargs...) ≈
        td.nth_derivative_from_series(pressure, 2) / (9T^2) rtol=2e-7
    md = Models.mass_derivatives(T, mu[1]; kwargs...)
    @test md.masses ≈ Models.calculate_mass_vec(model, SVector{3}(context.x0[1:3])) rtol=1e-12
    @test all(isfinite, md.dM_dT)
    @test all(isfinite, md.dM_dmu)
    @test_throws ArgumentError td.pressure_series_direction(T, mu, SVector(1.0, 0.0, 0.0);
        order=1, kwargs..., thermo_quadrature_kwargs=(; thermo_quadrature_policy=:tensor_gauss))
    @test_throws ArgumentError Models.solve_gap(model, 0.0, mu; quadrature...)
    @test_throws ArgumentError Models.chi2_B(0.0, 3mu[1]; thermo_quadrature_kwargs=quadrature)
    @test_throws ArgumentError Models.chi_BQS(T, 3mu[1], 0.0, 0.0;
        orders=(1, 1, 0), thermo_quadrature_kwargs=quadrature)
end
