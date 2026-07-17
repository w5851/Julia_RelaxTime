using Test
using QuadGK
using StaticArrays

const PROJECT_ROOT_PHASE_QUAD = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_PHASE_QUAD, "src", "models", "Models.jl"))
end

function _tensor_free_2d_oracle(masses, Φ, Φbar, mu_vec, T_fm, xi)
    total, _ = quadgk(0.0, 1.0; rtol=2e-8, atol=1e-10) do c
        flavor_total = 0.0
        for i in eachindex(masses)
            mass = masses[i]
            mu = mu_vec[i]
            p_fermi = abs(mu) > abs(mass) ?
                sqrt((mu * mu - mass * mass) / (1 + xi * c * c)) : 0.0
            integrand = p -> begin
                energy = sqrt(mass * mass + p * p * (1 + xi * c * c))
                return p * p * (-2 * T_fm) *
                    Models.PNJLIntegrals.calculate_log_term(energy, mu, T_fm, Φ, Φbar)
            end
            value, _ = if p_fermi > 0.0
                quadgk(integrand, 0.0, p_fermi, Inf; rtol=2e-9, atol=1e-11)
            else
                quadgk(integrand, 0.0, Inf; rtol=2e-9, atol=1e-11)
            end
            flavor_total += value
        end
        return 2 * flavor_total / (2 * π)^2
    end
    return total
end

@testset "PNJL phase RS-reduced quadrature validation" begin
    masses = SVector(0.3, 0.31, 0.5)
    Φ, Φbar = 0.3, 0.25
    xi_values = (-0.5, 0.0, 0.5)

    for T_MeV in (1.0, 5.0, 10.0, 20.0, 40.0, 60.0, 120.0, 200.0, 240.0)
        T_fm = T_MeV / 197.327
        mu_vec = SVector(0.8, 0.75, 0.9)
        for xi in xi_values
            adaptive = Models.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive_with_error(
                masses, Φ, Φbar, mu_vec, T_fm, xi;
                rtol=1e-9, atol=1e-11,
            )
            oracle = _tensor_free_2d_oracle(masses, Φ, Φbar, mu_vec, T_fm, xi)
            @test adaptive.value ≈ oracle rtol=3e-7 atol=2e-10
            @test adaptive.error <= 1e-7 * max(abs(adaptive.value), 1.0)
        end
    end

    T_fm = 20.0 / 197.327
    for mu_vec in (
        SVector(0.1, 0.1, 0.1),
        SVector(0.3, 0.31, 0.5),
        SVector(0.8, 0.75, 0.9),
    ), xi in xi_values
        adaptive = Models.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive(
            masses, Φ, Φbar, mu_vec, T_fm, xi;
            rtol=1e-9, atol=1e-11,
        )
        oracle = _tensor_free_2d_oracle(masses, Φ, Φbar, mu_vec, T_fm, xi)
        @test adaptive ≈ oracle rtol=3e-7 atol=2e-10
    end
end
