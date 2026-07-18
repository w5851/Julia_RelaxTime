using Test
using StaticArrays

const PROJECT_ROOT_PHASE_QUAD = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_PHASE_QUAD, "src", "models", "Models.jl"))
end

using Main.GaussLegendre: standard_nodes_weights

function _fixed_gauss_legendre(f, a, b, nodes, weights)
    b > a || return 0.0
    half = (b - a) / 2
    center = (a + b) / 2
    total = 0.0
    @inbounds for i in eachindex(nodes)
        total += weights[i] * f(muladd(half, nodes[i], center))
    end
    return half * total
end

function _tensor_fixed_2d_oracle(
    masses,
    Φ,
    Φbar,
    mu_vec,
    T_fm,
    xi;
    n_cos::Int=48,
    n_radial::Int=192,
    tail_exponent::Float64=40.0,
)
    cos_nodes, cos_weights = standard_nodes_weights(n_cos)
    radial_nodes, radial_weights = standard_nodes_weights(n_radial)
    total = 0.0

    @inbounds for j in eachindex(cos_nodes)
        c = (cos_nodes[j] + 1) / 2
        c_weight = cos_weights[j] / 2
        flavor_total = 0.0
        for i in eachindex(masses)
            mass = masses[i]
            mu = mu_vec[i]
            angular_scale = 1 + xi * c * c
            p_fermi = abs(mu) > abs(mass) ?
                sqrt((mu * mu - mass * mass) / angular_scale) : 0.0
            energy_cutoff = max(abs(mu), abs(mass)) + tail_exponent * T_fm
            p_max = sqrt(max(energy_cutoff * energy_cutoff - mass * mass, 0.0) / angular_scale)
            integrand = p -> begin
                energy = sqrt(mass * mass + p * p * angular_scale)
                return p * p * (-2 * T_fm) *
                    Models.PNJLIntegrals.calculate_log_term(energy, mu, T_fm, Φ, Φbar)
            end
            flavor_total += _fixed_gauss_legendre(
                integrand, 0.0, p_fermi, radial_nodes, radial_weights,
            )
            flavor_total += _fixed_gauss_legendre(
                integrand, p_fermi, p_max, radial_nodes, radial_weights,
            )
        end
        total += c_weight * flavor_total
    end
    return 2 * total / (2 * π)^2
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
            oracle = _tensor_fixed_2d_oracle(masses, Φ, Φbar, mu_vec, T_fm, xi)
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
        oracle = _tensor_fixed_2d_oracle(masses, Φ, Φbar, mu_vec, T_fm, xi)
        @test adaptive ≈ oracle rtol=3e-7 atol=2e-10
    end

    @testset "固定高节点 oracle 自收敛" begin
        for (T_MeV, xi) in ((1.0, -0.5), (60.0, 0.0), (240.0, 0.5))
            T_fm = T_MeV / 197.327
            mu_vec = SVector(0.8, 0.75, 0.9)
            medium = _tensor_fixed_2d_oracle(
                masses, Φ, Φbar, mu_vec, T_fm, xi;
                n_cos=48, n_radial=192, tail_exponent=40.0,
            )
            strict = _tensor_fixed_2d_oracle(
                masses, Φ, Φbar, mu_vec, T_fm, xi;
                n_cos=64, n_radial=256, tail_exponent=44.0,
            )
            @test medium ≈ strict rtol=5e-8 atol=2e-11
        end
    end
end
