using Test
using ForwardDiff

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _CONSTANTS_PATH = normpath(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PATH)
end

const _RELAXTIME_PATH = normpath(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH)
end

using Main.Constants_PNJL: G_fm2, K_fm5
using Main.RelaxTime.OneLoopIntegrals: B0, A
using Main.RelaxTime.GaussLegendre: gauleg
using Main.RelaxTime.PolarizationAniso: polarization_with_width
using Main.RelaxTime.EffectiveCouplings: calculate_effective_couplings
using Main.RelaxTime.MesonPropagator: meson_propagator_simple

@testset "meson phase-shift dual compatibility smoke" begin
    T0 = 0.18
    Φ = 0.2
    Φbar = 0.2
    m1 = 0.3
    m2 = 0.3
    μ1 = 0.0
    μ2 = 0.0
    γ = 0.05

    @testset "B0 accepts Dual on k=0 and k>0 branches" begin
        f_k0(x) = first(B0(x, 0.0, m1, μ1, m2, μ2, T0; Φ=Φ, Φbar=Φbar))
        f_kp(x) = first(B0(x, 0.4, m1, μ1, m2, μ2, T0; Φ=Φ, Φbar=Φbar))

        d_k0 = ForwardDiff.derivative(f_k0, 0.2)
        d_kp = ForwardDiff.derivative(f_kp, 0.2)

        @test isfinite(d_k0)
        @test isfinite(d_kp)
    end

    nodes_p, weights_p = gauleg(0.0, 20.0, 24)
    A1 = A(m1, μ1, T0, Φ, Φbar, nodes_p, weights_p)
    A2 = A(m2, μ2, T0, Φ, Φbar, nodes_p, weights_p)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, 1.2, 1.0)

    @testset "polarization_with_width and meson_propagator_simple accept Dual for xi=0" begin
        re_fun(ω) = begin
            Πr, Πi = polarization_with_width(:P, ω, γ, 0.4, m1, m2, μ1, μ2, T0, Φ, Φbar, 0.0, A1, A2, 0)
            real(meson_propagator_simple(:pi, K_coeffs, Complex(Πr, Πi)))
        end
        im_fun(ω) = begin
            Πr, Πi = polarization_with_width(:P, ω, γ, 0.4, m1, m2, μ1, μ2, T0, Φ, Φbar, 0.0, A1, A2, 0)
            imag(meson_propagator_simple(:pi, K_coeffs, Complex(Πr, Πi)))
        end

        d_re = ForwardDiff.derivative(re_fun, 0.2)
        d_im = ForwardDiff.derivative(im_fun, 0.2)

        @test isfinite(d_re)
        @test isfinite(d_im)
    end

    @testset "anisotropic correction also accepts Dual after minimal λ-generic upgrade" begin
        f_xi0(ω) = first(polarization_with_width(:P, ω, γ, 0.0, m1, m2, μ1, μ2, T0, Φ, Φbar, 0.3, A1, A2, 0))
        d_xi0 = ForwardDiff.derivative(f_xi0, 0.2)
        @test isfinite(d_xi0)

        re_fun_xi(ω) = begin
            Πr, Πi = polarization_with_width(:P, ω, γ, 0.4, m1, m2, μ1, μ2, T0, Φ, Φbar, 0.3, A1, A2, 0)
            real(meson_propagator_simple(:pi, K_coeffs, Complex(Πr, Πi)))
        end
        @test isfinite(ForwardDiff.derivative(re_fun_xi, 0.2))
    end
end
