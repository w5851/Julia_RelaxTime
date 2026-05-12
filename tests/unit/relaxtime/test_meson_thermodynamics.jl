using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _RELAXTIME_PATH_MT = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH_MT)
end

using Main.RelaxTime.MesonThermodynamics: bosonic_log_pressure_factor,
                                          phase_shift_meson_pressure,
                                          phase_shift_meson_pressure_summary,
                                          stable_meson_pressure,
                                          stable_meson_pressure_summary,
                                          strict_bw_meson_pressure,
                                          strict_bw_meson_pressure_summary
using Main.Constants_PNJL: Λ_inv_fm

@testset "MesonThermodynamics stable pressure" begin
    @test bosonic_log_pressure_factor(1.0, 0.0, 0.2) > 0.0
    @test bosonic_log_pressure_factor(1.0, 0.0, 0.0) == 0.0

    T = 0.3
    pressure = stable_meson_pressure(0.0, T; μ=-1e-6, degeneracy=2, qmax=12.0, num_q_nodes=512)
    expected = 2.0 * π^2 / 90.0 * T^4
    @test pressure ≈ expected rtol=3e-3

    summary = stable_meson_pressure_summary(0.14, 0.49, 0.20; num_q_nodes=128)
    @test summary.P_pi > 0.0
    @test summary.P_K > 0.0
    @test summary.P_pi > summary.P_K
    @test summary.P_meson ≈ summary.P_pi + summary.P_K rtol=1e-12
end

@testset "MesonThermodynamics supports sigma_pi channel" begin
    qp = (m=(u=0.098, d=0.098, s=0.42), μ=(u=0.0, d=0.0, s=0.0), A=(u=0.1, d=0.1, s=0.08))
    tp = (T=0.18, Φ=0.25, Φbar=0.25, ξ=0.0)

    sigma = phase_shift_meson_pressure(
        :sigma_pi,
        qp,
        tp;
        scheme=:current,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )
    @test isfinite(sigma.pressure)
    @test sigma.pressure ≈ sigma.pressure_qp + sigma.pressure_ld rtol=1e-12

    summary = phase_shift_meson_pressure_summary(
        qp,
        tp;
        pi_channel=:pi,
        k_channel=:sigma_pi,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )
    @test summary.k_channel == :sigma_pi
    @test summary.d_K == 1
    @test isfinite(summary.P_K)
    @test summary.P_K ≈ summary.P_K_qp + summary.P_K_ld rtol=1e-12
end

@testset "MesonThermodynamics strict BW pressure" begin
    T = 0.20
    stable = stable_meson_pressure(0.14, T; degeneracy=3, qmax=12.0, num_q_nodes=48)
    bw_zero = strict_bw_meson_pressure(
        0.14,
        0.0,
        T;
        degeneracy=3,
        qmax=12.0,
        q_nodes=48,
        omega_max=10.0,
        omega_nodes=24,
    )
    @test bw_zero.mode == :stable_limit
    @test bw_zero.pressure ≈ stable rtol=1e-10

    bw_finite = strict_bw_meson_pressure(
        0.14,
        0.08,
        T;
        degeneracy=3,
        qmax=12.0,
        q_nodes=24,
        omega_max=10.0,
        omega_nodes=24,
    )
    @test bw_finite.mode == :strict_bw_reduced
    @test isfinite(bw_finite.pressure)
    @test bw_finite.pressure > 0.0

    summary = strict_bw_meson_pressure_summary(
        0.14, 0.08,
        0.49, 0.12,
        T;
        qmax=12.0,
        q_nodes=24,
        omega_max=10.0,
        omega_nodes=24,
    )
    @test summary.P_pi > 0.0
    @test summary.P_K > 0.0
    @test summary.P_meson ≈ summary.P_pi + summary.P_K rtol=1e-12
end

@testset "MesonThermodynamics phase-shift pressure schemes" begin
    qp = (m=(u=0.098, d=0.098, s=0.42), μ=(u=0.0, d=0.0, s=0.0), A=(u=0.1, d=0.1, s=0.08))
    tp = (T=0.18, Φ=0.25, Φbar=0.25, ξ=0.0)

    current = phase_shift_meson_pressure(
        :pi_plus,
        qp,
        tp;
        scheme=:current,
        degeneracy=1,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )
    gbu = phase_shift_meson_pressure(
        :pi_plus,
        qp,
        tp;
        scheme=:gbu_reference,
        degeneracy=1,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )

    @test current.scheme == :phase_shift_current
    @test gbu.scheme == :phase_shift_gbu_reference
    @test isfinite(current.pressure)
    @test isfinite(gbu.pressure)
    @test current.pressure > 0.0
    @test gbu.pressure > 0.0
    @test current.pressure ≈ current.pressure_qp + current.pressure_ld rtol=1e-12
    @test gbu.pressure ≈ gbu.pressure_qp + gbu.pressure_ld rtol=1e-12
    @test isfinite(current.pressure_qp)
    @test current.pressure_ld ≥ 0.0
    @test current.ld_cutoff ≈ min(Λ_inv_fm, 4.0) rtol=1e-12
    @test current.ld_cutoff_mode == :match_model_lambda
    @test !isapprox(current.pressure, gbu.pressure; rtol=1e-6, atol=1e-10)

    ld_trimmed = phase_shift_meson_pressure(
        :pi_plus,
        qp,
        tp;
        scheme=:current,
        degeneracy=1,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
        ld_cutoff=1.0,
        ld_cutoff_mode=:explicit,
    )
    @test ld_trimmed.pressure ≈ ld_trimmed.pressure_qp + ld_trimmed.pressure_ld rtol=1e-12
    @test ld_trimmed.pressure_ld < current.pressure_ld
    @test ld_trimmed.pressure < current.pressure

    summary = phase_shift_meson_pressure_summary(
        qp,
        tp;
        scheme=:gbu_reference,
        qmax=4.0,
        q_nodes=6,
        omega_min=0.05,
        omega_max=3.0,
        omega_nodes=6,
    )
    @test summary.scheme == :phase_shift_gbu_reference
    @test summary.P_pi > 0.0
    @test summary.P_K > 0.0
    @test summary.P_pi ≈ summary.P_pi_qp + summary.P_pi_ld rtol=1e-12
    @test summary.P_K ≈ summary.P_K_qp + summary.P_K_ld rtol=1e-12
    @test summary.P_meson ≈ summary.P_meson_qp + summary.P_meson_ld rtol=1e-12
    @test summary.P_meson ≈ summary.P_pi + summary.P_K rtol=1e-12
    @test summary.ld_cutoff ≈ min(Λ_inv_fm, 4.0) rtol=1e-12
    @test summary.ld_cutoff_mode == :match_model_lambda
end
