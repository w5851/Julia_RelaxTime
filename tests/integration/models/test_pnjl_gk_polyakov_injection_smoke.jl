using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :solve_gap))
    Base.include(Main, _models_entry)
end

@inline function _polyakov_U_manual(consts, Φ, Φbar, T_fm)
    T0 = consts.T0_inv_fm
    T_ratio = T0 / T_fm
    Ta = consts.a0 + consts.a1 * T_ratio + consts.a2 * T_ratio^2
    Tb = consts.b3 * T_ratio^3
    value = 1 - 6 * Φbar * Φ + 4 * (Φbar^3 + Φ^3) - 3 * (Φbar * Φ)^2
    value <= 0 && (value = 1e-16)
    value < 1e-16 && (value = 1e-16)
    return T_fm^4 * (-0.5 * Ta * Φbar * Φ + Tb * log(value))
end

@testset "Models PNJL G/K + Polyakov params injection" begin
    m_default = Models.create_model(:PNJL; profile="default", physics_profile="default")
    m_custom = Models.create_model(:PNJL; profile="unittest_gk_polyakov", physics_profile="default")

    @test m_custom.consts.G_fm2 != m_default.consts.G_fm2
    @test m_custom.consts.K_fm5 != m_default.consts.K_fm5
    @test m_custom.consts.T0_inv_fm != m_default.consts.T0_inv_fm

    # Ensure b4 is injected (even if current log potential doesn't use it directly).
    @test isapprox(m_custom.consts.b4, 123.456; rtol=0, atol=0)

    φ = SVector{3, Float64}(0.010, 0.020, 0.030)

    masses = Models.calculate_mass_vec(m_custom, φ)
    chi = Models.calculate_chiral(m_custom, φ)

    # Manual checks: must match injected-consts formulas.
    c = m_custom.consts
    m_ud0 = c.m_ud0_inv_fm
    m_s0 = c.m_s0_inv_fm
    G = c.G_fm2
    K = c.K_fm5

    expected_masses = SVector{3, Float64}(
        m_ud0 - 4 * G * φ[1] + 2 * K * φ[2] * φ[3],
        m_ud0 - 4 * G * φ[2] + 2 * K * φ[1] * φ[3],
        m_s0 - 4 * G * φ[3] + 2 * K * φ[1] * φ[2],
    )
    expected_chi = 2 * G * sum(φ .^ 2) - 4 * K * prod(φ)

    @test all(isfinite, masses)
    @test isfinite(chi)
    @test isapprox(masses, expected_masses; rtol=0, atol=1e-12)
    @test isapprox(chi, expected_chi; rtol=0, atol=1e-12)

    # Should differ from default profile if injection is actually used.
    masses_default = Models.calculate_mass_vec(m_default, φ)
    chi_default = Models.calculate_chiral(m_default, φ)
    @test masses != masses_default
    @test chi != chi_default

    # Polyakov potential must reflect injected params.
    Φ = 0.12
    Φbar = 0.34
    T_fm = 0.20

    U_custom = Models.polyakov_potential(m_custom, Φ, Φbar, T_fm)
    U_manual = _polyakov_U_manual(m_custom.consts, Φ, Φbar, T_fm)
    U_default = Models.polyakov_potential(m_default, Φ, Φbar, T_fm)

    @test isfinite(U_custom)
    @test isfinite(U_default)
    @test isapprox(U_custom, U_manual; rtol=0, atol=1e-12)
    @test U_custom != U_default
end
