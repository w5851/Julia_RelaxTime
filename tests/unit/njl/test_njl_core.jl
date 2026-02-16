using Test
using StaticArrays

include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "njl", "NJLCore.jl"))

@testset "NJLCore config initialization" begin
    p = NJLCore.njl_params(profile="default")

    @test p.N_color == 3
    @test p.N_flavor == 3

    hbarc = 197.3269804
    @test isapprox(p.hbarc_MeV_fm, hbarc; atol=0.0, rtol=0.0)

    # Standard parameter table (docs/reference/formula/models/njl/NJL_core.md)
    Λ_MeV = 602.3
    GΛ2 = 1.835
    KΛ5 = 12.36
    m_ud0_MeV = 5.5
    m_s0_MeV = 140.7

    @test isapprox(p.Lambda_inv_fm, Λ_MeV / hbarc; rtol=1e-14)
    @test isapprox(p.m_ud0_inv_fm, m_ud0_MeV / hbarc; rtol=1e-14)
    @test isapprox(p.m_s0_inv_fm, m_s0_MeV / hbarc; rtol=1e-14)

    G_MeV_inv2 = GΛ2 / (Λ_MeV^2)
    K_MeV_inv5 = KΛ5 / (Λ_MeV^5)
    @test isapprox(p.G_fm2, G_MeV_inv2 * hbarc^2; rtol=1e-14)
    @test isapprox(p.K_fm5, K_MeV_inv5 * hbarc^5; rtol=1e-14)
end

@testset "NJLCore mass equation and chiral term" begin
    p = NJLCore.njl_params(profile="default")

    φ = SVector{3, Float64}(0.10, 0.20, 0.30)
    masses = NJLCore.calculate_mass_vec(p, φ)

    G = p.G_fm2
    K = p.K_fm5
    m_ud0 = p.m_ud0_inv_fm
    m_s0 = p.m_s0_inv_fm

    expected = SVector(
        m_ud0 - 4 * G * φ[1] + 2 * K * φ[2] * φ[3],
        m_ud0 - 4 * G * φ[2] + 2 * K * φ[1] * φ[3],
        m_s0 - 4 * G * φ[3] + 2 * K * φ[1] * φ[2],
    )
    @test masses == expected

    χ = NJLCore.chiral_potential(p, φ)
    χ_expected = 2 * G * sum(φ .^ 2) - 4 * K * prod(φ)
    @test χ == χ_expected
end
