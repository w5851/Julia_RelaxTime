using Test
using StaticArrays

include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "njl", "NJL2Core.jl"))

@testset "NJL2Core config initialization" begin
    p = NJL2Core.njl2_params(profile="default")

    @test p.N_color == 3
    @test p.N_flavor == 2

    hbarc = p.hbarc_MeV_fm
    @test isfinite(hbarc)

    Λ_MeV = 631.0
    GΛ2 = 2.14
    m_ud0_MeV = 5.5

    @test isapprox(p.Lambda_inv_fm, Λ_MeV / hbarc; rtol=1e-12)
    @test isapprox(p.m_ud0_inv_fm, m_ud0_MeV / hbarc; rtol=1e-12)

    G_MeV_inv2 = GΛ2 / (Λ_MeV^2)
    @test isapprox(p.G_fm2, G_MeV_inv2 * hbarc^2; rtol=1e-12)
end

@testset "NJL2Core mass equation and chiral term" begin
    p = NJL2Core.njl2_params(profile="default")

    φ = SVector{3, Float64}(0.10, 0.20, 0.30)
    masses = NJL2Core.calculate_mass_vec(p, φ)

    G = p.G_fm2
    m0 = p.m_ud0_inv_fm

    @test masses[1] == m0 - 4 * G * φ[1]
    @test masses[2] == m0 - 4 * G * φ[2]
    @test masses[3] == 0.0

    χ = NJL2Core.chiral_potential(p, φ)
    χ_expected = 2 * G * (φ[1]^2 + φ[2]^2)
    @test χ == χ_expected
end
