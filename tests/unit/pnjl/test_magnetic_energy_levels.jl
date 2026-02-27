using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.legacy_pnjl_module()
const PNJL = Models.legacy_pnjl_module()
const alpha_n = getproperty(PNJL, :alpha_n)
const energy_landau = getproperty(PNJL, :energy_landau)
const resolve_nmax_from_cutoff = getproperty(PNJL, :resolve_nmax_from_cutoff)
const smooth_cutoff = getproperty(PNJL, :smooth_cutoff)

@testset "magnetic energy levels" begin
    @test alpha_n(0) == 1.0
    @test alpha_n(1) == 2.0
    @test alpha_n(5) == 2.0

    m = 1.2
    pz = 0.8
    q = 2 / 3
    eB = 0.15

    E0 = energy_landau(m, pz, 0, q, eB)
    E1 = energy_landau(m, pz, 1, q, eB)
    @test isfinite(E0)
    @test isfinite(E1)
    @test E1 > E0

    nmax = resolve_nmax_from_cutoff(m, 1.0, q, eB)
    @test nmax >= 0

    val = smooth_cutoff(1.0; Λ=3.0, N=10)
    @test 0.0 < val <= 1.0
end
