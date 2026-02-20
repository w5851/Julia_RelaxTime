using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "core", "MagneticIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "core", "MagneticThermodynamics.jl"))

@testset "magnetic densities" begin
    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)
    mu_vec = SVector{3, Float64}(0.5, 0.5, 0.5)
    T_fm = 0.7
    conf = Main.MagneticThermodynamics.MagneticConfig(eB_fm2=0.08, p_num=24, pz_max=10.0)

    rho = Main.MagneticThermodynamics.calculate_magnetic_rho(x_state, mu_vec, T_fm, conf)
    @test all(isfinite.(rho))

    nd = Main.MagneticThermodynamics.calculate_magnetic_number_densities(x_state, mu_vec, T_fm, conf)
    @test all(isfinite.(nd.quark))
    @test isfinite(nd.baryon)
end
