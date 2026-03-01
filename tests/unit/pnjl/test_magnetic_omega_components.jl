using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

const PNJL = Models.pnjl_module()
const cached_nodes = getproperty(PNJL, :cached_nodes)
const calculate_omega = getproperty(PNJL, :calculate_omega)
const MagneticConfig = getproperty(PNJL, :MagneticConfig)
const calculate_magnetic_omega_components = getproperty(PNJL, :calculate_magnetic_omega_components)
const calculate_magnetic_omega = getproperty(PNJL, :calculate_magnetic_omega)

@testset "magnetic omega components" begin
    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)
    mu_vec = SVector{3, Float64}(0.4, 0.4, 0.4)
    T_fm = 0.8

    conf = MagneticConfig(eB_fm2=0.12, p_num=24, pz_max=10.0)
    comp = calculate_magnetic_omega_components(x_state, mu_vec, T_fm, conf)

    @test isfinite(comp.omega)
    @test isfinite(comp.vac)
    @test isfinite(comp.therm)
    @test comp.n_max >= 0
    @test comp.G_B > 0

    thermal_nodes = cached_nodes(24, 8)
    omega_legacy = calculate_omega(x_state, mu_vec, T_fm, thermal_nodes, 0.0)
    omega_b0 = calculate_magnetic_omega(x_state, mu_vec, T_fm, MagneticConfig(eB_fm2=0.0, p_num=24, pz_max=10.0))
    @test isapprox(omega_b0, omega_legacy; rtol=1e-7, atol=1e-8)
end
