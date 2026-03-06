# PNJL Solver 单元测试（使用新架构）

using Test
using StaticArrays
using Base.MathConstants: π

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL_MODEL = Models.create_model(:PNJL)
const PNJL = Models.pnjl_module()
const vacuum_integral = getproperty(PNJL, :vacuum_integral)
const calculate_energy_sum = getproperty(PNJL, :calculate_energy_sum)
const calculate_mass_vec = getproperty(PNJL, :calculate_mass_vec)
const Integrals = getproperty(PNJL, :Integrals)
const DEFAULT_MOMENTUM_NODES = getproperty(Integrals, :DEFAULT_MOMENTUM_NODES)
const DEFAULT_MOMENTUM_WEIGHTS = getproperty(Integrals, :DEFAULT_MOMENTUM_WEIGHTS)

const Λ = Constants_PNJL.Λ_inv_fm
const GAUSS_NODES = 256

# 使用 GaussLegendre 生成高精度节点
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
const GL_NODES, GL_WEIGHTS = GaussLegendre.gauleg(0.0, Λ, GAUSS_NODES)

NumericVacuum(m) = begin
    integrand = @. GL_NODES^2 * sqrt(GL_NODES^2 + m^2)
    sum(GL_WEIGHTS .* integrand) / (2π^2)
end

@testset "Vacuum integral matches numeric quadrature" begin
    for mass in (Constants_PNJL.m_ud0_inv_fm, Constants_PNJL.m_s0_inv_fm, 0.0, 0.25)
        analytic = vacuum_integral(mass)
        numeric = NumericVacuum(mass)
        @test isapprox(analytic, numeric; atol=1e-8, rtol=1e-6)
    end
end

@testset "calculate_energy_sum aggregates vacuum_integral" begin
    phi = @SVector [-0.02, -0.02, -0.25]
    masses = Models.calculate_mass_vec(PNJL_MODEL, phi)
    direct = calculate_energy_sum(masses)
    manual = -2 * Constants_PNJL.N_color * sum(vacuum_integral(masses[i]) for i in 1:3)
    @test isapprox(direct, manual; atol=1e-10, rtol=1e-10)
end
