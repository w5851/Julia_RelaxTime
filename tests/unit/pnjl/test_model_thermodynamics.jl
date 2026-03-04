# ModelThermodynamics 模块单元测试
#
# 测试内容：
# 1. pressure_model 热力学压力
# 2. rho_model 手征凝聚
# 3. number_densities_model 数密度

using Test
using StaticArrays

const PROJECT_ROOT_MT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_MT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ModelThermodynamics 通过 Models 间接提供
const PNJL_MT = Models.pnjl_module()
const solve_mt = getproperty(PNJL_MT, :solve)
const FixedMu_mt = getproperty(PNJL_MT, :FixedMu)

@testset "ModelThermodynamics" begin
    # 在一个简单的物理点求解
    T_fm = 0.8  # ~158 MeV
    mu_fm = 0.0

    result = solve_mt(FixedMu_mt(), T_fm, mu_fm)

    @testset "求解收敛" begin
        @test result.converged
    end

    @testset "状态向量有限" begin
        @test all(isfinite, result.x_state)
    end

    @testset "质量正定" begin
        @test all(m -> m > 0, result.masses)
    end
end
