# EquilibriumFacade 模块单元测试
#
# 测试内容：
# 1. pnjl_model_kind 返回正确模型类型
# 2. solve_equilibrium_backend 基本收敛

using Test
using StaticArrays

const PROJECT_ROOT_EF = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_EF, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# EquilibriumFacade 位于 pnjl_physics/core/ 下，通过 Models 间接加载
const _EQ_FACADE_PATH = normpath(joinpath(PROJECT_ROOT_EF, "src", "models", "pnjl_physics", "core", "EquilibriumFacade.jl"))
if !isdefined(Main, :EquilibriumFacade)
    Base.include(Main, _EQ_FACADE_PATH)
end

using Main.EquilibriumFacade: pnjl_model_kind

@testset "EquilibriumFacade" begin
    @testset "pnjl_model_kind 返回符号" begin
        kind = pnjl_model_kind(:legacy)
        @test kind isa Symbol
        @test kind == :PNJL
    end
end
