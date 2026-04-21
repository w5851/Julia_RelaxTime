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
using Main.EquilibriumFacade: solve_equilibrium_backend

@testset "EquilibriumFacade" begin
    @testset "pnjl_model_kind 返回符号" begin
        kind = pnjl_model_kind(:legacy)
        @test kind isa Symbol
        @test kind == :PNJL
    end

    @testset "models backend 默认走 FixedMu 治理链" begin
        res = solve_equilibrium_backend(
            0.15,
            0.0;
            solver_backend=:models,
            p_num=8,
            t_num=4,
            seed_state=Models.HADRON_SEED_5,
            models_residual_norm_max=1e-4,
        )

        @test res.converged === true
        @test !ismissing(res.iterations)
        @test !ismissing(res.residual_norm)
        @test isfinite(Float64(res.residual_norm))
    end

    @testset "models backend 显式 models_solver 仍走 solve_gap 路径" begin
        custom_solver = Models.NLsolveGapSolver(method=:trust_region, jacobian=:forward)
        res = solve_equilibrium_backend(
            0.15,
            0.0;
            solver_backend=:models,
            models_solver=custom_solver,
            p_num=8,
            t_num=4,
            seed_state=Models.HADRON_SEED_5,
            models_residual_norm_max=1e-4,
        )

        @test res.converged === true
        @test ismissing(res.iterations)
        @test ismissing(res.residual_norm)
    end
end
