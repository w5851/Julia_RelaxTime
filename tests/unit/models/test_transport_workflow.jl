# TransportWorkflow.jl 单元测试
#
# 测试内容：
# 1. 模块加载
# 2. solve_gap_and_transport 接口存在
# 3. TransportIntegrationConfig 可用
# 4. build_equilibrium_params 接口存在

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

const _TW = Models.TransportWorkflow

# ============================================================================

@testset "TransportWorkflow" begin

    @testset "solve_gap_and_transport 接口存在" begin
        @test isdefined(_TW, :solve_gap_and_transport)
        @test _TW.solve_gap_and_transport isa Function
    end

    @testset "build_equilibrium_params 接口存在" begin
        @test isdefined(_TW, :build_equilibrium_params)
        @test _TW.build_equilibrium_params isa Function
    end

    @testset "solve_transport_from_equilibrium 接口存在" begin
        @test isdefined(_TW, :solve_transport_from_equilibrium)
        @test _TW.solve_transport_from_equilibrium isa Function
    end

    @testset "reset_transport_workflow_config_cache! 接口存在" begin
        @test isdefined(_TW, :reset_transport_workflow_config_cache!)
        _TW.reset_transport_workflow_config_cache!()
        @test true  # 确保不抛异常
    end
end
