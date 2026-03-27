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

function _contains_mutable_container(x)
    if x isa Dict || x isa Vector
        return true
    elseif x isa NamedTuple
        for v in values(x)
            _contains_mutable_container(v) && return true
        end
    end
    return false
end

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

    @testset "WorkflowCache 注入与定向 reset" begin
        @test isdefined(_TW, :WorkflowCache)

        cache = _TW.WorkflowCache()
        profile = get(ENV, "PHYSICS_PARAM_PROFILE", "default")

        @test isempty(cache.model_cache)
        @test isempty(cache.prefer_energy_aniso_cache)
        @test isempty(cache.a_builder_config_cache)

        model = _TW._get_model(cache, :PNJL)
        prefer = _TW._default_prefer_energy_aniso_from_toml(cache)
        a_cfg = _TW._default_a_builder_config_from_toml(cache)

        @test model !== nothing
        @test haskey(cache.model_cache, :PNJL)
        @test haskey(cache.prefer_energy_aniso_cache, profile)
        @test haskey(cache.a_builder_config_cache, profile)
        @test prefer isa Bool
        @test a_cfg.p_nodes > 0

        _TW.reset_transport_workflow_config_cache!(cache)

        @test haskey(cache.model_cache, :PNJL)
        @test isempty(cache.prefer_energy_aniso_cache)
        @test isempty(cache.a_builder_config_cache)
    end

    @testset "fallback schema and local dict isolation" begin
        @test isdefined(_TW, :DEFAULT_PHYSICS_FALLBACK_NT)
        @test isdefined(_TW, :DEFAULT_A_BUILDER_FALLBACK)
        @test !_contains_mutable_container(getproperty(_TW, :DEFAULT_PHYSICS_FALLBACK_NT))
        @test !_contains_mutable_container(getproperty(_TW, :DEFAULT_A_BUILDER_FALLBACK))

        @test isdefined(_TW, :_default_physics_fallback_dict)
        d1 = _TW._default_physics_fallback_dict()
        d2 = _TW._default_physics_fallback_dict()
        d1["transport_workflow"]["prefer_energy_aniso"] = false
        @test d2["transport_workflow"]["prefer_energy_aniso"] == true
    end
end
