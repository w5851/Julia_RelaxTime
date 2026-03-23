# Models factory.jl 单元测试
#
# 测试内容：
# 1. create_model(::Symbol) 各模型种类
# 2. create_model(::Type) 泛型路径
# 3. 未知模型抛异常
# 4. 仅支持现有模型种类

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# 确保 PNJL 模块已加载
Models.pnjl_module()

# ============================================================================
# create_model(::Symbol) 各种类
# ============================================================================

@testset "factory create_model" begin

    @testset ":NJL" begin
        m = Models.create_model(:NJL)
        @test m isa Models.AbstractNJLModel
    end

    @testset ":NJL2" begin
        m = Models.create_model(:NJL2)
        @test m isa Models.AbstractNJLModel
    end

    @testset ":PNJL" begin
        m = Models.create_model(:PNJL)
        @test m isa Models.AbstractPNJLModel
    end

    @testset ":PNJLMagnetic" begin
        m = Models.create_model(:PNJLMagnetic; eB_fm2=0.0)
        @test m isa Models.AbstractPNJLModel
    end

    @testset ":RPNJL" begin
        m = Models.create_model(:RPNJL)
        @test m isa Models.AbstractPNJLModel
    end

    @testset ":Rotation" begin
        m = Models.create_model(:Rotation)
        @test m isa Models.AbstractQCDModel
    end

    # ============================================================================
    # 未知模型
    # ============================================================================

    @testset "未知 Symbol 抛异常" begin
        @test_throws ErrorException Models.create_model(:FooBar)
    end

    # ============================================================================
    # create_model(::Type)
    # ============================================================================

    @testset "create_model(Type)" begin
        m = Models.create_model(Models.NJLModel)
        @test m isa Models.NJLModel
    end

end
