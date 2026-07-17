# CrossoverLine.jl 单元测试
#
# 测试内容：
# 1. CrossoverResult 结构
# 2. detect_crossover / scan_crossover_line / build_crossover_line 接口存在
# 3. 基本可调用性

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================

@testset "CrossoverLine" begin

    @testset "CrossoverResult 结构" begin
        res = Models.CrossoverResult(method=:chiral_susceptibility)
        @test res.found == false
        @test res.method === :chiral_susceptibility
        @test res.T_crossover === nothing
    end

    @testset "detect_crossover 接口存在" begin
        @test isdefined(Models, :detect_crossover)
        @test Models.detect_crossover isa Function
    end

    @testset "scan_crossover_line 接口存在" begin
        @test isdefined(Models, :scan_crossover_line)
        @test Models.scan_crossover_line isa Function
    end

    @testset "build_crossover_line 接口存在" begin
        @test isdefined(Models, :build_crossover_line)
        @test Models.build_crossover_line isa Function
    end

    @testset "legacy solver backend is rejected" begin
        @test_throws ArgumentError Models.detect_crossover(
            0.0,
            (0.45, 0.55);
            solver_backend=:legacy,
            n_scan=3,
            max_iter=1,
            p_num=8,
            t_num=4,
        )
    end

    @testset "热积分策略在求解前校验" begin
        @test_throws ArgumentError Models.detect_crossover(
            0.0,
            (0.45, 0.55);
            thermo_quadrature_policy=:unknown,
        )
    end
end
