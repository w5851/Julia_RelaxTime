# CEPDetector.jl 单元测试
#
# 测试内容：
# 1. find_cep 接口存在
# 2. CEPResult 结构
# 3. 内部辅助函数可访问

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================

@testset "CEPDetector" begin

    @testset "find_cep 接口存在" begin
        @test isdefined(Models, :find_cep)
        @test Models.find_cep isa Function
    end

    @testset "CEPResult 结构" begin
        @test isdefined(Models, :CEPResult)
        res = Models.CEPResult(
            found=false,
            T_cep_MeV=NaN,
            mu_cep_MeV=NaN,
            uncertainty_T_MeV=NaN,
            eval_count=0,
            unknown_count=0,
            reason="test",
            method=:bisection,
        )
        @test res.found == false
        @test res.reason == "test"
    end

    @testset "find_cep 空 curves dict → 未找到" begin
        curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
        result = Models.find_cep(curves)
        @test result isa Models.CEPResult
        @test result.found == false
    end
end
