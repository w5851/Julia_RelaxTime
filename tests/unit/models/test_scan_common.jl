# ScanCommon 模块单元测试

using Test

const PROJECT_ROOT_SC = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_SC, "src", "models", "Models.jl"))
end

# ScanCommon 使用 ..SeedStrategies 等相对导入，只能通过 Models 子模块访问
using Main.Models.ScanCommon: fmt, clean_message, quote_csv, join_messages, key3, key2, attempt_with_candidates

@testset "ScanCommon" begin
    @testset "fmt 格式化数值" begin
        @test fmt(1.23456789) isa String
        @test contains(fmt(1.23456789), "1.234568")  # %.6f 格式
    end

    @testset "clean_message 清理消息字符串" begin
        msg = clean_message("hello\nworld\r\ntest")
        @test !contains(msg, "\n")
    end

    @testset "quote_csv CSV 引用" begin
        q = quote_csv("hello,world")
        @test contains(q, "\"")
    end

    @testset "key3 返回三元组" begin
        k = key3(0.1234, 0.5678, 0.9012; digits=4)
        @test k isa NTuple{3, Float64}
    end

    @testset "key2 返回二元组" begin
        k = key2(0.1234, 0.5678; digits=4)
        @test k isa NTuple{2, Float64}
    end

    @testset "attempt_with_candidates 错误消息包含统一错误语义" begin
        candidates = [(label="seed-a", state=[1.0, 2.0, 3.0])]
        result, message = attempt_with_candidates(candidates;
            solve_point=_ -> throw(ArgumentError("bad_seed")),
            refine=r -> (r, ""),
            promote=r -> (r, ""),
            is_success=_ -> false,
            stop_on_first_success=true,
            evaluate_all_attempts=false,
        )
        @test result === nothing
        @test occursin("error.kind=constraint_error", message)
        @test occursin("bad_seed", message)
    end
end
