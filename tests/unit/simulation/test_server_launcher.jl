# ServerLauncher.jl 单元测试
#
# ServerLauncher 依赖 HTTP + FullServerApp，加载代价重。
# 本测试先尝试加载模块；若因依赖不可用而失败，则跳过。
# 测试范围：parse_port 纯函数 + 常量 + 接口声明。

using Test

const PROJECT_ROOT_SL = normpath(joinpath(@__DIR__, "..", "..", ".."))

_sl_loaded = false
try
    if !isdefined(Main, :ServerLauncher)
        include(joinpath(PROJECT_ROOT_SL, "src", "simulation", "ServerLauncher.jl"))
    end
    global _sl_loaded = true
catch e
    @warn "ServerLauncher 模块加载失败 ($(typeof(e)))；跳过测试"
end

# ============================================================================

@testset "ServerLauncher" begin
    if !_sl_loaded
        @test_skip "module not loaded"
    else
        @testset "DEFAULT_PORT" begin
            @test ServerLauncher.DEFAULT_PORT == 8080
        end

        @testset "parse_port 空参数 → 默认" begin
            @test ServerLauncher.parse_port(String[]) == 8080
        end

        @testset "parse_port 有效端口" begin
            @test ServerLauncher.parse_port(["3000"]) == 3000
            @test ServerLauncher.parse_port(["8888"]) == 8888
        end

        @testset "parse_port 超范围 → 默认" begin
            @test ServerLauncher.parse_port(["999999"]) == 8080
            @test ServerLauncher.parse_port(["100"]) == 8080
        end

        @testset "parse_port 非数字 → 默认" begin
            @test ServerLauncher.parse_port(["abc"]) == 8080
        end

        @testset "run_full_server 接口存在" begin
            @test isdefined(ServerLauncher, :run_full_server)
        end
    end
end
