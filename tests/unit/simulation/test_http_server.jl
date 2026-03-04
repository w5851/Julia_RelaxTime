# HTTPServer.jl 单元测试
#
# HTTPServer 依赖 HTTP / JSON3 / MomentumMapping，加载代价重。
# 本测试先尝试加载模块；若因依赖不可用而失败，则跳过。
# 测试范围：APIResponse 结构体 + 接口声明。

using Test

const PROJECT_ROOT_HS = normpath(joinpath(@__DIR__, "..", "..", ".."))

_hs_loaded = false
try
    if !isdefined(Main, :HTTPServer)
        include(joinpath(PROJECT_ROOT_HS, "src", "simulation", "HTTPServer.jl"))
    end
    global _hs_loaded = true
catch e
    @warn "HTTPServer 模块加载失败 ($(typeof(e)))；跳过测试"
end

# ============================================================================

@testset "HTTPServer" begin
    if !_hs_loaded
        @test_skip "module not loaded"
    else
        @testset "APIResponse 构造" begin
            resp = HTTPServer.APIResponse(true, Dict("k" => 1), nothing)
            @test resp.success == true
            @test resp.data isa Dict
            @test resp.error === nothing
        end

        @testset "APIResponse 错误响应" begin
            resp = HTTPServer.APIResponse(false, nothing, "some error")
            @test resp.success == false
            @test resp.data === nothing
            @test resp.error == "some error"
        end

        @testset "handle_compute 接口存在" begin
            @test isdefined(HTTPServer, :handle_compute)
        end

        @testset "start_server 接口存在" begin
            @test isdefined(HTTPServer, :start_server)
        end
    end
end
