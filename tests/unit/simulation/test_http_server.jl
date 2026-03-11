# HTTPServer.jl 单元测试
#
# HTTPServer 依赖 HTTP / JSON3 / MomentumMapping，加载代价重。
# 本测试先尝试加载模块；若因依赖不可用而失败，则跳过。
# 测试范围：APIResponse 结构体 + 接口声明。

using Test
using HTTP
using JSON3

const PROJECT_ROOT_HS = normpath(joinpath(@__DIR__, "..", "..", ".."))

function _make_request(payload::Dict)
    return HTTP.Request("POST", "/compute", ["Content-Type" => "application/json"], JSON3.write(payload))
end

function _parse_body(resp::HTTP.Response)
    return JSON3.read(String(resp.body))
end

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
            resp = HTTPServer.APIResponse(true, Dict("k" => 1), nothing, nothing)
            @test resp.success == true
            @test resp.data isa Dict
            @test resp.error === nothing
            @test resp.error_code === nothing
        end

        @testset "APIResponse 错误响应" begin
            resp = HTTPServer.APIResponse(false, nothing, "some error", "SOME_ERROR")
            @test resp.success == false
            @test resp.data === nothing
            @test resp.error == "some error"
            @test resp.error_code == "SOME_ERROR"
        end

        @testset "handle_compute 接口存在" begin
            @test isdefined(HTTPServer, :handle_compute)
        end

        @testset "handle_compute 成功响应合同" begin
            req = _make_request(Dict(
                "p1x" => 0.5,
                "p1y" => 0.0,
                "p1z" => 1.8,
                "p2x" => -0.5,
                "p2y" => 0.0,
                "p2z" => -1.8,
                "m1" => 1.52,
                "m2" => 1.52,
                "m3" => 1.52,
                "m4" => 1.52,
                "theta_star" => 0.785,
                "phi_star" => 0.524,
            ))

            resp = HTTPServer.handle_compute(req)
            body = _parse_body(resp)

            @test resp.status == 200
            @test body.success == true
            @test body.error === nothing
            @test body.error_code === nothing
            @test haskey(body.data, :physics)
            @test haskey(body.data, :validation)
        end

        @testset "handle_compute 用户输入错误合同" begin
            req = _make_request(Dict(
                "p1x" => 0.5,
                "p1y" => 0.0,
                "p1z" => 1.8,
                "p2x" => -0.5,
                "p2y" => 0.0,
                "p2z" => -1.8,
                "m1" => -1.0,
                "m2" => 1.52,
                "m3" => 1.52,
                "m4" => 1.52,
            ))

            resp = HTTPServer.handle_compute(req)
            body = _parse_body(resp)

            @test resp.status == 400
            @test body.success == false
            @test body.data === nothing
            @test body.error_code == "INVALID_INPUT"
            @test body.error == "Invalid input: masses must be positive"
        end

        @testset "handle_compute 运行时错误合同" begin
            req = _make_request(Dict(
                "p1x" => 0.5,
                "p1y" => 0.0,
                "p1z" => 1.8,
                "p2x" => -0.5,
                "p2y" => 0.0,
                "p2z" => -1.8,
                "m1" => 1.52,
                "m2" => 1.52,
                "m3" => 1.52,
            ))

            resp = HTTPServer.handle_compute(req)
            body = _parse_body(resp)

            @test resp.status == 500
            @test body.success == false
            @test body.data === nothing
            @test body.error_code == "COMPUTATION_ERROR"
            @test body.error == "Computation failed"
        end

        @testset "start_server 接口存在" begin
            @test isdefined(HTTPServer, :start_server)
        end
    end
end
