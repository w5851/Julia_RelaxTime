#!/usr/bin/env julia

"""
完整HTTP服务器 - 同时提供API和静态文件服务

使用方法:
    julia scripts/server/server_full.jl [port]

默认端口: 8080

示例:
    julia scripts/server/server_full.jl          # 使用默认端口8080
    julia scripts/server/server_full.jl 8081     # 使用端口8081
"""

using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(REPO_ROOT)

using HTTP

include(joinpath(REPO_ROOT, "src", "simulation", "FullServerApp.jl"))
using .FullServerApp

const DEFAULT_PORT = 8080
port = DEFAULT_PORT
if length(ARGS) >= 1
    try
        global port = parse(Int, ARGS[1])
        if port < 1024 || port > 65535
            @warn "端口号应在1024-65535之间，使用默认端口8080"
            global port = DEFAULT_PORT
        end
    catch
        @warn "无效的端口号，使用默认端口8080"
    end
end

app = FullServerApp.build_app(REPO_ROOT)

println("\n" * "="^60)
println("🚀 散射计算服务器启动中...")
println("="^60)
println("📍 服务地址: http://localhost:$port")
println("📡 API端点:")
println("   POST http://localhost:$port/compute")
println("   GET  http://localhost:$port/health")
println("   GET  http://localhost:$port/api/modules")
println("   POST http://localhost:$port/api/modules/pnjl-gap/run")
println("\n📁 静态文件:")
println("   http://localhost:$port/")
println("   http://localhost:$port/index.html")
println("   http://localhost:$port/simple_test.html")
println("\n💡 提示:")
println("   • 在浏览器中打开: http://localhost:$port")
println("   • 按 Ctrl+C 停止服务器")
println("="^60 * "\n")

try
    HTTP.serve(app, "0.0.0.0", port; verbose=false)
catch e
    if e isa Base.IOError || e isa ArgumentError
        @error "无法启动服务器: 端口 $port 可能已被占用" exception=e
        println("\n💡 尝试使用其他端口:")
        println("   julia scripts/server/server_full.jl 8081")
        rethrow(e)
    else
        rethrow(e)
    end
end
