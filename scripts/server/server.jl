#!/usr/bin/env julia

"""
散射计算HTTP服务器启动脚本

使用方法:
    julia --project=. scripts/server/server.jl [port]

默认端口: 8080

示例:
    julia --project=. scripts/server/server.jl          # 使用默认端口8080
    julia --project=. scripts/server/server.jl 8081     # 使用端口8081
"""

# 激活项目环境
using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(REPO_ROOT)

# 加载HTTP服务器模块（使用绝对路径）
include(joinpath(REPO_ROOT, "src", "simulation", "HTTPServer.jl"))
using .HTTPServer

# 解析命令行参数
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

# 启动服务器
println("\n🎯 正在启动散射计算服务器...")

try
    server = start_server(port=port, verbose=true)
    
    # 保持服务器运行
    println("\n✨ 服务器运行中... (按 Ctrl+C 停止)\n")
    
    # 等待中断信号
    try
        wait(server)
    catch e
        if e isa InterruptException
            println("\n\n👋 收到停止信号，正在关闭服务器...")
        else
            rethrow(e)
        end
    end
    
finally
    println("✅ 服务器已停止")
end
