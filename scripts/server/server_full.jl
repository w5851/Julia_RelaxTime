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

include(joinpath(REPO_ROOT, "src", "simulation", "ServerLauncher.jl"))
using .ServerLauncher

ServerLauncher.run_full_server(REPO_ROOT, ARGS)
