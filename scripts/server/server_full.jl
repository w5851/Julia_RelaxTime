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

# 激活项目环境
using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(REPO_ROOT)

using HTTP
using JSON3

# 加载HTTP服务器模块（使用绝对路径）
include(joinpath(REPO_ROOT, "src", "simulation", "MomentumMapping.jl"))
using .MomentumMapping
using LinearAlgebra

# ==================== API处理函数 ====================

"""
处理 /compute 端点
"""
function handle_compute(req::HTTP.Request)
    try
        # 解析请求体
        body = JSON3.read(String(req.body))
        
        # 提取参数
        p1 = [Float64(body.p1x), Float64(body.p1y), Float64(body.p1z)]
        p2 = [Float64(body.p2x), Float64(body.p2y), Float64(body.p2z)]
        m1 = Float64(body.m1)
        m2 = Float64(body.m2)
        m3 = Float64(body.m3)
        m4 = Float64(body.m4)
        
        theta_star = haskey(body, :theta_star) ? Float64(body.theta_star) : π/4
        phi_star = haskey(body, :phi_star) ? Float64(body.phi_star) : π/6
        
        # 输入验证
        if any(isnan.([p1; p2; m1; m2; m3; m4; theta_star; phi_star]))
            return HTTP.Response(400, ["Content-Type" => "application/json"], 
                JSON3.write(Dict("success" => false, "error" => "Invalid input: NaN detected")))
        end
        
        # 计算散射运动学
        result = calculate_outgoing_momenta(p1, p2, m1, m2, m3, m4, theta_star, phi_star)
        
        # 验证物理约束
        is_valid, checks = validate_kinematics(result, m1, m2, m3, m4, tol=1e-9)
        
        # 构造响应数据
        response_data = Dict(
            "success" => true,
            "data" => Dict(
                "ellipsoid" => Dict(
                    "center" => result.ellipsoid.center,
                    "axes_directions" => [result.ellipsoid.axes_directions[:, i] for i in 1:3],
                    "half_lengths" => result.ellipsoid.half_lengths
                ),
                "momenta" => Dict(
                    "p1" => result.p1_lab,
                    "p2" => result.p2_lab,
                    "p3" => result.p3_lab,
                    "p4" => result.p4_lab,
                    "E1" => result.E1,
                    "E2" => result.E2,
                    "E3" => result.E3,
                    "E4" => result.E4
                ),
                "physics" => Dict(
                    "s" => result.s,
                    "sqrt_s" => sqrt(result.s),
                    "p_star" => result.p_star,
                    "beta" => norm(result.beta),
                    "beta_vector" => result.beta,
                    "gamma" => result.gamma,
                    "theta_star" => result.theta_star,
                    "phi_star" => result.phi_star
                ),
                "validation" => Dict(
                    "is_valid" => is_valid,
                    "energy_conservation" => checks["energy_conservation"][1],
                    "momentum_conservation" => checks["momentum_conservation"][1]
                )
            ),
            "error" => nothing
        )
        
        headers = [
            "Content-Type" => "application/json; charset=utf-8",
            "Access-Control-Allow-Origin" => "*"
        ]
        return HTTP.Response(200, headers, JSON3.write(response_data))
        
    catch e
        error_msg = sprint(showerror, e, catch_backtrace())
        @error "Computation error" exception=(e, catch_backtrace())
        
        response_data = Dict("success" => false, "data" => nothing, "error" => error_msg)
        headers = [
            "Content-Type" => "application/json; charset=utf-8",
            "Access-Control-Allow-Origin" => "*"
        ]
        return HTTP.Response(400, headers, JSON3.write(response_data))
    end
end

"""
CORS中间件
"""
function cors_middleware(handler)
    return function(req::HTTP.Request)
        if req.method == "OPTIONS"
            return HTTP.Response(200, [
                "Access-Control-Allow-Origin" => "*",
                "Access-Control-Allow-Methods" => "POST, GET, OPTIONS",
                "Access-Control-Allow-Headers" => "Content-Type",
                "Access-Control-Max-Age" => "86400"
            ])
        end
        
        resp = handler(req)
        HTTP.setheader(resp, "Access-Control-Allow-Origin" => "*")
        HTTP.setheader(resp, "Access-Control-Allow-Methods" => "POST, GET, OPTIONS")
        HTTP.setheader(resp, "Access-Control-Allow-Headers" => "Content-Type")
        
        return resp
    end
end

"""
静态文件服务
"""
function serve_static_file(path::String)
    # 安全检查：防止目录遍历
    if contains(path, "..") || contains(path, "\\")
        return HTTP.Response(403, "Forbidden")
    end
    
    # 确定文件路径
    # 移除路径开头的 /
    clean_path = startswith(path, "/") ? path[2:end] : path
    
    file_path = if path == "/" || path == ""
        joinpath(REPO_ROOT, "web", "index.html")
    else
        joinpath(REPO_ROOT, "web", clean_path)
    end
    
    # 检查文件是否存在
    if !isfile(file_path)
        return HTTP.Response(404, "Not Found: $path")
    end
    
    # 确定Content-Type
    content_type = if endswith(file_path, ".html")
        "text/html; charset=utf-8"
    elseif endswith(file_path, ".css")
        "text/css; charset=utf-8"
    elseif endswith(file_path, ".js")
        "application/javascript; charset=utf-8"
    elseif endswith(file_path, ".json")
        "application/json; charset=utf-8"
    else
        "application/octet-stream"
    end
    
    # 读取并返回文件
    try
        content = read(file_path)
        return HTTP.Response(200, ["Content-Type" => content_type], body=content)
    catch e
        @error "Error serving file" file=file_path exception=e
        return HTTP.Response(500, "Internal Server Error")
    end
end

"""
请求路由
"""
function route_request(req::HTTP.Request)
    path = String(HTTP.URIs.unescapeuri(req.target))
    
    # API端点
    if path == "/health"
        return HTTP.Response(200, ["Content-Type" => "text/plain"], "OK")
    elseif path == "/compute"
        return handle_compute(req)
    # 静态文件
    else
        # 移除查询参数并确保是String
        path = String(split(path, '?')[1])
        return serve_static_file(path)
    end
end

# ==================== 启动服务器 ====================

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

# 应用CORS中间件
app = cors_middleware(route_request)

# 打印启动信息
println("\n" * "="^60)
println("🚀 散射计算服务器启动中...")
println("="^60)
println("📍 服务地址: http://localhost:$port")
println("📡 API端点:")
println("   POST http://localhost:$port/compute")
println("   GET  http://localhost:$port/health")
println("\n📁 静态文件:")
println("   http://localhost:$port/")
println("   http://localhost:$port/index.html")
println("   http://localhost:$port/simple_test.html")
println("\n💡 提示:")
println("   • 在浏览器中打开: http://localhost:$port")
println("   • 按 Ctrl+C 停止服务器")
println("="^60 * "\n")

# 启动服务器
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
