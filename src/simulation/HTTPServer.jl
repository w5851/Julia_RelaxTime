"""
    HTTPServer

HTTP服务器模块，提供REST API用于前端调用散射计算。

功能：
- POST /compute: 接收输入参数，返回计算结果
- CORS支持
- 错误处理和JSON序列化
"""
module HTTPServer

using HTTP
using JSON3
using LinearAlgebra

# 导入计算模块
include(joinpath(@__DIR__, "MomentumMapping.jl"))
using .MomentumMapping

export start_server, handle_compute

"""
API响应结构体
"""
struct APIResponse
    success::Bool
    data::Union{Dict, Nothing}
    error::Union{String, Nothing}
    error_code::Union{String, Nothing}
end

# 定义JSON3序列化规则
JSON3.StructType(::Type{APIResponse}) = JSON3.Struct()

@inline function _json_response(status::Int, resp::APIResponse)
    headers = ["Content-Type" => "application/json; charset=utf-8"]
    return HTTP.Response(status, headers, JSON3.write(resp))
end

@inline function _error_response(status::Int, code::String, message::String)
    return _json_response(status, APIResponse(false, nothing, message, code))
end

"""
    handle_compute(req::HTTP.Request)

处理 /compute 端点的POST请求。

# 请求格式 (JSON)
```json
{
    "p1x": 0.5, "p1y": 0.0, "p1z": 1.8,
    "p2x": -0.5, "p2y": 0.0, "p2z": -1.8,
    "m1": 1.52, "m2": 1.52, "m3": 1.52, "m4": 1.52,
    "theta_star": 0.785,  // 可选，默认π/4
    "phi_star": 0.524     // 可选，默认π/6
}
```

# 响应格式 (JSON)
成功时:
```json
{
    "success": true,
    "data": {
        "ellipsoid": {...},
        "momenta": {...},
        "physics": {...}
    },
    "error": null
}
```

失败时:
```json
{
    "success": false,
    "data": null,
    "error": "错误信息"
}
```
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
        
        # 可选参数，提供默认值
        theta_star = haskey(body, :theta_star) ? Float64(body.theta_star) : π/4
        phi_star = haskey(body, :phi_star) ? Float64(body.phi_star) : π/6
        
        # 输入验证
        if any(isnan.([p1; p2; m1; m2; m3; m4; theta_star; phi_star]))
            return _error_response(400, "INVALID_INPUT", "Invalid input: NaN detected")
        end
        
        if any([m1, m2, m3, m4] .<= 0)
            return _error_response(400, "INVALID_INPUT", "Invalid input: masses must be positive")
        end
        
        # 计算散射运动学
        result = calculate_outgoing_momenta(p1, p2, m1, m2, m3, m4, 
                                           theta_star, phi_star)
        
        # 验证物理约束
        is_valid, checks = validate_kinematics(result, m1, m2, m3, m4, tol=1e-9)
        
        if !is_valid
            @warn "Physics constraints not satisfied" checks
        end
        
        # 构造响应数据
        response_data = Dict(
            "ellipsoid" => Dict(
                "center" => result.ellipsoid.center,
                "axes_directions" => [result.ellipsoid.axes_directions[:, i] 
                                      for i in 1:3],
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
        )
        
        return _json_response(200, APIResponse(true, response_data, nothing, nothing))
        
    catch e
        @error "Computation error" exception=(e, catch_backtrace())

        if e isa ArgumentError || e isa DomainError
            return _error_response(400, "INVALID_INPUT", "Invalid request parameters")
        end
        return _error_response(500, "COMPUTATION_ERROR", "Computation failed")
    end
end

"""
    cors_middleware(handler)

CORS中间件，允许跨域请求。
"""
function cors_middleware(handler)
    return function(req::HTTP.Request)
        # 处理OPTIONS预检请求
        if req.method == "OPTIONS"
            return HTTP.Response(200, [
                "Access-Control-Allow-Origin" => "*",
                "Access-Control-Allow-Methods" => "POST, GET, OPTIONS",
                "Access-Control-Allow-Headers" => "Content-Type",
                "Access-Control-Max-Age" => "86400"
            ])
        end
        
        # 处理实际请求
        resp = handler(req)
        
        # 添加CORS头
        HTTP.setheader(resp, "Access-Control-Allow-Origin" => "*")
        HTTP.setheader(resp, "Access-Control-Allow-Methods" => "POST, GET, OPTIONS")
        HTTP.setheader(resp, "Access-Control-Allow-Headers" => "Content-Type")
        
        return resp
    end
end

"""
    start_server(; port=8080, host="0.0.0.0", verbose=true)

启动HTTP服务器。

# 参数
- `port::Int`: 服务器端口（默认8080）
- `host::String`: 绑定地址（默认"0.0.0.0"，监听所有接口）
- `verbose::Bool`: 是否打印详细信息（默认true）

# 返回
- `HTTP.Server`: 服务器对象

# 使用
```julia
server = start_server(port=8080)
# 按 Ctrl+C 停止服务器
```
"""
function start_server(; port::Int=8080, host::String="0.0.0.0", verbose::Bool=true)
    # 创建路由
    router = HTTP.Router()
    
    # 注册路由
    HTTP.register!(router, "POST", "/compute", handle_compute)
    
    # 健康检查端点
    HTTP.register!(router, "GET", "/health", req -> HTTP.Response(200, "OK"))
    
    # 应用CORS中间件
    app = cors_middleware(router)
    
    # 打印启动信息
    if verbose
        println("\n" * "="^60)
        println("🚀 散射计算服务器启动中...")
        println("="^60)
        println("📍 服务地址: http://$host:$port")
        println("📡 API端点:")
        println("   POST http://localhost:$port/compute")
        println("   GET  http://localhost:$port/health")
        println("\n💡 提示:")
        println("   • 前端页面: 打开 web/index.html")
        println("   • 按 Ctrl+C 停止服务器")
        println("="^60 * "\n")
    end
    
    # 启动服务器
    try
        server = HTTP.serve!(app, host, port; verbose=verbose)
        return server
    catch e
        if e isa Base.IOError || e isa ArgumentError
            @error "无法启动服务器: 端口 $port 可能已被占用" exception=e
            println("\n💡 尝试使用其他端口:")
            println("   julia> start_server(port=8081)")
            rethrow(e)
        else
            rethrow(e)
        end
    end
end

"""
    test_api_endpoint(; port=8080)

测试API端点是否正常工作。

# 示例
```julia
test_api_endpoint()
```
"""
function test_api_endpoint(; port::Int=8080)
    url = "http://localhost:$port/compute"
    
    # 测试数据
    test_data = Dict(
        "p1x" => 0.5, "p1y" => 0.0, "p1z" => 1.8,
        "p2x" => -0.5, "p2y" => 0.0, "p2z" => -1.8,
        "m1" => 1.52, "m2" => 1.52, "m3" => 1.52, "m4" => 1.52
    )
    
    println("测试API端点: $url")
    println("发送数据: ", test_data)
    
    try
        response = HTTP.post(url, 
                           ["Content-Type" => "application/json"],
                           JSON3.write(test_data))
        
        if response.status == 200
            result = JSON3.read(String(response.body))
            println("\n✅ API测试成功!")
            println("响应: ", result.success ? "成功" : "失败")
            if result.success
                println("√s = ", result.data["physics"]["sqrt_s"], " fm⁻¹")
            end
        else
            println("\n❌ API测试失败")
            println("状态码: ", response.status)
        end
    catch e
        println("\n❌ 连接失败: $e")
        println("💡 请确保服务器已启动: julia server.jl")
    end
end

end # module
