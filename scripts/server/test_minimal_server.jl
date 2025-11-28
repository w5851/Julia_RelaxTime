#!/usr/bin/env julia

"""
最小化测试服务器 - 只提供静态文件服务
"""

using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(REPO_ROOT)

using HTTP

"""
静态文件服务
"""
function serve_static_file(path::String)
    println("  [静态文件] 请求: $path")
    
    # 安全检查：防止目录遍历
    if contains(path, "..") || contains(path, "\\")
        println("  [静态文件] 403 Forbidden")
        return HTTP.Response(403, "Forbidden")
    end
    
    # 移除路径开头的 /
    clean_path = startswith(path, "/") ? path[2:end] : path
    
    file_path = if path == "/" || path == ""
        joinpath(REPO_ROOT, "web", "index.html")
    else
        joinpath(REPO_ROOT, "web", clean_path)
    end
    
    println("  [静态文件] 文件路径: $file_path")
    
    # 检查文件是否存在
    if !isfile(file_path)
        println("  [静态文件] 404 Not Found")
        return HTTP.Response(404, "Not Found")
    end
    
    # 确定Content-Type
    ext = lowercase(splitext(file_path)[2])
    content_type = if ext == ".html"
        "text/html; charset=utf-8"
    elseif ext == ".css"
        "text/css; charset=utf-8"
    elseif ext == ".js"
        "application/javascript; charset=utf-8"
    elseif ext == ".json"
        "application/json; charset=utf-8"
    else
        "application/octet-stream"
    end
    
    # 读取并返回文件
    try
        content = read(file_path)
        println("  [静态文件] ✓ 200 OK ($(length(content)) bytes)")
        return HTTP.Response(200, 
            ["Content-Type" => content_type,
             "Access-Control-Allow-Origin" => "*"], 
            body=content)
    catch e
        println("  [静态文件] ✗ 500 Error: $e")
        @error "Error serving file" file=file_path exception=e
        return HTTP.Response(500, "Internal Server Error")
    end
end

"""
请求路由
"""
function route_request(req::HTTP.Request)
    println("\n[请求] $(req.method) $(req.target)")
    
    path = HTTP.URIs.unescapeuri(req.target)
    
    # 移除查询参数
    path = split(path, '?')[1]
    
    return serve_static_file(path)
end

# 启动服务器
const DEFAULT_PORT = 8080
port = DEFAULT_PORT

println("\n" * "="^60)
println("🧪 最小化测试服务器")
println("="^60)
println("📍 http://localhost:$port")
println("="^60 * "\n")

try
    HTTP.serve(route_request, "0.0.0.0", port; verbose=true)
catch e
    @error "Server error" exception=e
    rethrow(e)
end
