# 测试服务器的独立脚本
using Pkg
Pkg.activate(".")

using HTTP
using JSON3

println("Testing static file serving...")

# 模拟serve_static_file的逻辑
function test_serve(path::String)
    println("\n测试路径: '$path'")
    
    # 安全检查
    if contains(path, "..") || contains(path, "\\")
        println("  ✗ 403 Forbidden")
        return
    end
    
    # 移除路径开头的 /
    clean_path = startswith(path, "/") ? path[2:end] : path
    println("  清理后: '$clean_path'")
    
    file_path = if path == "/" || path == ""
        joinpath(@__DIR__, "web", "index.html")
    else
        joinpath(@__DIR__, "web", clean_path)
    end
    
    println("  文件路径: $file_path")
    println("  文件存在: $(isfile(file_path))")
    
    if isfile(file_path)
        try
            content = read(file_path)
            println("  ✓ 文件大小: $(length(content)) bytes")
        catch e
            println("  ✗ 读取错误: $e")
        end
    else
        println("  ✗ 文件不存在")
    end
end

# 测试各种路径
test_serve("/")
test_serve("/index.html")
test_serve("/js/api.js")
test_serve("/css/style.css")
test_serve("/simple_test.html")
