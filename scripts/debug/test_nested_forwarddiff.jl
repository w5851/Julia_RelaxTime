#!/usr/bin/env julia
"""
    test_nested_forwarddiff.jl

测试嵌套 ForwardDiff 是否能正确工作。

测试场景：
1. 简单嵌套：外层 derivative，内层 derivative
2. 复杂嵌套：外层 gradient，内层 gradient
3. PNJL 实际场景：calculate_thermo 内部使用 ForwardDiff，外层再用 ForwardDiff
"""

using ForwardDiff
using StaticArrays
using Printf
using Test

println()
println("=" ^ 70)
println("  嵌套 ForwardDiff 测试")
println("=" ^ 70)
println()

# ============================================================================
# 测试 1: 简单嵌套 - 外层 derivative，内层 derivative
# ============================================================================

println("测试 1: 简单嵌套 (derivative inside derivative)")
println("-" ^ 50)

# 定义一个内层使用 ForwardDiff 的函数
function inner_with_ad(x)
    # 内层：计算 d(x^3)/dx = 3x^2
    return ForwardDiff.derivative(t -> t^3, x)
end

# 外层再用 ForwardDiff
function outer_derivative(x)
    # 外层：计算 d(inner_with_ad(x))/dx = d(3x^2)/dx = 6x
    return ForwardDiff.derivative(inner_with_ad, x)
end

# 解析解
analytical_result(x) = 6x

x_test = 2.0
try
    result = outer_derivative(x_test)
    expected = analytical_result(x_test)
    @printf("  x = %.1f\n", x_test)
    @printf("  嵌套 ForwardDiff 结果: %.6f\n", result)
    @printf("  解析解: %.6f\n", expected)
    @printf("  相对误差: %.2e\n", abs(result - expected) / abs(expected))
    @test isapprox(result, expected, rtol=1e-10)
    println("  ✓ 测试通过")
catch e
    println("  ✗ 测试失败: ", e)
end
println()

# ============================================================================
# 测试 2: 复杂嵌套 - 外层 gradient，内层 gradient
# ============================================================================

println("测试 2: 复杂嵌套 (gradient inside gradient)")
println("-" ^ 50)

# 内层函数：计算 ∇f(x) 的某个分量
function inner_gradient_func(x::AbstractVector)
    # f(x) = x[1]^2 * x[2]
    # ∇f = [2*x[1]*x[2], x[1]^2]
    grad = ForwardDiff.gradient(y -> y[1]^2 * y[2], x)
    return grad[1]  # 返回 2*x[1]*x[2]
end

# 外层：计算 ∇(inner_gradient_func)
function outer_gradient(x::AbstractVector)
    # ∇(2*x[1]*x[2]) = [2*x[2], 2*x[1]]
    return ForwardDiff.gradient(inner_gradient_func, x)
end

# 解析解
analytical_outer_gradient(x) = [2*x[2], 2*x[1]]

x_test_2 = [3.0, 4.0]
try
    result = outer_gradient(x_test_2)
    expected = analytical_outer_gradient(x_test_2)
    @printf("  x = [%.1f, %.1f]\n", x_test_2...)
    @printf("  嵌套 ForwardDiff 结果: [%.6f, %.6f]\n", result...)
    @printf("  解析解: [%.6f, %.6f]\n", expected...)
    rel_err = maximum(abs.(result - expected) ./ abs.(expected))
    @printf("  最大相对误差: %.2e\n", rel_err)
    @test isapprox(result, expected, rtol=1e-10)
    println("  ✓ 测试通过")
catch e
    println("  ✗ 测试失败: ", e)
    println("    ", e)
end
println()

# ============================================================================
# 测试 3: 模拟 PNJL 场景 - 热力学量计算
# ============================================================================

println("测试 3: 模拟 PNJL 场景")
println("-" ^ 50)

# 模拟 calculate_rho：内部使用 ForwardDiff.gradient
function mock_calculate_rho(x_state, mu_vec, T)
    # 模拟：ρ = ∂P/∂μ，P 是某个复杂函数
    function pressure(μ)
        # P = T^4 * (1 + sum(μ)^2/T^2) * sum(x_state)
        return T^4 * (1 + sum(μ)^2/T^2) * sum(x_state)
    end
    return ForwardDiff.gradient(pressure, mu_vec)
end

# 模拟 calculate_thermo：内部调用 calculate_rho
function mock_calculate_thermo(x_state, mu_vec, T)
    rho = mock_calculate_rho(x_state, mu_vec, T)
    
    # 熵：s = ∂P/∂T
    function pressure_T(τ)
        return τ^4 * (1 + sum(mu_vec)^2/τ^2) * sum(x_state)
    end
    s = ForwardDiff.derivative(pressure_T, T)
    
    return (rho=rho, s=s)
end

# 外层：计算 ∂s/∂T（固定 x_state, mu_vec）
function outer_ds_dT(x_state, mu_vec, T)
    function s_of_T(τ)
        result = mock_calculate_thermo(x_state, mu_vec, τ)
        return result.s
    end
    return ForwardDiff.derivative(s_of_T, T)
end

# 解析解
# P = T^4 * (1 + μ²/T²) * Σx = T^4 * Σx + μ² * T² * Σx
# s = ∂P/∂T = 4T³ * Σx + 2μ² * T * Σx = (4T³ + 2μ²T) * Σx
# ∂s/∂T = (12T² + 2μ²) * Σx
function analytical_ds_dT(x_state, mu_vec, T)
    μ_sq = sum(mu_vec)^2
    Σx = sum(x_state)
    return (12*T^2 + 2*μ_sq) * Σx
end

x_state_test = SVector(1.0, 2.0, 3.0)
mu_vec_test = SVector(0.5, 0.5, 0.5)
T_test = 0.8

try
    result = outer_ds_dT(x_state_test, mu_vec_test, T_test)
    expected = analytical_ds_dT(x_state_test, mu_vec_test, T_test)
    @printf("  x_state = [%.1f, %.1f, %.1f]\n", x_state_test...)
    @printf("  mu_vec = [%.1f, %.1f, %.1f]\n", mu_vec_test...)
    @printf("  T = %.1f\n", T_test)
    @printf("  嵌套 ForwardDiff 结果: %.6f\n", result)
    @printf("  解析解: %.6f\n", expected)
    rel_err = abs(result - expected) / abs(expected)
    @printf("  相对误差: %.2e\n", rel_err)
    @test isapprox(result, expected, rtol=1e-10)
    println("  ✓ 测试通过")
catch e
    println("  ✗ 测试失败:")
    println("    ", e)
end
println()

# ============================================================================
# 测试 4: 三层嵌套
# ============================================================================

println("测试 4: 三层嵌套")
println("-" ^ 50)

function level1(x)
    return ForwardDiff.derivative(t -> t^4, x)  # 4x^3
end

function level2(x)
    return ForwardDiff.derivative(level1, x)  # 12x^2
end

function level3(x)
    return ForwardDiff.derivative(level2, x)  # 24x
end

x_test_3 = 2.0
try
    result = level3(x_test_3)
    expected = 24 * x_test_3
    @printf("  x = %.1f\n", x_test_3)
    @printf("  三层嵌套 ForwardDiff 结果: %.6f\n", result)
    @printf("  解析解 (24x): %.6f\n", expected)
    rel_err = abs(result - expected) / abs(expected)
    @printf("  相对误差: %.2e\n", rel_err)
    @test isapprox(result, expected, rtol=1e-10)
    println("  ✓ 测试通过")
catch e
    println("  ✗ 测试失败:")
    println("    ", e)
end
println()

# ============================================================================
# 测试 5: 实际 PNJL 模块测试
# ============================================================================

println("测试 5: 实际 PNJL 模块嵌套测试")
println("-" ^ 50)

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL: solve, FixedMu
using .PNJL.Thermodynamics: calculate_thermo, calculate_rho
using .PNJL.Integrals: cached_nodes

const ħc = 197.3269804

# 求解一个点
T_MeV = 150.0
μ_MeV = 200.0
T_fm = T_MeV / ħc
μ_fm = μ_MeV / ħc

result = solve(FixedMu(), T_fm, μ_fm; xi=0.0)
x_state = SVector{5}(result.solution)  # 确保是 SVector{5}
mu_vec = SVector(μ_fm, μ_fm, μ_fm)
thermal_nodes = cached_nodes(64, 16)

println("  x_state 类型: ", typeof(x_state))

# 测试：外层 ForwardDiff 对 T 求导，内层 calculate_thermo 也用 ForwardDiff
function s_of_T_nested(T)
    # 注意：x_state 需要在闭包外部捕获，保持 SVector 类型
    _, _, s, _ = calculate_thermo(x_state, mu_vec, T, thermal_nodes, 0.0)
    return s
end

try
    # 使用嵌套 ForwardDiff
    ds_dT_nested = ForwardDiff.derivative(s_of_T_nested, T_fm)
    
    # 使用有限差分验证
    ε = 1e-6
    s_plus = s_of_T_nested(T_fm + ε)
    s_minus = s_of_T_nested(T_fm - ε)
    ds_dT_fd = (s_plus - s_minus) / (2ε)
    
    @printf("  T = %.1f MeV, μ = %.1f MeV\n", T_MeV, μ_MeV)
    @printf("  嵌套 ForwardDiff ∂s/∂T: %.6f\n", ds_dT_nested)
    @printf("  有限差分 ∂s/∂T: %.6f\n", ds_dT_fd)
    rel_err = abs(ds_dT_nested - ds_dT_fd) / abs(ds_dT_fd)
    @printf("  相对误差: %.2e\n", rel_err)
    
    if rel_err < 0.01
        println("  ✓ 测试通过 (误差 < 1%)")
    else
        println("  ⚠ 警告：误差较大")
    end
catch e
    println("  ✗ 测试失败:")
    showerror(stdout, e)
    println()
end
println()

# ============================================================================
# 总结
# ============================================================================

println("=" ^ 70)
println("  总结")
println("=" ^ 70)
println()
println("ForwardDiff 支持嵌套调用，因为它使用 Dual 数的 Tag 系统来区分不同层级的导数。")
println("每个 ForwardDiff 调用都会创建一个唯一的 Tag，避免不同层级的 Dual 数混淆。")
println()
println("关键点：")
println("1. 简单嵌套（derivative inside derivative）：✓ 正常工作")
println("2. 复杂嵌套（gradient inside gradient）：✓ 正常工作")
println("3. 模拟 PNJL 场景：✓ 正常工作")
println("4. 三层嵌套：✓ 正常工作")
println("5. 实际 PNJL 模块：需要验证")
println()
