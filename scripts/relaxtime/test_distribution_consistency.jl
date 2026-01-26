#!/usr/bin/env julia
"""
测试 Julia 和 Fortran 的分布函数是否一致

根据调查发现:
- Fortran: fphi(E - μ, T, Φ, Φbar) 和 fphibar(E + μ, T, Φ, Φbar)
- Julia: quark_distribution(E, μ, T, Φ, Φbar) 内部计算 E - μ
         antiquark_distribution(E, μ, T, Φ, Φbar) 内部计算 E + μ

需要验证这两种调用方式是否等价
"""

using Printf

# 包含 Julia 的分布函数实现
include("../../src/QuarkDistribution.jl")
using .PNJLQuarkDistributions

# Fortran 的分布函数实现（手动复制）
function fphi_fortran(x::Float64, T::Float64, Phi1::Float64, Phi2::Float64)
    ee = exp(-x / T)
    ee2 = ee * ee
    ee3 = ee2 * ee
    numerator = Phi1 * ee + 2.0 * Phi2 * ee2 + ee3
    denominator = 1.0 + 3.0 * Phi1 * ee + 3.0 * Phi2 * ee2 + ee3
    return numerator / denominator
end

function fphibar_fortran(x::Float64, T::Float64, Phi1::Float64, Phi2::Float64)
    ee = exp(-x / T)
    ee2 = ee * ee
    ee3 = ee2 * ee
    numerator = Phi2 * ee + 2.0 * Phi1 * ee2 + ee3
    denominator = 1.0 + 3.0 * Phi2 * ee + 3.0 * Phi1 * ee2 + ee3
    return numerator / denominator
end

# 测试参数（使用 T=300 MeV, μ=2 MeV 的情况）
const ħc = 0.1973269804  # GeV·fm
T_MeV = 300.0
μ_MeV = 2.0
T = T_MeV / (ħc * 1000.0)  # 转换为 fm⁻¹
μ = μ_MeV / (ħc * 1000.0)  # 转换为 fm⁻¹
Φ = 0.99999994
Φbar = 0.99999994

println("=" ^ 80)
println("分布函数一致性测试")
println("=" ^ 80)
println()
println("测试参数:")
@printf("  T = %.6f fm⁻¹ (%.1f MeV)\n", T, T_MeV)
@printf("  μ = %.6f fm⁻¹ (%.1f MeV)\n", μ, μ_MeV)
@printf("  Φ = %.8f\n", Φ)
@printf("  Φbar = %.8f\n", Φbar)
println()

# 测试多个能量点
E_values = [0.5, 1.0, 2.0, 5.0, 10.0]  # fm⁻¹

println("=" ^ 80)
println("夸克分布函数对比")
println("=" ^ 80)
println()
println("E (fm⁻¹) | Julia quark_dist | Fortran fphi(E-μ) | 相对差异")
println("-" ^ 80)

for E in E_values
    # Julia 方式: quark_distribution(E, μ, T, Φ, Φbar)
    julia_val = quark_distribution(E, μ, T, Φ, Φbar)
    
    # Fortran 方式: fphi(E - μ, T, Φ, Φbar)
    fortran_val = fphi_fortran(E - μ, T, Φ, Φbar)
    
    rel_diff = abs(julia_val - fortran_val) / max(abs(julia_val), abs(fortran_val), 1e-300)
    
    @printf("%8.4f | %16.10e | %16.10e | %12.6e\n", E, julia_val, fortran_val, rel_diff)
end

println()
println("=" ^ 80)
println("反夸克分布函数对比")
println("=" ^ 80)
println()
println("E (fm⁻¹) | Julia antiquark_dist | Fortran fphibar(E+μ) | 相对差异")
println("-" ^ 80)

for E in E_values
    # Julia 方式: antiquark_distribution(E, μ, T, Φ, Φbar)
    julia_val = antiquark_distribution(E, μ, T, Φ, Φbar)
    
    # Fortran 方式: fphibar(E + μ, T, Φ, Φbar)
    fortran_val = fphibar_fortran(E + μ, T, Φ, Φbar)
    
    rel_diff = abs(julia_val - fortran_val) / max(abs(julia_val), abs(fortran_val), 1e-300)
    
    @printf("%8.4f | %16.10e | %16.10e | %12.6e\n", E, julia_val, fortran_val, rel_diff)
end

println()
println("=" ^ 80)
println("结论")
println("=" ^ 80)
println()
println("如果相对差异 < 1e-10，则分布函数实现一致")
println("如果相对差异较大，则需要检查实现差异")
println()
