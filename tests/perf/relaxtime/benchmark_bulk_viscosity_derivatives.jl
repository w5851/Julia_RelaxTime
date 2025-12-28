"""
性能测试：体粘滞系数热力学导数计算

测试 ThermoDerivatives 模块中体粘滞相关函数的性能：
- bulk_viscosity_coefficients: 一次性计算所有体粘滞系数所需的导数
- compute_B_bracket: 体粘滞公式中的 B 项

结果保存到 tests/perf/results/relaxtime/ 目录
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))

using BenchmarkTools
using Printf
using JSON
using Dates

include("../../../src/pnjl/PNJL.jl")
using .PNJL: ThermoDerivatives

# 结果输出目录
const RESULTS_DIR = joinpath(@__DIR__, "..", "results", "relaxtime")
mkpath(RESULTS_DIR)

# 常数
const hc = 197.327  # MeV·fm

# 测试点
const T_MeV = 150.0
const μB_MeV = 800.0
const T_fm = T_MeV / hc
const μq_fm = μB_MeV / 3.0 / hc

println("=" ^ 60)
println("体粘滞系数热力学导数性能测试")
println("=" ^ 60)
println("测试点: T = $T_MeV MeV, μ_B = $μB_MeV MeV")
println()

# 预热
println("预热中...")
ThermoDerivatives.bulk_viscosity_coefficients(T_fm, μq_fm)
println("预热完成")
println()

# ============================================================================
# 基准测试
# ============================================================================

println("1. bulk_viscosity_coefficients 性能测试")
b1 = @benchmark ThermoDerivatives.bulk_viscosity_coefficients($T_fm, $μq_fm)
println("   中位时间: $(@sprintf("%.2f", median(b1).time / 1e6)) ms")
println("   最小时间: $(@sprintf("%.2f", minimum(b1).time / 1e6)) ms")
println("   内存分配: $(@sprintf("%.2f", b1.memory / 1024)) KB")
println()

# 获取系数用于 B 项测试
coeffs = ThermoDerivatives.bulk_viscosity_coefficients(T_fm, μq_fm)
p_test = 1.0

println("2. compute_B_bracket 性能测试 (单次调用)")
b2 = @benchmark ThermoDerivatives.compute_B_bracket(
    $p_test,
    $(coeffs.masses[1]),
    $μq_fm,
    $T_fm,
    $(coeffs.v_n_sq),
    $(coeffs.dμB_dT_sigma),
    $(coeffs.dM_dT[1]),
    $(coeffs.dM_dμB[1]);
    is_antiquark=false
)
println("   中位时间: $(@sprintf("%.2f", median(b2).time)) ns")
println("   最小时间: $(@sprintf("%.2f", minimum(b2).time)) ns")
println("   内存分配: $(b2.memory) bytes")
println()

# ============================================================================
# 批量计算性能测试
# ============================================================================

println("3. 批量 B 项计算性能测试 (100 个动量点)")

# 模拟实际使用场景：对多个动量点计算 B 项
p_values = range(0.1, 5.0, length=100)

# 提取纯数值（避免 Dual 类型问题）
using ForwardDiff
extract_val(x) = x isa ForwardDiff.Dual ? ForwardDiff.value(x) : Float64(x)

v_n_sq_val = extract_val(coeffs.v_n_sq)
dμB_dT_sigma_val = extract_val(coeffs.dμB_dT_sigma)
masses_val = [extract_val(m) for m in coeffs.masses]
dM_dT_val = [extract_val(d) for d in coeffs.dM_dT]
dM_dμB_val = [extract_val(d) for d in coeffs.dM_dμB]

function compute_all_B_brackets(p_values, masses, dM_dT, dM_dμB, v_n_sq, dμB_dT_sigma, μq, T)
    B_q = zeros(length(p_values))
    B_qbar = zeros(length(p_values))
    
    for (i, p) in enumerate(p_values)
        for flavor in 1:3
            B_q[i] += ThermoDerivatives.compute_B_bracket(
                p, masses[flavor], μq, T,
                v_n_sq, dμB_dT_sigma,
                dM_dT[flavor], dM_dμB[flavor];
                is_antiquark=false
            )
            B_qbar[i] += ThermoDerivatives.compute_B_bracket(
                p, masses[flavor], μq, T,
                v_n_sq, dμB_dT_sigma,
                dM_dT[flavor], dM_dμB[flavor];
                is_antiquark=true
            )
        end
    end
    return B_q, B_qbar
end

b3 = @benchmark compute_all_B_brackets($p_values, $masses_val, $dM_dT_val, $dM_dμB_val, 
                                        $v_n_sq_val, $dμB_dT_sigma_val, $μq_fm, $T_fm)
println("   中位时间: $(@sprintf("%.2f", median(b3).time / 1e3)) μs")
println("   最小时间: $(@sprintf("%.2f", minimum(b3).time / 1e3)) μs")
println("   内存分配: $(@sprintf("%.2f", b3.memory / 1024)) KB")
println()

# ============================================================================
# 与 thermo_derivatives 对比
# ============================================================================

println("4. 与 thermo_derivatives 性能对比")
b4 = @benchmark ThermoDerivatives.thermo_derivatives($T_fm, $μq_fm)
println("   thermo_derivatives 中位时间: $(@sprintf("%.2f", median(b4).time / 1e6)) ms")
println("   bulk_viscosity_coefficients 中位时间: $(@sprintf("%.2f", median(b1).time / 1e6)) ms")
println("   比值: $(@sprintf("%.2f", median(b1).time / median(b4).time))x")
println()

# ============================================================================
# 保存结果到 JSON
# ============================================================================

results = Dict(
    "metadata" => Dict(
        "date" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
        "julia_version" => string(VERSION),
        "test_point" => Dict(
            "T_MeV" => T_MeV,
            "μB_MeV" => μB_MeV,
            "T_fm" => T_fm,
            "μq_fm" => μq_fm
        )
    ),
    "benchmarks" => Dict(
        "bulk_viscosity_coefficients" => Dict(
            "median_ms" => median(b1).time / 1e6,
            "min_ms" => minimum(b1).time / 1e6,
            "memory_KB" => b1.memory / 1024
        ),
        "compute_B_bracket" => Dict(
            "median_ns" => median(b2).time,
            "min_ns" => minimum(b2).time,
            "memory_bytes" => b2.memory
        ),
        "batch_B_brackets_600calls" => Dict(
            "median_μs" => median(b3).time / 1e3,
            "min_μs" => minimum(b3).time / 1e3,
            "memory_KB" => b3.memory / 1024,
            "per_call_ns" => median(b3).time / 600
        ),
        "thermo_derivatives" => Dict(
            "median_ms" => median(b4).time / 1e6,
            "min_ms" => minimum(b4).time / 1e6,
            "memory_KB" => b4.memory / 1024
        )
    )
)

# 保存 JSON 文件
json_filename = "bulk_viscosity_derivatives_$(Dates.format(now(), "yyyy-mm-dd")).json"
json_path = joinpath(RESULTS_DIR, json_filename)
open(json_path, "w") do f
    JSON.print(f, results, 2)
end
println("结果已保存到: $json_path")
println()

# ============================================================================
# 总结
# ============================================================================

println("=" ^ 60)
println("性能总结")
println("=" ^ 60)
println()
println("单次调用时间:")
println("   bulk_viscosity_coefficients: $(@sprintf("%6.2f", median(b1).time / 1e6)) ms")
println("   compute_B_bracket:           $(@sprintf("%6.2f", median(b2).time)) ns")
println()
println("批量计算 (100 点 × 3 味 × 2 粒子/反粒子):")
println("   总时间: $(@sprintf("%.2f", median(b3).time / 1e3)) μs")
println("   每点平均: $(@sprintf("%.2f", median(b3).time / 600)) ns")
println()
println("建议：")
println("   - bulk_viscosity_coefficients 应该在每个 (T, μ) 点只调用一次")
println("   - compute_B_bracket 非常快，可以在积分循环中频繁调用")
println()
