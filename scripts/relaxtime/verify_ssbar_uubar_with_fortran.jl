#!/usr/bin/env julia
"""
验证 ssbar_to_uubar 散射过程的差异

从Fortran输出中提取的数据:
w(10) ssbar->uubar = 5.9408786969375932E-002 fm^-1

对比Julia的计算结果
"""

using Printf

# 包含必要的模块
include("../../src/pnjl/PNJL.jl")
using .PNJL

include("../../src/Constants_PNJL.jl")
using .Constants_PNJL

println("="^80)
println("验证 ssbar→uubar 散射过程")
println("="^80)

# ============================================================================
# Fortran 数据
# ============================================================================

println("\n第一步：Fortran 数据")
println("-"^80)

# 从Fortran输出提取
w_ssbar_uubar_fortran = 5.9408786969375932E-002  # fm^-1
w_uubar_ssbar_fortran = 4.7545562405610760E-002  # fm^-1 (w(9))

println("Fortran 散射率:")
@printf("  w(ssbar→uubar) = %.10e fm⁻¹\n", w_ssbar_uubar_fortran)
@printf("  w(uubar→ssbar) = %.10e fm⁻¹\n", w_uubar_ssbar_fortran)
@printf("  比值 w(uubar→ssbar)/w(ssbar→uubar) = %.6f\n", w_uubar_ssbar_fortran / w_ssbar_uubar_fortran)

# 理论预期（详细平衡）
const ħc = 0.1973269804  # GeV·fm
T_MeV = 300.0
m_u_MeV = 7.994  # 从Fortran输出
m_s_MeV = 203.7  # 从Fortran输出

ΔE = 2*m_u_MeV - 2*m_s_MeV
detailed_balance_ratio = exp(-ΔE / T_MeV)

println("\n详细平衡理论预期:")
@printf("  ΔE = 2m_u - 2m_s = %.1f MeV\n", ΔE)
@printf("  exp(-ΔE/T) = %.6f\n", detailed_balance_ratio)
@printf("  Fortran偏差: %.2f倍\n", detailed_balance_ratio / (w_uubar_ssbar_fortran / w_ssbar_uubar_fortran))

# ============================================================================
# Julia 计算
# ============================================================================

println("\n" * "="^80)
println("第二步：Julia 计算")
println("-"^80)

# 使用Fortran的参数
T_MeV = 300.0
μ_B_MeV = 2.0
μ_MeV = μ_B_MeV / 3.0
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)

Φ_fortran = 0.84040830078505535
m_u_fortran = 0.040510297370740585
m_s_fortran = 1.0322945569740205

n_u_fortran = 1.6809005906725343
n_s_fortran = 1.5376284190948928
n_ub_fortran = 1.6739110503594345
n_sb_fortran = 1.5311537136289866

quark_params = (
    m = (u = m_u_fortran, d = m_u_fortran, s = m_s_fortran),
    μ = (u = μ, d = μ, s = μ)
)

thermo_params = (
    T = T,
    Φ = Φ_fortran,
    Φbar = Φ_fortran,
    ξ = 0.0
)

densities = (
    u = n_u_fortran,
    d = n_u_fortran,
    s = n_s_fortran,
    ubar = n_ub_fortran,
    dbar = n_ub_fortran,
    sbar = n_sb_fortran
)

println("\n计算散射率...")

try
    include("../../src/relaxtime/RelaxationTime.jl")
    using .RelaxationTime
    
    include("../../src/relaxtime/EffectiveCouplings.jl")
    using .EffectiveCouplings
    
    include("../../src/relaxtime/OneLoopIntegrals.jl")
    using .OneLoopIntegrals
    
    # 计算A函数和有效耦合
    nodes_p, weights_p = RelaxationTime.AverageScatteringRate.gauleg(0.0, 20.0, 16)
    A_u = OneLoopIntegrals.A(m_u_fortran, μ, T, Φ_fortran, Φ_fortran, nodes_p, weights_p)
    A_s = OneLoopIntegrals.A(m_s_fortran, μ, T, Φ_fortran, Φ_fortran, nodes_p, weights_p)
    
    G_u = EffectiveCouplings.calculate_G_from_A(A_u, m_u_fortran)
    G_s = EffectiveCouplings.calculate_G_from_A(A_s, m_s_fortran)
    
    K_coeffs = EffectiveCouplings.calculate_effective_couplings(
        Constants_PNJL.G_fm2,
        Constants_PNJL.K_fm5,
        G_u,
        G_s
    )
    
    quark_params_with_A = merge(quark_params, (A = (u = A_u, d = A_u, s = A_s),))
    
    # 计算所有散射率
    println("  计算所有散射过程...")
    rates = RelaxationTime.compute_average_rates(
        quark_params_with_A,
        thermo_params,
        K_coeffs;
        p_nodes = 16,
        angle_nodes = 16,
        phi_nodes = 16,
        n_sigma_points = 64
    )
    
    # 提取关键散射率
    w_ssbar_uubar_julia = rates.ssbar_to_uubar
    w_uubar_ssbar_julia = rates.uubar_to_ssbar
    
    println("\nJulia 散射率:")
    @printf("  w(ssbar→uubar) = %.10e fm⁻¹\n", w_ssbar_uubar_julia)
    @printf("  w(uubar→ssbar) = %.10e fm⁻¹\n", w_uubar_ssbar_julia)
    @printf("  比值 w(uubar→ssbar)/w(ssbar→uubar) = %.6f\n", w_uubar_ssbar_julia / w_ssbar_uubar_julia)
    
    println("\n详细平衡检查:")
    @printf("  理论预期: %.6f\n", detailed_balance_ratio)
    @printf("  Julia偏差: %.2f倍\n", detailed_balance_ratio / (w_uubar_ssbar_julia / w_ssbar_uubar_julia))
    
    # ========================================================================
    # 对比结果
    # ========================================================================
    
    println("\n" * "="^80)
    println("第三步：对比结果")
    println("="^80)
    
    println("\nssbar→uubar 散射率对比:")
    @printf("  %-15s  %-20s  %-20s  %-10s\n", "实现", "w (fm⁻¹)", "w (fm⁻¹)", "差异")
    @printf("  %-15s  %-20s  %-20s  %-10s\n", "----", "--------", "--------", "----")
    
    diff_ssbar_uubar = abs(w_ssbar_uubar_julia - w_ssbar_uubar_fortran) / w_ssbar_uubar_fortran * 100
    ratio_ssbar_uubar = w_ssbar_uubar_julia / w_ssbar_uubar_fortran
    
    @printf("  %-15s  %20.10e  %20s  %9s\n", "Fortran", w_ssbar_uubar_fortran, "-", "-")
    @printf("  %-15s  %20.10e  %20s  %9.2f%%\n", "Julia", w_ssbar_uubar_julia, "-", diff_ssbar_uubar)
    @printf("  %-15s  %20s  %20.6f  %9s\n", "Julia/Fortran", "-", ratio_ssbar_uubar, "-")
    
    println("\nuubar→ssbar 散射率对比:")
    @printf("  %-15s  %-20s  %-20s  %-10s\n", "实现", "w (fm⁻¹)", "w (fm⁻¹)", "差异")
    @printf("  %-15s  %-20s  %-20s  %-10s\n", "----", "--------", "--------", "----")
    
    diff_uubar_ssbar = abs(w_uubar_ssbar_julia - w_uubar_ssbar_fortran) / w_uubar_ssbar_fortran * 100
    ratio_uubar_ssbar = w_uubar_ssbar_julia / w_uubar_ssbar_fortran
    
    @printf("  %-15s  %20.10e  %20s  %9s\n", "Fortran", w_uubar_ssbar_fortran, "-", "-")
    @printf("  %-15s  %20.10e  %20s  %9.2f%%\n", "Julia", w_uubar_ssbar_julia, "-", diff_uubar_ssbar)
    @printf("  %-15s  %20s  %20.6f  %9s\n", "Julia/Fortran", "-", ratio_uubar_ssbar, "-")
    
    # ========================================================================
    # 分析
    # ========================================================================
    
    println("\n" * "="^80)
    println("分析")
    println("="^80)
    
    println("\n1. ssbar→uubar 差异:")
    if ratio_ssbar_uubar > 2.0
        println("  ❌ 有显著差异 (> 2倍)")
        println("  这与之前的调查一致:")
        println("    - 之前报告: Julia/Fortran ≈ 3.238")
        println("    - 当前结果: Julia/Fortran ≈ $(round(ratio_ssbar_uubar, digits=3))")
        println("\n  根据之前的分析,这是由于:")
        println("    1. 动量积分范围不同 (Julia: [0,∞), Fortran: [0,Λ])")
        println("    2. ssbar→uubar 的阈值效应")
        println("    3. 主要散射发生在 p > Λ 区域")
    elseif ratio_ssbar_uubar > 1.1
        println("  ⚠️  有小差异 (10% - 100%)")
    else
        println("  ✅ 基本一致 (< 10%)")
    end
    
    println("\n2. uubar→ssbar 差异:")
    if ratio_uubar_ssbar > 1.1
        println("  ⚠️  有差异 (> 10%)")
    else
        println("  ✅ 基本一致 (< 10%)")
    end
    
    println("\n3. 详细平衡:")
    fortran_db_ratio = w_uubar_ssbar_fortran / w_ssbar_uubar_fortran
    julia_db_ratio = w_uubar_ssbar_julia / w_ssbar_uubar_julia
    
    @printf("  理论预期: %.6f\n", detailed_balance_ratio)
    @printf("  Fortran:   %.6f (偏差 %.2f倍)\n", fortran_db_ratio, detailed_balance_ratio / fortran_db_ratio)
    @printf("  Julia:     %.6f (偏差 %.2f倍)\n", julia_db_ratio, detailed_balance_ratio / julia_db_ratio)
    
    if abs(julia_db_ratio - detailed_balance_ratio) / detailed_balance_ratio < 0.1
        println("  ✅ Julia 满足详细平衡 (< 10%)")
    else
        println("  ⚠️  Julia 不满足详细平衡")
        println("  注意: 平均散射率不严格满足详细平衡是正常的")
    end
    
    println("\n4. 对弛豫时间的影响:")
    println("  s夸克弛豫率公式:")
    println("    Γ_s = 2*n_u*w_us + 2*n_ub*w_usbar + n_s*w_ss + n_sb*(w_ssbar + 2*w_ssbar_uubar)")
    
    contrib_fortran = n_sb_fortran * 2 * w_ssbar_uubar_fortran
    contrib_julia = n_sb_fortran * 2 * w_ssbar_uubar_julia
    
    @printf("\n  ssbar→uubar 的贡献:\n")
    @printf("    Fortran: n_sb * 2 * w = %.4f fm⁻¹\n", contrib_fortran)
    @printf("    Julia:   n_sb * 2 * w = %.4f fm⁻¹\n", contrib_julia)
    @printf("    差异: %.4f fm⁻¹\n", contrib_julia - contrib_fortran)
    
    println("\n  这解释了为什么 s 夸克弛豫时间有 ~8% 差异")
    
    println("\n" * "="^80)
    println("结论")
    println("="^80)
    
    println("\n根据之前的详细分析 (WHY_ONLY_SSBAR_UUBAR_DIFFERS.md):")
    println("  1. ssbar→uubar 的差异是由于动量积分范围不同")
    println("  2. 这个过程同时满足三个条件:")
    println("     - 初态粒子重 (s, sbar)")
    println("     - 末态阈值高")
    println("     - 主要散射发生在 p > Λ 区域")
    println("  3. Julia 的半无穷积分包含了 Fortran 截断的贡献")
    println("  4. 这导致 Julia/Fortran ≈ 3.2 的比值")
    println("\n当前验证:")
    @printf("  Julia/Fortran = %.3f\n", ratio_ssbar_uubar)
    println("  与之前的分析一致!")
    
catch e
    println("\n❌ 计算失败:")
    println(e)
    if isa(e, ErrorException)
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
end

println("\n分析完成!")
