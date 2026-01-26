#!/usr/bin/env julia
"""
对比 Fortran 和 Julia 的弛豫时间计算结果

从 Fortran 输出中提取的数据:
T = 300 MeV, μ_B = 2 MeV
"""

using Printf

# 包含必要的模块
include("../../src/pnjl/PNJL.jl")
using .PNJL

include("../../src/Constants_PNJL.jl")
using .Constants_PNJL

include("../../src/pnjl/workflows/TransportWorkflow.jl")
using .TransportWorkflow

println("="^80)
println("Fortran vs Julia 弛豫时间对比")
println("="^80)

# ============================================================================
# Fortran 数据 (从输出中提取)
# ============================================================================

println("\n第一步：Fortran 数据")
println("-"^80)

# 温度和化学势
const ħc = 0.1973269804  # GeV·fm
T_MeV = 300.0
μ_B_MeV = 2.0
μ_MeV = μ_B_MeV / 3.0
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)

println("温度和化学势:")
@printf("  T = %.6f fm⁻¹ (%.1f MeV)\n", T, T_MeV)
@printf("  μ_B = %.6f fm⁻¹ (%.1f MeV)\n", 3*μ, μ_B_MeV)
@printf("  μ_q = %.6f fm⁻¹ (%.1f MeV)\n", μ, μ_MeV)

# Fortran 的能隙方程解
Φ_fortran = 0.84040830078505535
m_u_fortran = 0.040510297370740585
m_s_fortran = 1.0322945569740205

println("\n能隙方程解:")
@printf("  Φ = %.10f\n", Φ_fortran)
@printf("  m_u = %.10f fm⁻¹\n", m_u_fortran)
@printf("  m_s = %.10f fm⁻¹\n", m_s_fortran)

# Fortran 的数密度
n_u_fortran = 1.6809005906725343
n_s_fortran = 1.5376284189948928
n_ub_fortran = 1.6739110503594345
n_sb_fortran = 1.5311537136289866

println("\n数密度:")
@printf("  n_u  = %.10f fm⁻³\n", n_u_fortran)
@printf("  n_s  = %.10f fm⁻³\n", n_s_fortran)
@printf("  n_ub = %.10f fm⁻³\n", n_ub_fortran)
@printf("  n_sb = %.10f fm⁻³\n", n_sb_fortran)

# Fortran 的弛豫时间
τ_u_fortran = 0.58096796228598324
τ_s_fortran = 0.59203596086110832
τ_ub_fortran = 0.57915958660547406
τ_sb_fortran = 0.59022691868275623

Γ_u_fortran = 1.7212653105090621
Γ_s_fortran = 1.6890865861349258
Γ_ub_fortran = 1.7266398124584688
Γ_sb_fortran = 1.6942636269991858

println("\n弛豫时间:")
@printf("  τ_u  = %.10f fm\n", τ_u_fortran)
@printf("  τ_s  = %.10f fm\n", τ_s_fortran)
@printf("  τ_ub = %.10f fm\n", τ_ub_fortran)
@printf("  τ_sb = %.10f fm\n", τ_sb_fortran)

println("\n散射率 (Γ = 1/τ):")
@printf("  Γ_u  = %.10f fm⁻¹\n", Γ_u_fortran)
@printf("  Γ_s  = %.10f fm⁻¹\n", Γ_s_fortran)
@printf("  Γ_ub = %.10f fm⁻¹\n", Γ_ub_fortran)
@printf("  Γ_sb = %.10f fm⁻¹\n", Γ_sb_fortran)

# ============================================================================
# Julia 计算
# ============================================================================

println("\n" * "="^80)
println("第二步：Julia 计算")
println("-"^80)

println("\n使用 Fortran 的能隙方程解进行计算...")

# 使用 Fortran 的参数
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

println("\n计算弛豫时间...")
println("(这可能需要几分钟...)")

try
    # 直接使用 RelaxationTime 模块
    include("../../src/relaxtime/RelaxationTime.jl")
    using .RelaxationTime
    
    include("../../src/relaxtime/EffectiveCouplings.jl")
    using .EffectiveCouplings
    
    # 计算有效耦合系数
    println("\n计算 A 函数和有效耦合系数...")
    include("../../src/relaxtime/OneLoopIntegrals.jl")
    using .OneLoopIntegrals
    
    nodes_p, weights_p = RelaxationTime.AverageScatteringRate.gauleg(0.0, 20.0, 16)
    A_u = OneLoopIntegrals.A(m_u_fortran, μ, T, Φ_fortran, Φ_fortran, nodes_p, weights_p)
    A_s = OneLoopIntegrals.A(m_s_fortran, μ, T, Φ_fortran, Φ_fortran, nodes_p, weights_p)
    
    @printf("  A_u = %.10f fm⁻²\n", A_u)
    @printf("  A_s = %.10f fm⁻²\n", A_s)
    
    # 计算 G 函数
    G_u = EffectiveCouplings.calculate_G_from_A(A_u, m_u_fortran)
    G_s = EffectiveCouplings.calculate_G_from_A(A_s, m_s_fortran)
    
    @printf("  G_u = %.10f\n", G_u)
    @printf("  G_s = %.10f\n", G_s)
    
    # 计算有效耦合系数
    K_coeffs = EffectiveCouplings.calculate_effective_couplings(
        Constants_PNJL.G_fm2,
        Constants_PNJL.K_fm5,
        G_u,
        G_s
    )
    
    @printf("  K₀⁺ = %.10f fm²\n", K_coeffs.K0_plus)
    @printf("  K₁₂₃⁻ = %.10f fm² (π)\n", K_coeffs.K123_minus)
    @printf("  K₄₅₆₇⁻ = %.10f fm² (K)\n", K_coeffs.K4567_minus)
    
    quark_params_with_A = merge(quark_params, (A = (u = A_u, d = A_u, s = A_s),))
    
    println("\n计算弛豫时间...")
    result = RelaxationTime.relaxation_times(
        quark_params_with_A,
        thermo_params,
        K_coeffs;
        densities = densities,
        p_nodes = 16,
        angle_nodes = 16,
        phi_nodes = 16,
        n_sigma_points = 64
    )
    
    println("\nJulia 计算结果:")
    @printf("  τ_u  = %.10f fm\n", result.tau.u)
    @printf("  τ_s  = %.10f fm\n", result.tau.s)
    @printf("  τ_ub = %.10f fm\n", result.tau.ubar)
    @printf("  τ_sb = %.10f fm\n", result.tau.sbar)
    
    println("\n散射率 (Γ = 1/τ):")
    @printf("  Γ_u  = %.10f fm⁻¹\n", result.tau_inv.u)
    @printf("  Γ_s  = %.10f fm⁻¹\n", result.tau_inv.s)
    @printf("  Γ_ub = %.10f fm⁻¹\n", result.tau_inv.ubar)
    @printf("  Γ_sb = %.10f fm⁻¹\n", result.tau_inv.sbar)
    
    # ========================================================================
    # 对比结果
    # ========================================================================
    
    println("\n" * "="^80)
    println("第三步：对比结果")
    println("="^80)
    
    println("\n弛豫时间对比:")
    @printf("  %-10s  %-18s  %-18s  %-10s\n", "夸克", "Fortran (fm)", "Julia (fm)", "差异")
    @printf("  %-10s  %-18s  %-18s  %-10s\n", "----", "------------", "----------", "----")
    
    diff_u = abs(result.tau.u - τ_u_fortran) / τ_u_fortran * 100
    diff_s = abs(result.tau.s - τ_s_fortran) / τ_s_fortran * 100
    diff_ub = abs(result.tau.ubar - τ_ub_fortran) / τ_ub_fortran * 100
    diff_sb = abs(result.tau.sbar - τ_sb_fortran) / τ_sb_fortran * 100
    
    @printf("  %-10s  %18.10f  %18.10f  %9.2f%%\n", "u", τ_u_fortran, result.tau.u, diff_u)
    @printf("  %-10s  %18.10f  %18.10f  %9.2f%%\n", "s", τ_s_fortran, result.tau.s, diff_s)
    @printf("  %-10s  %18.10f  %18.10f  %9.2f%%\n", "ubar", τ_ub_fortran, result.tau.ubar, diff_ub)
    @printf("  %-10s  %18.10f  %18.10f  %9.2f%%\n", "sbar", τ_sb_fortran, result.tau.sbar, diff_sb)
    
    println("\n散射率对比 (Γ = 1/τ):")
    @printf("  %-10s  %-18s  %-18s  %-10s\n", "夸克", "Fortran (fm⁻¹)", "Julia (fm⁻¹)", "差异")
    @printf("  %-10s  %-18s  %-18s  %-10s\n", "----", "--------------", "------------", "----")
    
    diff_Γ_u = abs(result.tau_inv.u - Γ_u_fortran) / Γ_u_fortran * 100
    diff_Γ_s = abs(result.tau_inv.s - Γ_s_fortran) / Γ_s_fortran * 100
    diff_Γ_ub = abs(result.tau_inv.ubar - Γ_ub_fortran) / Γ_ub_fortran * 100
    diff_Γ_sb = abs(result.tau_inv.sbar - Γ_sb_fortran) / Γ_sb_fortran * 100
    
    @printf("  %-10s  %18.10f  %18.10f  %9.2f%%\n", "u", Γ_u_fortran, result.tau_inv.u, diff_Γ_u)
    @printf("  %-10s  %18.10f  %18.10f  %9.2f%%\n", "s", Γ_s_fortran, result.tau_inv.s, diff_Γ_s)
    @printf("  %-10s  %18.10f  %18.10f  %9.2f%%\n", "ubar", Γ_ub_fortran, result.tau_inv.ubar, diff_Γ_ub)
    @printf("  %-10s  %18.10f  %18.10f  %9.2f%%\n", "sbar", Γ_sb_fortran, result.tau_inv.sbar, diff_Γ_sb)
    
    # 判断结果
    max_diff = max(diff_u, diff_s, diff_ub, diff_sb)
    println("\n" * "="^80)
    println("结论")
    println("="^80)
    @printf("\n最大差异: %.2f%%\n", max_diff)
    
    if max_diff < 1.0
        println("\n✅ 结果一致! (差异 < 1%)")
        println("\n之前报告的 3 倍差异可能是:")
        println("  - 使用了不同的参数 (Φ, m, T, μ)")
        println("  - 单位换算错误")
        println("  - 数据提取错误")
    elseif max_diff < 10.0
        println("\n⚠️  有小差异 (1% - 10%)")
        println("\n这可能是由于:")
        println("  - 数值积分精度不同")
        println("  - 截断误差")
        println("  - 舍入误差累积")
    else
        println("\n❌ 有显著差异 (> 10%)")
        println("\n需要逐步检查计算链:")
        println("  1. 对比所有平均散射率 w_ij")
        println("  2. 对比总截面 σ(s)")
        println("  3. 对比散射振幅 M²")
        println("  4. 对比微分截面 dσ/dt")
    end
    println("="^80)
    
catch e
    println("\n❌ Julia 计算失败:")
    println(e)
    if isa(e, ErrorException)
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
end

println("\n分析完成!")
