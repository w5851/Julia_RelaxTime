#!/usr/bin/env julia
"""
提取并对比 Fortran 和 Julia 的弛豫时间计算

目标：
1. 从 Fortran 输出中提取弛豫时间值
2. 使用 Julia 计算相同条件下的弛豫时间
3. 逐步对比计算链的每个环节

计算链：
能隙方程 → (Φ, m) ✅
    ↓
A 函数 → A_f ✅
    ↓
G^f = -3/(4π²) × m_f × A_f ✅
    ↓
有效耦合系数 K
    ↓
散射振幅 M²
    ↓
微分截面 dσ/dt
    ↓
总截面 σ(s)
    ↓
散射率 Γ_ij
    ↓
平均散射率 Γ_i = Σ_j ρ_j × w_ij
    ↓
弛豫时间 τ_i = 1/Γ_i
"""

using Printf

# 添加项目路径到 LOAD_PATH
project_root = joinpath(@__DIR__, "..", "..")
src_path = joinpath(project_root, "src")
if !(src_path in LOAD_PATH)
    push!(LOAD_PATH, src_path)
end

# 导入必要的模块
include(joinpath(src_path, "relaxtime", "RelaxationTime.jl"))
include(joinpath(src_path, "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(src_path, "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(src_path, "QuarkDistribution.jl"))
include(joinpath(src_path, "Constants_PNJL.jl"))

using .RelaxationTime
using .EffectiveCouplings
using .OneLoopIntegrals
using .QuarkDistribution
using .Constants_PNJL

println("="^80)
println("Fortran vs Julia 弛豫时间对比")
println("="^80)

# 使用与 Fortran 相同的参数
T = 0.15  # fm⁻¹
μ_B = 0.4  # fm⁻¹
μ_u = μ_B / 3.0
μ_s = μ_B / 3.0

println("\n温度和化学势:")
println("  T = $T fm⁻¹ ($(T * 197.3) MeV)")
println("  μ_B = $μ_B fm⁻¹ ($(μ_B * 197.3) MeV)")
println("  μ_u = μ_d = $μ_u fm⁻¹")
println("  μ_s = $μ_s fm⁻¹")

# ============================================================================
# 第一步：从 Fortran 输出中提取弛豫时间
# ============================================================================

println("\n" * "="^80)
println("第一步：从 Fortran 输出中提取弛豫时间")
println("="^80)

# 从上面的 Fortran 输出中手动提取
τ_u_fortran = 0.58096796228598324
τ_s_fortran = 0.59203596086110832
τ_ub_fortran = 0.57915958660547406
τ_sb_fortran = 0.59022691868275623

Γ_u_fortran = 1.7212653105090621
Γ_s_fortran = 1.6890865861349258
Γ_ub_fortran = 1.7266398124584688
Γ_sb_fortran = 1.6942636269991858

println("\nFortran 弛豫时间:")
println("  τ_u  = $τ_u_fortran fm")
println("  τ_s  = $τ_s_fortran fm")
println("  τ_ub = $τ_ub_fortran fm")
println("  τ_sb = $τ_sb_fortran fm")

println("\nFortran 散射率 (Γ = 1/τ):")
println("  Γ_u  = $Γ_u_fortran fm⁻¹")
println("  Γ_s  = $Γ_s_fortran fm⁻¹")
println("  Γ_ub = $Γ_ub_fortran fm⁻¹")
println("  Γ_sb = $Γ_sb_fortran fm⁻¹")

# ============================================================================
# 第二步：使用 Julia 计算弛豫时间
# ============================================================================

println("\n" * "="^80)
println("第二步：使用 Julia 计算弛豫时间")
println("="^80)

# 从 Fortran 调试输出中提取的能隙方程解
Φ_fortran = 0.8404
m_u_fortran = 0.04051  # fm⁻¹
m_s_fortran = 1.0323   # fm⁻¹

println("\n使用 Fortran 的能隙方程解:")
println("  Φ = $Φ_fortran")
println("  m_u = $m_u_fortran fm⁻¹")
println("  m_s = $m_s_fortran fm⁻¹")

# 构造参数
quark_params = (
    m = (u = m_u_fortran, d = m_u_fortran, s = m_s_fortran),
    μ = (u = μ_u, d = μ_u, s = μ_s)
)

thermo_params = (
    T = T,
    Φ = Φ_fortran,
    Φbar = Φ_fortran,  # 假设 Φ = Φbar
    ξ = 0.0
)

# 计算 A 函数
println("\n计算 A 函数...")
nodes_p, weights_p = RelaxationTime.AverageScatteringRate.gauleg(0.0, 20.0, 16)
A_u = OneLoopIntegrals.A(m_u_fortran, μ_u, T, Φ_fortran, Φ_fortran, nodes_p, weights_p)
A_s = OneLoopIntegrals.A(m_s_fortran, μ_s, T, Φ_fortran, Φ_fortran, nodes_p, weights_p)

println("  A_u = $A_u fm⁻²")
println("  A_s = $A_s fm⁻²")

# 添加 A 到参数中
quark_params = merge(quark_params, (A = (u = A_u, d = A_u, s = A_s),))

# 计算有效耦合系数
println("\n计算有效耦合系数...")
K_coeffs = EffectiveCouplings.calculate_effective_couplings(quark_params)

println("  K_π = $(K_coeffs.K_pi) fm⁵")
println("  K_K = $(K_coeffs.K_K) fm⁵")
println("  K_η = $(K_coeffs.K_eta) fm⁵")

# 计算数密度
println("\n计算数密度...")
n_u = QuarkDistribution.number_density(m_u_fortran, μ_u, T, Φ_fortran, Φ_fortran)
n_s = QuarkDistribution.number_density(m_s_fortran, μ_s, T, Φ_fortran, Φ_fortran)

densities = (
    u = n_u,
    d = n_u,
    s = n_s,
    ubar = n_u,  # 假设对称
    dbar = n_u,
    sbar = n_s
)

println("  n_u = $n_u fm⁻³")
println("  n_s = $n_s fm⁻³")

# 计算弛豫时间
println("\n计算弛豫时间...")
println("(这可能需要几分钟...)")

try
    result = RelaxationTime.relaxation_times(
        quark_params,
        thermo_params,
        K_coeffs;
        densities = densities,
        p_nodes = 16,
        angle_nodes = 16,
        phi_nodes = 16,
        n_sigma_points = 64
    )
    
    println("\nJulia 计算结果:")
    println("  τ_u  = $(result.tau.u) fm")
    println("  τ_s  = $(result.tau.s) fm")
    println("  τ_ub = $(result.tau.ubar) fm")
    println("  τ_sb = $(result.tau.sbar) fm")
    
    println("\n平均散射率 (Γ = 1/τ):")
    println("  Γ_u  = $(result.tau_inv.u) fm⁻¹")
    println("  Γ_s  = $(result.tau_inv.s) fm⁻¹")
    println("  Γ_ub = $(result.tau_inv.ubar) fm⁻¹")
    println("  Γ_sb = $(result.tau_inv.sbar) fm⁻¹")
    
    # 输出部分散射率以供检查
    println("\n部分散射过程的平均散射率:")
    for (process, rate) in pairs(result.rates)
        if rate > 1e-10  # 只显示非零的
            @printf("  %-25s = %.6e fm⁻¹\n", process, rate)
        end
    end
    
    # 如果有 Fortran 数据，进行对比
    println("\n" * "="^80)
    println("第三步：对比 Fortran 和 Julia 的结果")
    println("="^80)
    
    println("\n弛豫时间对比:")
    @printf("  %-10s  %-15s  %-15s  %-10s\n", "夸克", "Fortran", "Julia", "差异")
    @printf("  %-10s  %-15s  %-15s  %-10s\n", "----", "-------", "-----", "----")
    
    diff_u = abs(result.tau.u - τ_u_fortran) / τ_u_fortran * 100
    diff_s = abs(result.tau.s - τ_s_fortran) / τ_s_fortran * 100
    diff_ub = abs(result.tau.ubar - τ_ub_fortran) / τ_ub_fortran * 100
    diff_sb = abs(result.tau.sbar - τ_sb_fortran) / τ_sb_fortran * 100
    
    @printf("  %-10s  %15.6e  %15.6e  %9.2f%%\n", "u", τ_u_fortran, result.tau.u, diff_u)
    @printf("  %-10s  %15.6e  %15.6e  %9.2f%%\n", "s", τ_s_fortran, result.tau.s, diff_s)
    @printf("  %-10s  %15.6e  %15.6e  %9.2f%%\n", "ubar", τ_ub_fortran, result.tau.ubar, diff_ub)
    @printf("  %-10s  %15.6e  %15.6e  %9.2f%%\n", "sbar", τ_sb_fortran, result.tau.sbar, diff_sb)
    
    println("\n散射率对比 (Γ = 1/τ):")
    @printf("  %-10s  %-15s  %-15s  %-10s\n", "夸克", "Fortran", "Julia", "差异")
    @printf("  %-10s  %-15s  %-15s  %-10s\n", "----", "-------", "-----", "----")
    
    diff_Γ_u = abs(result.tau_inv.u - Γ_u_fortran) / Γ_u_fortran * 100
    diff_Γ_s = abs(result.tau_inv.s - Γ_s_fortran) / Γ_s_fortran * 100
    diff_Γ_ub = abs(result.tau_inv.ubar - Γ_ub_fortran) / Γ_ub_fortran * 100
    diff_Γ_sb = abs(result.tau_inv.sbar - Γ_sb_fortran) / Γ_sb_fortran * 100
    
    @printf("  %-10s  %15.6e  %15.6e  %9.2f%%\n", "u", Γ_u_fortran, result.tau_inv.u, diff_Γ_u)
    @printf("  %-10s  %15.6e  %15.6e  %9.2f%%\n", "s", Γ_s_fortran, result.tau_inv.s, diff_Γ_s)
    @printf("  %-10s  %15.6e  %15.6e  %9.2f%%\n", "ubar", Γ_ub_fortran, result.tau_inv.ubar, diff_Γ_ub)
    @printf("  %-10s  %15.6e  %15.6e  %9.2f%%\n", "sbar", Γ_sb_fortran, result.tau_inv.sbar, diff_Γ_sb)
    
    # 判断结果
    max_diff = max(diff_u, diff_s, diff_ub, diff_sb)
    println("\n" * "="^80)
    if max_diff < 1.0
        println("✅ 结果一致! 最大差异 = $(round(max_diff, digits=2))% < 1%")
        println("之前报告的 3 倍差异可能是:")
        println("  - 使用了不同的参数 (Φ, m, T, μ)")
        println("  - 单位换算错误")
        println("  - 数据提取错误")
    elseif max_diff < 10.0
        println("⚠️  有小差异: 最大差异 = $(round(max_diff, digits=2))% (1% - 10%)")
        println("这可能是由于:")
        println("  - 数值积分精度不同")
        println("  - 截断误差")
        println("  - 舍入误差累积")
    else
        println("❌ 有显著差异: 最大差异 = $(round(max_diff, digits=2))% > 10%")
        println("需要逐步检查计算链:")
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

println("\n" * "="^80)
println("分析完成")
println("="^80)

println("\n下一步:")
println("1. 在 Fortran 代码中添加弛豫时间的调试输出")
println("2. 重新编译并运行 Fortran 代码")
println("3. 重新运行此脚本进行对比")
println("4. 如果弛豫时间仍有差异，逐步对比计算链的每个环节")
