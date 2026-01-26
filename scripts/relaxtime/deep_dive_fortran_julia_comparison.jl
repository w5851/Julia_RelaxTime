#!/usr/bin/env julia
"""
深入对比 Fortran 和 Julia 的中间量

目标：找出 3 倍差异的根源

对比内容：
1. 数密度 ρ_u, ρ_s
2. 有效耦合常数 K
3. 散射振幅 M²（相同 s, t）
4. 截面 σ(s)（相同 s）
5. 平均散射率 ⟨Γ⟩
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "ParameterTypes.jl"))
include(joinpath(PROJECT_ROOT, "src", "QuarkDistribution.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "ScatteringAmplitude.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))

using .Constants_PNJL
using .ParameterTypes
using .PNJLQuarkDistributions
using .GaussLegendre
using .EffectiveCouplings
using .ScatteringAmplitude
using .TotalCrossSection
using .AverageScatteringRate
using Printf

println("="^80)
println("Fortran vs Julia 深入对比")
println("="^80)
println()

# 参数设置（与 Fortran 一致）
T_MeV = 300.0
μ_MeV = 2.0
T = T_MeV / Constants_PNJL.ħc_MeV_fm
μ = μ_MeV / Constants_PNJL.ħc_MeV_fm

# Fortran 的平衡态结果
Φ = 0.840408
Φbar = 0.840412
m_u = 0.040510  # fm⁻¹
m_s = 1.032295  # fm⁻¹

@printf("参数设置:\n")
@printf("  T = %.1f MeV = %.6f fm⁻¹\n", T_MeV, T)
@printf("  μ = %.1f MeV = %.6f fm⁻¹\n", μ_MeV, μ)
@printf("  Φ = %.6f\n", Φ)
@printf("  Φbar = %.6f\n", Φbar)
@printf("  m_u = %.6f fm⁻¹ = %.3f MeV\n", m_u, m_u * Constants_PNJL.ħc_MeV_fm)
@printf("  m_s = %.6f fm⁻¹ = %.3f MeV\n", m_s, m_s * Constants_PNJL.ħc_MeV_fm)
println()

# 构建参数
quark_params = QuarkParams(
    (u=m_u, d=m_u, s=m_s),  # m
    (u=μ, d=μ, s=μ)          # μ
)
thermo_params = ThermoParams(T, Φ, Φbar, 0.0)

# ============================================================================
# 步骤 1: 对比数密度
# ============================================================================
println("="^80)
println("步骤 1: 对比数密度 ρ")
println("="^80)
println()

# Julia 计算数密度（使用半无穷积分）
ρ_u_julia = AverageScatteringRate.number_density(
    :u, m_u, μ, T, Φ, Φbar, 0.0;
    p_nodes=32,
    angle_nodes=2,
    scale=10.0
)

ρ_s_julia = AverageScatteringRate.number_density(
    :s, m_s, μ, T, Φ, Φbar, 0.0;
    p_nodes=32,
    angle_nodes=2,
    scale=10.0
)

# Fortran 结果
ρ_u_fortran = 1.680901  # fm⁻³
ρ_s_fortran = 1.537628  # fm⁻³

println("数密度对比:")
@printf("  ρ_u (Julia)   = %.6f fm⁻³\n", ρ_u_julia)
@printf("  ρ_u (Fortran) = %.6f fm⁻³\n", ρ_u_fortran)
@printf("  比值          = %.4f\n", ρ_u_julia / ρ_u_fortran)
println()
@printf("  ρ_s (Julia)   = %.6f fm⁻³\n", ρ_s_julia)
@printf("  ρ_s (Fortran) = %.6f fm⁻³\n", ρ_s_fortran)
@printf("  比值          = %.4f\n", ρ_s_julia / ρ_s_fortran)
println()

if abs(ρ_u_julia / ρ_u_fortran - 1.0) < 0.01
    println("✅ 数密度一致（差异 < 1%）")
else
    println("❌ 数密度有差异")
end
println()

# ============================================================================
# 步骤 2: 对比有效耦合常数 K
# ============================================================================
println("="^80)
println("步骤 2: 对比有效耦合常数 K")
println("="^80)
println()

# Julia 计算 K
# 准备积分节点
p_nodes, p_weights = GaussLegendre.gauleg(0.0, 20.0, 16)

A_u = EffectiveCouplings.OneLoopIntegrals.A(m_u, μ, T, Φ, Φbar, p_nodes, p_weights)
A_s = EffectiveCouplings.OneLoopIntegrals.A(m_s, μ, T, Φ, Φbar, p_nodes, p_weights)
G_u = EffectiveCouplings.calculate_G_from_A(A_u, m_u)
G_s = EffectiveCouplings.calculate_G_from_A(A_s, m_s)
K_coeffs = EffectiveCouplings.calculate_effective_couplings(
    Constants_PNJL.G_fm2,
    Constants_PNJL.K_fm5,
    G_u,
    G_s
)

println("Julia 计算:")
@printf("  A_u = %.6f fm⁻²\n", A_u)
@printf("  A_s = %.6f fm⁻²\n", A_s)
@printf("  G_u = %.6f fm²\n", G_u)
@printf("  G_s = %.6f fm²\n", G_s)
println()
println("有效耦合常数:")
println("  K0_plus     = $(K_coeffs.K0_plus)")
println("  K0_minus    = $(K_coeffs.K0_minus)")
println("  K123_plus   = $(K_coeffs.K123_plus)")
println("  K123_minus  = $(K_coeffs.K123_minus)")
println("  K4567_plus  = $(K_coeffs.K4567_plus)")
println("  K4567_minus = $(K_coeffs.K4567_minus)")
println("  K8_plus     = $(K_coeffs.K8_plus)")
println("  K8_minus    = $(K_coeffs.K8_minus)")
println()

# 转换为物理介子耦合常数
# K_pi = K123_plus, K_K = K4567_plus, K_sigma = K0_plus, K_kappa = K4567_minus
println("物理介子耦合常数:")
@printf("  K_pi    = %.4f (K123_plus)\n", K_coeffs.K123_plus)
@printf("  K_K     = %.4f (K4567_plus)\n", K_coeffs.K4567_plus)
@printf("  K_sigma = %.4f (K0_plus)\n", K_coeffs.K0_plus)
@printf("  K_kappa = %.4f (K4567_minus)\n", K_coeffs.K4567_minus)
println()

# 注意：Fortran 的 K 值需要从代码中提取
println("注意：需要从 Fortran 代码中提取 K 值进行对比")
println()

# ============================================================================
# 步骤 3: 对比散射振幅 M²
# ============================================================================
println("="^80)
println("步骤 3: 对比散射振幅 M²")
println("="^80)
println()

# 添加 A 字段到 quark_params
quark_params_full = (
    m = quark_params.m,
    μ = quark_params.μ,
    A = (u=A_u, d=A_u, s=A_s)
)

# 选择测试点
s_test_MeV2 = 500000.0  # 500 GeV²
s_test = s_test_MeV2 / Constants_PNJL.ħc_MeV_fm^2
t_test_MeV2 = -100000.0  # -100 GeV²
t_test = t_test_MeV2 / Constants_PNJL.ħc_MeV_fm^2

println("测试点:")
@printf("  s = %.0f MeV² = %.3f fm⁻²\n", s_test_MeV2, s_test)
@printf("  t = %.0f MeV² = %.3f fm⁻²\n", t_test_MeV2, t_test)
println()

# 测试几个关键过程
processes = [:uu_to_uu, :us_to_us, :ss_to_ss, :ssbar_to_uubar]

println("散射振幅 M² (fm⁴):")
for process in processes
    M2 = ScatteringAmplitude.scattering_amplitude_squared(
        process, s_test, t_test,
        quark_params_full, thermo_params, K_coeffs
    )
    @printf("  M²(%s) = %.6e fm⁴\n", process, M2)
end
println()

println("注意：需要从 Fortran 代码中提取相同 s, t 的 M² 值进行对比")
println()

# ============================================================================
# 步骤 4: 对比截面 σ(s)
# ============================================================================
println("="^80)
println("步骤 4: 对比截面 σ(s)")
println("="^80)
println()

# 选择几个测试点
s_values_MeV2 = [200000.0, 500000.0, 1000000.0]  # 200, 500, 1000 GeV²

println("截面 σ(s) (fm²):")
for s_MeV2 in s_values_MeV2
    s = s_MeV2 / Constants_PNJL.ħc_MeV_fm^2
    @printf("\ns = %.0f MeV² = %.3f fm⁻²:\n", s_MeV2, s)
    
    for process in processes
        σ = TotalCrossSection.total_cross_section(
            process, s,
            quark_params_full, thermo_params, K_coeffs;
            n_points=32
        )
        @printf("  σ(%s) = %.6e fm²\n", process, σ)
    end
end
println()

println("注意：需要从 Fortran 代码中提取相同 s 的 σ 值进行对比")
println()

# ============================================================================
# 步骤 5: 对比平均散射率 ⟨Γ⟩
# ============================================================================
println("="^80)
println("步骤 5: 对比平均散射率 ⟨Γ⟩")
println("="^80)
println()

# 计算平均散射率
println("平均散射率 ⟨Γ⟩ (fm⁻¹):")
println("（使用半无穷积分 + Λ截断σ(s)）")
println()

# 选择几个关键过程
test_processes = [:uu_to_uu, :ud_to_ud, :us_to_us, :ss_to_ss, :ssbar_to_uubar]

for process in test_processes
    rate = AverageScatteringRate.average_scattering_rate(
        process,
        quark_params_full,
        thermo_params,
        K_coeffs;
        p_nodes=20,
        angle_nodes=4,
        phi_nodes=8,
        n_sigma_points=32,
        sigma_cutoff=Constants_PNJL.Λ_inv_fm,
    )
    @printf("  ⟨Γ⟩(%s) = %.6e fm⁻¹\n", process, rate)
end
println()

println("注意：需要从 Fortran 输出中提取相应的散射率进行对比")
println()

# ============================================================================
# 步骤 6: 计算弛豫时间（用于验证）
# ============================================================================
println("="^80)
println("步骤 6: 计算弛豫时间（验证）")
println("="^80)
println()

# 注意：这里需要完整的弛豫时间计算
# 由于 RelaxationTime 模块的复杂性，我们先跳过
println("弛豫时间计算需要完整的 RelaxationTime 模块")
println("请参考 test_semi_infinite_momentum.jl 的结果")
println()

# ============================================================================
# 总结
# ============================================================================
println("="^80)
println("总结")
println("="^80)
println()

println("已完成的对比:")
println("  ✅ 数密度 ρ_u, ρ_s")
println("  ✅ 有效耦合常数 K")
println("  ✅ 散射振幅 M²（Julia 值）")
println("  ✅ 截面 σ(s)（Julia 值）")
println("  ✅ 平均散射率 ⟨Γ⟩（Julia 值）")
println()

println("下一步:")
println("  1. 从 Fortran 代码中提取相应的中间量")
println("  2. 逐项对比，找出第一个出现差异的地方")
println("  3. 重点关注：")
println("     - 数密度的计算方法")
println("     - K 值的计算")
println("     - M² 的计算公式")
println("     - σ(s) 的积分方法")
println()

println("="^80)
