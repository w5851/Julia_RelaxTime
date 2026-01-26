#!/usr/bin/env julia
"""
步骤 1: 对比 Fortran 和 Julia 在相同 s 值下的 σ(s)
"""

using Printf

# 包含必要的模块
include("../../src/Constants_PNJL.jl")
using .Constants_PNJL

include("../../src/relaxtime/TotalCrossSection.jl")
using .TotalCrossSection

include("../../src/relaxtime/EffectiveCouplings.jl")
using .EffectiveCouplings

include("../../src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegrals

println("="^80)
println("步骤 1: 对比 σ(s) 的计算")
println("="^80)

# 参数
const ħc = 0.1973269804  # GeV·fm
T_MeV = 300.0
μ_MeV = 0.67
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)

Φ = 0.8404
m_u = 0.04051  # fm⁻¹
m_s = 1.0323   # fm⁻¹

println("\n参数:")
@printf("  T = %.4f fm⁻¹ (%.1f MeV)\n", T, T_MeV)
@printf("  μ = %.6f fm⁻¹ (%.2f MeV)\n", μ, μ_MeV)
@printf("  Φ = %.4f\n", Φ)
@printf("  m_u = %.5f fm⁻¹\n", m_u)
@printf("  m_s = %.5f fm⁻¹\n", m_s)

# 准备参数
quark_params = (
    m = (u = m_u, d = m_u, s = m_s),
    μ = (u = μ, d = μ, s = μ)
)

thermo_params = (
    T = T,
    Φ = Φ,
    Φbar = Φ,
    ξ = 0.0
)

# 计算 A 和 K
include("../../src/integration/GaussLegendre.jl")
using .GaussLegendre: gauleg

nodes_p, weights_p = gauleg(0.0, 20.0, 16)
A_u = OneLoopIntegrals.A(m_u, μ, T, Φ, Φ, nodes_p, weights_p)
A_s = OneLoopIntegrals.A(m_s, μ, T, Φ, Φ, nodes_p, weights_p)

G_u = EffectiveCouplings.calculate_G_from_A(A_u, m_u)
G_s = EffectiveCouplings.calculate_G_from_A(A_s, m_s)

K_coeffs = EffectiveCouplings.calculate_effective_couplings(
    Constants_PNJL.G_fm2,
    Constants_PNJL.K_fm5,
    G_u,
    G_s
)

quark_params_with_A = merge(quark_params, (A = (u = A_u, d = A_u, s = A_s),))

println("\n" * "="^80)
println("Fortran 输出的 σ(s) 值")
println("="^80)

# Fortran 输出的数据
fortran_data = [
    (1, 4.2667907376421894, 0.57078927045445838),
    (13, 7.3853387287831733, 4.7682136329242494E-002),
    (25, 10.503886719924161, 5.2929204092785384E-002),
    (37, 13.622434711065145, 5.8581038999614317E-002),
    (49, 16.740982702206132, 6.3742718970422246E-002),
    (61, 19.859530693347118, 6.8466726042647461E-002),
    (73, 22.978078684488104, 7.2901592560502382E-002),
    (85, 26.096626675629086, 7.7252956389068220E-002),
    (97, 29.215174666770075, 8.1590667972460118E-002),
    (109, 32.333722657911053, 8.6384351238286181E-002),
    (121, 35.452270649052039, 9.1913025488137470E-002),
]

println("\n| i | s (fm⁻²) | σ_Fortran (fm²) |")
println("|---|----------|-----------------|")
for (i, s, σ) in fortran_data
    @printf("| %3d | %8.4f | %.6e |\n", i, s, σ)
end

println("\n" * "="^80)
println("Julia 计算相同 s 值的 σ(s)")
println("="^80)

println("\n| i | s (fm⁻²) | σ_Julia (fm²) | σ_Fortran (fm²) | 差异 |")
println("|---|----------|---------------|-----------------|------|")

results = []
for (i, s, σ_fortran) in fortran_data
    σ_julia = total_cross_section(
        :ssbar_to_uubar,
        s,
        quark_params_with_A,
        thermo_params,
        K_coeffs;
        n_points = 32  # 使用足够多的点
    )
    
    diff = abs(σ_julia - σ_fortran) / σ_fortran * 100
    push!(results, (i, s, σ_julia, σ_fortran, diff))
    
    @printf("| %3d | %8.4f | %.6e | %.6e | %.2f%% |\n", 
            i, s, σ_julia, σ_fortran, diff)
end

max_diff = maximum(r[5] for r in results)

println("\n" * "="^80)
println("结论")
println("="^80)

@printf("\n最大差异: %.2f%%\n", max_diff)

if max_diff < 1.0
    println("\n✅ σ(s) 计算一致! (差异 < 1%)")
    println("\n结论: 散射截面的计算不是问题的根源")
    println("下一步: 检查插值方法的影响")
elseif max_diff < 5.0
    println("\n⚠️  σ(s) 有小差异 (1-5%)")
    println("\n可能的原因:")
    println("  1. t 积分的节点数不同")
    println("  2. 数值积分方法的细微差异")
    println("\n下一步: 增加 t 积分节点数,看是否能减小差异")
else
    println("\n❌ σ(s) 有显著差异! (> 5%)")
    println("\n结论: 问题在散射截面的计算")
    println("下一步: 检查微分截面 dσ/dt 的计算")
end

println("\n" * "="^80)
println("额外分析: 阈值附近的行为")
println("="^80)

s_th = (2 * m_s)^2
println("\n阈值 s_th = $(s_th) fm⁻²")

# 检查第一个点 (最接近阈值)
s1, σ1_fortran = fortran_data[1][2], fortran_data[1][3]
@printf("\n第一个点: s = %.4f fm⁻² (s/s_th = %.4f)\n", s1, s1/s_th)
@printf("  σ_Fortran = %.6e fm²\n", σ1_fortran)

# 这个点的 σ 异常大,可能是阈值奇点
if σ1_fortran > 0.1
    println("\n⚠️  注意: 第一个点的 σ 异常大!")
    println("  这可能是因为太接近阈值,导致数值不稳定")
    println("  建议: 检查 Fortran 的 s_bo 设置")
end

println("\n测试完成!")
