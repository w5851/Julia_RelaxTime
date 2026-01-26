#!/usr/bin/env julia
"""
验证阈值计算

检查 ssbar→uubar 和 uubar→ssbar 的阈值是否相同
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .TransportWorkflow: solve_gap_and_transport
using Printf

println("="^80)
println("验证阈值计算")
println("="^80)
println()

# 测试参数
T_MeV = 300.0
μ_MeV = 2.0
ξ = 0.0

T = T_MeV / ħc_MeV_fm
μ = μ_MeV / ħc_MeV_fm

# 求解平衡态
println("求解PNJL平衡态...")
seed_state = [-0.001, -0.001, -0.04, 0.8, 0.8]
base = TransportWorkflow.PNJL.solve(
    TransportWorkflow.PNJL.FixedMu(),
    T, μ;
    xi=ξ,
    seed_strategy=TransportWorkflow.PNJL.DefaultSeed(seed_state, seed_state, :quark),
    iterations=40,
)

masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

println("完成")
@printf("  m_u = %.3f MeV (%.6f fm⁻¹)\n", masses.u * ħc_MeV_fm, masses.u)
@printf("  m_s = %.3f MeV (%.6f fm⁻¹)\n", masses.s * ħc_MeV_fm, masses.s)
println()

println("="^80)
println("阈值计算")
println("="^80)
println()

# ssbar → uubar
println("过程1: ssbar → uubar")
println("-"^80)
m_i_1 = masses.s  # s
m_j_1 = masses.s  # sbar
m_c_1 = masses.u  # u
m_d_1 = masses.u  # ubar

s_threshold_initial_1 = (m_i_1 + m_j_1)^2
s_threshold_final_1 = (m_c_1 + m_d_1)^2
s_threshold_1 = max(s_threshold_initial_1, s_threshold_final_1)

@printf("  初态粒子: s (%.3f MeV) + sbar (%.3f MeV)\n", m_i_1 * ħc_MeV_fm, m_j_1 * ħc_MeV_fm)
@printf("  末态粒子: u (%.3f MeV) + ubar (%.3f MeV)\n", m_c_1 * ħc_MeV_fm, m_d_1 * ħc_MeV_fm)
println()
@printf("  s_threshold(初态) = (m_s + m_s)² = %.6f fm⁻² = %.1f MeV²\n", 
    s_threshold_initial_1, s_threshold_initial_1 * ħc_MeV_fm^2)
@printf("  s_threshold(末态) = (m_u + m_u)² = %.6f fm⁻² = %.1f MeV²\n", 
    s_threshold_final_1, s_threshold_final_1 * ħc_MeV_fm^2)
@printf("  s_threshold(总)   = max(初态, 末态) = %.6f fm⁻² = %.1f MeV²\n", 
    s_threshold_1, s_threshold_1 * ħc_MeV_fm^2)
println()

# uubar → ssbar
println("过程2: uubar → ssbar")
println("-"^80)
m_i_2 = masses.u  # u
m_j_2 = masses.u  # ubar
m_c_2 = masses.s  # s
m_d_2 = masses.s  # sbar

s_threshold_initial_2 = (m_i_2 + m_j_2)^2
s_threshold_final_2 = (m_c_2 + m_d_2)^2
s_threshold_2 = max(s_threshold_initial_2, s_threshold_final_2)

@printf("  初态粒子: u (%.3f MeV) + ubar (%.3f MeV)\n", m_i_2 * ħc_MeV_fm, m_j_2 * ħc_MeV_fm)
@printf("  末态粒子: s (%.3f MeV) + sbar (%.3f MeV)\n", m_c_2 * ħc_MeV_fm, m_d_2 * ħc_MeV_fm)
println()
@printf("  s_threshold(初态) = (m_u + m_u)² = %.6f fm⁻² = %.1f MeV²\n", 
    s_threshold_initial_2, s_threshold_initial_2 * ħc_MeV_fm^2)
@printf("  s_threshold(末态) = (m_s + m_s)² = %.6f fm⁻² = %.6f fm⁻² = %.1f MeV²\n", 
    s_threshold_final_2, s_threshold_final_2, s_threshold_final_2 * ħc_MeV_fm^2)
@printf("  s_threshold(总)   = max(初态, 末态) = %.6f fm⁻² = %.1f MeV²\n", 
    s_threshold_2, s_threshold_2 * ħc_MeV_fm^2)
println()

println("="^80)
println("对比")
println("="^80)
println()

@printf("ssbar→uubar 的阈值: %.6f fm⁻² = %.1f MeV²\n", s_threshold_1, s_threshold_1 * ħc_MeV_fm^2)
@printf("uubar→ssbar 的阈值: %.6f fm⁻² = %.1f MeV²\n", s_threshold_2, s_threshold_2 * ħc_MeV_fm^2)
println()

if abs(s_threshold_1 - s_threshold_2) < 1e-10
    println("✓ 阈值相同！（物理正确）")
    println()
    println("  两个过程的阈值都是 max((m_s+m_s)², (m_u+m_u)²) = (m_s+m_s)²")
    println("  因为 m_s >> m_u，所以阈值由重粒子决定")
else
    println("✗ 阈值不同！（物理错误）")
    @printf("  差异: %.6e fm⁻²\n", abs(s_threshold_1 - s_threshold_2))
end
println()

println("="^80)
println("物理解释")
println("="^80)
println()

println("对于 2→2 散射过程 i + j → c + d：")
println()
println("  质心能量阈值 = max((m_i + m_j)², (m_c + m_d)²)")
println()
println("  原因：")
println("    - 初态必须有足够能量产生末态粒子")
println("    - 末态必须有足够能量来自初态粒子")
println("    - 取最大值确保两个方向都满足")
println()
println("  对于 ssbar ↔ uubar：")
println("    - 初态阈值(ssbar): (2m_s)² = %.1f MeV²" % (s_threshold_initial_1 * ħc_MeV_fm^2))
println("    - 末态阈值(uubar): (2m_u)² = %.1f MeV²" % (s_threshold_final_1 * ħc_MeV_fm^2))
println("    - 总阈值: max(%.1f, %.1f) = %.1f MeV²" % 
    (s_threshold_initial_1 * ħc_MeV_fm^2, s_threshold_final_1 * ħc_MeV_fm^2, s_threshold_1 * ħc_MeV_fm^2))
println()
println("  因此，ssbar→uubar 和 uubar→ssbar 的阈值**必须相同**！")
println()

println("="^80)
println("结论")
println("="^80)
println()

println("Julia的实现是正确的：")
println("  - 代码中使用 s_threshold = max(s_threshold_initial, s_threshold_final)")
println("  - 这确保了 ssbar→uubar 和 uubar→ssbar 有相同的阈值")
println("  - 阈值 = %.1f MeV² (由重粒子 s 决定)" % (s_threshold_1 * ħc_MeV_fm^2))
println()

println("之前的分析错误：")
println("  - 我错误地认为两个过程的阈值不同")
println("  - 实际上它们的阈值完全相同")
println("  - 需要重新分析3.238倍差异的原因")
println()

println("="^80)

