#!/usr/bin/env julia
"""
在相同s值下详细对比 ssbar→uubar 和 uubar→ssbar

目标：理解为什么在相同阈值下，两个过程的截面差异巨大
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: Λ_inv_fm, ħc_MeV_fm, G_fm2, K_fm5
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .TotalCrossSection: total_cross_section, calculate_t_bounds
using Printf

println("="^80)
println("在相同s值下详细对比 ssbar→uubar 和 uubar→ssbar")
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

Φ = Float64(base.x_state[4])
Φbar = Float64(base.x_state[5])
masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

println("完成")
@printf("  m_u = %.3f MeV\n", masses.u * ħc_MeV_fm)
@printf("  m_s = %.3f MeV\n", masses.s * ħc_MeV_fm)
println()

# 计算K系数
nodes = DEFAULT_MOMENTUM_NODES
weights = DEFAULT_MOMENTUM_WEIGHTS
A_u = A(masses.u, μ, T, Φ, Φbar, nodes, weights)
A_s = A(masses.s, μ, T, Φ, Φbar, nodes, weights)
G_u = calculate_G_from_A(A_u, masses.u)
G_s = calculate_G_from_A(A_s, masses.s)
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

quark_params = (
    m = masses,
    μ = (u=μ, d=μ, s=μ),
    A = (u=A_u, d=A_u, s=A_s)
)

thermo_params = (
    T = T,
    Φ = Φ,
    Φbar = Φbar,
    ξ = ξ
)

# 计算阈值
m_u = masses.u
m_s = masses.s

s_threshold = max((2*m_s)^2, (2*m_u)^2)

@printf("阈值: s_threshold = %.6f fm⁻² = %.1f MeV²\n", s_threshold, s_threshold * ħc_MeV_fm^2)
println()

# 测试不同的s值
s_factors = [1.01, 1.05, 1.1, 1.2, 1.5, 2.0]

println("="^80)
println("在不同s值下的对比")
println("="^80)
println()

for factor in s_factors
    s = s_threshold * factor
    
    println("-"^80)
    @printf("s = %.2f × s_threshold = %.1f MeV²\n", factor, s * ħc_MeV_fm^2)
    println("-"^80)
    println()
    
    # ssbar → uubar
    println("过程1: ssbar → uubar")
    mi_1, mj_1 = m_s, m_s
    mc_1, md_1 = m_u, m_u
    
    # 计算质心动量
    p_cm_in_1 = sqrt(max(0.0, (s - (mi_1+mj_1)^2) * (s - (mi_1-mj_1)^2))) / (2*sqrt(s))
    p_cm_out_1 = sqrt(max(0.0, (s - (mc_1+md_1)^2) * (s - (mc_1-md_1)^2))) / (2*sqrt(s))
    
    @printf("  初态: s + sbar, m_i = m_j = %.3f MeV\n", mi_1 * ħc_MeV_fm)
    @printf("  末态: u + ubar, m_c = m_d = %.3f MeV\n", mc_1 * ħc_MeV_fm)
    @printf("  p_cm(初态) = %.6f fm⁻¹ = %.3f MeV\n", p_cm_in_1, p_cm_in_1 * ħc_MeV_fm)
    @printf("  p_cm(末态) = %.6f fm⁻¹ = %.3f MeV\n", p_cm_out_1, p_cm_out_1 * ħc_MeV_fm)
    
    # 计算t积分范围
    t_bounds_1 = calculate_t_bounds(s, mi_1, mj_1, mc_1, md_1)
    @printf("  t_min = %.6f fm⁻², t_max = %.6f fm⁻²\n", t_bounds_1.t_min, t_bounds_1.t_max)
    @printf("  Δt = %.6f fm⁻²\n", t_bounds_1.t_max - t_bounds_1.t_min)
    
    # 计算截面
    σ_1 = total_cross_section(:ssbar_to_uubar, s, quark_params, thermo_params, K_coeffs)
    @printf("  σ(s) = %.6e fm²\n", σ_1)
    println()
    
    # uubar → ssbar
    println("过程2: uubar → ssbar")
    mi_2, mj_2 = m_u, m_u
    mc_2, md_2 = m_s, m_s
    
    # 计算质心动量
    p_cm_in_2 = sqrt(max(0.0, (s - (mi_2+mj_2)^2) * (s - (mi_2-mj_2)^2))) / (2*sqrt(s))
    p_cm_out_2 = sqrt(max(0.0, (s - (mc_2+md_2)^2) * (s - (mc_2-md_2)^2))) / (2*sqrt(s))
    
    @printf("  初态: u + ubar, m_i = m_j = %.3f MeV\n", mi_2 * ħc_MeV_fm)
    @printf("  末态: s + sbar, m_c = m_d = %.3f MeV\n", mc_2 * ħc_MeV_fm)
    @printf("  p_cm(初态) = %.6f fm⁻¹ = %.3f MeV\n", p_cm_in_2, p_cm_in_2 * ħc_MeV_fm)
    @printf("  p_cm(末态) = %.6f fm⁻¹ = %.3f MeV\n", p_cm_out_2, p_cm_out_2 * ħc_MeV_fm)
    
    # 计算t积分范围
    t_bounds_2 = calculate_t_bounds(s, mi_2, mj_2, mc_2, md_2)
    @printf("  t_min = %.6f fm⁻², t_max = %.6f fm⁻²\n", t_bounds_2.t_min, t_bounds_2.t_max)
    @printf("  Δt = %.6f fm⁻²\n", t_bounds_2.t_max - t_bounds_2.t_min)
    
    # 计算截面
    σ_2 = total_cross_section(:uubar_to_ssbar, s, quark_params, thermo_params, K_coeffs)
    @printf("  σ(s) = %.6e fm²\n", σ_2)
    println()
    
    # 对比
    println("对比:")
    @printf("  p_cm_in(ssbar) / p_cm_in(uubar) = %.3f\n", p_cm_in_1 / p_cm_in_2)
    @printf("  p_cm_out(ssbar) / p_cm_out(uubar) = %.3f\n", p_cm_out_1 / p_cm_out_2)
    @printf("  Δt(ssbar) / Δt(uubar) = %.3f\n", 
        (t_bounds_1.t_max - t_bounds_1.t_min) / (t_bounds_2.t_max - t_bounds_2.t_min))
    if σ_2 > 0.0
        @printf("  σ(ssbar→uubar) / σ(uubar→ssbar) = %.3f\n", σ_1 / σ_2)
    else
        @printf("  σ(ssbar→uubar) / σ(uubar→ssbar) = Inf (σ_2 ≈ 0)\n")
    end
    println()
    
    # 分析
    println("分析:")
    println("  - 截面 ∝ 1/p_cm_in²")
    @printf("  - 预期比值(仅考虑p_cm): (p_cm_in_2/p_cm_in_1)² = %.3f\n", (p_cm_in_2/p_cm_in_1)^2)
    if σ_2 > 0.0
        @printf("  - 实际比值: %.3f\n", σ_1 / σ_2)
        @printf("  - 差异因子: %.3f (可能来自M²或其他因素)\n", (σ_1/σ_2) / (p_cm_in_2/p_cm_in_1)^2)
    end
    println()
end

println("="^80)
println("总结")
println("="^80)
println()

println("关键发现:")
println("  1. 阈值相同：两个过程的s_threshold都是 %.1f MeV²" % (s_threshold * ħc_MeV_fm^2))
println("  2. 但质心动量不同：")
println("     - ssbar→uubar: p_cm_in小（接近阈值），p_cm_out大")
println("     - uubar→ssbar: p_cm_in大，p_cm_out小（接近阈值）")
println("  3. 截面 ∝ 1/p_cm_in²，所以：")
println("     - ssbar→uubar: σ大（因为p_cm_in小）")
println("     - uubar→ssbar: σ小（因为p_cm_in大）")
println()

println("="^80)

