#!/usr/bin/env julia
"""
对比Julia和Fortran的单个散射率 w̄_ij

目标：找出绝对值差异的根源
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxationTime.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))

using .Constants_PNJL: Λ_inv_fm, ħc_MeV_fm, G_fm2, K_fm5
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .RelaxationTime: compute_average_rates
using Printf

println("="^80)
println("Julia vs Fortran 散射率详细对比")
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

# 计算Julia的散射率
println("计算Julia散射率...")
@time rates_julia = compute_average_rates(
    quark_params,
    thermo_params,
    K_coeffs;
    p_nodes=20,
    angle_nodes=4,
    phi_nodes=8,
    sigma_cutoff=Λ_inv_fm
)

println()
println("="^80)
println("散射率对比")
println("="^80)
println()

# Fortran的散射率（从w_ij.dat）
# 格式：rho_B/rho_0, mu, T, w(1:11)
# w(1)=ud->ud, w(2)=uu->uu, w(3)=us->us, w(4)=ss->ss,
# w(5)=udbar->udbar, w(6)=uubar->uubar, w(7)=uubar->ddbar,
# w(8)=usbar->usbar, w(9)=uubar->ssbar, w(10)=ssbar->uubar, w(11)=ssbar->ssbar

fortran_w = [
    0.5776136603e-01,  # w(1) = ud->ud
    0.5759193845e-01,  # w(2) = uu->uu
    0.5899150396e-01,  # w(3) = us->us
    0.5334502064e-01,  # w(4) = ss->ss
    0.2406399323e+00,  # w(5) = udbar->udbar
    0.2970841963e+00,  # w(6) = uubar->uubar
    0.5892755992e-01,  # w(7) = uubar->ddbar
    0.2291186797e+00,  # w(8) = usbar->usbar
    0.4692409224e-01,  # w(9) = uubar->ssbar
    0.6074054988e-01,  # w(10) = ssbar->uubar (注意：这是ssbar->uubar，不是uubar->ssbar)
    0.2955235064e+00,  # w(11) = ssbar->ssbar
]

# Fortran的反夸克散射率（从w_ija.dat）
# wa(1)=ubdb->ubdb, wa(2)=ubub->ubub, wa(3)=ubsb->ubsb, wa(4)=sbsb->sbsb,
# wa(5)=dubar->dubar, wa(6)=subar->subar

fortran_wa = [
    0.5786655731e-01,  # wa(1) = ubdb->ubdb
    0.5769676643e-01,  # wa(2) = ubub->ubub
    0.5909173088e-01,  # wa(3) = ubsb->ubsb
    0.5342869800e-01,  # wa(4) = sbsb->sbsb
    0.2406399323e+00,  # wa(5) = dubar->dubar
    0.2291622652e+00,  # wa(6) = subar->subar
]

# 映射到Julia的命名
process_map = [
    (:ud_to_ud, fortran_w[1], "ud→ud (qq)"),
    (:uu_to_uu, fortran_w[2], "uu→uu (qq)"),
    (:us_to_us, fortran_w[3], "us→us (qq)"),
    (:ss_to_ss, fortran_w[4], "ss→ss (qq)"),
    (:udbar_to_udbar, fortran_w[5], "uđ̄→uđ̄ (qq̄)"),
    (:uubar_to_uubar, fortran_w[6], "uū→uū (qq̄)"),
    (:uubar_to_ddbar, fortran_w[7], "uū→đđ̄ (qq̄)"),
    (:usbar_to_usbar, fortran_w[8], "us̄→us̄ (qq̄)"),
    (:uubar_to_ssbar, fortran_w[9], "uū→ss̄ (qq̄)"),
    (:ssbar_to_uubar, fortran_w[10], "ss̄→uū (qq̄)"),
    (:ssbar_to_ssbar, fortran_w[11], "ss̄→ss̄ (qq̄)"),
    (:ubardbar_to_ubardbar, fortran_wa[1], "ūđ̄→ūđ̄ (q̄q̄)"),
    (:ubarubar_to_ubarubar, fortran_wa[2], "ūū→ūū (q̄q̄)"),
    (:ubarsbar_to_ubarsbar, fortran_wa[3], "ūs̄→ūs̄ (q̄q̄)"),
    (:sbarsbar_to_sbarsbar, fortran_wa[4], "s̄s̄→s̄s̄ (q̄q̄)"),
    (:dubar_to_dubar, fortran_wa[5], "đū→đū (qq̄)"),
    (:subar_to_subar, fortran_wa[6], "sū→sū (qq̄)"),
]

println("单个散射过程对比 (fm⁻³):")
println("-"^80)
@printf("%-25s %15s %15s %12s\n", "过程", "Julia", "Fortran", "Julia/Fortran")
println("-"^80)

ratios = Float64[]
for (proc, fortran_val, desc) in process_map
    julia_val = rates_julia[proc]
    ratio = julia_val / fortran_val
    push!(ratios, ratio)
    @printf("%-25s %15.6e %15.6e %12.3f\n", desc, julia_val, fortran_val, ratio)
end

println("-"^80)
@printf("%-25s %15s %15s %12.3f\n", "平均比值", "", "", sum(ratios)/length(ratios))
@printf("%-25s %15s %15s %12.3f\n", "最小比值", "", "", minimum(ratios))
@printf("%-25s %15s %15s %12.3f\n", "最大比值", "", "", maximum(ratios))
println()

println("="^80)
println("关键发现")
println("="^80)
println()

avg_ratio = sum(ratios) / length(ratios)
@printf("1. 平均比值: Julia/Fortran = %.3f\n", avg_ratio)
println()

if avg_ratio < 0.3
    println("2. Julia的散射率比Fortran小约3-5倍")
    println("   这解释了为什么弛豫时间大约3倍")
    println("   因为 τ = ρ / (Σ ρ_j w̄_ij)")
    println("   如果 w̄_ij 小3倍，则 τ 大3倍")
    println()
elseif avg_ratio > 3.0
    println("2. Julia的散射率比Fortran大约3-5倍")
    println("   这与弛豫时间的差异方向不一致")
    println("   需要检查公式")
    println()
else
    println("2. Julia和Fortran的散射率在同一数量级")
    println("   差异可能来自其他因素")
    println()
end

println("3. 可能的原因:")
println("   a) 归一化约定不同")
println("      - Fortran可能有额外的因子（如2π, ħ等）")
println("      - 需要检查averaged_rate()的定义")
println()
println("   b) 积分方法不同")
println("      - Julia: 半无穷积分 [0, ∞)")
println("      - Fortran: 有限积分 [0, Λ]")
println()
println("   c) 散射截面σ(s)的计算不同")
println("      - 可能使用不同的近似或公式")
println()

println("="^80)
println("下一步")
println("="^80)
println()

println("需要检查:")
println("  1. Fortran的averaged_rate()函数定义")
println("  2. 是否有额外的归一化因子")
println("  3. 积分方法的具体实现")
println("  4. 散射截面σ(s)的计算公式")
println()

println("="^80)
