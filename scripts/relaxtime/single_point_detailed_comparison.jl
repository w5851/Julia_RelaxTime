"""
单点详细对比：μ=2 MeV, T=300 MeV

与Fortran版本精确对比，输出所有中间计算细节
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm, ρ0_inv_fm3
using .TransportWorkflow: solve_gap_and_transport, FixedMu
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using Printf

println("="^80)
println("Julia vs Fortran 单点详细对比")
println("="^80)

# 精确匹配Fortran的参数
T_MeV = 300.0
μ_MeV = 2.0  # Fortran使用2 MeV（代码中的最小值保护）

T_fm = T_MeV / ħc_MeV_fm
μ_fm = μ_MeV / ħc_MeV_fm

println("\n输入参数:")
println("-"^80)
@printf("温度 T = %.1f MeV = %.6f fm⁻¹\n", T_MeV, T_fm)
@printf("化学势 μ = %.1f MeV = %.6f fm⁻¹\n", μ_MeV, μ_fm)
@printf("各向异性参数 ξ = %.1f\n", 0.0)

# 使用高温初值（与Fortran一致）
seed_state = [-0.001, -0.001, -0.04, 0.8, 0.8]

println("\n" * "="^80)
println("步骤 1: 求解PNJL平衡态")
println("="^80)

# 求解平衡态
base = TransportWorkflow.PNJL.solve(
    TransportWorkflow.PNJL.FixedMu(),
    T_fm,
    μ_fm;
    xi=0.0,
    p_num=12,
    t_num=6,
    seed_strategy=TransportWorkflow.PNJL.DefaultSeed(seed_state, seed_state, :quark),
    iterations=40,
)

println("\n平衡态求解结果:")
println("-"^80)
@printf("收敛状态: %s\n", base.converged ? "✓ 收敛" : "✗ 未收敛")
@printf("迭代次数: %d\n", base.iterations)
@printf("残差范数: %.6e\n", base.residual_norm)

# 提取Polyakov loop和有效质量
Φ = Float64(base.x_state[4])
Φbar = Float64(base.x_state[5])
masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

println("\nPolyakov loop:")
@printf("  Φ     = %.6f\n", Φ)
@printf("  Φ̄     = %.6f\n", Φbar)

println("\n有效夸克质量 (fm⁻¹):")
@printf("  m_u   = %.6f fm⁻¹ = %.2f MeV\n", masses.u, masses.u * ħc_MeV_fm)
@printf("  m_d   = %.6f fm⁻¹ = %.2f MeV\n", masses.d, masses.d * ħc_MeV_fm)
@printf("  m_s   = %.6f fm⁻¹ = %.2f MeV\n", masses.s, masses.s * ħc_MeV_fm)

println("\n热力学量:")
@printf("  Ω     = %.6e fm⁻⁴ = %.6e MeV/fm³\n", base.omega, base.omega * ħc_MeV_fm)
@printf("  P     = %.6e fm⁻⁴ = %.6e MeV/fm³\n", base.pressure, base.pressure * ħc_MeV_fm)
@printf("  ε     = %.6e fm⁻⁴ = %.6e MeV/fm³\n", base.energy, base.energy * ħc_MeV_fm)
@printf("  s     = %.6e fm⁻³\n", base.entropy)

println("\n" * "="^80)
println("步骤 2: 计算有效耦合常数 K")
println("="^80)

# 计算K系数
nodes = DEFAULT_MOMENTUM_NODES
weights = DEFAULT_MOMENTUM_WEIGHTS
A_u = A(masses.u, μ_fm, T_fm, Φ, Φbar, nodes, weights)
A_s = A(masses.s, μ_fm, T_fm, Φ, Φbar, nodes, weights)
G_u = calculate_G_from_A(A_u, masses.u)
G_s = calculate_G_from_A(A_s, masses.s)
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

println("\n单圈积分 A:")
@printf("  A_u   = %.6e\n", A_u)
@printf("  A_s   = %.6e\n", A_s)

println("\n有效耦合 G:")
@printf("  G_u   = %.6e fm²\n", G_u)
@printf("  G_s   = %.6e fm²\n", G_s)

println("\n有效耦合常数 K:")
for (key, val) in pairs(K_coeffs)
    @printf("  %-15s = %.6e fm⁵\n", string(key), val)
end

println("\n" * "="^80)
println("步骤 3: 计算弛豫时间")
println("="^80)

# 计算弛豫时间（使用与Fortran相同的积分节点数）
result = solve_gap_and_transport(
    T_fm,
    μ_fm;
    xi=0.0,
    equilibrium=base,
    compute_tau=true,
    K_coeffs=K_coeffs,
    p_num=12,
    t_num=6,
    solver_kwargs=(iterations=40,),
    tau_kwargs=(
        p_nodes=20,      # 动量积分节点
        angle_nodes=4,   # 角度积分节点
        phi_nodes=8,     # φ积分节点
        n_sigma_points=64,  # σ(s)积分点数
    ),
)

eq = result.equilibrium
dens = result.densities
tau = result.tau
tau_inv = result.tau_inv

println("\n粒子数密度 (fm⁻³):")
@printf("  n_u    = %.6e\n", dens.u)
@printf("  n_d    = %.6e\n", dens.d)
@printf("  n_s    = %.6e\n", dens.s)
@printf("  n_ūbar = %.6e\n", dens.ubar)
@printf("  n_d̄bar = %.6e\n", dens.dbar)
@printf("  n_s̄bar = %.6e\n", dens.sbar)

# 计算重子数密度
rho_quark_net = (dens.u - dens.ubar) + (dens.d - dens.dbar) + (dens.s - dens.sbar)
rho_baryon = rho_quark_net / 3.0
rho_norm = rho_baryon / ρ0_inv_fm3

println("\n重子数密度:")
@printf("  ρ_B    = %.6e fm⁻³\n", rho_baryon)
@printf("  ρ_B/ρ₀ = %.6f\n", rho_norm)

println("\n弛豫时间 τ (fm):")
println("-"^80)
@printf("  τ_u    = %.6f fm\n", tau.u)
@printf("  τ_d    = %.6f fm\n", tau.d)
@printf("  τ_s    = %.6f fm\n", tau.s)
@printf("  τ_ūbar = %.6f fm\n", tau.ubar)
@printf("  τ_d̄bar = %.6f fm\n", tau.dbar)
@printf("  τ_s̄bar = %.6f fm\n", tau.sbar)

println("\n散射率 τ⁻¹ (fm⁻¹):")
println("-"^80)
@printf("  τ_u⁻¹  = %.6e fm⁻¹\n", tau_inv.u)
@printf("  τ_d⁻¹  = %.6e fm⁻¹\n", tau_inv.d)
@printf("  τ_s⁻¹  = %.6e fm⁻¹\n", tau_inv.s)
@printf("  τ_ū⁻¹  = %.6e fm⁻¹\n", tau_inv.ubar)
@printf("  τ_d̄⁻¹  = %.6e fm⁻¹\n", tau_inv.dbar)
@printf("  τ_s̄⁻¹  = %.6e fm⁻¹\n", tau_inv.sbar)

println("\n弛豫时间比值:")
println("-"^80)
@printf("  τ_u/τ_s = %.6f\n", tau.u / tau.s)
@printf("  τ_ū/τ_s̄ = %.6f\n", tau.ubar / tau.sbar)

println("\n" * "="^80)
println("步骤 4: 计算输运系数")
println("="^80)

tr = result.transport

println("\n输运系数:")
@printf("  η      = %.6e fm⁻³\n", tr.eta)
@printf("  σ      = %.6e fm⁻³\n", tr.sigma)
@printf("  ζ      = %.6e fm⁻³\n", tr.zeta)

println("\n无量纲比值:")
@printf("  η/s    = %.6f\n", tr.eta / eq.entropy)
@printf("  ζ/s    = %.6f\n", tr.zeta / eq.entropy)

println("\n" * "="^80)
println("与Fortran对比")
println("="^80)

println("\nFortran结果 (T=300 MeV, μ=2 MeV):")
println("  τ_u    = 0.5839 fm")
println("  τ_d    = 0.5839 fm")
println("  τ_s    = 0.5932 fm")
println("  τ_ū    = 0.5820 fm")
println("  τ_d̄    = 0.5820 fm")
println("  τ_s̄    = 0.5913 fm")
println("  τ_u/τ_s = 0.9843")

println("\nJulia结果 (T=300 MeV, μ=2 MeV):")
@printf("  τ_u    = %.4f fm\n", tau.u)
@printf("  τ_d    = %.4f fm\n", tau.d)
@printf("  τ_s    = %.4f fm\n", tau.s)
@printf("  τ_ū    = %.4f fm\n", tau.ubar)
@printf("  τ_d̄    = %.4f fm\n", tau.dbar)
@printf("  τ_s̄    = %.4f fm\n", tau.sbar)
@printf("  τ_u/τ_s = %.4f\n", tau.u / tau.s)

println("\n差异分析:")
fortran_tau_u = 0.5839
fortran_tau_s = 0.5932
julia_tau_u = tau.u
julia_tau_s = tau.s

diff_u = abs(julia_tau_u - fortran_tau_u) / fortran_tau_u * 100
diff_s = abs(julia_tau_s - fortran_tau_s) / fortran_tau_s * 100
diff_ratio = abs((julia_tau_u/julia_tau_s) - (fortran_tau_u/fortran_tau_s)) / (fortran_tau_u/fortran_tau_s) * 100

@printf("  Δτ_u   = %.1f%%\n", diff_u)
@printf("  Δτ_s   = %.1f%%\n", diff_s)
@printf("  Δ(τ_u/τ_s) = %.1f%%\n", diff_ratio)

println("\n" * "="^80)
println("可能的差异来源:")
println("="^80)
println("1. 积分精度设置（节点数、截断值）")
println("2. 散射截面计算方法的细节差异")
println("3. 数值积分算法（Gauss-Legendre vs 其他）")
println("4. 有效质量和Polyakov loop的微小差异")
println("5. 单圈积分A的计算方法")
println("="^80)
