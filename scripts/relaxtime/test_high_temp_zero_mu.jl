"""
快速测试：高温零化学势下 τ_u 和 τ_s 是否相等

在各向同性PNJL模型中，当化学势为0且温度足够高时，
由于手征对称性恢复，u和s夸克的有效质量应该接近，
因此它们的弛豫时间也应该接近。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm
using .TransportWorkflow: solve_gap_and_transport, FixedMu
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using Printf

println("="^70)
println("测试：高温零化学势下 τ_u 和 τ_s 的比较")
println("="^70)

# 测试温度范围：300-400 MeV
T_values_MeV = [300.0, 320.0, 340.0, 360.0, 380.0, 400.0]
μ_MeV = 1.0  # 接近零化学势

println("\n温度(MeV)  μ(MeV)  m_u(fm⁻¹)  m_s(fm⁻¹)  τ_u(fm)  τ_s(fm)  τ_u/τ_s  收敛")
println("-"^70)

for T_MeV in T_values_MeV
    T_fm = T_MeV / ħc_MeV_fm
    μ_fm = μ_MeV / ħc_MeV_fm
    
    try
        # 使用高温初值
        seed_state = [-0.001, -0.001, -0.04, 0.8, 0.8]
        
        # 首先求解平衡态
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
        
        if !base.converged
            println("$T_MeV MeV: 平衡态未收敛")
            continue
        end
        
        # 提取质量和Polyakov loop
        Φ = Float64(base.x_state[4])
        Φbar = Float64(base.x_state[5])
        masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))
        
        # 计算K系数
        nodes = DEFAULT_MOMENTUM_NODES
        weights = DEFAULT_MOMENTUM_WEIGHTS
        A_u = A(masses.u, μ_fm, T_fm, Φ, Φbar, nodes, weights)
        A_s = A(masses.s, μ_fm, T_fm, Φ, Φbar, nodes, weights)
        G_u = calculate_G_from_A(A_u, masses.u)
        G_s = calculate_G_from_A(A_s, masses.s)
        K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
        
        # 现在求解能隙和输运系数
        result = solve_gap_and_transport(
            T_fm,
            μ_fm;
            xi=0.0,  # 各向同性
            equilibrium=base,
            compute_tau=true,
            K_coeffs=K_coeffs,
            p_num=12,
            t_num=6,
            solver_kwargs=(iterations=40,),
            tau_kwargs=(
                p_nodes=20,
                angle_nodes=4,
                phi_nodes=8,
                n_sigma_points=64,
            ),
        )
        
        eq = result.equilibrium
        tau = result.tau
        
        # 计算比值
        ratio = tau.u / tau.s
        
        # 输出结果
        @printf("%8.1f  %6.1f  %9.4f  %9.4f  %7.3f  %7.3f  %7.4f  %s\n",
            T_MeV, μ_MeV, eq.masses[1], eq.masses[3],
            tau.u, tau.s, ratio,
            eq.converged ? "✓" : "✗")
        
    catch e
        println("$T_MeV MeV: 计算失败 - $e")
    end
end

println("\n" * "="^70)
println("说明：")
println("- 在高温极限下，手征对称性恢复，m_u ≈ m_s ≈ 0")
println("- 因此散射截面相似，τ_u/τ_s 应该接近 1.0")
println("- 如果比值显著偏离1，可能表明Julia实现与C++参考有差异")
println("="^70)
