#!/usr/bin/env julia
"""
使用正确的参数重新分析弛豫时间差异

现在我们知道:
1. A 函数实现一致 (差异 < 0.1%)
2. 能隙方程解基本一致 (Φ 差异 < 0.3%)

需要检查:
1. G^f 的计算是否一致
2. 有效耦合系数 K 的计算是否一致
3. 弛豫时间的计算是否一致
"""

using Printf

# 包含必要的模块
include("../../src/pnjl/PNJL.jl")
using .PNJL

include("../../src/Constants_PNJL.jl")
using .Constants_PNJL

include("../../src/pnjl/workflows/TransportWorkflow.jl")
using .TransportWorkflow

include("../../src/relaxtime/EffectiveCouplings.jl")
using .EffectiveCouplings

println("=" ^ 80)
println("重新分析弛豫时间差异 (使用正确的参数)")
println("=" ^ 80)
println()

# 测试参数 (与 Fortran 一致)
const ħc = 0.1973269804  # GeV·fm
T_MeV = 300.0
μ_B_MeV = 6.0
μ_MeV = μ_B_MeV / 3.0
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)

println("测试条件:")
@printf("  T = %.6f fm⁻¹ (%.1f MeV)\n", T, T_MeV)
@printf("  μ_B = %.6f fm⁻¹ (%.1f MeV)\n", 3*μ, μ_B_MeV)
@printf("  μ_q = %.6f fm⁻¹ (%.1f MeV)\n", μ, μ_MeV)
println()

println("=" ^ 80)
println("步骤 1: 求解能隙方程")
println("=" ^ 80)
println()

try
    # 求解能隙方程
    base = PNJL.solve(
        PNJL.FixedMu(),
        T,
        μ;
        xi=0.0,
        p_num=64,
        t_num=16
    )
    
    println("能隙方程求解成功!")
    @printf("  收敛: %s\n", base.converged)
    @printf("  迭代次数: %d\n", base.iterations)
    println()
    
    # 提取参数
    x_state = base.x_state
    φ_u = x_state[1]
    φ_s = x_state[3]
    Φ = x_state[4]
    Φbar = x_state[5]
    
    masses = base.masses
    m_u = masses[1]
    m_s = masses[3]
    
    @printf("序参量:\n")
    @printf("  φ_u = %.6e fm⁻³\n", φ_u)
    @printf("  φ_s = %.6e fm⁻³\n", φ_s)
    @printf("  Φ = %.8f\n", Φ)
    @printf("  Φbar = %.8f\n", Φbar)
    println()
    
    @printf("有效质量:\n")
    @printf("  m_u = %.5f fm⁻¹ (%.2f MeV)\n", m_u, m_u * ħc * 1000.0)
    @printf("  m_s = %.5f fm⁻¹ (%.2f MeV)\n", m_s, m_s * ħc * 1000.0)
    println()
    
    println("=" ^ 80)
    println("步骤 2: 计算 A 函数")
    println("=" ^ 80)
    println()
    
    # 手动计算 A 函数
    include("../../src/relaxtime/OneLoopIntegrals.jl")
    using .OneLoopIntegrals
    include("../../src/integration/GaussLegendre.jl")
    using .GaussLegendre
    
    nodes_p, weights_p = gauleg(0.0, 20.0, 64)
    A_u = A(m_u, μ, T, Φ, Φbar, nodes_p, weights_p)
    A_s = A(m_s, μ, T, Φ, Φbar, nodes_p, weights_p)
    
    @printf("A 函数:\n")
    @printf("  A_u = %.6f fm⁻²\n", A_u)
    @printf("  A_s = %.6f fm⁻²\n", A_s)
    println()
    
    # 对比 Fortran
    A_u_fortran = -4.9721620988593660
    A_s_fortran = -5.1655414336119714
    
    @printf("对比 Fortran:\n")
    @printf("  ΔA_u = %.6f fm⁻² (%.2f%%)\n", 
            A_u - A_u_fortran,
            100.0 * abs(A_u - A_u_fortran) / abs(A_u_fortran))
    @printf("  ΔA_s = %.6f fm⁻² (%.2f%%)\n",
            A_s - A_s_fortran,
            100.0 * abs(A_s - A_s_fortran) / abs(A_s_fortran))
    println()
    
    println("=" ^ 80)
    println("步骤 3: 计算 G^f")
    println("=" ^ 80)
    println()
    
    # 计算 G^f (使用 Julia 的公式)
    # G^f = -3 / (4π²) × m_f × A_f
    G_u_julia = -3.0 / (4.0 * π^2) * m_u * A_u
    G_s_julia = -3.0 / (4.0 * π^2) * m_s * A_s
    
    @printf("G^f (Julia):\n")
    @printf("  G_u = %.6e fm⁻³\n", G_u_julia)
    @printf("  G_s = %.6e fm⁻³\n", G_s_julia)
    println()
    
    # 对比 Fortran
    # 从 Fortran 调试输出: arrG(u) = 0.000714 fm²
    # 如果 arrG = K_f × G^f, 则:
    # G^f (Fortran) = arrG / K_f
    
    K_f = K_fm5
    arrG_u_fortran = 0.000714  # fm²
    arrG_s_fortran = 0.01891   # fm²
    G_u_fortran_implied = arrG_u_fortran / K_f
    G_s_fortran_implied = arrG_s_fortran / K_f
    
    @printf("G^f (从 Fortran arrG 推算,假设 arrG = K_f × G^f):\n")
    @printf("  arrG_u = %.6e fm²\n", arrG_u_fortran)
    @printf("  arrG_s = %.6e fm²\n", arrG_s_fortran)
    @printf("  K_f = %.6e fm⁵\n", K_f)
    @printf("  G_u = arrG_u / K_f = %.6e fm⁻³\n", G_u_fortran_implied)
    @printf("  G_s = arrG_s / K_f = %.6e fm⁻³\n", G_s_fortran_implied)
    println()
    
    @printf("差异:\n")
    @printf("  G_u (Julia) / G_u (Fortran) = %.2f\n", G_u_julia / G_u_fortran_implied)
    @printf("  G_s (Julia) / G_s (Fortran) = %.2f\n", G_s_julia / G_s_fortran_implied)
    println()
    
    println("=" ^ 80)
    println("步骤 4: 检查 Fortran 的 arrG 计算")
    println("=" ^ 80)
    println()
    
    # Fortran 的 arrG 计算公式:
    # arrG(:) = - arrMass(:) * arrA(:) * (3d0*K_f) / (4d0*pi**2)
    
    # 使用 Fortran 的参数重新计算
    m_u_fortran = 0.04051
    A_u_fortran_val = -4.9721620988593660
    K_f_fortran = 0.046657314097986152
    
    arrG_u_recalc = - m_u_fortran * A_u_fortran_val * (3.0 * K_f_fortran) / (4.0 * π^2)
    
    @printf("重新计算 Fortran 的 arrG:\n")
    @printf("  m_u = %.5f fm⁻¹\n", m_u_fortran)
    @printf("  A_u = %.10f fm⁻²\n", A_u_fortran_val)
    @printf("  K_f = %.10f fm⁵\n", K_f_fortran)
    @printf("  arrG_u = -m_u × A_u × (3K_f)/(4π²)\n")
    @printf("         = %.6e fm²\n", arrG_u_recalc)
    println()
    
    @printf("对比 Fortran 输出:\n")
    @printf("  arrG_u (输出) = %.6e fm²\n", arrG_u_fortran)
    @printf("  arrG_u (重算) = %.6e fm²\n", arrG_u_recalc)
    @printf("  差异: %.2f%%\n", 100.0 * abs(arrG_u_recalc - arrG_u_fortran) / abs(arrG_u_fortran))
    println()
    
    println("=" ^ 80)
    println("步骤 5: 理解 G^f 的定义")
    println("=" ^ 80)
    println()
    
    println("Julia 的定义:")
    println("  G^f = -3 / (4π²) × m_f × A_f")
    println("  单位: fm⁻³")
    println()
    
    println("Fortran 的计算:")
    println("  arrG = -m × A × (3K_f) / (4π²)")
    println("       = [-3 / (4π²) × m × A] × K_f")
    println("       = G^f × K_f")
    println("  单位: fm⁻³ × fm⁵ = fm²")
    println()
    
    println("结论:")
    println("  Fortran 的 arrG 确实等于 K_f × G^f")
    println("  这是一个预先计算的组合量,用于后续计算")
    println()
    
    println("=" ^ 80)
    println("分析")
    println("=" ^ 80)
    println()
    
    println("关键发现:")
    println()
    
    if abs(A_u - A_u_fortran) / abs(A_u_fortran) < 0.02
        println("✅ A 函数一致 (差异 < 2%)")
        println("   这证明 A 函数的实现没有问题")
        println()
    end
    
    if abs(arrG_u_recalc - arrG_u_fortran) / abs(arrG_u_fortran) < 0.02
        println("✅ Fortran 的 arrG 计算公式验证正确")
        println("   arrG = -m × A × (3K_f) / (4π²) = G^f × K_f")
        println()
    end
    
    if abs(G_u_julia / G_u_fortran_implied - 1.0) > 0.1
        println("❌ G^f 有显著差异!")
        println()
        println("可能的原因:")
        println("1. A 函数虽然数值接近,但可能有微小的归一化差异")
        println("2. 质量 m 的微小差异被放大")
        println("3. 需要检查 G^f 在后续计算中的使用方式")
        println()
        
        # 详细分析
        println("详细分析:")
        @printf("  Julia: G_u = -3/(4π²) × %.5f × %.6f = %.6e fm⁻³\n", 
                m_u, A_u, G_u_julia)
        @printf("  Fortran: G_u = arrG_u / K_f = %.6e / %.6e = %.6e fm⁻³\n",
                arrG_u_fortran, K_f, G_u_fortran_implied)
        println()
        
        # 检查是否是质量差异导致的
        G_u_julia_with_fortran_m = -3.0 / (4.0 * π^2) * m_u_fortran * A_u
        @printf("  如果使用 Fortran 的质量:\n")
        @printf("    G_u = -3/(4π²) × %.5f × %.6f = %.6e fm⁻³\n",
                m_u_fortran, A_u, G_u_julia_with_fortran_m)
        @printf("    差异: %.2f%%\n",
                100.0 * abs(G_u_julia_with_fortran_m - G_u_fortran_implied) / abs(G_u_fortran_implied))
        println()
    else
        println("✅ G^f 基本一致")
        println()
    end
    
    println("=" ^ 80)
    println("下一步")
    println("=" ^ 80)
    println()
    
    println("需要:")
    println("1. 检查 Fortran 代码中 arrG 如何在传播子中使用")
    println("2. 检查 Julia 代码中 G^f 如何在有效耦合系数中使用")
    println("3. 对比有效耦合系数 K 的计算")
    println("4. 对比散射振幅的计算")
    println("5. 最终对比弛豫时间")
    println()
    
catch e
    println("❌ 计算失败!")
    println("错误信息: ", e)
    println()
    showerror(stdout, e, catch_backtrace())
    println()
end
