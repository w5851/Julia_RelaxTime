using Test
using Logging
using Printf

include("../../src/QuarkDistribution_Aniso.jl")

using .PNJLQuarkDistributions_Aniso
using .PNJLQuarkDistributions_Aniso: PNJLQuarkDistributions_Aniso as Aniso
using .PNJLQuarkDistributions_Aniso.PNJLQuarkDistributions: PNJLQuarkDistributions as Dist

@testset "各向异性分布函数测试" begin
    # 测试参数设置
    p_inv_fm = 0.5      # 动量
    m_inv_fm = 0.3      # 质量
    μ_inv_fm = 0.1      # 化学势
    T_inv_fm = 0.15     # 温度
    Φ = 0.5             # Polyakov loop
    Φbar = 0.5          # Polyakov loop conjugate
    
    println("\n" * "="^80)
    println("测试1: 不同各向异性参数ξ下distribution_aniso与distribution的差距")
    println("="^80)
    println("参数: p=$p_inv_fm fm⁻¹, m=$m_inv_fm fm⁻¹, μ=$μ_inv_fm fm⁻¹, T=$T_inv_fm fm⁻¹")
    println("      Φ=$Φ, Φ̄=$Φbar")
    println("-"^80)
    
    # 测试不同ξ值（包括负值）
    ξ_values = [-1.0, -0.5, -0.2, -0.1, -0.05, -0.01, 0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
    cosθ = 0.5  # 固定角度
    
    println("夸克分布函数 (cosθ=$cosθ):")
    println(@sprintf("%-10s %-20s %-20s %-20s", "ξ", "distribution", "distribution_aniso", "相对差距(%)"))
    println("-"^80)
    
    E_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2)
    f0 = Dist.quark_distribution(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
    
    for ξ in ξ_values
        f_aniso = Aniso.quark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
        rel_diff = abs(f_aniso - f0) / max(abs(f0), 1e-10) * 100
        println(@sprintf("%-10.2f %-20.10e %-20.10e %-20.10e", ξ, f0, f_aniso, rel_diff))
        
        # 当ξ=0时，两者应该完全相等
        if ξ == 0.0
            @test isapprox(f_aniso, f0, atol=1e-12)
        end
    end
    
    println("\n反夸克分布函数 (cosθ=$cosθ):")
    println(@sprintf("%-10s %-20s %-20s %-20s", "ξ", "distribution", "distribution_aniso", "相对差距(%)"))
    println("-"^80)
    
    fbar0 = Dist.antiquark_distribution(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
    
    for ξ in ξ_values
        fbar_aniso = Aniso.antiquark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
        rel_diff = abs(fbar_aniso - fbar0) / max(abs(fbar0), 1e-10) * 100
        println(@sprintf("%-10.2f %-20.10e %-20.10e %-20.10e", ξ, fbar0, fbar_aniso, rel_diff))
        
        # 当ξ=0时，两者应该完全相等
        if ξ == 0.0
            @test isapprox(fbar_aniso, fbar0, atol=1e-12)
        end
    end
    
    println("\n" * "="^80)
    println("测试2: 不同角度cosθ下的各向异性效应")
    println("="^80)
    println("参数: ξ=0.1")
    println("-"^80)
    
    ξ = 0.1
    cosθ_values = [-1.0, -0.5, 0.0, 0.5, 1.0]
    
    println("夸克分布函数:")
    println(@sprintf("%-10s %-20s %-20s %-20s", "cosθ", "distribution", "distribution_aniso", "相对差距(%)"))
    println("-"^80)
    
    for cosθ in cosθ_values
        f_aniso = Aniso.quark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
        rel_diff = abs(f_aniso - f0) / max(abs(f0), 1e-10) * 100
        println(@sprintf("%-10.2f %-20.10e %-20.10e %-20.10e", cosθ, f0, f_aniso, rel_diff))
    end
    
    println("\n" * "="^80)
    println("测试3: 不同物理参数下的各向异性效应")
    println("="^80)
    
    # 测试不同温度
    println("\n不同温度 (ξ=0.1, cosθ=0.5):")
    println(@sprintf("%-15s %-20s %-20s %-20s", "T (fm⁻¹)", "distribution", "distribution_aniso", "相对差距(%)"))
    println("-"^80)
    
    ξ = 0.1
    cosθ = 0.5
    T_values = [0.05, 0.1, 0.15, 0.2, 0.3]
    
    for T in T_values
        E_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2)
        f0_T = Dist.quark_distribution(E_inv_fm, μ_inv_fm, T, Φ, Φbar)
        f_aniso_T = Aniso.quark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T, Φ, Φbar, ξ, cosθ)
        rel_diff = abs(f_aniso_T - f0_T) / max(abs(f0_T), 1e-10) * 100
        println(@sprintf("%-15.3f %-20.10e %-20.10e %-20.10e", T, f0_T, f_aniso_T, rel_diff))
    end
    
    # 测试不同动量
    println("\n不同动量 (ξ=0.1, cosθ=0.5, T=0.15 fm⁻¹):")
    println(@sprintf("%-15s %-20s %-20s %-20s", "p (fm⁻¹)", "distribution", "distribution_aniso", "相对差距(%)"))
    println("-"^80)
    
    p_values = [0.1, 0.3, 0.5, 1.0, 2.0]
    
    for p in p_values
        E_inv_fm = sqrt(p^2 + m_inv_fm^2)
        f0_p = Dist.quark_distribution(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
        f_aniso_p = Aniso.quark_distribution_aniso(p, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
        rel_diff = abs(f_aniso_p - f0_p) / max(abs(f0_p), 1e-10) * 100
        println(@sprintf("%-15.3f %-20.10e %-20.10e %-20.10e", p, f0_p, f_aniso_p, rel_diff))
    end
    
    println("\n" * "="^80)
    println("测试4: distribution_aniso与线性近似的比较")
    println("="^80)
    println("比较 (distribution_aniso - distribution) 与 correction 的差距")
    println("在小|ξ|情况下，两者应该接近")
    println("-"^80)
    
    p_inv_fm = 0.5
    cosθ = 0.5
    ξ_values_all = [-1.0, -0.8, -0.5, -0.2, -0.1, -0.05, -0.01, -0.001, 0.0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1.0]
    
    println("夸克分布函数 (包含大ξ值):")
    println(@sprintf("%-10s %-20s %-20s %-20s %-20s", "ξ", "Δf(完整)", "correction", "比例(Δf/corr)", "相对差距(%)"))
    println("-"^80)
    
    E_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2)
    f0 = Dist.quark_distribution(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
    
    for ξ in ξ_values_all
        # 完整形式的差值：f_aniso - f0
        f_aniso = Aniso.quark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
        delta_f_full = f_aniso - f0
        
        # 线性近似（一阶修正项）
        correction = Aniso.quark_distribution_aniso_correction(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
        
        # 计算比例和相对差距
        ratio = delta_f_full / (abs(correction) > 1e-15 ? correction : 1e-15)
        abs_diff = abs(delta_f_full - correction)
        # 相对差距：||比例|-1|，表示偏离±1的百分比
        rel_diff = abs(abs(ratio) - 1.0) * 100
        
        println(@sprintf("%-10.3f %-20.10e %-20.10e %-20.10e %-20.10e", ξ, delta_f_full, correction, ratio, rel_diff))
    end
    
    println("\n反夸克分布函数 (包含大ξ值):")
    println(@sprintf("%-10s %-20s %-20s %-20s %-20s", "ξ", "Δf̄(完整)", "correction", "比例(Δf̄/corr)", "相对差距(%)"))
    println("-"^80)
    
    fbar0 = Dist.antiquark_distribution(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
    
    for ξ in ξ_values_all
        # 完整形式的差值：f_aniso - f0
        fbar_aniso = Aniso.antiquark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
        delta_fbar_full = fbar_aniso - fbar0
        
        # 线性近似（一阶修正项）
        correction = Aniso.antiquark_distribution_aniso_correction(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
        
        # 计算比例和相对差距
        ratio = delta_fbar_full / (abs(correction) > 1e-15 ? correction : 1e-15)
        abs_diff = abs(delta_fbar_full - correction)
        # 相对差距：||比例|-1|，表示偏离±1的百分比
        rel_diff = abs(abs(ratio) - 1.0) * 100
        
        println(@sprintf("%-10.3f %-20.10e %-20.10e %-20.10e %-20.10e", ξ, delta_fbar_full, correction, ratio, rel_diff))
    end
    
    println("\n" * "="^80)
    println("测试5: 使用统一接口distribution_aniso进行测试")
    println("="^80)
    
    p_inv_fm = 0.5
    ξ = 0.1
    cosθ = 0.5
    
    println("\n测试:quark符号:")
    f_quark = distribution_aniso(:quark, p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    f_quark_direct = Aniso.quark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    @test f_quark == f_quark_direct
    println(@sprintf("distribution_aniso(:quark, ...) = %.10e", f_quark))
    println(@sprintf("quark_distribution_aniso(...)  = %.10e", f_quark_direct))
    
    println("\n测试:antiquark符号:")
    f_antiquark = distribution_aniso(:antiquark, p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    f_antiquark_direct = Aniso.antiquark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    @test f_antiquark == f_antiquark_direct
    println(@sprintf("distribution_aniso(:antiquark, ...) = %.10e", f_antiquark))
    println(@sprintf("antiquark_distribution_aniso(...)  = %.10e", f_antiquark_direct))
    
    println("\n测试:correction统一接口:")
    corr_quark = distribution_aniso_correction(:quark, p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    corr_quark_direct = Aniso.quark_distribution_aniso_correction(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    @test corr_quark == corr_quark_direct
    println(@sprintf("distribution_aniso_correction(:quark, ...) = %.10e", corr_quark))
    println(@sprintf("quark_distribution_aniso_correction(...)  = %.10e", corr_quark_direct))
    
    println("\n" * "="^80)
    println("测试6: 极端情况测试")
    println("="^80)
    
    # cosθ=0时，各向异性效应应该消失
    println("\n测试cosθ=0 (垂直于各向异性方向):")
    ξ = 0.5
    cosθ = 0.0
    E_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2)
    f0 = Dist.quark_distribution(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
    f_aniso = Aniso.quark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    correction = Aniso.quark_distribution_aniso_correction(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    
    println(@sprintf("f0 = %.10e", f0))
    println(@sprintf("f_aniso = %.10e", f_aniso))
    println(@sprintf("correction = %.10e", correction))
    println(@sprintf("差距 = %.10e", abs(f_aniso - f0)))
    @test isapprox(f_aniso, f0, atol=1e-12)
    @test isapprox(correction, 0.0, atol=1e-12)
    
    println("\n" * "="^80)
    println("所有测试完成!")
    println("="^80)
end

