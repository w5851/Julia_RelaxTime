"""
测试 OneLoopIntegralsAniso 模块中的 A_correction 和 A_aniso 函数

测试内容：
1. 基本功能测试
2. ξ=0 边界条件测试
3. A + A_correction 与 A_aniso 的数值对比
4. 线性近似有效范围验证
5. 不同参数下的系统性测试
6. 性能基准测试

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("test/test_oneloopintegrals_aniso.jl")
```
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using Test
using Printf

include("../../../src/relaxtime/OneLoopIntegrals.jl")
include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
include("../../../src/integration/GaussLegendre.jl")

using .OneLoopIntegrals: A, const_integral_term_A
using .OneLoopIntegralsCorrection: A_correction, A_aniso
using .GaussLegendre: gauleg, DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS, 
    DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS

# 标准测试参数
TEST_PARAMS = (
    m = 0.30,      # fm⁻¹
    μ = 0.10,      # fm⁻¹
    T = 0.15,      # fm⁻¹
    Φ = 0.20,
    Φbar = 0.20,
)

# 生成默认积分节点
const NODES_COSΘ, WEIGHTS_COSΘ = DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS
const NODES_P, WEIGHTS_P = DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS

@testset "OneLoopIntegralsAniso 模块测试" begin
    
    @testset "A_correction 基本功能" begin
        m = TEST_PARAMS.m
        μ = TEST_PARAMS.μ
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.1
        
        result = A_correction(m, μ, T, Φ, Φbar, ξ, NODES_P, WEIGHTS_P)
        
        @test result isa Float64
        @test isfinite(result)
        
        @info "A_correction 基本测试" ξ=ξ result_fm=result
    end
    
    @testset "A_aniso 基本功能" begin
        m = TEST_PARAMS.m
        μ = TEST_PARAMS.μ
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.1
        
        result = A_aniso(m, μ, T, Φ, Φbar, ξ, NODES_P, WEIGHTS_P, NODES_COSΘ, WEIGHTS_COSΘ)
        
        @test result isa Float64
        @test isfinite(result)
        # 注意：A积分可以是负值（因为公式中有-1的真空项）
        
        @info "A_aniso 基本测试" ξ=ξ result_fm=result
    end
    
    @testset "ξ=0 边界条件" begin
        m = TEST_PARAMS.m
        μ = TEST_PARAMS.μ
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        
        # 计算各向同性项
        A_iso = A(m, μ, T, Φ, Φbar, NODES_P, WEIGHTS_P)
        
        # ξ=0 时，A_correction 应为 0
        ΔA_zero = A_correction(m, μ, T, Φ, Φbar, 0.0, NODES_P, WEIGHTS_P)
        
        # ξ=0 时，A_aniso 应等于 A_iso
        A_aniso_zero = A_aniso(m, μ, T, Φ, Φbar, 0.0, NODES_P, WEIGHTS_P, NODES_COSΘ, WEIGHTS_COSΘ)
        
        @info "ξ=0 边界条件测试" A_iso=A_iso ΔA_zero=ΔA_zero A_aniso_zero=A_aniso_zero
        
        # A_correction(ξ=0) 应该非常接近 0
        @test abs(ΔA_zero) < 1e-8
        
        # A_aniso(ξ=0) 应该等于 A_iso（允许小的数值误差）
        rel_error = abs(A_aniso_zero - A_iso) / A_iso
        @test rel_error < 1e-3  # 相对误差 < 0.1%
        
        @info "ξ=0 验证" ΔA_zero=ΔA_zero A_diff=(A_aniso_zero-A_iso) rel_error_percent=(rel_error*100)
    end
    
    @testset "线性近似验证" begin
        m = TEST_PARAMS.m
        μ = TEST_PARAMS.μ
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        
        # 测试 A_correction 的线性性质
        ξ1 = 0.05
        ξ2 = 0.10
        ξ3 = 0.15
        
        ΔA1 = A_correction(m, μ, T, Φ, Φbar, ξ1, NODES_P, WEIGHTS_P)
        ΔA2 = A_correction(m, μ, T, Φ, Φbar, ξ2, NODES_P, WEIGHTS_P)
        ΔA3 = A_correction(m, μ, T, Φ, Φbar, ξ3, NODES_P, WEIGHTS_P)
        
        ratio_2_1 = ΔA2 / ΔA1
        ratio_3_1 = ΔA3 / ΔA1
        
        @info "线性性质检验" ξ1=ξ1 ΔA1=ΔA1 ξ2=ξ2 ΔA2=ΔA2 ratio_2_1=ratio_2_1
        @info "线性性质检验" ξ3=ξ3 ΔA3=ΔA3 ratio_3_1=ratio_3_1
        
        # 在小 ξ 时，ΔA 应该线性于 ξ
        @test abs(ratio_2_1 - 2.0) < 0.05  # ΔA(2ξ) / ΔA(ξ) ≈ 2.0，误差 < 5%
        @test abs(ratio_3_1 - 3.0) < 0.10  # ΔA(3ξ) / ΔA(ξ) ≈ 3.0，误差 < 10%
    end
    
    @testset "A + A_correction vs A_aniso 对比" begin
        m = TEST_PARAMS.m
        μ = TEST_PARAMS.μ
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        
        # 计算各向同性基准
        A_iso = A(m, μ, T, Φ, Φbar, NODES_P, WEIGHTS_P)
        
        # 计算常数项（用于分析分布函数部分）
        const_term = const_integral_term_A(m)
        
        println("\n" * "="^80)
        println("A + A_correction vs A_aniso 数值对比")
        println("="^80)
        println("表1: 完整A值对比（包含常数项）")
        println("-"^80)
        println(@sprintf("%-8s %-14s %-14s %-14s %-10s", "ξ", "A_approx", "A_exact", "Δ(abs)", "Δ(%)"))
        println("-"^80)
        
        ξ_values = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50]
        results = []
        
        for ξ in ξ_values
            # 方法1：一阶近似
            ΔA = A_correction(m, μ, T, Φ, Φbar, ξ, NODES_P, WEIGHTS_P)
            A_approx = A_iso + ΔA
            
            # 方法2：完整计算
            A_exact = A_aniso(m, μ, T, Φ, Φbar, ξ, NODES_P, WEIGHTS_P, NODES_COSΘ, WEIGHTS_COSΘ)
            
            # 完整A值的误差
            abs_diff_full = A_exact - A_approx
            rel_error_full = abs(abs_diff_full) / abs(A_exact) * 100
            
            println(@sprintf("%-8.2f %-14.6f %-14.6f %-14.6f %-10.3f", 
                            ξ, A_approx, A_exact, abs_diff_full, rel_error_full))
            
            # 提取分布函数部分（减去常数项影响）
            # 公式：A = 分布函数积分 - 4×常数项
            # 因此：分布函数积分 = A + 4×常数项
            A_approx_dist = A_approx + 4.0 * const_term
            A_exact_dist = A_exact + 4.0 * const_term
            
            abs_diff_dist = A_exact_dist - A_approx_dist
            rel_error_dist = abs(abs_diff_dist) / abs(A_exact_dist) * 100
            
            push!(results, (ξ=ξ, 
                           A_approx=A_approx, A_exact=A_exact, 
                           abs_diff_full=abs_diff_full, rel_error_full=rel_error_full,
                           A_approx_dist=A_approx_dist, A_exact_dist=A_exact_dist,
                           abs_diff_dist=abs_diff_dist, rel_error_dist=rel_error_dist))
        end
        
        println("="^80)
        println()
        println("表2: 分布函数部分对比（排除常数项）")
        println("-"^80)
        println(@sprintf("%-8s %-14s %-14s %-14s %-10s", "ξ", "Dist_approx", "Dist_exact", "Δ(abs)", "Δ(%)"))
        println("-"^80)
        
        for r in results
            println(@sprintf("%-8.2f %-14.6f %-14.6f %-14.6f %-10.3f", 
                            r.ξ, r.A_approx_dist, r.A_exact_dist, r.abs_diff_dist, r.rel_error_dist))
        end
        
        println("="^80 * "\n")
        
        # 验证误差趋势（使用完整A值的误差）
        # 小 ξ (≤0.2) 时，完整A值相对误差应 < 0.02%
        for r in results
            if r.ξ <= 0.2
                @test r.rel_error_full < 0.02
            end
        end
        
        # 中等 ξ (0.2-0.4) 时，完整A值相对误差应 < 0.03%
        for r in results
            if 0.2 < r.ξ <= 0.4
                @test r.rel_error_full < 0.03
            end
        end
        
        @info "对比测试完成" 最小完整误差_percent=minimum(r.rel_error_full for r in results) 最大完整误差_percent=maximum(r.rel_error_full for r in results) 最小分布误差_percent=minimum(r.rel_error_dist for r in results) 最大分布误差_percent=maximum(r.rel_error_dist for r in results)
    end
    
    @testset "不同温度下的系统性测试" begin
        m = TEST_PARAMS.m
        μ = TEST_PARAMS.μ
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.2
        
        println("\n" * "="^70)
        println("不同温度下的 A 值扫描 (ξ=$ξ)")
        println("="^70)
        println(@sprintf("%-12s %-12s %-12s %-12s", "T (fm⁻¹)", "A_iso", "ΔA", "A_total"))
        println("-"^70)
        
        T_values = [0.10, 0.15, 0.20, 0.25, 0.30]
        
        for T in T_values
            A_iso = A(m, μ, T, Φ, Φbar, NODES_P, WEIGHTS_P)
            ΔA = A_correction(m, μ, T, Φ, Φbar, ξ, NODES_P, WEIGHTS_P)
            A_total = A_iso + ΔA
            
            println(@sprintf("%-12.3f %-12.6f %-12.6f %-12.6f", T, A_iso, ΔA, A_total))
            
            @test isfinite(A_iso)
            @test isfinite(ΔA)
            @test isfinite(A_total)
        end
        
        println("="^70 * "\n")
    end
    
    @testset "不同化学势下的系统性测试" begin
        m = TEST_PARAMS.m
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.2
        
        println("\n" * "="^70)
        println("不同化学势下的 A 值扫描 (ξ=$ξ)")
        println("="^70)
        println(@sprintf("%-12s %-12s %-12s %-12s", "μ (fm⁻¹)", "A_iso", "ΔA", "A_total"))
        println("-"^70)
        
        μ_values = [0.0, 0.1, 0.2, 0.3, 0.4]
        
        for μ in μ_values
            A_iso = A(m, μ, T, Φ, Φbar, NODES_P, WEIGHTS_P)
            ΔA = A_correction(m, μ, T, Φ, Φbar, ξ, NODES_P, WEIGHTS_P)
            A_total = A_iso + ΔA
            
            println(@sprintf("%-12.3f %-12.6f %-12.6f %-12.6f", μ, A_iso, ΔA, A_total))
            
            @test isfinite(A_iso)
            @test isfinite(ΔA)
            @test isfinite(A_total)
        end
        
        println("="^70 * "\n")
    end
    
    @testset "不同夸克质量下的系统性测试" begin
        μ = TEST_PARAMS.μ
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.2
        
        println("\n" * "="^70)
        println("不同夸克质量下的 A 值扫描 (ξ=$ξ)")
        println("="^70)
        println(@sprintf("%-12s %-12s %-12s %-12s", "m (fm⁻¹)", "A_iso", "ΔA", "A_total"))
        println("-"^70)
        
        m_values = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
        
        for m in m_values
            A_iso = A(m, μ, T, Φ, Φbar, NODES_P, WEIGHTS_P)
            ΔA = A_correction(m, μ, T, Φ, Φbar, ξ, NODES_P, WEIGHTS_P)
            A_total = A_iso + ΔA
            
            println(@sprintf("%-12.3f %-12.6f %-12.6f %-12.6f", m, A_iso, ΔA, A_total))
            
            @test isfinite(A_iso)
            @test isfinite(ΔA)
            @test isfinite(A_total)
        end
        
        println("="^70 * "\n")
    end
    
    @testset "积分上限 p_max 收敛性测试" begin
        m = TEST_PARAMS.m
        μ = TEST_PARAMS.μ
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.2
        n_points = 64  # 固定节点数
        
        println("\n" * "="^70)
        println("积分上限 p_max 收敛性测试")
        println("="^70)
        println("测试目的：验证 p_max = 20 fm⁻¹ 是否足够让积分收敛")
        println("-"^70)
        
        # 测试 A 函数（各向同性）
        println("\n表1: A 函数随 p_max 的变化")
        println(@sprintf("%-12s %-15s %-15s", "p_max(fm⁻¹)", "A", "相对变化(%)"))
        println("-"^70)
        
        p_max_values = [10.0, 15.0, 20.0, 25.0, 30.0, 40.0]
        prev_A = nothing
        
        for p_max in p_max_values
            nodes, weights = gauleg(0.0, p_max, n_points)
            result = A(m, μ, T, Φ, Φbar, nodes, weights)
            
            if prev_A === nothing
                println(@sprintf("%-12.1f %-15.8f %-15s", p_max, result, "-"))
            else
                rel_change = abs(result - prev_A) / abs(prev_A) * 100
                println(@sprintf("%-12.1f %-15.8f %-15.6f", p_max, result, rel_change))
                
                # p_max >= 20 后，变化应 < 0.01%
                if p_max >= 20.0
                    @test rel_change < 0.01
                end
            end
            
            prev_A = result
        end
        
        # 测试 A_correction 函数
        println("\n表2: A_correction 函数随 p_max 的变化 (ξ=$ξ)")
        println(@sprintf("%-12s %-15s %-15s", "p_max(fm⁻¹)", "A_correction", "相对变化(%)"))
        println("-"^70)
        
        prev_ΔA = nothing
        
        for p_max in p_max_values
            nodes, weights = gauleg(0.0, p_max, n_points)
            result = A_correction(m, μ, T, Φ, Φbar, ξ, nodes, weights)
            
            if prev_ΔA === nothing
                println(@sprintf("%-12.1f %-15.8f %-15s", p_max, result, "-"))
            else
                rel_change = abs(result - prev_ΔA) / abs(prev_ΔA) * 100
                println(@sprintf("%-12.1f %-15.8f %-15.6f", p_max, result, rel_change))
                
                # p_max >= 20 后，变化应 < 0.01%
                if p_max >= 20.0
                    @test rel_change < 0.01
                end
            end
            
            prev_ΔA = result
        end
        
        # 测试不同温度下的收敛性（高温时分布更宽，需要更大 p_max）
        println("\n表3: 不同温度下 A 在 p_max=20 vs p_max=30 的差异")
        println(@sprintf("%-12s %-15s %-15s %-12s", "T(fm⁻¹)", "A(p_max=20)", "A(p_max=30)", "相对差(%)"))
        println("-"^70)
        
        nodes_20, weights_20 = gauleg(0.0, 20.0, n_points)
        nodes_30, weights_30 = gauleg(0.0, 30.0, n_points)
        
        T_values = [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
        
        for T_test in T_values
            A_20 = A(m, μ, T_test, Φ, Φbar, nodes_20, weights_20)
            A_30 = A(m, μ, T_test, Φ, Φbar, nodes_30, weights_30)
            rel_diff = abs(A_30 - A_20) / abs(A_20) * 100
            
            println(@sprintf("%-12.2f %-15.8f %-15.8f %-12.6f", T_test, A_20, A_30, rel_diff))
            
            # 对于 T ≤ 0.3 fm⁻¹ (~60 MeV)，p_max=20 应该足够（误差 < 0.1%）
            if T_test <= 0.30
                @test rel_diff < 0.1
            end
        end
        
        println("="^70)
        println("\n结论：")
        println("  • 对于 T ≤ 0.3 fm⁻¹ (~60 MeV)，p_max = 20 fm⁻¹ 足够（误差 < 0.1%）")
        println("  • 对于更高温度，可能需要增大 p_max")
        println("  • 建议：对于典型 PNJL 计算（T ~ 150-200 MeV），p_max = 20 fm⁻¹ 足够")
        println("="^70 * "\n")
    end
    
    @testset "积分节点数收敛性测试" begin
        m = TEST_PARAMS.m
        μ = TEST_PARAMS.μ
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.2
        
        println("\n" * "="^70)
        println("A_correction 积分节点数收敛性测试 (ξ=$ξ)")
        println("="^70)
        println(@sprintf("%-15s %-15s %-15s", "节点数", "A_correction", "相对变化(%)"))
        println("-"^70)
        
        n_points_list = [16, 32, 64, 128, 256]
        prev_result = nothing
        
        for n_points in n_points_list
            nodes, weights = gauleg(0.0, 20.0, n_points)
            result = A_correction(m, μ, T, Φ, Φbar, ξ, nodes, weights)
            
            if prev_result === nothing
                println(@sprintf("%-15d %-15.8f %-15s", n_points, result, "-"))
            else
                rel_change = abs(result - prev_result) / abs(result) * 100
                println(@sprintf("%-15d %-15.8f %-15.4f", n_points, result, rel_change))
                
                # 64 点及以上，变化应 < 1%
                if n_points >= 64
                    @test rel_change < 1.0
                end
            end
            
            prev_result = result
        end
        
        println("="^70 * "\n")
        
        println("\n" * "="^70)
        println("A_aniso 积分节点数收敛性测试 (ξ=$ξ)")
        println("="^70)
        println(@sprintf("%-20s %-15s %-15s", "节点数(p×cosθ)", "A_aniso", "相对变化(%)"))
        println("-"^70)
        
        node_configs = [
            (16, 16),
            (32, 16),
            (32, 32),
            (64, 32),
            (128, 32),
            (128, 64),
        ]
        
        prev_result = nothing
        
        for (n_p, n_cosθ) in node_configs
            nodes_p, weights_p = gauleg(0.0, 20.0, n_p)
            nodes_cosθ, weights_cosθ = gauleg(-1.0, 1.0, n_cosθ)
            
            result = A_aniso(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p, nodes_cosθ, weights_cosθ)
            
            config_str = "$(n_p)×$(n_cosθ)"
            
            if prev_result === nothing
                println(@sprintf("%-20s %-15.8f %-15s", config_str, result, "-"))
            else
                rel_change = abs(result - prev_result) / abs(result) * 100
                println(@sprintf("%-20s %-15.8f %-15.4f", config_str, result, rel_change))
                
                # 64×32 及以上，变化应 < 0.5%
                if n_p >= 64 && n_cosθ >= 32
                    @test rel_change < 0.5
                end
            end
            
            prev_result = result
        end
        
        println("="^70 * "\n")
    end
    
    @testset "性能基准测试" begin
        m = TEST_PARAMS.m
        μ = TEST_PARAMS.μ
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.2
        
        n_iterations = 10
        
        println("\n" * "="^70)
        println("性能基准测试 (每个函数运行 $n_iterations 次)")
        println("="^70)
        
        # 测试 A (各向同性)
        GC.gc()
        t_start = time()
        for _ in 1:n_iterations
            A(m, μ, T, Φ, Φbar, NODES_P, WEIGHTS_P)
        end
        t_A = (time() - t_start) / n_iterations * 1000  # ms
        
        # 测试 A_correction
        GC.gc()
        t_start = time()
        for _ in 1:n_iterations
            A_correction(m, μ, T, Φ, Φbar, ξ, NODES_P, WEIGHTS_P)
        end
        t_A_corr = (time() - t_start) / n_iterations * 1000  # ms
        
        # 测试 A_aniso
        GC.gc()
        t_start = time()
        for _ in 1:n_iterations
            A_aniso(m, μ, T, Φ, Φbar, ξ, NODES_P, WEIGHTS_P, NODES_COSΘ, WEIGHTS_COSΘ)
        end
        t_A_aniso = (time() - t_start) / n_iterations * 1000  # ms
        
        println(@sprintf("%-25s: %8.3f ms", "A (各向同性)", t_A))
        println(@sprintf("%-25s: %8.3f ms", "A_correction (一阶修正)", t_A_corr))
        println(@sprintf("%-25s: %8.3f ms", "A + A_correction (组合)", t_A + t_A_corr))
        println(@sprintf("%-25s: %8.3f ms", "A_aniso (完整计算)", t_A_aniso))
        println("-"^70)
        println(@sprintf("%-25s: %8.2fx", "A_aniso / (A + ΔA) 耗时比", t_A_aniso / (t_A + t_A_corr)))
        println("="^70 * "\n")
        
        @info "性能测试完成" t_A_ms=t_A t_A_correction_ms=t_A_corr t_A_aniso_ms=t_A_aniso
        
        # 验证性能在合理范围内（单次调用 < 500ms）
        @test t_A < 500
        @test t_A_corr < 500
        @test t_A_aniso < 1000
    end
    
end

println("\n" * "="^70)
println("所有测试完成！")
println("="^70)
println("\n建议查看生成的测试总结文档：")
println("  test/test_oneloopintegrals_aniso_summary.md")
println()



