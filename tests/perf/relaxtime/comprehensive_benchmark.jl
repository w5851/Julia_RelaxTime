"""
综合基准测试：精度 + 实际用时（多参数组合）

对应系统流程步骤：
- 单圈积分修正：`OneLoopIntegralsCorrection.tilde_B0_correction_k_positive`
- 参考值：QuadGK (`quadgk`) 作为高精度对照

测试内容：
- 对一组“有根/无根”参数，比较不同策略（QUADGK/CLUSTER_GL/HYBRID）的：
    - 平均耗时（多次循环取均值）
    - 相对误差（相对 QuadGK 参考值）

运行方式：
- `julia --project=. tests/perf/relaxtime/comprehensive_benchmark.jl`

备注：
- 该脚本兼顾精度与性能，用于选择默认策略/节点数的经验依据。
"""
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
using QuadGK: quadgk

println("=" ^ 70)
println("Comprehensive Benchmark: Accuracy + Timing")
println("=" ^ 70)

# 更多测试参数，特别是无根情况
test_cases = [
    # 有根情况
    (name="有根1", λ=-1.0, k=0.01, m=0.3, m_prime=0.3, has_roots=true),
    (name="有根2", λ=-1.0, k=0.1, m=0.3, m_prime=0.3, has_roots=true),
    (name="有根3", λ=-0.8, k=0.02, m=0.25, m_prime=0.25, has_roots=true),
    # 无根情况 - 多种参数组合
    (name="无根1", λ=-0.5, k=0.05, m=0.3, m_prime=0.3, has_roots=false),
    (name="无根2", λ=0.5, k=0.1, m=0.3, m_prime=0.35, has_roots=false),
    (name="无根3", λ=-0.3, k=0.08, m=0.35, m_prime=0.35, has_roots=false),
    (name="无根4", λ=0.8, k=0.15, m=0.25, m_prime=0.3, has_roots=false),
    (name="无根5", λ=-0.2, k=0.03, m=0.4, m_prime=0.4, has_roots=false),
]

ξ = -0.2; T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

strategies = [
    (OneLoopIntegralsCorrection.STRATEGY_QUADGK, "QUADGK"),
    (OneLoopIntegralsCorrection.STRATEGY_CLUSTER_GL, "CLUSTER_GL"),
    (OneLoopIntegralsCorrection.STRATEGY_HYBRID, "HYBRID"),
]

# 预热
for tc in test_cases[1:2]
    for (strat, _) in strategies
        OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
            :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
            strategy=strat, cluster_n=32
        )
    end
end

println("\n1. 实际用时测试 (n=32, 100次平均):")
println("-" ^ 70)

timing_results = Dict()
n_runs = 100

for tc in test_cases
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
    
    timing_results[tc.name] = Dict()
    
    print(@sprintf("%-8s: ", tc.name))
    for (strat, sname) in strategies
        times = Float64[]
        err = 0.0
        for _ in 1:n_runs
            t0 = time_ns()
            result = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
                :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
                strategy=strat, cluster_n=32
            )
            push!(times, (time_ns() - t0) / 1e6)  # ms
            err = abs((result[1] - ref) / ref)
        end
        avg_time = sum(times) / length(times)
        timing_results[tc.name][sname] = (time_ms=avg_time, err=err)
        print(@sprintf("%s=%.3fms ", sname, avg_time))
    end
    println()
end

# 汇总平均用时
println("\n平均用时汇总:")
for (_, sname) in strategies
    avg = sum(timing_results[tc.name][sname].time_ms for tc in test_cases) / length(test_cases)
    println(@sprintf("  %-12s: %.3f ms", sname, avg))
end

println("\n" * "=" ^ 70)
println("2. 无根情况下不同聚簇方法对比:")
println("=" ^ 70)

no_root_cases = filter(tc -> !tc.has_roots, test_cases)

println("\n被积函数特性分析:")
for tc in no_root_cases
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    
    # 采样找最大变化位置
    Es = range(Emin + 0.001, Emax - 0.001, length=200)
    vals = integrand.(Es)
    diffs = abs.(diff(vals))
    max_diff_idx = argmax(diffs)
    max_change_pos = Es[max_diff_idx]
    
    # 判断变化剧烈位置是在左端还是右端
    rel_pos = (max_change_pos - Emin) / (Emax - Emin)
    side = rel_pos < 0.5 ? "左端" : "右端"
    
    println(@sprintf("  %-8s: 最大变化在 E=%.3f (%s, 相对位置=%.1f%%)", 
        tc.name, max_change_pos, side, rel_pos * 100))
end

println("\n不同聚簇方法精度对比 (n=32):")
println(@sprintf("%-8s %12s %12s %12s %12s", "Case", "标准GL", "power_left", "power_right", "tanh(β=8)"))
println("-" ^ 60)

for tc in no_root_cases
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
    
    n = 32
    
    # 标准 GL
    nodes, weights = OneLoopIntegralsCorrection.GaussLegendre.gauleg(Emin, Emax, n)
    err_gl = abs((sum(weights .* integrand.(nodes)) - ref) / ref)
    
    # power_left
    xs_l, wx_l = OneLoopIntegralsCorrection.power_left_nodes(Emin, Emax, n)
    err_left = abs((sum(wx_l .* integrand.(xs_l)) - ref) / ref)
    
    # power_right
    xs_r, wx_r = OneLoopIntegralsCorrection.power_right_nodes(Emin, Emax, n)
    err_right = abs((sum(wx_r .* integrand.(xs_r)) - ref) / ref)
    
    # tanh
    xs_t, wx_t = OneLoopIntegralsCorrection.clustered_gl_nodes(Emin, Emax, n; beta=8.0)
    err_tanh = abs((sum(wx_t .* integrand.(xs_t)) - ref) / ref)
    
    println(@sprintf("%-8s %12.2e %12.2e %12.2e %12.2e", 
        tc.name, err_gl, err_left, err_right, err_tanh))
end

println("\n" * "=" ^ 70)
println("3. 结论:")
println("=" ^ 70)
