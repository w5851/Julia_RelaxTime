"""
QuadGK（高精度参考）vs HYBRID（固定节点）性能对比

对应系统流程步骤：
- 单圈积分修正：`OneLoopIntegralsCorrection.tilde_B0_correction_k_positive`

测试内容：
- 对一组“有根/无根”参数：
    - 计算 QuadGK 参考值（较宽松 rtol，仅作对比基线）
    - 计算 HYBRID(n=32)
    - 比较用时与相对误差

运行方式：
- `julia --project=. tests/perf/relaxtime/compare_quadgk_vs_hybrid.jl`
"""
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
using QuadGK: quadgk

println("=" ^ 70)
println("QuadGK (无限制) vs HYBRID (n=32) 性能对比")
println("=" ^ 70)

# 测试参数
test_cases = [
    # 有根情况
    (name="有根1", λ=-1.0, k=0.01, m=0.3, m_prime=0.3, has_roots=true),
    (name="有根2", λ=-1.0, k=0.1, m=0.3, m_prime=0.3, has_roots=true),
    (name="有根3", λ=-0.8, k=0.02, m=0.25, m_prime=0.25, has_roots=true),
    # 无根情况
    (name="无根1", λ=-0.5, k=0.05, m=0.3, m_prime=0.3, has_roots=false),
    (name="无根2", λ=0.5, k=0.1, m=0.3, m_prime=0.35, has_roots=false),
    (name="无根3", λ=-0.3, k=0.08, m=0.35, m_prime=0.35, has_roots=false),
    (name="无根4", λ=0.8, k=0.15, m=0.25, m_prime=0.3, has_roots=false),
    (name="无根5", λ=-0.2, k=0.03, m=0.4, m_prime=0.4, has_roots=false),
]

ξ = -0.2; T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

# 预热
for tc in test_cases[1:2]
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    quadgk(integrand, Emin, Emax; rtol=1e-3)
    OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
        strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=32)
end

println("\n性能测试 (100次平均):")
println("-" ^ 70)
println(@sprintf("%-8s %12s %12s %12s %12s %12s", 
    "Case", "QuadGK(ms)", "HYBRID(ms)", "加速比", "QuadGK误差", "HYBRID误差"))
println("-" ^ 80)

n_runs = 100
total_quadgk_time = 0.0
total_hybrid_time = 0.0

for tc in test_cases
    local Emin = tc.m
    local Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    
    # 高精度参考值
    ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
    
    # 测试 QuadGK (默认 rtol=1e-3)
    times_quadgk = Float64[]
    result_quadgk = 0.0
    for _ in 1:n_runs
        t0 = time_ns()
        result_quadgk, _ = quadgk(integrand, Emin, Emax; rtol=1e-3)
        push!(times_quadgk, (time_ns() - t0) / 1e6)
    end
    avg_quadgk = sum(times_quadgk) / length(times_quadgk)
    err_quadgk = abs((result_quadgk - ref) / ref)
    
    # 测试 HYBRID (n=32)
    times_hybrid = Float64[]
    result_hybrid = 0.0
    for _ in 1:n_runs
        t0 = time_ns()
        result_hybrid, _ = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
            :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
            strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=32)
        push!(times_hybrid, (time_ns() - t0) / 1e6)
    end
    avg_hybrid = sum(times_hybrid) / length(times_hybrid)
    err_hybrid = abs((result_hybrid - ref) / ref)
    
    speedup = avg_quadgk / avg_hybrid
    
    global total_quadgk_time += avg_quadgk
    global total_hybrid_time += avg_hybrid
    
    println(@sprintf("%-8s %12.3f %12.3f %12.2fx %12.2e %12.2e", 
        tc.name, avg_quadgk, avg_hybrid, speedup, err_quadgk, err_hybrid))
end

println("-" ^ 80)
avg_quadgk_all = total_quadgk_time / length(test_cases)
avg_hybrid_all = total_hybrid_time / length(test_cases)
println(@sprintf("%-8s %12.3f %12.3f %12.2fx", 
    "平均", avg_quadgk_all, avg_hybrid_all, avg_quadgk_all / avg_hybrid_all))

println("\n" * "=" ^ 70)
println("结论:")
println("=" ^ 70)
println("QuadGK: 自适应积分，节点数不固定，根据 rtol 自动调整")
println("HYBRID: 固定 n=32 节点/区间，使用根分割 + 自适应变换")
