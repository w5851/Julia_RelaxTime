# 分析 HYBRID 策略的时间开销分布 (v2)
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))

println("=" ^ 70)
println("HYBRID 策略时间开销分析 (v2)")
println("=" ^ 70)

tc = (λ=-1.0, k=0.01, m=0.3, m_prime=0.3)
ξ = -0.2; T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0
Emin = tc.m
Emax = OneLoopIntegrals.energy_cutoff(tc.m)

n_runs = 1000  # 更多次数以获得稳定结果

# 预热
for _ in 1:10
    OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
    OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
        strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=32)
end

# 1. 单独测量根查找
t0 = time_ns()
for _ in 1:n_runs
    OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
end
avg_roots = (time_ns() - t0) / n_runs / 1e6

# 2. 单独测量节点生成 (假设已有根和区间)
roots = OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
intervals = OneLoopIntegralsCorrection.build_intervals_from_roots(roots, Emin, Emax)
t0 = time_ns()
for _ in 1:n_runs
    for (idx, (a, b)) in enumerate(intervals)
        sing_pos = if idx == 1
            OneLoopIntegralsCorrection.SING_RIGHT
        elseif idx == length(intervals)
            OneLoopIntegralsCorrection.SING_LEFT
        else
            OneLoopIntegralsCorrection.SING_BOTH
        end
        xs, wx = OneLoopIntegralsCorrection.hybrid_nodes(a, b, 32, sing_pos)
    end
end
avg_nodes = (time_ns() - t0) / n_runs / 1e6

# 3. 单独测量被积函数求值 (假设已有节点)
integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
    :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
all_xs = Float64[]
for (idx, (a, b)) in enumerate(intervals)
    sing_pos = if idx == 1
        OneLoopIntegralsCorrection.SING_RIGHT
    elseif idx == length(intervals)
        OneLoopIntegralsCorrection.SING_LEFT
    else
        OneLoopIntegralsCorrection.SING_BOTH
    end
    xs, _ = OneLoopIntegralsCorrection.hybrid_nodes(a, b, 32, sing_pos)
    append!(all_xs, xs)
end
t0 = time_ns()
for _ in 1:n_runs
    s = 0.0
    for x in all_xs
        s += integrand(x)
    end
end
avg_eval = (time_ns() - t0) / n_runs / 1e6

# 4. 完整 HYBRID 调用
t0 = time_ns()
for _ in 1:n_runs
    OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
        strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=32)
end
avg_total = (time_ns() - t0) / n_runs / 1e6

println("\n时间开销分解 ($(n_runs)次平均, 有根情况):")
println("-" ^ 50)
println(@sprintf("根查找 (find_roots_AB):     %.4f ms", avg_roots))
println(@sprintf("节点生成 (3区间):           %.4f ms", avg_nodes))
println(@sprintf("被积函数求值 (%d点):        %.4f ms", length(all_xs), avg_eval))
println(@sprintf("其他开销:                   %.4f ms", avg_total - avg_roots - avg_nodes - avg_eval))
println("-" ^ 50)
println(@sprintf("总计:                       %.4f ms", avg_total))
println(@sprintf("各部分占比: 根查找=%.1f%%, 节点=%.1f%%, 求值=%.1f%%", 
    100*avg_roots/avg_total, 100*avg_nodes/avg_total, 100*avg_eval/avg_total))

# 对比无根情况
println("\n" * "=" ^ 70)
println("无根情况对比:")
println("=" ^ 70)

tc_noroot = (λ=-0.5, k=0.05, m=0.3, m_prime=0.3)
Emin_nr = tc_noroot.m
Emax_nr = OneLoopIntegrals.energy_cutoff(tc_noroot.m)

# 根查找（无根）
t0 = time_ns()
for _ in 1:n_runs
    OneLoopIntegralsCorrection.find_roots_AB(tc_noroot.λ, tc_noroot.k, tc_noroot.m, tc_noroot.m_prime, Emin_nr, Emax_nr)
end
avg_roots_nr = (time_ns() - t0) / n_runs / 1e6

# 完整调用（无根）
t0 = time_ns()
for _ in 1:n_runs
    OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
        :quark, tc_noroot.λ, tc_noroot.k, tc_noroot.m, tc_noroot.m_prime, μ, T, Φ, Φbar, ξ;
        strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=32)
end
avg_total_nr = (time_ns() - t0) / n_runs / 1e6

println(@sprintf("无根情况: 根查找=%.4f ms, 总计=%.4f ms", avg_roots_nr, avg_total_nr))
println(@sprintf("根查找占比: %.1f%%", 100*avg_roots_nr/avg_total_nr))

println("\n" * "=" ^ 70)
println("优化建议:")
println("=" ^ 70)
println("1. 根查找使用 N=2000 点扫描，占用约 %.1f%% 时间" % (100*avg_roots/avg_total))
println("   - 可考虑减少到 N=500 或使用更智能的根查找算法")
println("2. 如果参数不变，可缓存根位置避免重复计算")
