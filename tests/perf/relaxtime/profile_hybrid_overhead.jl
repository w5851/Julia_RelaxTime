# 分析 HYBRID 策略的时间开销分布
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))

println("=" ^ 70)
println("HYBRID 策略时间开销分析")
println("=" ^ 70)

tc = (λ=-1.0, k=0.01, m=0.3, m_prime=0.3)
ξ = -0.2; T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0
Emin = tc.m
Emax = OneLoopIntegrals.energy_cutoff(tc.m)

n_runs = 100

# 1. 根查找时间
times_roots = Float64[]
for _ in 1:n_runs
    t0 = time_ns()
    roots = OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
    push!(times_roots, (time_ns() - t0) / 1e6)
end
avg_roots = sum(times_roots) / length(times_roots)

# 2. 区间构建时间
roots = OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
times_intervals = Float64[]
for _ in 1:n_runs
    t0 = time_ns()
    intervals = OneLoopIntegralsCorrection.build_intervals_from_roots(roots, Emin, Emax)
    push!(times_intervals, (time_ns() - t0) / 1e6)
end
avg_intervals = sum(times_intervals) / length(times_intervals)

# 3. 节点生成时间
intervals = OneLoopIntegralsCorrection.build_intervals_from_roots(roots, Emin, Emax)
times_nodes = Float64[]
for _ in 1:n_runs
    t0 = time_ns()
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
    push!(times_nodes, (time_ns() - t0) / 1e6)
end
avg_nodes = sum(times_nodes) / length(times_nodes)

# 4. 纯被积函数求值时间（不含节点生成）
integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
    :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)

# 预先生成所有节点
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

times_eval = Float64[]
for _ in 1:n_runs
    t0 = time_ns()
    for x in all_xs
        integrand(x)
    end
    push!(times_eval, (time_ns() - t0) / 1e6)
end
avg_eval = sum(times_eval) / length(times_eval)
println("总节点数: ", length(all_xs))

# 5. 完整 HYBRID 调用时间
times_total = Float64[]
for _ in 1:n_runs
    t0 = time_ns()
    OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
        strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=32)
    push!(times_total, (time_ns() - t0) / 1e6)
end
avg_total = sum(times_total) / length(times_total)

println("\n时间开销分解 (100次平均):")
println("-" ^ 50)
println(@sprintf("根查找 (find_roots_AB):     %.3f ms (%.1f%%)", avg_roots, 100*avg_roots/avg_total))
println(@sprintf("区间构建:                   %.3f ms (%.1f%%)", avg_intervals, 100*avg_intervals/avg_total))
println(@sprintf("节点生成 (hybrid_nodes):    %.3f ms (%.1f%%)", avg_nodes, 100*avg_nodes/avg_total))
println(@sprintf("被积函数求值 (%d点):        %.3f ms (%.1f%%)", length(all_xs), avg_eval, 100*avg_eval/avg_total))
println("-" ^ 50)
println(@sprintf("总计:                       %.3f ms", avg_total))

println("\n" * "=" ^ 70)
println("优化建议:")
println("=" ^ 70)
if avg_roots / avg_total > 0.5
    println("- 根查找占用 >50% 时间，可考虑减少扫描点数或缓存根位置")
end
if avg_nodes / avg_total > 0.3
    println("- 节点生成占用 >30% 时间，可考虑预计算/缓存节点")
end
