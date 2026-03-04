# 详细性能对比：策略 × 节点数（误差与用时）
#
# 对应系统流程步骤：
# - `OneLoopIntegralsCorrection.tilde_B0_correction_k_positive`
# - QuadGK (`quadgk`) 作为参考值
#
# 测试内容：
# - 固定一组典型参数，对比不同策略在不同 `cluster_n` 下的误差与耗时
#
# 运行方式：
# - `julia --project=. tests/perf/relaxtime/performance_comparison.jl`
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
using QuadGK: quadgk

println("=" ^ 70)
println("Performance Comparison: Strategies × Node Counts")
println("=" ^ 70)

test_cases = [
    (name="有根(标准)", λ=-1.0, k=0.01, m=0.3, m_prime=0.3),
    (name="有根(大k)", λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
    (name="无根", λ=-0.5, k=0.05, m=0.3, m_prime=0.3),
    (name="无根(正λ)", λ=0.5, k=0.1, m=0.3, m_prime=0.35),
]

ξ = -0.2; T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

strategies = [
    (OneLoopIntegralsCorrection.STRATEGY_QUADGK, "QUADGK"),
    (OneLoopIntegralsCorrection.STRATEGY_INTERVAL_GL, "INTERVAL_GL"),
    (OneLoopIntegralsCorrection.STRATEGY_CLUSTER_GL, "CLUSTER_GL"),
    (OneLoopIntegralsCorrection.STRATEGY_HYBRID, "HYBRID"),
]

ns = [16, 24, 32, 48, 64]

# 收集所有结果
results = Dict()

for tc in test_cases
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
    
    results[tc.name] = Dict()
    
    for (strat, sname) in strategies
        results[tc.name][sname] = Dict()
        for n in ns
            # 多次运行取平均时间
            times = Float64[]
            err = 0.0
            for _ in 1:5
                t0 = time()
                result = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
                    :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
                    strategy=strat, cluster_n=n
                )
                push!(times, (time() - t0) * 1000)
                err = abs((result[1] - ref) / ref)
            end
            results[tc.name][sname][n] = (err=err, time_ms=sum(times)/length(times))
        end
    end
end

# 打印结果表格
println("\n相对误差 (relerr):")
println("-" ^ 70)

for tc in test_cases
    println("\n$(tc.name):")
    print(@sprintf("%-12s", "n"))
    for n in ns
        print(@sprintf("%10d", n))
    end
    println()
    
    for (_, sname) in strategies
        print(@sprintf("%-12s", sname))
        for n in ns
            err = results[tc.name][sname][n].err
            if err < 1e-4
                print(@sprintf("%10.1e", err))
            else
                print(@sprintf("%9.1e*", err))  # 标记超过阈值
            end
        end
        println()
    end
end

# 找出 n=32 时各策略的表现
println("\n" * "=" ^ 70)
println("n=32 时各策略对比 (目标: relerr < 1e-4)")
println("=" ^ 70)

for tc in test_cases
    println("\n$(tc.name):")
    for (_, sname) in strategies
        r = results[tc.name][sname][32]
        status = r.err < 1e-4 ? "✓" : "✗"
        println(@sprintf("  %s %-12s: relerr=%.2e", status, sname, r.err))
    end
end

# 统计 n=32 时 HYBRID 的通过率
println("\n" * "=" ^ 70)
println("Summary for n=32:")
println("=" ^ 70)

for (_, sname) in strategies
    pass_count = 0
    max_err = 0.0
    for tc in test_cases
        err = results[tc.name][sname][32].err
        if err < 1e-4
            pass_count += 1
        end
        max_err = max(max_err, err)
    end
    println(@sprintf("%-12s: %d/4 pass, max_err=%.2e", sname, pass_count, max_err))
end

# 如果 n=32 不够，找最小满足的 n
println("\n" * "=" ^ 70)
println("Minimum n for all cases < 1e-4:")
println("=" ^ 70)

for (_, sname) in strategies
    min_n = 0
    for n in ns
        all_pass = true
        for tc in test_cases
            if results[tc.name][sname][n].err >= 1e-4
                all_pass = false
                break
            end
        end
        if all_pass
            min_n = n
            break
        end
    end
    if min_n > 0
        println(@sprintf("%-12s: n=%d", sname, min_n))
    else
        println(@sprintf("%-12s: n>64 needed", sname))
    end
end
