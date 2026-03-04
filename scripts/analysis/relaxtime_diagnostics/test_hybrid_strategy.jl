using Test

# 注意：该文件历史上是诊断/对比脚本（包含 println 与耗时统计），不属于单元测试。
# 为保持 UNIT_PROFILE=full 可用，这里默认只做一个轻量占位。
# 如需运行诊断脚本，请手动执行：
#   julia --project=. tests/analysis/relaxtime_diagnostics/test_hybrid_strategy.jl

if get(ENV, "RUN_HYBRID_DIAGNOSTIC", "0") != "1"
    @testset "Hybrid strategy diagnostic (skipped in unit)" begin
        @test true
    end
else
    # 原诊断脚本（保留，仅在显式开启时执行）
    using Printf
    include(joinpath(@__DIR__, "../../src/relaxtime/OneLoopIntegralsAniso.jl"))
    include(joinpath(@__DIR__, "../../src/relaxtime/OneLoopIntegrals.jl"))
    using QuadGK: quadgk

println("=" ^ 70)
println("Testing STRATEGY_HYBRID vs other strategies")
println("=" ^ 70)

# 测试参数
λ = -1.0; k = 0.01; m = 0.3; m_prime = 0.3; ξ = -0.2
T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

Emin = m
Emax = OneLoopIntegrals.energy_cutoff(m)

# 参考值
integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, λ, k, m, m_prime, E, ξ, T, μ, Φ, Φbar)
val_ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-14)
println(@sprintf("Reference: %.15g\n", val_ref))

strategies = [
    (OneLoopIntegralsCorrection.STRATEGY_QUADGK, "QUADGK"),
    (OneLoopIntegralsCorrection.STRATEGY_INTERVAL_GL, "INTERVAL_GL"),
    (OneLoopIntegralsCorrection.STRATEGY_CLUSTER_GL, "CLUSTER_GL"),
    (OneLoopIntegralsCorrection.STRATEGY_HYBRID, "HYBRID"),
]

println("Comparison across node counts:")
println("-" ^ 70)

for n in [32, 64, 128]
    println(@sprintf("\nn = %d:", n))
    for (strat, name) in strategies
        result = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
            :quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=strat, cluster_n=n, diagnostics=true
        )
        real_part, _, diag = result
        err = abs((real_part - val_ref) / val_ref)
        println(@sprintf("  %-12s: relerr=%.2e  time=%.2fms", name, err, diag.elapsed_ms))
    end
end

# 详细诊断
println("\n" * "=" ^ 70)
println("Detailed diagnostics for STRATEGY_HYBRID (n=64)")
println("=" ^ 70)

result = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
    :quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
    strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID,
    cluster_n=64,
    diagnostics=true
)
real_part, imag_part, diag = result

println(@sprintf("Strategy: %s", diag.strategy))
println(@sprintf("Roots found: %d", diag.n_roots))
println(@sprintf("Root positions: %s", diag.roots))
println(@sprintf("Intervals: %d", diag.n_intervals))
for (i, interval) in enumerate(diag.intervals)
    println(@sprintf("  [%d] (%.8f, %.8f)", i, interval[1], interval[2]))
end
println(@sprintf("Real part: %.15g", diag.real_part))
println(@sprintf("Imag part: %.15g", diag.imag_part))
println(@sprintf("Reference: %.15g", val_ref))
println(@sprintf("Rel error: %.2e", abs((real_part - val_ref) / val_ref)))

# 多参数测试
println("\n" * "=" ^ 70)
println("Multi-parameter test")
println("=" ^ 70)

test_cases = [
    (λ=-1.0, k=0.01, m=0.3, m_prime=0.3),
    (λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
    (λ=-0.5, k=0.05, m=0.3, m_prime=0.3),
    (λ=0.5, k=0.1, m=0.3, m_prime=0.35),
]

for tc in test_cases
    Emin_tc = tc.m
    Emax_tc = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand_tc(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    ref_tc, _ = quadgk(integrand_tc, Emin_tc, Emax_tc; rtol=1e-12)
    
    # CLUSTER_GL
    r_cluster = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
        strategy=OneLoopIntegralsCorrection.STRATEGY_CLUSTER_GL, cluster_n=64
    )
    err_cluster = abs((r_cluster[1] - ref_tc) / ref_tc)
    
    # HYBRID
    r_hybrid = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
        strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=64
    )
    err_hybrid = abs((r_hybrid[1] - ref_tc) / ref_tc)
    
    improvement = err_cluster / max(err_hybrid, 1e-16)
    println(@sprintf("λ=%.1f k=%.2f: CLUSTER=%.2e HYBRID=%.2e improvement=%.1fx",
        tc.λ, tc.k, err_cluster, err_hybrid, improvement))
end

println("\n" * "=" ^ 70)
println("Test completed!")
println("=" ^ 70)

end
