# 详细用时分析
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))

println("=" ^ 60)
println("Timing Analysis")
println("=" ^ 60)

λ = -1.0; k = 0.01; m = 0.3; m_prime = 0.3; ξ = -0.2
T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

# 预热所有代码路径
for _ in 1:3
    OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
        :quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=32
    )
end

println("\n1. 单次调用用时 (预热后, 1000次平均):")
n_runs = 1000

for (strat, sname) in [
    (OneLoopIntegralsCorrection.STRATEGY_QUADGK, "QUADGK"),
    (OneLoopIntegralsCorrection.STRATEGY_CLUSTER_GL, "CLUSTER_GL"),
    (OneLoopIntegralsCorrection.STRATEGY_HYBRID, "HYBRID"),
]
    times = Float64[]
    for _ in 1:n_runs
        t0 = time_ns()
        OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
            :quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=strat, cluster_n=32
        )
        push!(times, (time_ns() - t0) / 1e6)
    end
    avg = sum(times) / length(times)
    std = sqrt(sum((t - avg)^2 for t in times) / length(times))
    println(@sprintf("  %-12s: %.3f ± %.3f ms", sname, avg, std))
end

println("\n2. 不同节点数的用时 (HYBRID):")
for n in [16, 24, 32, 48, 64]
    times = Float64[]
    for _ in 1:500
        t0 = time_ns()
        OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
            :quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=n
        )
        push!(times, (time_ns() - t0) / 1e6)
    end
    avg = sum(times) / length(times)
    println(@sprintf("  n=%2d: %.3f ms", n, avg))
end

println("\n3. 节点变换开销分析:")
Emin = m
Emax = OneLoopIntegrals.energy_cutoff(m)
n = 32

# 标准 GL
t0 = time_ns()
for _ in 1:10000
    OneLoopIntegralsCorrection.GaussLegendre.gauleg(Emin, Emax, n)
end
t_gl = (time_ns() - t0) / 1e6 / 10000
println(@sprintf("  标准 GL:     %.4f ms/call", t_gl))

# power_left
t0 = time_ns()
for _ in 1:10000
    OneLoopIntegralsCorrection.power_left_nodes(Emin, Emax, n)
end
t_power = (time_ns() - t0) / 1e6 / 10000
println(@sprintf("  power_left:  %.4f ms/call", t_power))

# clustered (tanh)
t0 = time_ns()
for _ in 1:10000
    OneLoopIntegralsCorrection.clustered_gl_nodes(Emin, Emax, n)
end
t_tanh = (time_ns() - t0) / 1e6 / 10000
println(@sprintf("  tanh:        %.4f ms/call", t_tanh))

# DE
t0 = time_ns()
for _ in 1:10000
    OneLoopIntegralsCorrection.de_nodes(Emin, Emax, n)
end
t_de = (time_ns() - t0) / 1e6 / 10000
println(@sprintf("  DE:          %.4f ms/call", t_de))

println("\n结论：节点变换开销很小，主要时间在被积函数求值")
