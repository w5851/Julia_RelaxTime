# 分析无根情况的被积函数特性
using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
using QuadGK: quadgk

# 无根情况参数
λ = -0.5; k = 0.05; m = 0.3; m_prime = 0.3; ξ = -0.1
T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

Emin = m
Emax = OneLoopIntegrals.energy_cutoff(m)

integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(:quark, λ, k, m, m_prime, E, ξ, T, μ, Φ, Φbar)

# 采样被积函数
Es = range(Emin + 0.001, Emax - 0.001, length=100)
vals = integrand.(Es)

println("无根情况被积函数分析:")
println(@sprintf("  区间: [%.3f, %.3f]", Emin, Emax))
println(@sprintf("  最大值: %.6e at E=%.3f", maximum(vals), Es[argmax(vals)]))
println(@sprintf("  最小值: %.6e at E=%.3f", minimum(vals), Es[argmin(vals)]))
println(@sprintf("  变化范围: %.2fx", maximum(abs.(vals)) / minimum(abs.(vals))))

# 检查是否有尖峰
diffs = diff(vals)
max_diff_idx = argmax(abs.(diffs))
println(@sprintf("  最大变化率位置: E≈%.3f", Es[max_diff_idx]))

# 对比不同方法
ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
println(@sprintf("\n参考值: %.10e", ref))

println("\n不同节点数的标准 GL 误差:")
for n in [16, 24, 32, 48, 64]
    nodes, weights = OneLoopIntegralsCorrection.GaussLegendre.gauleg(Emin, Emax, n)
    val = sum(weights .* integrand.(nodes))
    err = abs((val - ref) / ref)
    println(@sprintf("  n=%2d: relerr=%.2e", n, err))
end

# 结论：无根情况下被积函数是否光滑？
println("\n结论：无根情况下被积函数相对光滑，标准 GL 收敛较慢是因为区间较大。")
println("可以考虑使用 tanh 聚簇来加速收敛。")
