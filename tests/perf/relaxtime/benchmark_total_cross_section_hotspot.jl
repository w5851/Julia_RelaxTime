# benchmark_total_cross_section_hotspot.jl
#
# 目的：定位“散射总截面 / 散射振幅”相关的性能热点，并提供可重复的基准入口。
#
# 对应系统流程：
# - 总截面：src/relaxtime/TotalCrossSection.jl::total_cross_section
# - 振幅：src/relaxtime/ScatteringAmplitude.jl::scattering_amplitude_squared
# - 传播子/极化：src/relaxtime/TotalPropagator.jl / src/relaxtime/PolarizationCache.jl
#
# 运行：
# - julia --project=. tests/perf/relaxtime/benchmark_total_cross_section_hotspot.jl

using BenchmarkTools
using Statistics
using Profile

include(joinpath(@__DIR__, "../../../src/Constants_PNJL.jl"))
include(joinpath(@__DIR__, "../../../src/integration/GaussLegendre.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/EffectiveCouplings.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/ScatteringAmplitude.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/TotalCrossSection.jl"))

using .Constants_PNJL
using .GaussLegendre: gausslegendre
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_effective_couplings, calculate_G_from_A
using .ScatteringAmplitude: scattering_amplitude_squared
using .TotalCrossSection: total_cross_section

# -----------------------------
# 1) 准备一组“尽量接近真实工作流”的参数
# -----------------------------
T = 0.15
Φ = 0.5
Φbar = 0.5
ξ = 0.0

m_u = 1.52
m_s = 2.50
μ_u = 0.3
μ_s = 0.0

# A 函数需要动量高斯求积
nodes_p, weights_p = gausslegendre(64)
p_max = 20.0
nodes_p = @. (nodes_p + 1.0) / 2.0 * p_max
weights_p = @. weights_p * p_max / 2.0

println("预计算 A_u/A_s（可能较慢）...")
A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)

G_u = calculate_G_from_A(A_u, m_u)
G_s = calculate_G_from_A(A_s, m_s)

# 注意：这里沿用项目内部的 fm 单位；耦合常数的具体换算以 src/Constants_PNJL.jl 为准
G_fm2 = Constants_PNJL.G_fm2
K_fm5 = Constants_PNJL.K_fm5
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

quark_params = (
    m = (u=m_u, d=m_u, s=m_s),
    μ = (u=μ_u, d=μ_u, s=μ_s),
    A = (u=A_u, d=A_u, s=A_s),
)
thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)

process = :uu_to_uu
s = 31.0
# 取一个典型的 t（避免极端阈值）
t = -2.0

println("\n预热（编译/缓存）...")
scattering_amplitude_squared(process, s, t, quark_params, thermo_params, K_coeffs)
total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=8)

# -----------------------------
# 2) 基准：散射振幅单点评估
# -----------------------------
println("\nBenchmark: scattering_amplitude_squared (single point)")
bench_M = @benchmark scattering_amplitude_squared($process, $s, $t, $quark_params, $thermo_params, $K_coeffs) samples=10 evals=1
show(stdout, MIME("text/plain"), bench_M)
println()

# -----------------------------
# 3) 基准：总截面（t 积分）
# -----------------------------
println("\nBenchmark: total_cross_section (t integration)")
bench_sigma = @benchmark total_cross_section($process, $s, $quark_params, $thermo_params, $K_coeffs; n_points=16) samples=10 evals=1
show(stdout, MIME("text/plain"), bench_sigma)
println()

# -----------------------------
# 4) Profile：看调用栈热点（输出可能较长）
# -----------------------------
println("\nProfile: total_cross_section (n_points=64, looped)")
Profile.clear()
@profile begin
    for _ in 1:200
        total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=64)
    end
end
Profile.print(maxdepth=18)

println("\nProfile: scattering_amplitude_squared (looped)")
Profile.clear()
@profile begin
    for _ in 1:10_000
        scattering_amplitude_squared(process, s, t, quark_params, thermo_params, K_coeffs)
    end
end
Profile.print(maxdepth=18)
