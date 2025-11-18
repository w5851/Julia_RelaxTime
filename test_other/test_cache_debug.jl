"""
调试缓存为什么不工作
"""

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../src/relaxtime"))

using Printf

# 加载模块
include("../src/Constants_PNJL.jl")
include("../src/integration/GaussLegendre.jl")
include("../src/relaxtime/EffectiveCouplings.jl")
include("../src/relaxtime/OneLoopIntegrals.jl")
include("../src/relaxtime/TotalPropagator.jl")

using .Constants_PNJL: G_fm2, K_fm5, ħc_MeV_fm
using .GaussLegendre: gauleg
using .EffectiveCouplings: calculate_effective_couplings, calculate_G_from_A
using .OneLoopIntegrals: A
using .TotalPropagator: calculate_all_propagators, reset_cache!, get_cache_stats  # 从TotalPropagator导入缓存函数

println("="^70)
println("缓存调试测试")
println("="^70)

# 设置物理参数
T_MeV = 150.0
m_u_MeV = 300.0
m_s_MeV = 500.0

T = T_MeV / ħc_MeV_fm
m_u = m_u_MeV / ħc_MeV_fm
m_s = m_s_MeV / ħc_MeV_fm
μ_u = 0.0
μ_s = 0.0
Φ = 0.5
Φbar = 0.5
ξ = 0.0

k0_MeV = 100.0
k_norm_MeV = 50.0
k0 = k0_MeV / ħc_MeV_fm
k_norm = k_norm_MeV / ħc_MeV_fm

# 计算A函数和K系数（与test_total_propagator_performance.jl保持一致）
nodes_p, weights_p = gauleg(0.0, 20.0, 64)
A_u = A(m_u, μ_u, T, Φ, Φbar, nodes_p, weights_p)
A_s = A(m_s, μ_s, T, Φ, Φbar, nodes_p, weights_p)
G_u = calculate_G_from_A(A_u)
G_s = calculate_G_from_A(A_s)
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

quark_params = (
    m = (u = m_u, d = m_u, s = m_s),
    μ = (u = μ_u, d = μ_u, s = μ_s),
    A = (u = A_u, d = A_u, s = A_s)
)

thermo_params = (T = T, Φ = Φ, Φbar = Φbar, ξ = ξ)

println("\n重置缓存...")
reset_cache!()

println("初始缓存统计:")
stats0 = get_cache_stats()
println("  total_calls: $(stats0.total_calls)")
println("  cache_hits: $(stats0.cache_hits)")
println("  cache_size: $(stats0.cache_size)")

println("\n执行第一次计算: uu_to_uu")
println("调用参数:")
println("  k0=$(k0), k_norm=$(k_norm)")
println("  T=$(thermo_params.T), ξ=$(thermo_params.ξ)")
println("  m_u=$(quark_params.m.u), A_u=$(quark_params.A.u)")

result1 = calculate_all_propagators(:uu_to_uu, k0, k_norm, quark_params, thermo_params, K_coeffs)
println("  result: t=$(result1.t), u=$(result1.u)")

println("\n第一次计算后缓存统计:")
stats1 = get_cache_stats()
println("  total_calls: $(stats1.total_calls)")
println("  cache_hits: $(stats1.cache_hits)")
println("  cache_misses: $(stats1.cache_misses)")
println("  hit_rate: $(stats1.hit_rate * 100)%")
println("  cache_size: $(stats1.cache_size)")

println("\n执行第二次计算: uu_to_uu (相同参数)")
result2 = calculate_all_propagators(:uu_to_uu, k0, k_norm, quark_params, thermo_params, K_coeffs)

println("第二次计算后缓存统计:")
stats2 = get_cache_stats()
println("  total_calls: $(stats2.total_calls)")
println("  cache_hits: $(stats2.cache_hits)")
println("  cache_misses: $(stats2.cache_misses)")
println("  hit_rate: $(stats2.hit_rate * 100)%")
println("  cache_size: $(stats2.cache_size)")

println("\n执行第三次计算: ud_to_ud")
result3 = calculate_all_propagators(:ud_to_ud, k0, k_norm, quark_params, thermo_params, K_coeffs)

println("第三次计算后缓存统计:")
stats3 = get_cache_stats()
println("  total_calls: $(stats3.total_calls)")
println("  cache_hits: $(stats3.cache_hits)")
println("  cache_misses: $(stats3.cache_misses)")
println("  hit_rate: $(stats3.hit_rate * 100)%")
println("  cache_size: $(stats3.cache_size)")

println("\n="^70)
if stats1.total_calls == 0
    println("❌ 问题确认：polarization_aniso_cached 没有被调用！")
    println("   需要检查 TotalPropagator.jl 的实现")
elseif stats2.cache_hits == 0
    println("❌ 缓存未命中：可能是参数差异或缓存机制问题")
else
    println("✅ 缓存正常工作")
    println("   命中率: $(stats3.hit_rate * 100)%")
end
println("="^70)
