# 详细追踪ss→ss散射中的警告来源

include("../src/Constants_PNJL.jl")
include("../src/relaxtime/ScatteringAmplitude.jl")

using .ScatteringAmplitude
using .ScatteringAmplitude.TotalPropagator

println("="^70)
println("详细追踪ss→ss散射中k0²-u<0警告")
println("="^70)

# 参数设置
T = 0.15
m_u = 0.3
m_s = 0.5
μ = 0.25
A_u = 0.05
A_s = 0.08

quark_params = (
    m = (u=m_u, d=m_u, s=m_s),
    μ = (u=μ, d=μ, s=μ),
    A = (u=A_u, d=A_u, s=A_s)
)

thermo_params = (T=T, Φ=0.5, Φbar=0.5, ξ=1.0)

K_coeffs = (
    K0_plus=10.0, K0_minus=10.0,
    K123_plus=5.0, K123_minus=5.0,
    K4567_plus=5.0, K4567_minus=5.0,
    K8_plus=10.0, K8_minus=10.0,
    K08_plus=0.0, K08_minus=0.0,
    det_K_plus=500.0, det_K_minus=500.0
)

s = 20.0
t = -0.2
process = :ss_to_ss

println("\n散射过程: $process")
println("s = $s fm⁻²")
println("t = $t fm⁻²")

# 获取质量
m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
u = m1^2 + m2^2 + m3^2 + m4^2 - s - t

println("质量: m1=$m1, m2=$m2, m3=$m3, m4=$m4 fm⁻¹")
println("u = $(round(u, digits=2)) fm⁻²")
println()

# 手动计算t道和u道的CMS动量
println("="^70)
println("手动计算CMS动量")
println("="^70)

println("\n## t道 CMS动量")
println("-"^70)
try
    cms_t = calculate_cms_momentum(process, s, t, :t, quark_params; u=u)
    println("✓ t道: k0 = $(cms_t.k0), k = $(cms_t.k) fm⁻¹")
catch e
    println("✗ t道计算失败: ", sprint(showerror, e))
end

println("\n## u道 CMS动量")
println("-"^70)
try
    cms_u = calculate_cms_momentum(process, s, t, :u, quark_params; u=u)
    println("✓ u道: k0 = $(cms_u.k0), k = $(cms_u.k) fm⁻¹")
catch e
    println("✗ u道计算失败: ", sprint(showerror, e))
end

# 现在调用calculate_all_propagators_by_channel，观察警告
println("\n" * "="^70)
println("调用calculate_all_propagators_by_channel（观察警告）")
println("="^70)

println("\n## 计算t道传播子")
println("-"^70)
try
    cms_t = calculate_cms_momentum(process, s, t, :t, quark_params; u=u)
    println("t道 CMS: k0 = $(cms_t.k0), k = $(cms_t.k)")
    
    println("调用 calculate_all_propagators_by_channel...")
    prop_t = calculate_all_propagators_by_channel(
        process, cms_t.k0, cms_t.k, quark_params, thermo_params, K_coeffs
    )
    println("✓ t道传播子计算成功")
    println("  D_t_S = $(prop_t.t_S)")
    println("  D_t_P = $(prop_t.t_P)")
catch e
    println("✗ t道传播子计算失败")
    println("错误: ", sprint(showerror, e)[1:min(500, end)])
end

println("\n## 计算u道传播子")
println("-"^70)
try
    cms_u = calculate_cms_momentum(process, s, t, :u, quark_params; u=u)
    println("u道 CMS: k0 = $(cms_u.k0), k = $(cms_u.k)")
    
    println("调用 calculate_all_propagators_by_channel...")
    prop_u = calculate_all_propagators_by_channel(
        process, cms_u.k0, cms_u.k, quark_params, thermo_params, K_coeffs
    )
    println("✓ u道传播子计算成功")
    println("  D_u_S = $(prop_u.u_S)")
    println("  D_u_P = $(prop_u.u_P)")
catch e
    println("✗ u道传播子计算失败")
    println("错误: ", sprint(showerror, e)[1:min(500, end)])
end

println("\n" * "="^70)
println("追踪完成")
println("="^70)
println("\n说明：警告可能来自内部极化函数计算，而不是外部CMS动量计算")
