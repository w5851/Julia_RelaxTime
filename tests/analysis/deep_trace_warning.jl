# 详细追踪警告来源 - 添加更多调试信息

include("../../src/Constants_PNJL.jl")
include("../../src/relaxtime/ScatteringAmplitude.jl")

using .ScatteringAmplitude

println("="^70)
println("追踪k0²-u<0警告的真实来源")
println("="^70)

# 参数
T = 0.15
quark_params = (
    m = (u=0.3, d=0.3, s=0.5),
    μ = (u=0.25, d=0.25, s=0.25),
    A = (u=0.05, d=0.05, s=0.08)
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

s_qq = 20.0
t_qq = -0.2

println("\n## 测试1: 单独调用 uu_to_uu")
println("-"^70)
println("uu_to_uu: s=$s_qq, t=$t_qq")
m_u = 0.3
u_uu = 4*m_u^2 - s_qq - t_qq
println("u(uu) = $(round(u_uu, digits=2)) fm⁻²")

M_sq = scattering_amplitude_squared(
    :uu_to_uu, s_qq, t_qq, quark_params, thermo_params, K_coeffs
)
println("结果: |M|² = $(round(M_sq, sigdigits=4)) fm⁻⁴")

println("\n## 测试2: 单独调用 ss_to_ss")
println("-"^70)
println("ss_to_ss: s=$s_qq, t=$t_qq")
m_s = 0.5
u_ss = 4*m_s^2 - s_qq - t_qq
println("u(ss) = $(round(u_ss, digits=2)) fm⁻²")

try
    M_sq = scattering_amplitude_squared(
        :ss_to_ss, s_qq, t_qq, quark_params, thermo_params, K_coeffs
    )
    println("结果: |M|² = $(round(M_sq, sigdigits=4)) fm⁻⁴")
catch e
    println("✗ 失败: ", sprint(showerror, e)[1:200])
end

println("\n## 测试3: 连续调用 uu → ss")
println("-"^70)
println("先调用 uu_to_uu...")
M_sq_uu = scattering_amplitude_squared(
    :uu_to_uu, s_qq, t_qq, quark_params, thermo_params, K_coeffs
)
println("  |M|²(uu) = $(round(M_sq_uu, sigdigits=4)) fm⁻⁴")

println("\n再调用 ss_to_ss...")
try
    M_sq_ss = scattering_amplitude_squared(
        :ss_to_ss, s_qq, t_qq, quark_params, thermo_params, K_coeffs
    )
    println("  |M|²(ss) = $(round(M_sq_ss, sigdigits=4)) fm⁻⁴")
catch e
    println("  ✗ 失败: ", sprint(showerror, e)[1:200])
end

# 关键分析：检查u道的k0²-u值
println("\n## 关键分析: u道的k0²-u")
println("-"^70)

using .ScatteringAmplitude.TotalPropagator

for (name, process, m_quark) in [("uu→uu", :uu_to_uu, m_u), ("ss→ss", :ss_to_ss, m_s)]
    println("\n$name:")
    u = 4*m_quark^2 - s_qq - t_qq
    println("  u = $u fm⁻²")
    
    m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
    
    # 计算u道CMS
    cms_u = calculate_cms_momentum(process, s_qq, t_qq, :u, quark_params; u=u)
    k0_u = cms_u.k0
    k_u = cms_u.k
    
    println("  k0(u道) = $k0_u fm⁻¹")
    println("  k(u道) = $k_u fm⁻¹")
    println("  k0² = $(k0_u^2) fm⁻²")
    println("  k0² - u = $(k0_u^2 - u) fm⁻²")
    
    if k0_u^2 - u < 0
        println("  ⚠️  k0² - u < 0 !")
    else
        println("  ✓ k0² - u > 0")
    end
end

println("\n" * "="^70)
println("分析完成")
println("="^70)

