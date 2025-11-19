# 验证ss→ss散射的推荐参数

include("../src/Constants_PNJL.jl")
include("../src/relaxtime/ScatteringAmplitude.jl")

using .ScatteringAmplitude

println("="^70)
println("验证ss→ss散射的推荐参数")
println("="^70)

# 基本参数
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

# 简化K系数
K_coeffs = (
    K0_plus=10.0, K0_minus=10.0,
    K123_plus=5.0, K123_minus=5.0,
    K4567_plus=5.0, K4567_minus=5.0,
    K8_plus=10.0, K8_minus=10.0,
    K08_plus=0.0, K08_minus=0.0,
    det_K_plus=500.0, det_K_minus=500.0
)

println("\n## 测试不同参数组合")
println("-"^70)

test_cases = [
    ("失败参数", 20.0, -0.2),   # 当前测试中失败的
    ("推荐1", 5.0, -1.0),        # 推荐参数
    ("推荐2", 8.0, -0.5),        # 推荐参数
    ("推荐3", 10.0, -0.8),       # 推荐参数
]

for (name, s, t) in test_cases
    u = 4*m_s^2 - s - t
    println("\n[$name] s=$s, t=$t, u=$(round(u, digits=2)) fm⁻²")
    
    try
        M_sq = scattering_amplitude_squared(
            :ss_to_ss, s, t, quark_params, thermo_params, K_coeffs
        )
        println("  ✓ 成功！ |M|² = $(round(M_sq, digits=2)) fm⁻⁴")
    catch e
        println("  ✗ 失败: ", sprint(showerror, e)[1:min(100, end)])
    end
end

println("\n" * "="^70)
println("验证完成")
println("="^70)
