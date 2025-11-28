# 模拟测试中的连续调用

include("../../src/Constants_PNJL.jl")
include("../../src/relaxtime/ScatteringAmplitude.jl")

using .ScatteringAmplitude

println("="^70)
println("模拟测试中的连续调用（观察警告出现时机）")
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

qq_processes = [:uu_to_uu, :ss_to_ss, :ud_to_ud, :us_to_us]

println("\n参数: s=$s_qq, t=$t_qq")
println("\n开始连续调用...")
println("-"^70)

for (i, process) in enumerate(qq_processes)
    println("\n[$i] 调用 $process...")
    try
        M_sq = scattering_amplitude_squared(
            process, s_qq, t_qq, quark_params, thermo_params, K_coeffs
        )
        println("  ✓ |M|² = $(round(M_sq, sigdigits=4)) fm⁻⁴")
    catch e
        println("  ✗ 失败: ", split(sprint(showerror, e), '\n')[1])
    end
end

println("\n" * "="^70)
println("连续调用完成")
println("="^70)

