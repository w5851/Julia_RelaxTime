# 调试脚本：逐点调用 total_cross_section，打印失败的 s 和回溯
# 只读，不写任何数据文件

# 加载模块（通过 include 定义模块），使用相对路径，不激活环境以保持只读行为
# 注意：脚本通过命令行的 `--project=.` 使用工程环境
base_dir = joinpath(@__DIR__, "..", "..", "..")
include(joinpath(base_dir, "src", "relaxtime", "TotalCrossSection.jl"))
include(joinpath(base_dir, "src", "relaxtime", "ScatteringAmplitude.jl"))
include(joinpath(base_dir, "src", "relaxtime", "TotalPropagator.jl"))
include(joinpath(base_dir, "src", "relaxtime", "PolarizationCache.jl"))
include(joinpath(base_dir, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(base_dir, "src", "Constants_PNJL.jl"))
include(joinpath(base_dir, "src", "relaxtime", "EffectiveCouplings.jl"))

using .EffectiveCouplings

using .TotalCrossSection
using .ScatteringAmplitude
using .TotalPropagator
using .PolarizationCache
using .OneLoopIntegrals
using .Constants_PNJL

# 尝试从已有脚本里复用参数构造方式，若失败则使用保守的占位参数
# 这里尽量构造与主脚本相同的 quark_params / thermo_params / K_coeffs

function build_minimal_params()
    # 物理常数转换: 例子里 T 用 MeV->fm^{-1}
    T = 150.0 / 197.327
    Φ = 0.5
    Φbar = 0.5
    ξ = 0.0

    # 质量和化学势：使用示例值（仅用于触发相同代码路径）
    m_u = 300.0 / 197.327
    m_s = 500.0 / 197.327
    μ_u = 800.0 / 197.327 / 3.0  # 粗略分配（若主脚本使用不同，请替换）
    μ_s = μ_u

    # 计算 A 函数（OneLoopIntegrals.A）如果存在
    A_u = try
        OneLoopIntegrals.A(T, μ_u, m_u, Φ, Φbar)
    catch
        -1.0
    end
    A_s = try
        OneLoopIntegrals.A(T, μ_s, m_s, Φ, Φbar)
    catch
        -1.0
    end

    quark_params = (
        m = (u = m_u, d = m_u, s = m_s),
        μ = (u = μ_u, d = μ_u, s = μ_s),
        A = (u = A_u, d = A_u, s = A_s)
    )

    thermo_params = (T = T, Φ = Φ, Φbar = Φbar, ξ = ξ)

    # K_coeffs 构造：若模块存在计算函数则调用，否则使用简单占位
    K_coeffs = try
        G_fm2 = Constants_PNJL.G_GeV_inv2 / (197.327^2)
        K_fm5 = Constants_PNJL.K_GeV_inv5 / (197.327^5)
        G_u = EffectiveCouplings.calculate_G_from_A(A_u, m_u)
        G_s = EffectiveCouplings.calculate_G_from_A(A_s, m_s)
        EffectiveCouplings.calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
    catch
        # 提供包含 MesonPropagator 需要字段的占位 K_coeffs，避免字段缺失导致错误
        (K123_plus=1.0, K4567_plus=1.0, K123_minus=1.0, K4567_minus=1.0,
         K0_plus=1.0, K8_plus=1.0, K08_plus=1.0, det_K_plus=1.0,
         K0_minus=1.0, K8_minus=1.0, K08_minus=1.0, det_K_minus=1.0)
    end

    return quark_params, thermo_params, K_coeffs
end

# 可配置参数
process = :uu_to_uu
s_min = 0.1
s_max = 100.0
s_step = 0.5
n_points = 6  # 低点数以加速

quark_params, thermo_params, K_coeffs = build_minimal_params()

println("Starting debug scan for process=", process)
println("Will test s in [", s_min, ",", s_max, "] step=", s_step)

for s in s_min:s_step:s_max
    try
        σ = TotalCrossSection.total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
        if !isfinite(σ)
            println("s=", s, " -> non-finite σ: ", σ)
        end
    catch err
        println("--- FAILURE at s=", s, " ---")
        println("Error: ", err)
        println("Backtrace:")
        Base.show_backtrace(stdout, catch_backtrace())
        println("-------------------------------")
    end
end

println("Debug scan finished.")
