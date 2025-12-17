"""
# OneLoopIntegrals.jl

有限温度-密度单圈积分实现，包含 B0 (two-line integral) 的数值计算。

公式参考 `doc/formula/B0.md`。
"""
module OneLoopIntegrals

include("../integration/GaussLegendre.jl")
include("../Constants_PNJL.jl")
include("../QuarkDistribution.jl")
using .GaussLegendre: gauleg
using .PNJLQuarkDistributions: quark_distribution, antiquark_distribution,
    quark_distribution_integral, antiquark_distribution_integral
using .Constants_PNJL: Λ_inv_fm
using QuadGK: quadgk

export B0, A

# ----------------------------------------------------------------------------
# 基础工具函数
const EPS_K = 1.0e-9 # 三动量k大小为零的判定阈值
const EPS_SEGMENT = 1.0e-12 # 分母的最小值判定阈值
const DEFAULT_RTOL = 1.0e-3 # 积分相对误差默认值
const DEFAULT_ATOL = 0.0 # 积分绝对误差默认值
"""计算给定质量下的能量截断值"""
@inline @fastmath function energy_cutoff(m::Float64)
    m_pos = max(m, 0.0)
    return sqrt(m_pos * m_pos + Λ_inv_fm * Λ_inv_fm)
end

"""计算给定能量和质量对应的动量"""
@inline @fastmath function internal_momentum(E::Float64, m::Float64)
    m_pos = max(m, 0.0)
    return sqrt(E * E - m_pos * m_pos)
end

"""计算夸克有效分布函数的值（供 A 积分使用：分别用夸克/反夸克分布）"""
@inline function distribution_value(mode::Symbol, sign_flag::Symbol, E::Float64, μ::Float64,
    T::Float64, Φ::Float64, Φbar::Float64)
    @assert sign_flag === :plus || sign_flag === :minus "sign_flag must be :plus or :minus"
    if mode === :pnjl
        if sign_flag === :plus
            return quark_distribution(E, μ, T, Φ, Φbar)
        else
            return antiquark_distribution(E, μ, T, Φ, Φbar)
        end
    else
        throw(ArgumentError("暂不支持的 distribution: $mode"))
    end
end

"""B0 中使用的“±”分布（仅在正能量 E 上评估）。

等价目标：与 C++ `fd(pm*E, mu, ...)` 一致，其中
- `pm=+1` 对应 `f^+(E,mu)`
- `pm=-1` 对应 `f^+(-E,mu) = 1 - f^-(E,mu)`

同时保持 C++ 对 `mu<0` 的约定：`fd(E,mu<0)=f^-(E,|mu|)`。
"""
@inline function distribution_value_b0(sign_flag::Symbol, E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
    @assert sign_flag === :plus || sign_flag === :minus "sign_flag must be :plus or :minus"

    if sign_flag === :plus
        # fd(+E, μ)
        if μ >= 0.0
            return quark_distribution(E, μ, T, Φ, Φbar)
        else
            return antiquark_distribution(E, -μ, T, Φ, Φbar)
        end
    else
        # fd(-E, μ)
        if μ >= 0.0
            # f^+(-E,μ) = 1 - f^-(E,μ)
            return 1.0 - antiquark_distribution(E, μ, T, Φ, Φbar)
        else
            # fd(-E,μ<0)=f^-( -E, |μ| ) = 1 - f^+(E,|μ|)
            return 1.0 - quark_distribution(E, -μ, T, Φ, Φbar)
        end
    end
end

@inline function distribution_integral_b0(sign_flag::Symbol, E_min::Float64, E_max::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
    @assert sign_flag === :plus || sign_flag === :minus "sign_flag must be :plus or :minus"

    if sign_flag === :plus
        if μ >= 0.0
            return quark_distribution_integral(E_min, E_max, μ, T, Φ, Φbar)
        else
            return antiquark_distribution_integral(E_min, E_max, -μ, T, Φ, Φbar)
        end
    else
        interval_len = E_max - E_min
        if μ >= 0.0
            return interval_len - antiquark_distribution_integral(E_min, E_max, μ, T, Φ, Φbar)
        else
            return interval_len - quark_distribution_integral(E_min, E_max, -μ, T, Φ, Φbar)
        end
    end
end

"""计算夸克有效分布函数在给定能量区间的积分"""
@inline function distribution_integral(mode::Symbol, sign_flag::Symbol, E_min::Float64   , E_max::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
    @assert sign_flag === :plus || sign_flag === :minus "sign_flag must be :plus or :minus"
    if mode === :pnjl
        if sign_flag === :plus
            return quark_distribution_integral(E_min, E_max, μ, T, Φ, Φbar)
        else
            return antiquark_distribution_integral(E_min, E_max, μ, T, Φ, Φbar)
        end
    else
        throw(ArgumentError("暂不支持的 distribution: $mode"))
    end
end
# -----------------------------------------------------------------------------
# k=0 时的积分计算相关函数
"""k=0时的积分实部被积函数"""
@inline function real_integrand_k_zero(sign_flag::Symbol, λ::Float64, m::Float64, denominator_term::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, E::Float64)
    p = internal_momentum(E, m)
    dist = distribution_value_b0(sign_flag, E, μ, T, Φ, Φbar)
    denominator = λ * E + denominator_term
    return p * dist / denominator
end

"""k=0时的奇点计算函数
返回位于 (Emin, Emax) 内的 E 奇点值列表；若无奇点则返回空列表"""
@inline function singularity_k_zero(λ::Float64, Emin::Float64, Emax::Float64, denominator_term::Float64)::Vector{Float64}
    E0 = -denominator_term / λ
    if !isfinite(E0)
        return Float64[]
    end

    # 允许 E0 贴近端点（数值误差内）也被识别出来
    tol = 1e-8 * max(1.0, abs(E0))
    if (E0 > Emin - tol) && (E0 < Emax + tol)
        return [E0]
    else
        return Float64[]
    end
end

"""三动量大小k=0(小于EPS_K)时的 B0分量 积分计算"""
@inline function tilde_B0_k_zero(sign_flag::Symbol, λ::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64; rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)
    m_pos = max(m, 0.0)
    m_prime_pos = max(m_prime, 0.0)
    Emin = m_pos
    Emax = energy_cutoff(m_pos)
    denominator_term = (λ ^ 2 + m_pos ^ 2 - m_prime_pos ^ 2) / 2.0
    singularity = singularity_k_zero(λ, Emin, Emax, denominator_term)
    integrand_fun(E) = real_integrand_k_zero(sign_flag, λ, m_pos, denominator_term,
        μ, T, Φ, Φbar, E) # 闭包被积函数

    imag_part = 0.0
    if isempty(singularity) # 无奇点
        real_part, _ = quadgk(integrand_fun, Emin, Emax; rtol=rtol, atol=atol)
    else # 有奇点：在 E0 附近留出对称的微小间隙，避免直接采样到 1/0
        E0 = singularity[1]
        eps = max(1e-6, 1e-6 * abs(E0))
        lo = max(Emin, E0 - eps)
        hi = min(Emax, E0 + eps)

        if hi <= lo
            # 极端情况下间隙退化，回退到直接积分（可能再次触发奇点，但覆盖面更小）
            real_part, _ = quadgk(integrand_fun, Emin, Emax; rtol=rtol, atol=atol)
        else
            left, _ = quadgk(integrand_fun, Emin, lo; rtol=rtol, atol=atol)
            right, _ = quadgk(integrand_fun, hi, Emax; rtol=rtol, atol=atol)
            real_part = left + right
        end

        p0 = internal_momentum(E0, m_pos)
        imag_part = 2.0 * π * p0 * distribution_value_b0(sign_flag, E0, μ, T, Φ, Φbar)
    end

    return real_part * 2.0, imag_part / λ
end
# ----------------------------------------------------------------------------
# k>0 时的积分计算相关函数
"""k>0 时的积分实部被积函数"""
@inline @fastmath function real_integrand_k_positive(sign_flag::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, E::Float64)
    p = internal_momentum(E, m)
    dist = distribution_value_b0(sign_flag, E, μ, T, Φ, Φbar)
    common_part = (λ + E)^2 - m_prime^2
    numerator = common_part - (p - k)^2
    denominator = common_part - (p + k)^2
    ratio = abs(numerator / denominator)
    log_term = log(ratio)
    return dist * log_term
end

"""
k>0 时的奇点计算函数
根据二次项系数 a = λ² - k² 的符号，返回虚部积分区间的元组 (intervals, sign_type)
- intervals: Vector{Tuple{Float64, Float64}} 虚部积分区间列表
- sign_type: Symbol, 表示区间的类型 (:between, :outside, :none)

根据 B0_extra.md 文档:
- a > 0 (λ² > k²): 虚部积分区间为 [E1, E2]
- a < 0 (λ² < k²): 虚部积分区间为 [Emin, E1] ∪ [E2, Emax]
- a ≈ 0 (λ² ≈ k²): 退化为线性情况，特殊处理
"""
@inline function singularity_k_positive(λ::Float64, k::Float64, m::Float64, m_prime::Float64,
    Emin::Float64, Emax::Float64)::Tuple{Vector{Tuple{Float64, Float64}}, Symbol}
    # 解析解：
    # E_{1,2} = [λ(λ^2 - k^2 + m^2 - m'^2) ± k * sqrt((λ^2 - k^2 + m^2 - m'^2)^2 - 4 m^2 (k^2 - λ^2))] / [2 (k^2 - λ^2)]
    
    intervals = Tuple{Float64, Float64}[]
    
    a = λ ^ 2 - k ^ 2                # 二次项系数
    d0 = k ^ 2 - λ ^ 2               # 分母的一半系数 (未乘 2) = -a
    denom = 2.0 * d0                 # 实际分母 2(k^2 - λ^2) = -2a
    
    # 情况3: a ≈ 0 (λ² ≈ k²)，退化为线性情况
    if abs(a) < EPS_SEGMENT
        # 原方程退化为线性: |2λE + m² - m'²| < 2k√(E² - m²)
        # 根据 B0_extra.md 文档:
        # 有解条件: λ > 0 且 m' > m
        # 解的范围: E > (m'² - m²)/(4λ) + λm²/(m'² - m²)
        
        if λ > 0.0 && m_prime > m
            # 计算临界能量
            m2_diff = m_prime ^ 2 - m ^ 2
            E_crit = m2_diff / (4.0 * λ) + λ * m ^ 2 / m2_diff
            
            # 虚部积分区间为 [max(E_crit, Emin), Emax]
            if E_crit < Emax
                E_start = max(E_crit, Emin)
                push!(intervals, (E_start, Emax))
            end
        end
        # 如果 λ ≤ 0 或 m' ≤ m，则无解，返回空区间
        return (intervals, :none)
    end
    
    A = λ ^ 2 - k ^ 2 + m ^ 2 - m_prime ^ 2
    disc = A ^ 2 + 4.0 * m ^ 2 * d0    # 判别式
    
    # 判别式为负无实根
    if disc < 0.0
        return (intervals, :none)
    end
    
    sqrt_disc = sqrt(disc)
    E1 = (λ * A - k * sqrt_disc) / denom
    E2 = (λ * A + k * sqrt_disc) / denom
    
    # 保证 E1 <= E2
    if E1 > E2
        E1, E2 = E2, E1
    end
    
    # 两根过于接近，视为无奇点
    if E2 - E1 < EPS_SEGMENT
        return (intervals, :none)
    end
    
    # 情况1: a > 0 (λ² > k², 即 d0 < 0)
    # 虚部积分区间为 [E1, E2]
    if a > 0.0
        # 检查 [E1, E2] 与 [Emin, Emax] 的交集
        if E2 > Emin && E1 < Emax
            E1_clipped = max(E1, Emin)
            E2_clipped = min(E2, Emax)
            if E2_clipped > E1_clipped
                push!(intervals, (E1_clipped, E2_clipped))
            end
        end
        return (intervals, :between)
    else
        # 情况2: a < 0 (λ² < k², 即 d0 > 0)
        # 虚部积分区间为 [Emin, E1] ∪ [E2, Emax]
        # 区间1: [Emin, E1]
        if E1 > Emin
            E1_clipped = min(E1, Emax)
            if E1_clipped > Emin
                push!(intervals, (Emin, E1_clipped))
            end
        end
        # 区间2: [E2, Emax]
        if E2 < Emax
            E2_clipped = max(E2, Emin)
            if Emax > E2_clipped
                push!(intervals, (E2_clipped, Emax))
            end
        end
        return (intervals, :outside)
    end
end

"""k>0 时的 B0分量 积分计算"""
function tilde_B0_k_positive(sign_flag::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64; rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)
    m_pos = max(m, 0.0)
    m_prime_pos = max(m_prime, 0.0)
    Emin = m_pos
    Emax = energy_cutoff(m_pos)
    integrand_fun(E) = real_integrand_k_positive(sign_flag, λ, k, m_pos, m_prime_pos, μ, T, Φ, Φbar, E) # 闭包被积函数
    intervals, sign_type = singularity_k_positive(λ, k, m_pos, m_prime_pos, Emin, Emax)
    
    imag_part = 0.0
    real_part, _ = quadgk(integrand_fun, Emin, Emax; rtol=rtol, atol=atol)
    
    # 根据区间类型计算虚部
    if !isempty(intervals)
        for (E1, E2) in intervals
            imag_part += distribution_integral_b0(sign_flag, E1, E2, μ, T, Φ, Φbar)
        end
        imag_part *= π * sign(λ)
    end
    
    return real_part / k, imag_part / k
end

"""计算单个 ̃B0 分量的函数"""
@inline function tilde_B0(sign_flag::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64; rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)
    if abs(k) < EPS_K
        return tilde_B0_k_zero(sign_flag, λ, m, m_prime, μ, T, Φ, Φbar; rtol=rtol, atol=atol)
    else
        return tilde_B0_k_positive(sign_flag, λ, k, m, m_prime, μ, T, Φ, Φbar; rtol=rtol, atol=atol)
    end
end

"""
    B0(λ, k, m1, μ1, m2, μ2, T; Φ=0.0, Φbar=0.0, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL)
根据文档公式组合四个 ̃B0 项得到完整的 B₀ 积分。
λ = k0 + μ1 - μ2, 其中k0是传播子能量
"""
function B0(λ::Float64, k::Float64, m1::Float64, μ1::Float64, m2::Float64, μ2::Float64, T::Float64;
    Φ::Float64=0.0, Φbar::Float64=0.0, rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)
    tol_kwargs = (; rtol=rtol, atol=atol)
    term1 = tilde_B0(:plus, -λ, k, m1, m2, μ1, T, Φ, Φbar; tol_kwargs...)
    term2 = tilde_B0(:minus, λ, k, m1, m2, μ1, T, Φ, Φbar; tol_kwargs...)
    term3 = tilde_B0(:plus, λ, k, m2, m1, μ2, T, Φ, Φbar; tol_kwargs...)
    term4 = tilde_B0(:minus, -λ, k, m2, m1, μ2, T, Φ, Φbar; tol_kwargs...)

    real_part = term1[1] - term2[1] + term3[1] - term4[1]
    imag_part = term1[2] - term2[2] + term3[2] - term4[2]
    return real_part, imag_part
end

# ---------------------------------------------------------------------------
"""计算A的辅助函数,处理对常数项的积分"""
function const_integral_term_A(m::Float64)
    m_pos = max(m, 0.0)
    # m -> 0 极限：term2 -> 0，term1 = Λ^2
    if m_pos < 1e-14
        return (Λ_inv_fm^2) / 2.0
    end

    term1 = Λ_inv_fm * sqrt(Λ_inv_fm^2 + m_pos^2)
    term2 = m_pos^2 * log((Λ_inv_fm + sqrt(Λ_inv_fm^2 + m_pos^2)) / m_pos)
    return (term1 - term2) / 2.0
end

"""
    A(m, μ, T, Φ, Φbar, nodes_p, weights_p)
计算单线积分函数A,需要传入预生成的动量的积分节点与权重
计算中的常数1需要截断,而分布函数的部分不需要截断
常数项的积分可以直接计算，见 const_integral_term_A 函数
因此只需计算分布函数的积分，传入的节点积分上限设为 10.0 fm⁻¹ 即可收敛

推荐使用 GaussLegendre 模块中的 DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS (p ∈ [0, 10] fm⁻¹)
"""
function A(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
    nodes_p::Vector{Float64}, weights_p::Vector{Float64})
    integral = -const_integral_term_A(m) # 计算常数项积分部分
    @inbounds @simd for i in eachindex(nodes_p)
        node_p = nodes_p[i]
        weight_p = weights_p[i]
        E = sqrt(node_p^2 + m^2)
        dist_quark = distribution_value(:pnjl, :plus, E, μ, T, Φ, Φbar)
        dist_antiquark = distribution_value(:pnjl, :minus, E, μ, T, Φ, Φbar)
        integral += weight_p * node_p^2 / E * (dist_quark + dist_antiquark)
    end
    return 4.0 * integral
end

end # module OneLoopIntegrals
