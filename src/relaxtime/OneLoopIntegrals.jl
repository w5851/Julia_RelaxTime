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
const DEFAULT_RTOL = 1.0e-2 # 积分相对误差默认值
const DEFAULT_ATOL = 0.0 # 积分绝对误差默认值
"""计算给定质量下的能量截断值"""
@inline @fastmath function energy_cutoff(m::Float64)
    return sqrt(m * m + Λ_inv_fm * Λ_inv_fm)
end

"""计算给定能量和质量对应的动量"""
@inline @fastmath function internal_momentum(E::Float64, m::Float64)
    return sqrt(E * E - m * m)
end

"""计算夸克有效分布函数的值"""
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
    dist = distribution_value(:pnjl, sign_flag, E, μ, T, Φ, Φbar)
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
    # 仅在严格位于积分区间内部且不靠近端点时视为奇点
    if (E0 > Emin) && (E0 < Emax)
        return [E0]
    else
        return Float64[]
    end
end

"""三动量大小k=0(小于EPS_K)时的 B0分量 积分计算"""
@inline function tilde_B0_k_zero(sign_flag::Symbol, λ::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64; rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)
    Emin = m
    Emax = energy_cutoff(m)
    denominator_term = (λ ^ 2 + m ^ 2 - m_prime ^ 2) / 2.0
    singularity = singularity_k_zero(λ, Emin, Emax, denominator_term)
    integrand_fun(E) = real_integrand_k_zero(sign_flag, λ, m, denominator_term,
        μ, T, Φ, Φbar, E) # 闭包被积函数

    imag_part = 0.0
    if isempty(singularity) # 无奇点
        real_part, _ = quadgk(integrand_fun, Emin, Emax; rtol=rtol, atol=atol)
    else # 有奇点
        real_part, _ = quadgk(integrand_fun, Emin, singularity..., Emax; rtol=rtol, atol=atol)
        p0 = internal_momentum(singularity[1], m)
        imag_part = 2.0 * π * p0 * distribution_value(:pnjl, sign_flag, singularity[1], μ, T, Φ, Φbar)
    end

    return real_part * 2.0, imag_part / λ
end
# ----------------------------------------------------------------------------
# k>0 时的积分计算相关函数
"""k>0 时的积分实部被积函数"""
@inline @fastmath function real_integrand_k_positive(sign_flag::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, E::Float64)
    p = internal_momentum(E, m)
    dist = distribution_value(:pnjl, sign_flag, E, μ, T, Φ, Φbar)
    common_part = (λ + E)^2 - m_prime^2
    numerator = common_part - (p - k)^2
    denominator = common_part - (p + k)^2
    ratio = abs(numerator / denominator)
    log_term = log(ratio)
    return dist * log_term
end

"""k>0 时的奇点计算函数
返回位于 (Emin, Emax) 内的 E 奇点值列表；若无奇点则返回空列表"""
@inline function singularity_k_positive(λ::Float64, k::Float64, m::Float64, m_prime::Float64,
    Emin::Float64, Emax::Float64)::Vector{Float64}
    # 解析解：
    # E_{1,2} = [λ(λ^2 - k^2 + m^2 - m'^2) ± k * sqrt((λ^2 - k^2 + m^2 - m'^2)^2 - 4 m^2 (k^2 - λ^2))] / [2 (k^2 - λ^2)]
    # 仅当判别式非负且分母不为 0 时有实根；并只返回落在 (Emin, Emax) 内的点。
    singularities = Float64[]

    d0 = k ^ 2 - λ ^ 2               # 分母的一半系数 (未乘 2)
    denom = 2.0 * d0                   # 实际分母 2(k^2 - λ^2)
    # 分母过小或为 0，解析式不稳定，直接返回空
    if abs(denom) < EPS_SEGMENT
        return singularities
    end
    A = λ ^ 2 - k ^ 2 + m ^ 2 - m_prime ^ 2
    disc = A ^ 2 - 4.0 * m ^ 2 * d0    # 判别式
    # 判别式为负无实根
    if disc < 0.0
        return singularities
    end
    sqrt_disc = sqrt(disc)
    E1 = (λ * A - k * sqrt_disc) / denom
    E2 = (λ * A + k * sqrt_disc) / denom
    # 保证 E1 <= E2
    if E1 > E2
        E1, E2 = E2, E1
    end
    tol = EPS_SEGMENT
    # 仅收集严格位于积分区间内部的点，避免端点奇异
    if isfinite(E1) && (E1 > Emin + tol) && (E1 < Emax - tol)
        push!(singularities, E1)
    end
    # 判别式为 0 时两根重合，仅添加一次
    if (disc > 0.0) && isfinite(E2) && (E2 > Emin + tol) && (E2 < Emax - tol)
        push!(singularities, E2)
    end
    #sort!(singularities)
    return singularities
end

"""k>0 时的 B0分量 积分计算"""
function tilde_B0_k_positive(sign_flag::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64; rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)
    Emin = m
    Emax = energy_cutoff(m)
    integrand_fun(E) = real_integrand_k_positive(sign_flag, λ, k, m, m_prime, μ, T, Φ, Φbar, E) # 闭包被积函数
    singularities = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)
    imag_part = 0.0
    if isempty(singularities) # 无奇点
        real_part, _ = quadgk(integrand_fun, Emin, Emax; rtol=rtol, atol=atol)
    else # 有奇点
        real_part, _ = quadgk(integrand_fun, Emin, singularities..., Emax; rtol=rtol, atol=atol)
        if length(singularities) == 2 # 有两个奇点
            imag_part = π * sign(λ) * distribution_integral(:pnjl, sign_flag, singularities[1], singularities[2],
                μ, T, Φ, Φbar)
        end
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
    term1 = Λ_inv_fm * sqrt(Λ_inv_fm^2 + m^2)
    term2 = m^2 * log((Λ_inv_fm + sqrt(Λ_inv_fm^2 + m^2)) / m)
    return (term1 - term2) / 2.0
end

"""
    A(m, μ, T, Φ, Φbar, nodes_p, weights_p)
计算单线积分函数A,需要传入预生成的动量的积分节点与权重
计算中的常数1需要截断,而分布函数的部分不需要截断
常数项的积分可以直接计算，见 const_integral_term_A 函数
因此只需计算分布函数的积分，传入的节点不需要截断(取0-足够大,如20)
"""
function A(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
    nodes_p::Vector{Float64}, weights_p::Vector{Float64})
    integral = -const_integral_term_A(m) # 计算常数项积分部分
    @inbounds for i in eachindex(nodes_p)
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
