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

export B0

# ----------------------------------------------------------------------------
# 基础工具函数
const EPS_K = 1.0e-9 # 三动量为零的判定阈值
@inline @fastmath function energy_cutoff(m::Float64)
    return sqrt(m * m + Λ_inv_fm * Λ_inv_fm)
end

@inline @fastmath function internal_momentum(E::Float64, m::Float64)
    return sqrt(E * E - m * m)
end

@inline function distribution_value(mode::Symbol, sign::Symbol, E::Float64, μ::Float64,
    T::Float64, Φ::Float64, Φbar::Float64)
    @assert sign === :plus || sign === :minus "sign must be :plus or :minus"
    if mode === :pnjl
        if sign === :plus
            return quark_distribution(E, μ, T, Φ, Φbar)
        else
            return antiquark_distribution(E, μ, T, Φ, Φbar)
        end
    else
        throw(ArgumentError("暂不支持的 distribution: $mode"))
    end
end
# -----------------------------------------------------------------------------
@inline function real_integrand_k_zero(sign::Symbol, λ::Float64, m::Float64, mp::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, E::Float64)
    p = internal_momentum(E, m)
    dist = distribution_value(distribution, sign, E, μ, T, Φ, Φbar)
    denom = λ * λ + 2.0 * λ * E + m * m - mp * mp
    return 4.0 * p * dist / denom
end

"""三动量大小为0(小于EPS_K)时的 B0分量 积分计算"""
@inline function tilde_B0_k_zero(sign::Symbol, λ::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64; cutoff::Float64=DEFAULT_CUTOFF, distribution::Symbol=:pnjl,
    Φ::Float64=0.0, Φbar::Float64=0.0, nE::Int=128, nCpv::Int=256)
    Emin = m
    Emax = energy_cutoff(m, cutoff)
    real_part = real_part_k_zero(sign, λ, m, m_prime, μ, T, Φ, Φbar,
        distribution, Emin, Emax, nE, nCpv)
    imag_part = imag_part_k_zero(sign, λ, m, mp, μ, T, Φ, Φbar,
        distribution, Emin, Emax)
    return complex(real_part, imag_part)
end

"""三动量大小大于0(大于EPS_K)时的 B0分量 积分计算"""
@inline function tilde_B0_k_positive(sign::Symbol, λ::Float64, k::Float64, m::Float64, mp::Float64, μ::Float64, T::Float64;
    cutoff::Float64=DEFAULT_CUTOFF, distribution::Symbol=:pnjl,
    Φ::Float64=0.0, Φbar::Float64=0.0, nE::Int=128, nCpv::Int=256,
    nscan::Int=DEFAULT_SCAN_POINTS, tol::Float64=EPS_SEGMENT)
    Emin = m
    Emax = energy_cutoff(m, cutoff)
    real_part = real_part_k_positive(sign, λ, k, m, mp, μ, T, Φ, Φbar,
        distribution, Emin, Emax, cutoff, nscan, tol)
    imag_part = imag_part_k_positive(sign, λ, k, m, mp, μ, T, Φ, Φbar,
        distribution, Emin, Emax, nE, nscan, tol, cutoff)
    return complex(real_part, imag_part)
end


@inline function sign_or_zero(x::Float64)
    return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0)
end

@inline function clamp_log_argument(x::Float64)
    if !isfinite(x) || x <= 0.0
        return MIN_LOG_ARG
    elseif x > MAX_LOG_ARG
        return MAX_LOG_ARG
    else
        return x
    end
end

@inline function fermi_dirac(E::Float64, μ_eff::Float64, invT::Float64)
    x = (E - μ_eff) * invT
    if x > 50.0
        return exp(-x)
    elseif x < -50.0
        return 1.0 - exp(x)
    else
        return 1.0 / (exp(x) + 1.0)
    end
end



@inline function condition_value(E::Float64, λ::Float64, k::Float64, m::Float64, mp::Float64)
    p = internal_momentum(E, m)
    expr = λ * λ + 2.0 * λ * E - k * k + m * m - mp * mp
    return 2.0 * p * k - abs(expr)
end

function refine_boundary(a::Float64, b::Float64, cond, tol::Float64)
    fa = cond(a)
    fb = cond(b)
    if !isfinite(fa) || !isfinite(fb)
        return 0.5 * (a + b)
    end
    if abs(fa) < tol
        return a
    end
    if abs(fb) < tol
        return b
    end
    if fa * fb > 0.0
        return 0.5 * (a + b)
    end
    left, right = a, b
    f_left, f_right = fa, fb
    for _ in 1:MAX_BISECTION_STEPS
        mid = 0.5 * (left + right)
        f_mid = cond(mid)
        if !isfinite(f_mid)
            return mid
        end
        if abs(f_mid) < tol || abs(right - left) < tol
            return mid
        end
        if f_left * f_mid <= 0.0
            right = mid
            f_right = f_mid
        else
            left = mid
            f_left = f_mid
        end
    end
    return 0.5 * (left + right)
end

function find_active_segments(λ::Float64, k::Float64, m::Float64, mp::Float64,
    cutoff::Float64; nscan::Int=DEFAULT_SCAN_POINTS,
    tol::Float64=EPS_SEGMENT)
    segments = Tuple{Float64,Float64}[]
    if abs(k) < EPS_K
        return segments
    end
    Emin = m
    Emax = energy_cutoff(m, cutoff)
    nscan = max(nscan, 2)
    cond(E) = condition_value(E, λ, k, m, mp)
    grid = collect(range(Emin, Emax; length=nscan))
    values = Float64[cond(x) for x in grid]
    i = 1
    while i <= length(grid)
        if values[i] >= 0.0
            start_idx = i
            i += 1
            while i <= length(grid) && values[i] >= 0.0
                i += 1
            end
            end_idx = min(i - 1, length(grid))
            left = grid[start_idx]
            if start_idx > 1
                left = refine_boundary(grid[start_idx-1], grid[start_idx], cond, tol)
            end
            right = grid[end_idx]
            if i <= length(grid)
                right = refine_boundary(grid[end_idx], grid[i], cond, tol)
            end
            left = max(left, Emin)
            right = min(right, Emax)
            if right - left > tol
                push!(segments, (left, right))
            end
        else
            i += 1
        end
    end
    return segments
end

function unique_sorted_points(points::Vector{Float64}; tol::Float64=EPS_SEGMENT)
    if isempty(points)
        return Float64[]
    end
    sort!(points)
    filtered = Float64[]
    for value in points
        if isempty(filtered) || abs(value - filtered[end]) > tol
            push!(filtered, value)
        end
    end
    return filtered
end

function singular_points(λ::Float64, k::Float64, m::Float64, mp::Float64, cutoff::Float64;
    nscan::Int=DEFAULT_SCAN_POINTS, tol::Float64=EPS_SEGMENT)
    if abs(k) < EPS_K
        return Float64[]
    end
    segments = find_active_segments(λ, k, m, mp, cutoff; nscan=nscan, tol=tol)
    isempty(segments) && return Float64[]
    Emin = m
    Emax = energy_cutoff(m, cutoff)
    points = Float64[]
    for (left, right) in segments
        if left > Emin + tol && left < Emax - tol
            push!(points, left)
        end
        if right > Emin + tol && right < Emax - tol
            push!(points, right)
        end
    end
    return unique_sorted_points(points; tol=tol)
end

# -----------------------------------------------------------------------------
# 实部计算

function real_part_k_positive(sign::Symbol, λ::Float64, k::Float64, m::Float64, mp::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
    distribution::Symbol, Emin::Float64, Emax::Float64,
    cutoff::Float64, nscan::Int, tol::Float64)
    invk = 1.0 / k
    integrand(E) = begin
        p = internal_momentum(E, m)
        dist = distribution_value(distribution, sign, E, μ, T, Φ, Φbar)
        num = (λ + E)^2 - (p - k)^2 - mp^2
        den = (λ + E)^2 - (p + k)^2 - mp^2
        ratio = clamp_log_argument(abs(num / den))
        return dist * log(ratio) * invk
    end
    breaks = singular_points(λ, k, m, mp, cutoff; nscan=nscan, tol=tol)
    if isempty(breaks)
        val, _ = quadgk(integrand, Emin, Emax; rtol=B0_QUADGK_RTOL,
            atol=B0_QUADGK_ATOL, maxevals=B0_QUADGK_MAXEVALS)
        return val
    else
        val, _ = quadgk(integrand, Emin, breaks..., Emax; rtol=B0_QUADGK_RTOL,
            atol=B0_QUADGK_ATOL, maxevals=B0_QUADGK_MAXEVALS)
        return val
    end
end

function real_part_k_zero(sign::Symbol, λ::Float64, m::Float64, mp::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
    distribution::Symbol, Emin::Float64, Emax::Float64,
    nE::Int, nCpv::Int)
    if abs(λ) < EPS_LAMBDA
        δ = m * m - mp * mp
        if abs(δ) < EPS_DEN
            throw(ArgumentError("lambda -> 0 且 m ~ m' 时 B0 实部发散 (分母趋于 0)"))
        end
        nodes, weights = gauleg(Emin, Emax, nE)
        total = 0.0
        @inbounds for idx in eachindex(nodes)
            E = nodes[idx]
            p = internal_momentum(E, m)
            dist = distribution_value(distribution, sign, E, μ, T, Φ, Φbar)
            total += weights[idx] * p * dist
        end
        return 4.0 * total / δ
    end

    E0 = -(λ * λ + m * m - mp * mp) / (2.0 * λ)
    if (E0 > Emin + EPS_SEGMENT) && (E0 < Emax - EPS_SEGMENT)
        nodes, weights, n_actual = cpvi(Emin, Emax, E0, nCpv)
        sum_val = 0.0
        @inbounds for idx in 1:n_actual
            E = nodes[idx]
            p = internal_momentum(E, m)
            dist = distribution_value(distribution, sign, E, μ, T, Φ, Φbar)
            sum_val += weights[idx] * (p * dist) / (E - E0)
        end
        return (2.0 / λ) * sum_val
    else
        nodes, weights = gauleg(Emin, Emax, nE)
        total = 0.0
        @inbounds for idx in eachindex(nodes)
            E = nodes[idx]
            p = internal_momentum(E, m)
            dist = distribution_value(distribution, sign, E, μ, T, Φ, Φbar)
            denom = λ * λ + 2.0 * λ * E + m * m - mp * mp
            total += weights[idx] * 4.0 * p * dist / denom
        end
        return total
    end
end

# -----------------------------------------------------------------------------
# 虚部计算

function imag_part_k_zero(sign::Symbol, λ::Float64, m::Float64, mp::Float64, μ::Float64,
    T::Float64, Φ::Float64, Φbar::Float64, distribution::Symbol,
    Emin::Float64, Emax::Float64)
    if abs(λ) < EPS_LAMBDA
        return 0.0
    end
    E0 = -(λ * λ + m * m - mp * mp) / (2.0 * λ)
    if E0 <= Emin + EPS_SEGMENT || E0 >= Emax - EPS_SEGMENT
        return 0.0
    end
    p0_sq = E0 * E0 - m * m
    if p0_sq <= 0.0
        return 0.0
    end
    p0 = sqrt(p0_sq)
    dist = distribution_value(distribution, sign, E0, μ, T, Φ, Φbar)
    return (2.0 * pi / λ) * p0 * dist
end

function imag_part_k_positive(sign::Symbol, λ::Float64, k::Float64, m::Float64, mp::Float64,
    μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
    distribution::Symbol, Emin::Float64, Emax::Float64,
    nE::Int, nscan::Int, tol::Float64, cutoff::Float64)
    segments = find_active_segments(λ, k, m, mp, cutoff; nscan=nscan, tol=tol)
    if isempty(segments)
        return 0.0
    end
    total_length = Emax - Emin
    total = 0.0
    @inbounds for (a, b) in segments
        seg_length = max(b - a, tol)
        base_nodes = max(16, Int(round(nE * seg_length / max(total_length, tol))))
        n_nodes = min(MAX_SEGMENT_NODES, base_nodes)
        nodes, weights = gauleg(a, b, n_nodes)
        for idx in eachindex(nodes)
            E = nodes[idx]
            dist = distribution_value(distribution, sign, E, μ, T, Φ, Φbar)
            total += weights[idx] * dist
        end
    end
    return (pi * sign_or_zero(λ) / k) * total
end

# -----------------------------------------------------------------------------
# 组合计算

function tilde_B0(sign::Symbol, λ::Float64, k::Float64, m::Float64, mp::Float64,
    μ::Float64, T::Float64; cutoff::Float64=DEFAULT_CUTOFF,
    distribution::Symbol=:pnjl, Φ::Float64=0.0, Φbar::Float64=0.0,
    nE::Int=128, nCpv::Int=256, nscan::Int=DEFAULT_SCAN_POINTS,
    tol::Float64=EPS_SEGMENT)
    Emin = m
    Emax = energy_cutoff(m, cutoff)
    real_part = if abs(k) < EPS_K
        real_part_k_zero(sign, λ, m, mp, μ, T, Φ, Φbar,
            distribution, Emin, Emax, nE, nCpv)
    else
        real_part_k_positive(sign, λ, k, m, mp, μ, T, Φ, Φbar,
            distribution, Emin, Emax, cutoff, nscan, tol)
    end

    imag_part = if abs(k) < EPS_K
        imag_part_k_zero(sign, λ, m, mp, μ, T, Φ, Φbar,
            distribution, Emin, Emax)
    else
        imag_part_k_positive(sign, λ, k, m, mp, μ, T, Φ, Φbar,
            distribution, Emin, Emax, nE, nscan, tol, cutoff)
    end
    return complex(real_part, imag_part)
end

"""
	B0(λ, k, m1, μ1, m2, μ2, T; cutoff=3.05, distribution=:pnjl,
	   Φ=0.0, Φbar=0.0, nE=128, nCpv=256, nscan=1024, tol=1e-8)

计算有限温度、有限密度场论中的两条传播子单圈积分 B0。

# 参数
- `λ`: 组合能量参数 λ = p₀ + μ₁ - μ₂ （单位 fm⁻¹）
- `k`: 外部三动量模 |k|（单位 fm⁻¹）
- `m1`, `m2`: 两条传播线的质量（单位 fm⁻¹）
- `μ1`, `μ2`: 对应化学势（单位 fm⁻¹）
- `T`: 温度（单位 fm⁻¹，必须 > 0）

# 关键字参数
- `cutoff`: 三动量截断 Λ（默认 3.05 fm⁻¹）
- `distribution`: `:pnjl` 使用 PNJL 有效分布；`:fermi` 使用费米-狄拉克分布
- `Φ`, `Φbar`: PNJL 模型的 Polyakov 环参数
- `nE`: 实部积分的高斯节点数
- `nCpv`: k = 0 时柯西主值积分的节点数
- `nscan`: 寻找虚部积分区间的扫描点数
- `tol`: 数值容差，控制分段与求根精度

# 返回
- `ComplexF64`: B0 的数值结果，实部对应质谱平移，虚部对应衰减宽度
"""
function B0(λ::Float64, k::Float64, m1::Float64, μ1::Float64,
    m2::Float64, μ2::Float64, T::Float64;
    cutoff::Float64=DEFAULT_CUTOFF, distribution::Symbol=:pnjl,
    Φ::Float64=0.0, Φbar::Float64=0.0, nE::Int=128, nCpv::Int=256,
    nscan::Int=DEFAULT_SCAN_POINTS, tol::Float64=EPS_SEGMENT)
    if T <= 0.0
        throw(ArgumentError("温度 T 必须为正"))
    end
    if m1 <= 0.0 || m2 <= 0.0
        throw(ArgumentError("质量 m₁、m₂ 必须大于 0"))
    end
    if cutoff <= 0.0
        throw(ArgumentError("截断 cutoff 必须大于 0"))
    end
    if nE < 16
        throw(ArgumentError("nE 至少为 16，以保证积分精度"))
    end
    if nCpv < 32
        throw(ArgumentError("nCpv 至少为 32"))
    end
    if nscan < 16
        throw(ArgumentError("nscan 至少为 16"))
    end
    if distribution ≠ :pnjl && distribution ≠ :fermi
        throw(ArgumentError("distribution 仅支持 :pnjl 或 :fermi"))
    end

    kwargs = (; cutoff=cutoff, distribution=distribution, Φ=Φ, Φbar=Φbar,
        nE=nE, nCpv=nCpv, nscan=nscan, tol=tol)

    term1 = tilde_B0(:plus, -λ, k, m1, m2, μ1, T; kwargs...)
    term2 = tilde_B0(:minus, λ, k, m1, m2, μ1, T; kwargs...)
    term3 = tilde_B0(:plus, λ, k, m2, m1, μ2, T; kwargs...)
    term4 = tilde_B0(:minus, -λ, k, m2, m1, μ2, T; kwargs...)

    return term1 - term2 + term3 - term4
end

end # module OneLoopIntegrals
