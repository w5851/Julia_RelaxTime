"""
PNJL 相变分析模块

提供 S 形曲线检测和 Maxwell 等面积构造功能。

主要功能：
- `detect_s_shape`: 检测 μ(ρ) 曲线的 S 形特征，定位 spinodal 点
- `maxwell_construction`: Maxwell 等面积构造，计算相变点

物理背景：
- 一阶相变区域的 μ(ρ) 曲线呈 S 形
- spinodal 点是 dμ/dρ = 0 的位置（亚稳态边界）
- Maxwell 构造确定两相共存的化学势
"""
module PhaseTransition

export SShapeResult, detect_s_shape
export MaxwellResult, maxwell_construction
export group_curves_by_temperature

# ============================================================================
# 常量
# ============================================================================

const EPS_SLOPE = 0.0
const DEFAULT_AREA_TOL = 1e-4
const DEFAULT_MAX_ITER = 60
const DEFAULT_CANDIDATE_STEPS = 64
const MAX_CANDIDATE_STEPS = 1024
const BRACKET_SHRINK_REL = 1e-3
const BRACKET_SHRINK_ABS = 1e-3

# ============================================================================
# S 形检测
# ============================================================================

"""
S 形检测结果

字段：
- `has_s_shape`: 是否检测到 S 形
- `mu_spinodal_hadron`: μ(ρ) 曲线的局部极大值（强子相侧 spinodal 的 μ）
- `mu_spinodal_quark`: μ(ρ) 曲线的局部极小值（夸克相侧 spinodal 的 μ）
- `rho_spinodal_hadron`: 强子相侧 spinodal 的 ρ（极大值点）
- `rho_spinodal_quark`: 夸克相侧 spinodal 的 ρ（极小值点）
- `derivative_sign_changes`: 导数符号变化次数
"""
struct SShapeResult
    has_s_shape::Bool
    mu_spinodal_hadron::Union{Nothing, Float64}
    mu_spinodal_quark::Union{Nothing, Float64}
    rho_spinodal_hadron::Union{Nothing, Float64}
    rho_spinodal_quark::Union{Nothing, Float64}
    derivative_sign_changes::Int
end

SShapeResult() = SShapeResult(false, nothing, nothing, nothing, nothing, 0)

"""按 ρ 升序排列曲线数据"""
function _sort_curve_by_rho(mu_vals::AbstractVector, rho_vals::AbstractVector)
    n = min(length(mu_vals), length(rho_vals))
    order = sortperm(rho_vals[1:n])
    return Float64.(mu_vals[order]), Float64.(rho_vals[order])
end

@inline _slope_sign(value) = abs(value) <= EPS_SLOPE ? 0 : (value > 0 ? 1 : -1)

"""
    _refine_extremum(mu_sorted, rho_sorted, idx; is_maximum=true)

使用二次插值在极值点附近细化，提高 spinodal 点的精度。
"""
function _refine_extremum(mu_sorted::Vector{Float64}, rho_sorted::Vector{Float64}, 
                          idx::Int; is_maximum::Bool=true)
    n = length(rho_sorted)
    
    if idx < 2 || idx >= n
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end
    
    i1, i2, i3 = idx - 1, idx, idx + 1
    if i3 > n
        i1, i2, i3 = idx - 2, idx - 1, idx
    end
    if i1 < 1
        i1, i2, i3 = 1, 2, 3
    end
    
    r1, r2, r3 = rho_sorted[i1], rho_sorted[i2], rho_sorted[i3]
    m1, m2, m3 = mu_sorted[i1], mu_sorted[i2], mu_sorted[i3]
    
    denom = (r1 - r2) * (r1 - r3) * (r2 - r3)
    if abs(denom) < 1e-15
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end
    
    a = (r3 * (m2 - m1) + r2 * (m1 - m3) + r1 * (m3 - m2)) / denom
    b = (r3^2 * (m1 - m2) + r2^2 * (m3 - m1) + r1^2 * (m2 - m3)) / denom
    c = (r2 * r3 * (r2 - r3) * m1 + r3 * r1 * (r3 - r1) * m2 + r1 * r2 * (r1 - r2) * m3) / denom
    
    if (is_maximum && a >= 0) || (!is_maximum && a <= 0)
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end
    
    rho_ext = -b / (2 * a)
    
    rho_min = min(r1, r2, r3)
    rho_max = max(r1, r2, r3)
    if rho_ext < rho_min || rho_ext > rho_max
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end
    
    mu_ext = a * rho_ext^2 + b * rho_ext + c
    
    return rho_ext, mu_ext
end

"""
    detect_s_shape(mu_vals, rho_vals; eps=EPS_SLOPE, min_points=5)

检测 μ(ρ) 曲线是否呈 S 形（导数 dμ/dρ 符号变化：正 → 负 → 正）。

# 参数
- `mu_vals`: 化学势数组
- `rho_vals`: 密度数组
- `eps`: 斜率判断阈值
- `min_points`: 最小数据点数

# 返回
`SShapeResult` 结构，包含 spinodal 点信息。

# 物理意义
- S 形曲线表示一阶相变区域
- `spinodal_hadron`: μ(ρ) 的局部极大值点，对应强子相的亚稳态边界
- `spinodal_quark`: μ(ρ) 的局部极小值点，对应夸克相的亚稳态边界
"""
function detect_s_shape(mu_vals::AbstractVector, rho_vals::AbstractVector; 
                        eps::Real=EPS_SLOPE, min_points::Int=5)
    n = min(length(mu_vals), length(rho_vals))
    n >= min_points || return SShapeResult()
    
    mu_sorted, rho_sorted = _sort_curve_by_rho(mu_vals, rho_vals)
    
    slopes = Float64[]
    slope_indices = Int[]
    for i in 1:(length(rho_sorted) - 1)
        drho = rho_sorted[i + 1] - rho_sorted[i]
        abs(drho) <= eps && continue
        push!(slopes, (mu_sorted[i + 1] - mu_sorted[i]) / drho)
        push!(slope_indices, i)
    end
    isempty(slopes) && return SShapeResult()

    signs = Int[]
    sign_to_slope_idx = Int[]
    last_sign = 0
    for (i, slope) in enumerate(slopes)
        sign = _slope_sign(slope)
        sign == 0 && continue
        if sign != last_sign
            push!(signs, sign)
            push!(sign_to_slope_idx, i)
            last_sign = sign
        end
    end
    
    length(signs) < 3 && return SShapeResult(false, nothing, nothing, nothing, nothing, max(0, length(signs) - 1))

    max_idx = nothing
    min_idx = nothing
    
    for i in 1:(length(signs) - 1)
        if signs[i] == 1 && signs[i + 1] == -1 && max_idx === nothing
            max_idx = i + 1
        elseif signs[i] == -1 && signs[i + 1] == 1 && max_idx !== nothing && min_idx === nothing
            min_idx = i + 1
        end
    end
    
    if max_idx === nothing || min_idx === nothing
        return SShapeResult(false, nothing, nothing, nothing, nothing, length(signs) - 1)
    end

    slope_idx_max = sign_to_slope_idx[max_idx]
    slope_idx_min = sign_to_slope_idx[min_idx]
    
    data_idx_max = slope_indices[slope_idx_max]
    data_idx_min = slope_indices[slope_idx_min]
    
    rho_hadron, mu_hadron = _refine_extremum(mu_sorted, rho_sorted, data_idx_max; is_maximum=true)
    rho_quark, mu_quark = _refine_extremum(mu_sorted, rho_sorted, data_idx_min; is_maximum=false)

    return SShapeResult(true, mu_hadron, mu_quark, rho_hadron, rho_quark, length(signs) - 1)
end

# ============================================================================
# Maxwell 构造
# ============================================================================

"""
Maxwell 等面积构造结果

字段：
- `converged`: 是否收敛
- `mu_transition`: 相变化学势 (MeV)
- `rho_hadron`: 强子相密度 (ρ/ρ₀)
- `rho_quark`: 夸克相密度 (ρ/ρ₀)
- `area_residual`: 面积残差
- `iterations`: 迭代次数
- `details`: 详细信息字典
"""
struct MaxwellResult
    converged::Bool
    mu_transition::Union{Nothing, Float64}
    rho_hadron::Union{Nothing, Float64}
    rho_quark::Union{Nothing, Float64}
    area_residual::Union{Nothing, Float64}
    iterations::Int
    details::Dict{Symbol, Any}
end

MaxwellResult() = MaxwellResult(false, nothing, nothing, nothing, nothing, 0, Dict{Symbol, Any}())

"""按 ρ 排序并准备曲线数据"""
function _prepare_curve(mu_vals::AbstractVector, rho_vals::AbstractVector)
    n = min(length(mu_vals), length(rho_vals))
    pairs = Vector{Tuple{Float64, Float64}}()
    sizehint!(pairs, n)
    for i in 1:n
        mu = Float64(mu_vals[i])
        rho = Float64(rho_vals[i])
        (isfinite(mu) && isfinite(rho)) || continue
        push!(pairs, (rho, mu))
    end
    sort!(pairs; by = first)
    rho_sorted = Float64[first(p) for p in pairs]
    mu_sorted = Float64[last(p) for p in pairs]
    return rho_sorted, mu_sorted
end

"""从 S 形结果获取 μ 区间"""
function _mu_bracket(hint::SShapeResult)
    mu_max = hint.mu_spinodal_hadron
    mu_min = hint.mu_spinodal_quark
    if mu_max === nothing || mu_min === nothing
        return nothing
    end
    if mu_min < mu_max
        return (mu_min, mu_max)
    else
        return (mu_max, mu_min)
    end
end

"""收缩区间边界"""
function _shrink_bracket(mu_lo::Float64, mu_hi::Float64)
    width = mu_hi - mu_lo
    width > 0 || return nothing
    δ = max(width * BRACKET_SHRINK_REL, BRACKET_SHRINK_ABS)
    if mu_lo + δ >= mu_hi - δ
        δ = width / 4
        mu_lo + δ < mu_hi - δ || return nothing
    end
    return (mu_lo + δ, mu_hi - δ)
end

"""计算面积差（等面积法则）"""
function _area_difference(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64})
    rho_left, rho_right = _find_outer_intersections(mu0, rho_vals, mu_vals)
    if rho_left === nothing || rho_right === nothing || rho_right - rho_left <= 1e-9
        return nothing
    end
    return _integrate_difference(rho_vals, mu_vals, rho_left, rho_right, mu0)
end

"""找到 μ = μ0 与曲线的最外侧交点"""
function _find_outer_intersections(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64}; atol::Real=1e-9)
    n = length(rho_vals)
    n >= 2 || return nothing, nothing
    left = nothing
    right = nothing
    for i in 1:(n - 1)
        r1 = rho_vals[i]
        r2 = rho_vals[i + 1]
        f1 = mu_vals[i] - mu0
        f2 = mu_vals[i + 1] - mu0
        if abs(f1) < atol
            left === nothing && (left = r1)
            right = r1
        end
        if f1 == f2
            continue
        end
        if f1 * f2 < 0
            α = f1 / (f1 - f2)
            crossing = r1 + α * (r2 - r1)
            left === nothing && (left = crossing)
            right = crossing
        elseif abs(f2) < atol
            left === nothing && (left = r2)
            right = r2
        end
    end
    return left, right
end

"""梯形积分计算面积差"""
function _integrate_difference(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        rho_left::Float64, rho_right::Float64, mu0::Float64)
    total = 0.0
    for i in 1:(length(rho_vals) - 1)
        r1 = rho_vals[i]
        r2 = rho_vals[i + 1]
        if r2 <= rho_left || r1 >= rho_right
            continue
        end
        left = max(r1, rho_left)
        right = min(r2, rho_right)
        if right <= left
            continue
        end
        m1 = mu_vals[i]
        m2 = mu_vals[i + 1]
        mu_left = _interp(r1, r2, m1, m2, left)
        mu_right = _interp(r1, r2, m1, m2, right)
        f_left = mu_left - mu0
        f_right = mu_right - mu0
        total += 0.5 * (f_left + f_right) * (right - left)
    end
    return total
end

@inline function _interp(r1::Float64, r2::Float64, m1::Float64, m2::Float64, target::Float64)
    r2 == r1 && return m1
    t = (target - r1) / (r2 - r1)
    return m1 + t * (m2 - m1)
end

"""在区间内搜索符号变化点"""
function _find_mu_bracket(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        mu_lo::Float64, mu_hi::Float64, steps::Int, tol_area::Real)
    attempt = max(steps, 3)
    while attempt <= MAX_CANDIDATE_STEPS
        prev_mu = nothing
        prev_area = nothing
        for mu0 in range(mu_lo, mu_hi; length=attempt)
            area = _area_difference(mu0, rho_vals, mu_vals)
            area === nothing && continue
            if abs(area) <= tol_area
                return (mu0, mu0, area, area)
            end
            if prev_area !== nothing && area * prev_area < 0
                return (prev_mu, mu0, prev_area, area)
            end
            prev_mu = mu0
            prev_area = area
        end
        attempt *= 2
    end
    return nothing
end

"""二分法求解等面积点"""
function _bisection_solve(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        mu_a::Float64, mu_b::Float64, area_a::Real, area_b::Real,
        tol_area::Real, max_iter::Int)

    if mu_a == mu_b
        return mu_a, area_a, 0
    end

    a, b = mu_a, mu_b
    fa, fb = area_a, area_b
    fa * fb <= 0 || return nothing, nothing, 0

    for iter in 1:max_iter
        mid = 0.5 * (a + b)
        area_mid = _area_difference(mid, rho_vals, mu_vals)
        area_mid === nothing && return nothing, nothing, iter
        if abs(area_mid) <= tol_area
            return mid, area_mid, iter
        end
        if fa * area_mid < 0
            b = mid
            fb = area_mid
        else
            a = mid
            fa = area_mid
        end
    end
    mid = 0.5 * (a + b)
    area_mid = _area_difference(mid, rho_vals, mu_vals)
    area_mid === nothing && return nothing, nothing, max_iter
    return mid, area_mid, max_iter
end

"""创建失败结果"""
function _failure(reason::AbstractString; kwargs...)
    details = Dict{Symbol, Any}(:reason => reason)
    for (k, v) in kwargs
        details[k] = v
    end
    return MaxwellResult(false, nothing, nothing, nothing, nothing, 0, details)
end

"""
    maxwell_construction(mu_vals, rho_vals; kwargs...)

Maxwell 等面积构造，计算一阶相变的共存化学势和共存密度。

# 参数
- `mu_vals`: 化学势数组
- `rho_vals`: 密度数组
- `min_samples`: 最小样本数 (默认 12)
- `detect_min_points`: S 形检测最小点数 (默认 6)
- `detect_eps`: S 形检测斜率阈值 (默认 1e-6)
- `candidate_steps`: 初始搜索步数 (默认 64)
- `max_iter`: 最大迭代次数 (默认 60)
- `tol_area`: 面积容差 (默认 1e-4)
- `spinodal_hint`: 预计算的 S 形结果 (可选)

# 返回
`MaxwellResult` 结构，包含相变点信息。

# 算法
1. 检测 S 形曲线，获取 spinodal 点
2. 从 spinodal 点估计 μ 搜索区间
3. 在区间内搜索面积差符号变化
4. 二分法精确求解等面积点
"""
function maxwell_construction(mu_vals::AbstractVector, rho_vals::AbstractVector;
        min_samples::Int=12, detect_min_points::Int=6, detect_eps::Real=1e-6,
        candidate_steps::Int=DEFAULT_CANDIDATE_STEPS,
        max_iter::Int=DEFAULT_MAX_ITER, tol_area::Real=DEFAULT_AREA_TOL,
        spinodal_hint::Union{Nothing, SShapeResult}=nothing)

    rho_sorted, mu_sorted = _prepare_curve(mu_vals, rho_vals)
    if length(rho_sorted) < min_samples
        return _failure("insufficient_points"; count=length(rho_sorted))
    end

    hint = isnothing(spinodal_hint) ?
        detect_s_shape(mu_vals, rho_vals; eps=detect_eps, min_points=detect_min_points) :
        spinodal_hint
    hint.has_s_shape || return _failure("no_s_shape")

    mu_bracket = _mu_bracket(hint)
    mu_bracket === nothing && return _failure("invalid_mu_bracket")
    mu_lo, mu_hi = mu_bracket
    
    tightened = _shrink_bracket(mu_lo, mu_hi)
    tightened === nothing && return _failure("degenerate_bracket"; bracket=(mu_lo, mu_hi))
    mu_lo, mu_hi = tightened

    bracket = _find_mu_bracket(rho_sorted, mu_sorted, mu_lo, mu_hi, candidate_steps, tol_area)
    bracket === nothing && return _failure("no_sign_change"; bracket=(mu_lo, mu_hi))
    mu_a, mu_b, area_a, area_b = bracket

    mu_root, area_root, iterations = _bisection_solve(rho_sorted, mu_sorted,
        mu_a, mu_b, area_a, area_b, tol_area, max_iter)
    mu_root === nothing && return _failure("bisection_failed"; bracket=(mu_a, mu_b))

    rho_left, rho_right = _find_outer_intersections(mu_root, rho_sorted, mu_sorted)
    if rho_left === nothing || rho_right === nothing || !(rho_left < rho_right)
        return _failure("no_crossings"; mu_transition=mu_root, bracket=(mu_a, mu_b))
    end

    details = Dict(
        :mu_bracket => (mu_a, mu_b),
        :rho_interval => (rho_left, rho_right),
        :spinodal_hint => (hint.mu_spinodal_hadron, hint.mu_spinodal_quark),
    )
    return MaxwellResult(true, mu_root, rho_left, rho_right, abs(area_root), iterations, details)
end

# ============================================================================
# 工具函数
# ============================================================================

"""
    group_curves_by_temperature(rows; xi=0.0, tol=1e-6)

按温度分组 CSV 数据行，用于下游分析。

# 参数
- `rows`: CSV 行的迭代器（每行是 Dict）
- `xi`: 要筛选的各向异性参数
- `tol`: xi 匹配容差

# 返回
`Dict{Float64, Vector{Tuple{Float64, Float64}}}`: 温度 → (μ, ρ) 对列表
"""
function group_curves_by_temperature(rows; xi::Real=0.0, tol::Real=1e-6)
    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    for row in rows
        T = try
            parse(Float64, row["T_MeV"])
        catch
            continue
        end
        xi_val = haskey(row, "xi") ? try
                parse(Float64, row["xi"])
            catch
                NaN
            end : NaN
        if isnan(xi_val) || abs(xi_val - xi) > tol
            continue
        end
        rho = try
            parse(Float64, row["rho"])
        catch
            continue
        end
        mu = haskey(row, "mu_avg_MeV") ? try
                parse(Float64, row["mu_avg_MeV"])
            catch
                NaN
            end : NaN
        if isnan(mu)
            mu = haskey(row, "mu_MeV") ? try
                    parse(Float64, row["mu_MeV"])
                catch
                    NaN
                end : NaN
        end
        if isnan(mu) || !isfinite(rho)
            continue
        end
        bucket = get!(grouped, T) do
            Vector{Tuple{Float64, Float64}}()
        end
        push!(bucket, (mu, rho))
    end
    return grouped
end

end # module
