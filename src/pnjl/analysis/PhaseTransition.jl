"""
PNJL 相变分析模块

提供 S 形曲线检测、Maxwell 等面积构造和 Crossover 检测功能。

主要功能：
- `detect_s_shape`: 检测 μ(ρ) 曲线的 S 形特征，定位 spinodal 点
- `maxwell_construction`: Maxwell 等面积构造，计算相变点
- `detect_crossover`: 检测 crossover 温度（峰值法或拐点法）
- `scan_crossover_line`: 扫描 crossover 线

物理背景：
- 一阶相变区域的 μ(ρ) 曲线呈 S 形
- spinodal 点是 dμ/dρ = 0 的位置（亚稳态边界）
- Maxwell 构造确定两相共存的化学势
- Crossover 区域（T > T_CEP）通过序参量导数的峰值或拐点定义
"""
module PhaseTransition

export SShapeResult, detect_s_shape
export MaxwellResult, maxwell_construction
export group_curves_by_temperature
export CrossoverResult, detect_crossover, scan_crossover_line

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

"""使用二次插值在极值点附近细化"""
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

"""检测 μ(ρ) 曲线是否呈 S 形"""
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
    return Float64[first(p) for p in pairs], Float64[last(p) for p in pairs]
end

function _mu_bracket(hint::SShapeResult)
    mu_max = hint.mu_spinodal_hadron
    mu_min = hint.mu_spinodal_quark
    (mu_max === nothing || mu_min === nothing) && return nothing
    return mu_min < mu_max ? (mu_min, mu_max) : (mu_max, mu_min)
end

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

function _find_outer_intersections(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64}; atol::Real=1e-9)
    n = length(rho_vals)
    n >= 2 || return nothing, nothing
    left, right = nothing, nothing
    for i in 1:(n - 1)
        r1, r2 = rho_vals[i], rho_vals[i + 1]
        f1, f2 = mu_vals[i] - mu0, mu_vals[i + 1] - mu0
        if abs(f1) < atol
            left === nothing && (left = r1)
            right = r1
        end
        f1 == f2 && continue
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

@inline function _interp(r1::Float64, r2::Float64, m1::Float64, m2::Float64, target::Float64)
    r2 == r1 && return m1
    return m1 + (target - r1) / (r2 - r1) * (m2 - m1)
end

function _integrate_difference(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        rho_left::Float64, rho_right::Float64, mu0::Float64)
    total = 0.0
    for i in 1:(length(rho_vals) - 1)
        r1, r2 = rho_vals[i], rho_vals[i + 1]
        (r2 <= rho_left || r1 >= rho_right) && continue
        left, right = max(r1, rho_left), min(r2, rho_right)
        right <= left && continue
        mu_left = _interp(r1, r2, mu_vals[i], mu_vals[i + 1], left)
        mu_right = _interp(r1, r2, mu_vals[i], mu_vals[i + 1], right)
        total += 0.5 * ((mu_left - mu0) + (mu_right - mu0)) * (right - left)
    end
    return total
end

function _area_difference(mu0::Float64, rho_vals::Vector{Float64}, mu_vals::Vector{Float64})
    rho_left, rho_right = _find_outer_intersections(mu0, rho_vals, mu_vals)
    (rho_left === nothing || rho_right === nothing || rho_right - rho_left <= 1e-9) && return nothing
    return _integrate_difference(rho_vals, mu_vals, rho_left, rho_right, mu0)
end

function _find_mu_bracket(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        mu_lo::Float64, mu_hi::Float64, steps::Int, tol_area::Real)
    attempt = max(steps, 3)
    while attempt <= MAX_CANDIDATE_STEPS
        prev_mu, prev_area = nothing, nothing
        for mu0 in range(mu_lo, mu_hi; length=attempt)
            area = _area_difference(mu0, rho_vals, mu_vals)
            area === nothing && continue
            abs(area) <= tol_area && return (mu0, mu0, area, area)
            if prev_area !== nothing && area * prev_area < 0
                return (prev_mu, mu0, prev_area, area)
            end
            prev_mu, prev_area = mu0, area
        end
        attempt *= 2
    end
    return nothing
end

function _bisection_solve(rho_vals::Vector{Float64}, mu_vals::Vector{Float64},
        mu_a::Float64, mu_b::Float64, area_a::Real, area_b::Real,
        tol_area::Real, max_iter::Int)
    mu_a == mu_b && return mu_a, area_a, 0
    a, b, fa, fb = mu_a, mu_b, area_a, area_b
    fa * fb <= 0 || return nothing, nothing, 0

    for iter in 1:max_iter
        mid = 0.5 * (a + b)
        area_mid = _area_difference(mid, rho_vals, mu_vals)
        area_mid === nothing && return nothing, nothing, iter
        abs(area_mid) <= tol_area && return mid, area_mid, iter
        if fa * area_mid < 0
            b, fb = mid, area_mid
        else
            a, fa = mid, area_mid
        end
    end
    mid = 0.5 * (a + b)
    area_mid = _area_difference(mid, rho_vals, mu_vals)
    area_mid === nothing && return nothing, nothing, max_iter
    return mid, area_mid, max_iter
end

function _failure(reason::AbstractString; kwargs...)
    details = Dict{Symbol, Any}(:reason => reason)
    for (k, v) in kwargs
        details[k] = v
    end
    return MaxwellResult(false, nothing, nothing, nothing, nothing, 0, details)
end

function maxwell_construction(mu_vals::AbstractVector, rho_vals::AbstractVector;
        min_samples::Int=12, detect_min_points::Int=6, detect_eps::Real=1e-6,
        candidate_steps::Int=DEFAULT_CANDIDATE_STEPS,
        max_iter::Int=DEFAULT_MAX_ITER, tol_area::Real=DEFAULT_AREA_TOL,
        spinodal_hint::Union{Nothing, SShapeResult}=nothing)

    rho_sorted, mu_sorted = _prepare_curve(mu_vals, rho_vals)
    length(rho_sorted) < min_samples && return _failure("insufficient_points"; count=length(rho_sorted))

    hint = isnothing(spinodal_hint) ?
        detect_s_shape(mu_vals, rho_vals; eps=detect_eps, min_points=detect_min_points) : spinodal_hint
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
    (rho_left === nothing || rho_right === nothing || !(rho_left < rho_right)) &&
        return _failure("no_crossings"; mu_transition=mu_root, bracket=(mu_a, mu_b))

    details = Dict(:mu_bracket => (mu_a, mu_b), :rho_interval => (rho_left, rho_right),
                   :spinodal_hint => (hint.mu_spinodal_hadron, hint.mu_spinodal_quark))
    return MaxwellResult(true, mu_root, rho_left, rho_right, abs(area_root), iterations, details)
end

# ============================================================================
# 工具函数
# ============================================================================

function group_curves_by_temperature(rows; xi::Real=0.0, tol::Real=1e-6)
    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    for row in rows
        T = try parse(Float64, row["T_MeV"]) catch; continue end
        xi_val = haskey(row, "xi") ? (try parse(Float64, row["xi"]) catch; NaN end) : NaN
        (isnan(xi_val) || abs(xi_val - xi) > tol) && continue
        rho = try parse(Float64, row["rho"]) catch; continue end
        mu = haskey(row, "mu_avg_MeV") ? (try parse(Float64, row["mu_avg_MeV"]) catch; NaN end) : NaN
        isnan(mu) && (mu = haskey(row, "mu_MeV") ? (try parse(Float64, row["mu_MeV"]) catch; NaN end) : NaN)
        (isnan(mu) || !isfinite(rho)) && continue
        bucket = get!(grouped, T) do; Vector{Tuple{Float64, Float64}}() end
        push!(bucket, (mu, rho))
    end
    return grouped
end

# ============================================================================
# Crossover 检测
# ============================================================================

"""
Crossover 检测结果

字段：
- `found`: 是否找到 crossover
- `T_crossover`: crossover 温度 (fm⁻¹)
- `rho`: crossover 点的密度 (ρ/ρ₀)
- `method`: 使用的方法 (:peak 或 :inflection)
- `derivative_value`: 在 crossover 点的导数值
- `iterations`: 迭代次数
- `details`: 详细信息
"""
struct CrossoverResult
    found::Bool
    T_crossover::Union{Nothing, Float64}
    rho::Union{Nothing, Float64}
    method::Symbol
    derivative_value::Union{Nothing, Float64}
    iterations::Int
    details::Dict{Symbol, Any}
end

CrossoverResult(; method::Symbol=:peak) = CrossoverResult(false, nothing, nothing, method, nothing, 0, Dict{Symbol, Any}())

# 导入求解器模块
const _SOLVER_PATH = normpath(joinpath(@__DIR__, "..", "solver", "Solver.jl"))
if !isdefined(@__MODULE__, :Solver)
    include(_SOLVER_PATH)
end

"""
    detect_crossover(μ_fm, T_range; method=:peak, variable=:phi_u, kwargs...)

检测 crossover 温度。

# 参数
- `μ_fm`: 化学势 (fm⁻¹)
- `T_range`: 温度范围 (T_min, T_max) (fm⁻¹)
- `method`: 检测方法 (:peak 找峰值, :inflection 找拐点)
- `variable`: 检测变量 (:phi_u 手征凝聚, :Phi Polyakov loop)
- `xi`: 各向异性参数
- `n_scan`: 初始扫描点数
- `tol`: 收敛容差 (fm⁻¹)
- `max_iter`: 最大迭代次数

# 返回
`CrossoverResult` 结构
"""
function detect_crossover(μ_fm::Real, T_range::Tuple{Real, Real};
                          method::Symbol=:peak, variable::Symbol=:phi_u,
                          xi::Real=0.0, n_scan::Int=20, tol::Real=1e-4,
                          max_iter::Int=20, p_num::Int=24, t_num::Int=12)
    
    T_min, T_max = Float64.(T_range)
    T_min < T_max || return CrossoverResult(method=method)
    
    var_idx = variable == :phi_u ? 1 : (variable == :Phi ? 4 : 1)
    
    function compute_derivatives(T::Float64)
        result = Solver.solve_with_derivatives(T, Float64(μ_fm); order=2, xi=xi, p_num=p_num, t_num=t_num)
        return result.dx_dT[var_idx], result.d2x_dT2[var_idx]
    end
    
    # 检测 crossover
    if method == :peak
        result = _detect_crossover_peak(T_min, T_max, compute_derivatives, n_scan, tol, max_iter, var_idx)
    elseif method == :inflection
        result = _detect_crossover_inflection(T_min, T_max, compute_derivatives, n_scan, tol, max_iter, var_idx)
    else
        error("Unknown method: $method. Use :peak or :inflection")
    end
    
    # 如果找到 crossover，计算对应的 ρ 值
    if result.found && result.T_crossover !== nothing
        rho = _compute_rho_at_crossover(result.T_crossover, Float64(μ_fm), xi, p_num, t_num)
        return CrossoverResult(result.found, result.T_crossover, rho, result.method, 
                               result.derivative_value, result.iterations, result.details)
    end
    
    return result
end

"""计算 crossover 点的密度"""
function _compute_rho_at_crossover(T_fm::Float64, μ_fm::Float64, xi::Real, p_num::Int, t_num::Int)
    try
        sol = Solver.solve(Solver.FixedMu(), T_fm, μ_fm; xi=xi, p_num=p_num, t_num=t_num)
        if sol.converged
            return sol.rho_norm
        end
    catch e
        # 忽略错误，返回 nothing
    end
    return nothing
end

"""峰值法检测 crossover - 找最大峰（手征 crossover）

从高温向低温扫描，找到最大的局部极大值。
使用黄金分割法在峰值附近细化。
"""
function _detect_crossover_peak(T_min::Float64, T_max::Float64, 
                                compute_derivatives::Function,
                                n_scan::Int, tol::Real, max_iter::Int, var_idx::Int)
    # 从高温到低温扫描
    T_vals = collect(range(T_max, T_min; length=n_scan))
    abs_derivs = Float64[]
    derivs = Float64[]
    
    for T in T_vals
        try
            dφ_dT, _ = compute_derivatives(T)
            push!(derivs, dφ_dT)
            push!(abs_derivs, abs(dφ_dT))
        catch
            push!(derivs, NaN)
            push!(abs_derivs, NaN)
        end
    end
    
    valid_mask = .!isnan.(abs_derivs)
    !any(valid_mask) && return CrossoverResult(method=:peak)
    
    # 找所有局部极大值
    local_maxima = Int[]
    for i in 2:(length(abs_derivs)-1)
        if valid_mask[i] && valid_mask[i-1] && valid_mask[i+1]
            if abs_derivs[i] > abs_derivs[i-1] && abs_derivs[i] > abs_derivs[i+1]
                push!(local_maxima, i)
            end
        end
    end
    
    # 如果没有找到局部极大值，使用全局最大值
    if isempty(local_maxima)
        valid_indices = findall(valid_mask)
        peak_idx = valid_indices[argmax(abs_derivs[valid_mask])]
    else
        # 选择最大的局部极大值（手征 crossover 对应最大峰）
        peak_idx = local_maxima[argmax(abs_derivs[local_maxima])]
    end
    
    # 确定细化区间
    left_idx = max(1, peak_idx - 1)
    right_idx = min(length(T_vals), peak_idx + 1)
    # 注意：T_vals 是降序的，所以 a > b
    a, b = T_vals[left_idx], T_vals[right_idx]
    if a < b
        a, b = b, a
    end
    
    # 黄金分割法细化（在局部极大值附近）
    φ = (sqrt(5) - 1) / 2
    iterations = 0
    
    for iter in 1:max_iter
        iterations = iter
        a - b < tol && break
        
        c = b + φ * (a - b)
        d = a - φ * (a - b)
        
        fc = try dφ_dT, _ = compute_derivatives(c); abs(dφ_dT) catch; NaN end
        fd = try dφ_dT, _ = compute_derivatives(d); abs(dφ_dT) catch; NaN end
        
        (isnan(fc) || isnan(fd)) && break
        
        if fc > fd
            b = d
        else
            a = c
        end
    end
    
    T_crossover = (a + b) / 2
    final_deriv, _ = try compute_derivatives(T_crossover) catch; (NaN, NaN) end
    
    # 记录所有峰值信息（转换回升序显示）
    all_peaks = isempty(local_maxima) ? Tuple{Float64, Float64}[] : 
                [(T_vals[i], abs_derivs[i]) for i in local_maxima]
    
    details = Dict{Symbol, Any}(:T_range => (T_min, T_max), :n_scan => n_scan,
                                :variable_index => var_idx,
                                :n_peaks_found => length(local_maxima),
                                :all_peaks => all_peaks,
                                :scan_data => collect(zip(T_vals, derivs)))
    
    return CrossoverResult(true, T_crossover, nothing, :peak, final_deriv, iterations, details)
end

"""拐点法检测 crossover - 找手征/退禁闭 crossover 的拐点

从高温向低温扫描，找到 ∂²φ/∂T² 从负变正的位置（峰值右侧的拐点）。
使用二分法细化拐点位置。

注意：会过滤掉 |dφ/dT| 过小的虚假拐点（高温饱和区域）。
"""
function _detect_crossover_inflection(T_min::Float64, T_max::Float64,
                                      compute_derivatives::Function,
                                      n_scan::Int, tol::Real, max_iter::Int, var_idx::Int)
    # 从高温到低温扫描
    T_vals = collect(range(T_max, T_min; length=n_scan))
    d1_vals, d2_vals = Float64[], Float64[]
    
    for T in T_vals
        try
            dφ_dT, d2φ_dT2 = compute_derivatives(T)
            push!(d1_vals, dφ_dT)
            push!(d2_vals, d2φ_dT2)
        catch
            push!(d1_vals, NaN)
            push!(d2_vals, NaN)
        end
    end
    
    # 计算 |dφ/dT| 的最大值，用于过滤虚假拐点
    valid_d1 = filter(!isnan, d1_vals)
    max_abs_d1 = isempty(valid_d1) ? 1.0 : maximum(abs, valid_d1)
    # 阈值：|dφ/dT| 至少要达到最大值的 15%（提高阈值以过滤更多虚假拐点）
    d1_threshold = max_abs_d1 * 0.15
    
    # 收集所有有效的拐点
    valid_inflections = Tuple{Int, Float64}[]  # (index, d1_at_inflection)
    for i in 1:(length(d2_vals) - 1)
        if !isnan(d2_vals[i]) && !isnan(d2_vals[i+1])
            if d2_vals[i] < 0 && d2_vals[i+1] > 0  # - → +（高温到低温方向）
                # 检查拐点附近的 |dφ/dT| 是否足够大
                d1_at_inflection = (abs(d1_vals[i]) + abs(d1_vals[i+1])) / 2
                if d1_at_inflection >= d1_threshold
                    push!(valid_inflections, (i, d1_at_inflection))
                end
            end
        end
    end
    
    # 如果有多个有效拐点，选择 |dφ/dT| 最大的那个（对应真正的 crossover）
    sign_change_idx = nothing
    if !isempty(valid_inflections)
        # 选择 |dφ/dT| 最大的拐点
        _, best_pos = findmax(x -> x[2], valid_inflections)
        sign_change_idx = valid_inflections[best_pos][1]
    end
    
    if sign_change_idx === nothing
        details = Dict{Symbol, Any}(:T_range => (T_min, T_max), :n_scan => n_scan,
                                    :variable_index => var_idx, :reason => "no_valid_sign_change",
                                    :d1_threshold => d1_threshold,
                                    :scan_data => collect(zip(T_vals, d1_vals, d2_vals)))
        return CrossoverResult(false, nothing, nothing, :inflection, nothing, 0, details)
    end
    
    # 二分法细化
    # T_vals[sign_change_idx] > T_vals[sign_change_idx + 1]（降序）
    a, b = T_vals[sign_change_idx], T_vals[sign_change_idx + 1]
    if a < b
        a, b = b, a
    end
    fa, fb = d2_vals[sign_change_idx], d2_vals[sign_change_idx + 1]
    
    iterations = 0
    for iter in 1:max_iter
        iterations = iter
        a - b < tol && break
        
        mid = (a + b) / 2
        _, f_mid = try compute_derivatives(mid) catch; (NaN, NaN) end
        isnan(f_mid) && break
        
        # fa < 0, fb > 0，找零点
        if fa * f_mid < 0
            b, fb = mid, f_mid
        else
            a, fa = mid, f_mid
        end
    end
    
    T_crossover = (a + b) / 2
    final_d1, final_d2 = try compute_derivatives(T_crossover) catch; (NaN, NaN) end
    
    # 记录所有拐点信息（包括被过滤的）
    all_inflections = Tuple{Float64, Float64, Bool}[]  # (T1, T2, is_valid)
    for i in 1:(length(d2_vals) - 1)
        if !isnan(d2_vals[i]) && !isnan(d2_vals[i+1])
            if d2_vals[i] < 0 && d2_vals[i+1] > 0  # - → +
                d1_at_inflection = (abs(d1_vals[i]) + abs(d1_vals[i+1])) / 2
                is_valid = d1_at_inflection >= d1_threshold
                push!(all_inflections, (T_vals[i], T_vals[i+1], is_valid))
            end
        end
    end
    
    details = Dict{Symbol, Any}(:T_range => (T_min, T_max), :n_scan => n_scan,
                                :variable_index => var_idx, :final_d2 => final_d2,
                                :d1_threshold => d1_threshold,
                                :n_inflections_found => length(all_inflections),
                                :all_inflections => all_inflections,
                                :scan_data => collect(zip(T_vals, d1_vals, d2_vals)))
    
    return CrossoverResult(true, T_crossover, nothing, :inflection, final_d1, iterations, details)
end

"""
    scan_crossover_line(mu_range, T_range; kwargs...)

扫描 crossover 线。

# 参数
- `mu_range`: 化学势范围 (μ_min, μ_max, n_points) (fm⁻¹)
- `T_range`: 温度搜索范围 (T_min, T_max) (fm⁻¹)
- `method`: 检测方法 (:peak 或 :inflection)
- `variable`: 检测变量 (:phi_u 或 :Phi)
- `xi`: 各向异性参数

# 返回
Vector{NamedTuple} 包含 (mu_fm, T_crossover_fm, rho, converged, derivative)
"""
function scan_crossover_line(mu_range::Tuple{Real, Real, Int}, T_range::Tuple{Real, Real};
                             method::Symbol=:peak, variable::Symbol=:phi_u,
                             xi::Real=0.0, kwargs...)
    μ_min, μ_max, n_mu = mu_range
    μ_vals = range(Float64(μ_min), Float64(μ_max); length=n_mu)
    
    results = NamedTuple{(:mu_fm, :T_crossover_fm, :rho, :converged, :derivative), 
                         Tuple{Float64, Union{Nothing, Float64}, Union{Nothing, Float64}, Bool, Union{Nothing, Float64}}}[]
    
    for μ in μ_vals
        result = detect_crossover(μ, T_range; method=method, variable=variable, xi=xi, kwargs...)
        push!(results, (mu_fm=μ, T_crossover_fm=result.T_crossover, rho=result.rho,
                       converged=result.found, derivative=result.derivative_value))
    end
    
    return results
end

end # module
