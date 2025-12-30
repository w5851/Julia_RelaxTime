#!/usr/bin/env julia
"""
PNJL 相结构计算脚本

计算并输出：
1. CEP (临界终点)
2. 一阶相变线 (T-ρ 扫描 + Maxwell 构造)
3. Spinodal 数据（亚稳态边界）

输出文件保存到 data/reference/pnjl/

用法：
    julia scripts/pnjl/calculate_phase_structure.jl [options]

选项：
    --xi=0.0          各向异性参数
    --T_min=50        最低温度 (MeV)
    --T_max=200       最高温度 (MeV)
    --T_step=10       温度步长 (MeV)
    --rho_max=4.0     最大密度 (ρ/ρ₀)
    --rho_step=0.05   密度步长
    --output_dir=...  输出目录
    --skip_trho       跳过 T-ρ 扫描（使用已有数据）
    --verbose         详细输出
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf
using Dates

# 加载模块
include(joinpath(@__DIR__, "..", "..", "src", "Constants_PNJL.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "pnjl", "PNJL.jl"))

using .PNJL
using .PNJL.TrhoScan
using .PNJL.PhaseTransition  # 使用模块化的相变分析功能

# ============================================================================
# 配置
# ============================================================================

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "data", "reference", "pnjl")

struct PhaseStructureConfig
    xi::Float64
    T_min::Float64
    T_max::Float64
    T_step::Float64
    rho_max::Float64
    rho_step::Float64
    output_dir::String
    skip_trho::Bool
    verbose::Bool
end

function parse_args(args)
    xi = 0.0
    T_min = 50.0
    T_max = 200.0
    T_step = 10.0
    rho_max = 4.0
    rho_step = 0.05
    output_dir = DEFAULT_OUTPUT_DIR
    skip_trho = false
    verbose = false
    
    for arg in args
        if startswith(arg, "--xi=")
            xi = parse(Float64, arg[6:end])
        elseif startswith(arg, "--T_min=")
            T_min = parse(Float64, arg[9:end])
        elseif startswith(arg, "--T_max=")
            T_max = parse(Float64, arg[9:end])
        elseif startswith(arg, "--T_step=")
            T_step = parse(Float64, arg[10:end])
        elseif startswith(arg, "--rho_max=")
            rho_max = parse(Float64, arg[11:end])
        elseif startswith(arg, "--rho_step=")
            rho_step = parse(Float64, arg[12:end])
        elseif startswith(arg, "--output_dir=")
            output_dir = arg[14:end]
        elseif arg == "--skip_trho"
            skip_trho = true
        elseif arg == "--verbose"
            verbose = true
        elseif arg == "--help" || arg == "-h"
            println("用法: julia calculate_phase_structure.jl [options]")
            println("选项:")
            println("  --xi=0.0          各向异性参数")
            println("  --T_min=50        最低温度 (MeV)")
            println("  --T_max=200       最高温度 (MeV)")
            println("  --T_step=10       温度步长 (MeV)")
            println("  --rho_max=4.0     最大密度 (ρ/ρ₀)")
            println("  --rho_step=0.05   密度步长")
            println("  --output_dir=...  输出目录")
            println("  --skip_trho       跳过 T-ρ 扫描")
            println("  --verbose         详细输出")
            exit(0)
        end
    end
    
    return PhaseStructureConfig(xi, T_min, T_max, T_step, rho_max, rho_step, 
                                 output_dir, skip_trho, verbose)
end

# ============================================================================
# S 形检测和 Maxwell 构造（基于 src/pnjl/analysis/ 的完善实现）
# ============================================================================

const EPS_SLOPE = 0.0
const DEFAULT_AREA_TOL = 1e-4
const DEFAULT_MAX_ITER = 60
const DEFAULT_CANDIDATE_STEPS = 64
const MAX_CANDIDATE_STEPS = 1024
const BRACKET_SHRINK_REL = 1e-3
const BRACKET_SHRINK_ABS = 1e-3

"""S 形检测结果"""
struct SShapeResult
    has_s_shape::Bool
    mu_spinodal_hadron::Union{Nothing, Float64}  # μ(ρ) 曲线的局部极大值（强子相侧）
    mu_spinodal_quark::Union{Nothing, Float64}   # μ(ρ) 曲线的局部极小值（夸克相侧）
    rho_spinodal_hadron::Union{Nothing, Float64} # 强子相侧 spinodal 的 ρ（极大值点）
    rho_spinodal_quark::Union{Nothing, Float64}  # 夸克相侧 spinodal 的 ρ（极小值点）
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
idx 是检测到符号变化的位置（斜率从正变负或从负变正）。
"""
function _refine_extremum(mu_sorted::Vector{Float64}, rho_sorted::Vector{Float64}, 
                          idx::Int; is_maximum::Bool=true)
    n = length(rho_sorted)
    
    # 确保有足够的点进行二次插值
    if idx < 2 || idx >= n
        # 回退到简单中点
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end
    
    # 使用三点二次插值找极值
    # 选择 idx-1, idx, idx+1 三个点
    i1, i2, i3 = idx - 1, idx, idx + 1
    if i3 > n
        i1, i2, i3 = idx - 2, idx - 1, idx
    end
    if i1 < 1
        i1, i2, i3 = 1, 2, 3
    end
    
    r1, r2, r3 = rho_sorted[i1], rho_sorted[i2], rho_sorted[i3]
    m1, m2, m3 = mu_sorted[i1], mu_sorted[i2], mu_sorted[i3]
    
    # 二次插值: μ(ρ) = a*ρ² + b*ρ + c
    # 极值点: ρ_ext = -b/(2a)
    denom = (r1 - r2) * (r1 - r3) * (r2 - r3)
    if abs(denom) < 1e-15
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end
    
    a = (r3 * (m2 - m1) + r2 * (m1 - m3) + r1 * (m3 - m2)) / denom
    b = (r3^2 * (m1 - m2) + r2^2 * (m3 - m1) + r1^2 * (m2 - m3)) / denom
    c = (r2 * r3 * (r2 - r3) * m1 + r3 * r1 * (r3 - r1) * m2 + r1 * r2 * (r1 - r2) * m3) / denom
    
    # 检查二次项系数符号是否符合预期
    # 极大值: a < 0, 极小值: a > 0
    if (is_maximum && a >= 0) || (!is_maximum && a <= 0)
        # 二次拟合不符合预期，回退到中点
        return (rho_sorted[idx] + rho_sorted[idx + 1]) / 2,
               (mu_sorted[idx] + mu_sorted[idx + 1]) / 2
    end
    
    rho_ext = -b / (2 * a)
    
    # 确保极值点在合理范围内
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
返回 SShapeResult，包含 spinodal 点（dμ/dρ = 0 的位置）。

物理意义：
- S 形曲线表示一阶相变区域
- spinodal_hadron: μ(ρ) 的局部极大值点，对应强子相的亚稳态边界
- spinodal_quark: μ(ρ) 的局部极小值点，对应夸克相的亚稳态边界
"""
function detect_s_shape(mu_vals::AbstractVector, rho_vals::AbstractVector; 
                        eps::Real=EPS_SLOPE, min_points::Int=5)
    n = min(length(mu_vals), length(rho_vals))
    n >= min_points || return SShapeResult()
    
    # 按 ρ 升序排列，计算 dμ/dρ
    mu_sorted, rho_sorted = _sort_curve_by_rho(mu_vals, rho_vals)
    
    # 计算斜率 dμ/dρ
    slopes = Float64[]
    slope_indices = Int[]  # 记录每个斜率对应的起始索引
    for i in 1:(length(rho_sorted) - 1)
        drho = rho_sorted[i + 1] - rho_sorted[i]
        abs(drho) <= eps && continue
        push!(slopes, (mu_sorted[i + 1] - mu_sorted[i]) / drho)
        push!(slope_indices, i)
    end
    isempty(slopes) && return SShapeResult()

    # 检测符号变化序列
    signs = Int[]
    sign_to_slope_idx = Int[]  # 记录每个符号变化对应的斜率索引
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

    # 查找 +1 → -1 → +1 模式（S 形特征）
    # 找到第一个 +1 → -1 的转折点（极大值，spinodal_hadron）
    # 找到第一个 -1 → +1 的转折点（极小值，spinodal_quark）
    max_idx = nothing  # 极大值点的符号索引
    min_idx = nothing  # 极小值点的符号索引
    
    for i in 1:(length(signs) - 1)
        if signs[i] == 1 && signs[i + 1] == -1 && max_idx === nothing
            max_idx = i + 1  # 转折发生在 i+1 位置
        elseif signs[i] == -1 && signs[i + 1] == 1 && max_idx !== nothing && min_idx === nothing
            min_idx = i + 1
        end
    end
    
    if max_idx === nothing || min_idx === nothing
        return SShapeResult(false, nothing, nothing, nothing, nothing, length(signs) - 1)
    end

    # 计算 spinodal 点位置（使用二次插值细化）
    slope_idx_max = sign_to_slope_idx[max_idx]
    slope_idx_min = sign_to_slope_idx[min_idx]
    
    data_idx_max = slope_indices[slope_idx_max]
    data_idx_min = slope_indices[slope_idx_min]
    
    # 使用二次插值细化极值点
    rho_hadron, mu_hadron = _refine_extremum(mu_sorted, rho_sorted, data_idx_max; is_maximum=true)
    rho_quark, mu_quark = _refine_extremum(mu_sorted, rho_sorted, data_idx_min; is_maximum=false)

    return SShapeResult(true, mu_hadron, mu_quark, rho_hadron, rho_quark, length(signs) - 1)
end

"""Maxwell 构造结果"""
struct MaxwellResult
    converged::Bool
    mu_coex_MeV::Union{Nothing, Float64}
    rho_gas::Union{Nothing, Float64}
    rho_liquid::Union{Nothing, Float64}
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
    # mu_spinodal_hadron 是 μ(ρ) 的极大值，mu_spinodal_quark 是极小值
    # Maxwell 构造的 μ_coex 在两者之间
    mu_max = hint.mu_spinodal_hadron  # μ 极大值
    mu_min = hint.mu_spinodal_quark   # μ 极小值
    if mu_max === nothing || mu_min === nothing
        return nothing
    end
    # 确保返回 (较小值, 较大值)
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
        if r2 <= rho_left
            continue
        end
        if r1 >= rho_right
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

"""创建失败结果"""
function _failure(reason::AbstractString; kwargs...)
    details = Dict{Symbol, Any}(:reason => reason)
    for (k, v) in kwargs
        details[k] = v
    end
    return MaxwellResult(false, nothing, nothing, nothing, nothing, 0, details)
end

"""
    maxwell_rho_mu(mu_vals, rho_vals; kwargs...)

Maxwell 等面积构造，使用完善的算法：
1. 检测 S 形曲线
2. 从 spinodal 点估计 μ 搜索区间
3. 在区间内搜索符号变化
4. 二分法精确求解等面积点
"""
function maxwell_rho_mu(mu_vals::AbstractVector, rho_vals::AbstractVector;
        min_samples::Int=12, detect_min_points::Int=6, detect_eps::Real=1e-6,
        candidate_steps::Int=DEFAULT_CANDIDATE_STEPS,
        max_iter::Int=DEFAULT_MAX_ITER, tol_area::Real=DEFAULT_AREA_TOL,
        spinodal_hint::Union{Nothing, SShapeResult}=nothing)

    rho_sorted, mu_sorted = _prepare_curve(mu_vals, rho_vals)
    if length(rho_sorted) < min_samples
        return _failure("insufficient_points"; count=length(rho_sorted))
    end

    # 检测 S 形
    hint = isnothing(spinodal_hint) ?
        detect_s_shape(mu_vals, rho_vals; eps=detect_eps, min_points=detect_min_points) :
        spinodal_hint
    hint.has_s_shape || return _failure("no_s_shape")

    # 获取 μ 搜索区间
    mu_bracket = _mu_bracket(hint)
    mu_bracket === nothing && return _failure("invalid_mu_bracket")
    mu_lo, mu_hi = mu_bracket
    
    # 收缩区间边界
    tightened = _shrink_bracket(mu_lo, mu_hi)
    tightened === nothing && return _failure("degenerate_bracket"; bracket=(mu_lo, mu_hi))
    mu_lo, mu_hi = tightened

    # 搜索符号变化点
    bracket = _find_mu_bracket(rho_sorted, mu_sorted, mu_lo, mu_hi, candidate_steps, tol_area)
    bracket === nothing && return _failure("no_sign_change"; bracket=(mu_lo, mu_hi))
    mu_a, mu_b, area_a, area_b = bracket

    # 二分法求解
    mu_root, area_root, iterations = _bisection_solve(rho_sorted, mu_sorted,
        mu_a, mu_b, area_a, area_b, tol_area, max_iter)
    mu_root === nothing && return _failure("bisection_failed"; bracket=(mu_a, mu_b))

    # 找到共存密度
    rho_left, rho_right = _find_outer_intersections(mu_root, rho_sorted, mu_sorted)
    if rho_left === nothing || rho_right === nothing || !(rho_left < rho_right)
        return _failure("no_crossings"; mu_coex=mu_root, bracket=(mu_a, mu_b))
    end

    details = Dict(
        :mu_bracket => (mu_a, mu_b),
        :rho_interval => (rho_left, rho_right),
        :spinodal_hint => (hint.mu_spinodal_hadron, hint.mu_spinodal_quark),
    )
    return MaxwellResult(true, mu_root, rho_left, rho_right, abs(area_root), iterations, details)
end

# ============================================================================
# 主函数
# ============================================================================

function main(args=ARGS)
    config = parse_args(args)
    
    println("=" ^ 60)
    println("PNJL 相结构计算")
    println("=" ^ 60)
    println("时间: $(now())")
    println("参数:")
    println("  xi = $(config.xi)")
    println("  T 范围: $(config.T_min) - $(config.T_max) MeV (步长 $(config.T_step))")
    println("  ρ 范围: 0 - $(config.rho_max) ρ₀ (步长 $(config.rho_step))")
    println("  输出目录: $(config.output_dir)")
    println()
    
    mkpath(config.output_dir)
    
    # Step 1: T-ρ 扫描
    trho_path = joinpath(config.output_dir, "trho_scan_xi$(config.xi).csv")
    curves = step1_trho_scan(config, trho_path)
    
    # Step 2: CEP 搜索
    cep_result = step2_find_cep(config, curves)
    
    # Step 3: Maxwell 构造
    boundary_results = step3_maxwell_construction(config, curves)
    
    # Step 4: 保存结果（包括 spinodal）
    step4_save_results(config, cep_result, boundary_results, curves)
    
    println("\n" * "=" ^ 60)
    println("计算完成!")
    println("=" ^ 60)
end

# ============================================================================
# Step 1: T-ρ 扫描
# ============================================================================

function step1_trho_scan(config::PhaseStructureConfig, output_path::String)
    println("\n[Step 1] T-ρ 扫描")
    println("-" ^ 40)
    
    if config.skip_trho && isfile(output_path)
        println("跳过扫描，加载已有数据: $output_path")
        return load_curves_from_csv(output_path, config.xi)
    end
    
    T_values = collect(config.T_min:config.T_step:config.T_max)
    rho_values = collect(0.0:config.rho_step:config.rho_max)
    
    println("温度点数: $(length(T_values))")
    println("密度点数: $(length(rho_values))")
    println("总点数: $(length(T_values) * length(rho_values))")
    
    # 进度回调
    progress_cb = nothing
    if config.verbose
        progress_cb = (point, result) -> begin
            if result !== nothing && result.converged
                print(".")
            else
                print("x")
            end
        end
    end
    
    println("\n开始扫描...")
    t_start = time()
    
    result = run_trho_scan(
        T_values = T_values,
        rho_values = rho_values,
        xi_values = [config.xi],
        output_path = output_path,
        overwrite = true,
        reverse_rho = true,  # 反向扫描避免 ρ=0 奇异点问题
        progress_cb = progress_cb
    )
    
    t_elapsed = time() - t_start
    
    println("\n扫描完成:")
    println("  成功: $(result.success) / $(result.total)")
    println("  失败: $(result.failure)")
    println("  耗时: $(round(t_elapsed, digits=1)) 秒")
    
    # 加载曲线数据
    return load_curves_from_csv(output_path, config.xi)
end

function load_curves_from_csv(path::String, xi::Float64)
    curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    
    isfile(path) || return curves
    
    # 按温度分组
    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    
    for line in eachline(path)
        startswith(line, "T_MeV") && continue
        isempty(strip(line)) && continue
        
        cols = split(line, ',')
        length(cols) < 21 && continue
        
        # TrhoScan 输出格式:
        # T_MeV(1), rho(2), xi(3), mu_u(4), mu_d(5), mu_s(6), mu_avg(7), ...
        # converged(21)
        T = tryparse(Float64, cols[1])
        rho = tryparse(Float64, cols[2])
        xi_val = tryparse(Float64, cols[3])
        mu = tryparse(Float64, cols[7])  # mu_avg_MeV
        converged = lowercase(strip(cols[21])) in ("true", "1")
        
        T === nothing && continue
        rho === nothing && continue
        mu === nothing && continue
        !converged && continue
        xi_val === nothing && continue
        abs(xi_val - xi) > 1e-6 && continue
        
        if !haskey(grouped, T)
            grouped[T] = Tuple{Float64, Float64}[]
        end
        push!(grouped[T], (mu, rho))
    end
    
    # 转换为 curves 格式
    for (T, samples) in grouped
        length(samples) < 3 && continue
        mu_vals = Float64[s[1] for s in samples]
        rho_vals = Float64[s[2] for s in samples]
        curves[T] = (mu_vals, rho_vals)
    end
    
    println("  加载了 $(length(curves)) 个温度点的曲线数据")
    
    return curves
end

# ============================================================================
# Step 2: CEP 搜索
# ============================================================================

"""
CEP 搜索：使用二分法细化 CEP 位置

算法：
1. 找到初始区间 [T_low, T_high]，其中 T_low 有 S 形，T_high 无 S 形
2. 二分法细化：计算 T_mid 的曲线，判断是否有 S 形
3. 重复直到区间宽度 < tol
"""
function step2_find_cep(config::PhaseStructureConfig, curves;
                        tol::Float64=0.01,  # 温度精度 (MeV)
                        max_bisect_iter::Int=20)
    println("\n[Step 2] CEP 搜索")
    println("-" ^ 40)
    
    if isempty(curves)
        println("警告: 无曲线数据，跳过 CEP 搜索")
        return (has_cep=false, T_cep=NaN, mu_cep=NaN)
    end
    
    println("可用温度点: $(length(curves))")
    
    # 按温度排序
    temperatures = sort(collect(keys(curves)))
    
    # 找到初始区间：最后一个有 S 形的温度和第一个没有 S 形的温度
    last_with_s = nothing
    first_without_s = nothing
    s_shape_cache = Dict{Float64, SShapeResult}()
    
    for T in temperatures
        mu_vals, rho_vals = curves[T]
        result = detect_s_shape(mu_vals, rho_vals)
        s_shape_cache[T] = result
        
        if result.has_s_shape
            last_with_s = (T, result)
        elseif last_with_s !== nothing && first_without_s === nothing
            first_without_s = T
        end
    end
    
    if last_with_s === nothing
        println("未找到 S 形曲线 (可能全为 crossover)")
        return (has_cep=false, T_cep=NaN, mu_cep=NaN)
    end
    
    T_low, res_low = last_with_s
    T_high = first_without_s !== nothing ? first_without_s : T_low + config.T_step
    
    println("初始区间: [$T_low, $T_high] MeV")
    
    # 二分法细化
    bisect_count = 0
    while (T_high - T_low) > tol && bisect_count < max_bisect_iter
        T_mid = (T_low + T_high) / 2
        
        # 检查是否已有该温度的曲线
        if haskey(curves, T_mid)
            mu_vals, rho_vals = curves[T_mid]
        else
            # 需要计算新的曲线（通过线性插值估计，或者实际计算）
            # 这里使用线性插值作为近似
            mu_vals, rho_vals = _interpolate_curve(curves, T_mid)
            if mu_vals === nothing
                println("  无法插值 T=$T_mid MeV 的曲线，停止二分")
                break
            end
            curves[T_mid] = (mu_vals, rho_vals)
        end
        
        result = detect_s_shape(mu_vals, rho_vals)
        s_shape_cache[T_mid] = result
        
        if result.has_s_shape
            T_low = T_mid
            res_low = result
            if config.verbose
                println("  T=$T_mid MeV: 有 S 形 → 更新下界")
            end
        else
            T_high = T_mid
            if config.verbose
                println("  T=$T_mid MeV: 无 S 形 → 更新上界")
            end
        end
        
        bisect_count += 1
    end
    
    # CEP 估计为区间中点
    T_cep = (T_low + T_high) / 2
    
    # μ_CEP 估计为 spinodal 中点
    mu_hadron = something(res_low.mu_spinodal_hadron, NaN)
    mu_quark = something(res_low.mu_spinodal_quark, NaN)
    mu_cep = (mu_hadron + mu_quark) / 2
    
    println("找到 CEP ($(bisect_count) 次二分):")
    println("  T_CEP ≈ $(round(T_cep, digits=2)) MeV (区间: $(round(T_low, digits=2)) - $(round(T_high, digits=2)))")
    println("  μ_CEP ≈ $(round(mu_cep, digits=2)) MeV")
    
    return (has_cep=true, T_cep=T_cep, mu_cep=mu_cep)
end

"""线性插值获取中间温度的曲线"""
function _interpolate_curve(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}, T_target::Float64)
    temps = sort(collect(keys(curves)))
    
    # 找到相邻的温度点
    T_below = nothing
    T_above = nothing
    for T in temps
        if T < T_target
            T_below = T
        elseif T > T_target && T_above === nothing
            T_above = T
        end
    end
    
    (T_below === nothing || T_above === nothing) && return nothing, nothing
    
    mu_below, rho_below = curves[T_below]
    mu_above, rho_above = curves[T_above]
    
    # 需要相同的 ρ 采样点
    length(rho_below) == length(rho_above) || return nothing, nothing
    
    # 线性插值 μ
    α = (T_target - T_below) / (T_above - T_below)
    mu_interp = mu_below .+ α .* (mu_above .- mu_below)
    
    return mu_interp, rho_below
end

# ============================================================================
# Step 3: Maxwell 构造
# ============================================================================

function step3_maxwell_construction(config::PhaseStructureConfig, curves)
    println("\n[Step 3] Maxwell 构造")
    println("-" ^ 40)
    
    if isempty(curves)
        println("警告: 无曲线数据，跳过 Maxwell 构造")
        return Dict{Float64, MaxwellResult}()
    end
    
    results = Dict{Float64, MaxwellResult}()
    
    for T in sort(collect(keys(curves)))
        mu_vals, rho_vals = curves[T]
        # 使用完善的 maxwell_rho_mu 算法
        result = maxwell_rho_mu(mu_vals, rho_vals;
            min_samples=8,           # 降低最小样本要求
            detect_min_points=5,     # S 形检测最小点数
            detect_eps=1e-6,
            candidate_steps=64,
            max_iter=60,
            tol_area=1e-4
        )
        results[T] = result
    end
    
    success_count = count(r -> r.converged, values(results))
    println("Maxwell 构造结果: $(success_count) / $(length(results)) 成功")
    
    if config.verbose
        println("\n详细结果:")
        for T in sort(collect(keys(results)))
            r = results[T]
            if r.converged
                area_str = r.area_residual !== nothing ? @sprintf("%.2e", r.area_residual) : "N/A"
                println("  T=$(T) MeV: μ_c=$(round(r.mu_coex_MeV, digits=2)) MeV, " *
                       "ρ_gas=$(round(r.rho_gas, digits=3)), ρ_liquid=$(round(r.rho_liquid, digits=3)), " *
                       "残差=$(area_str), 迭代=$(r.iterations)")
            else
                reason = get(r.details, :reason, "unknown")
                println("  T=$(T) MeV: 失败 ($(reason))")
            end
        end
    end
    
    return results
end

# ============================================================================
# Step 4: 保存结果
# ============================================================================

function step4_save_results(config::PhaseStructureConfig, cep_result, boundary_results, curves)
    println("\n[Step 4] 保存结果")
    println("-" ^ 40)
    
    # 保存 CEP
    cep_path = joinpath(config.output_dir, "cep.csv")
    save_cep(cep_path, config.xi, cep_result)
    println("CEP 数据: $cep_path")
    
    # 保存相变线
    boundary_path = joinpath(config.output_dir, "boundary.csv")
    save_boundary(boundary_path, config.xi, boundary_results)
    println("相变线数据: $boundary_path")
    
    # 保存 spinodal（亚稳态边界）
    spinodal_path = joinpath(config.output_dir, "spinodals.csv")
    save_spinodals(spinodal_path, config.xi, curves)
    println("Spinodal 数据: $spinodal_path")
end

function save_cep(path::String, xi::Float64, result)
    # 读取现有数据（如果存在）
    existing = Dict{Float64, Tuple{Float64, Float64}}()
    if isfile(path)
        for line in eachline(path)
            startswith(line, "xi") && continue
            cols = split(line, ',')
            length(cols) >= 3 || continue
            xi_val = tryparse(Float64, cols[1])
            T_cep = tryparse(Float64, cols[2])
            mu_cep = tryparse(Float64, cols[3])
            xi_val === nothing && continue
            existing[xi_val] = (T_cep, mu_cep)
        end
    end
    
    # 更新当前 xi 的数据
    if result.has_cep
        existing[xi] = (result.T_cep, result.mu_cep)
    end
    
    # 写入文件
    open(path, "w") do io
        println(io, "xi,T_CEP_MeV,mu_CEP_MeV")
        for xi_val in sort(collect(keys(existing)))
            T_cep, mu_cep = existing[xi_val]
            println(io, "$xi_val,$T_cep,$mu_cep")
        end
    end
end

function save_boundary(path::String, xi::Float64, results::Dict{Float64, MaxwellResult})
    # 读取现有数据（如果存在）
    existing = Dict{Tuple{Float64, Float64}, Tuple{Float64, Float64, Float64}}()
    if isfile(path)
        for line in eachline(path)
            startswith(line, "xi") && continue
            cols = split(line, ',')
            length(cols) >= 5 || continue
            xi_val = tryparse(Float64, cols[1])
            T = tryparse(Float64, cols[2])
            mu = tryparse(Float64, cols[3])
            rho_gas = tryparse(Float64, cols[4])
            rho_liquid = tryparse(Float64, cols[5])
            xi_val === nothing && continue
            existing[(xi_val, T)] = (mu, rho_gas, rho_liquid)
        end
    end
    
    # 更新当前 xi 的数据
    for (T, r) in results
        if r.converged && r.mu_coex_MeV !== nothing && r.rho_gas !== nothing && r.rho_liquid !== nothing
            existing[(xi, T)] = (r.mu_coex_MeV, r.rho_gas, r.rho_liquid)
        end
    end
    
    # 写入文件
    open(path, "w") do io
        println(io, "xi,T_MeV,mu_transition_MeV,rho_hadron,rho_quark")
        for key in sort(collect(keys(existing)))
            xi_val, T = key
            mu, rho_hadron, rho_quark = existing[key]
            println(io, "$xi_val,$T,$mu,$rho_hadron,$rho_quark")
        end
    end
end

"""保存 spinodal（亚稳态边界）数据"""
function save_spinodals(path::String, xi::Float64, curves)
    # 读取现有数据（如果存在）
    existing = Dict{Tuple{Float64, Float64}, Tuple{Float64, Float64, Float64, Float64}}()
    if isfile(path)
        for line in eachline(path)
            startswith(line, "xi") && continue
            cols = split(line, ',')
            length(cols) >= 6 || continue
            xi_val = tryparse(Float64, cols[1])
            T = tryparse(Float64, cols[2])
            mu_low = tryparse(Float64, cols[3])
            mu_high = tryparse(Float64, cols[4])
            rho_low = tryparse(Float64, cols[5])
            rho_high = tryparse(Float64, cols[6])
            xi_val === nothing && continue
            existing[(xi_val, T)] = (mu_low, mu_high, rho_low, rho_high)
        end
    end
    
    # 从曲线数据提取 spinodal
    for (T, (mu_vals, rho_vals)) in curves
        result = detect_s_shape(mu_vals, rho_vals)
        if result.has_s_shape && 
           result.mu_spinodal_hadron !== nothing && 
           result.mu_spinodal_quark !== nothing &&
           result.rho_spinodal_hadron !== nothing &&
           result.rho_spinodal_quark !== nothing
            existing[(xi, T)] = (
                result.mu_spinodal_hadron,
                result.mu_spinodal_quark,
                result.rho_spinodal_hadron,
                result.rho_spinodal_quark
            )
        end
    end
    
    # 写入文件
    open(path, "w") do io
        println(io, "xi,T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,rho_spinodal_hadron,rho_spinodal_quark")
        for key in sort(collect(keys(existing)))
            xi_val, T = key
            mu_low, mu_high, rho_hadron, rho_quark = existing[key]
            println(io, "$xi_val,$T,$mu_low,$mu_high,$rho_hadron,$rho_quark")
        end
    end
end

# ============================================================================
# 入口
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
