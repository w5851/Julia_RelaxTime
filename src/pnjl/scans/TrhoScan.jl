"""
    TrhoScan

T-ρ 参数空间扫描模块（使用新求解器架构）。

## 功能
- 在 (T, ρ, ξ) 参数空间进行网格扫描
- 支持断点续扫（resume）
- 连续性跟踪初值策略
- 自适应密度网格

## 使用示例
```julia
using PNJL.TrhoScan

# 基本扫描
result = run_trho_scan()

# 自定义参数
result = run_trho_scan(
    T_values = 50.0:10.0:200.0,
    rho_values = 0.0:0.1:3.0,
    xi_values = [0.0],
    output_path = "my_scan.csv"
)
```
"""
module TrhoScan

using Printf
using StaticArrays

# 导入新架构模块
using ..Constants_PNJL: ħc_MeV_fm
using ..ConstraintModes: FixedRho, ConstraintMode
using ..SeedStrategies: SeedStrategy, DefaultSeed, ContinuitySeed
using ..SeedStrategies: get_seed, update!, reset!, extend_seed
using ..SeedStrategies: HADRON_SEED_5, QUARK_SEED_5, MEDIUM_SEED_5, HIGH_DENSITY_SEED_5
using ..SeedStrategies: HADRON_SEED_8, MEDIUM_SEED_8, HIGH_DENSITY_SEED_8
using ..ImplicitSolver: solve, SolverResult

export run_trho_scan, DEFAULT_T_VALUES, DEFAULT_RHO_VALUES, DEFAULT_OUTPUT_PATH
export build_default_rho_grid

# ============================================================================
# 默认配置
# ============================================================================

const DEFAULT_T_VALUES = collect(50.0:10.0:200.0)

"""构建默认密度网格（多分辨率）"""
function build_default_rho_grid(;
    rho_max::Float64=3.0,
    coarse_step::Float64=0.05,
    medium_switch::Float64=1.0,
    medium_step::Float64=0.02,
    fine_switch::Float64=0.3,
    fine_step::Float64=0.01,
    ultra_fine_switch::Float64=0.15,
    ultra_fine_step::Float64=0.005
)
    rho_max > 0 || error("rho_max must be positive")
    values = Float64[]
    
    # 粗网格
    append!(values, collect(0.0:coarse_step:rho_max))
    # 中等网格
    append!(values, collect(0.0:medium_step:medium_switch))
    # 细网格
    append!(values, collect(0.0:fine_step:fine_switch))
    # 超细网格
    append!(values, collect(0.0:ultra_fine_step:ultra_fine_switch))
    
    unique_vals = unique(round.(values; digits=6))
    sort!(unique_vals)
    return unique_vals
end

const DEFAULT_RHO_VALUES = build_default_rho_grid()
const DEFAULT_OUTPUT_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "outputs", "results", "pnjl", "trho_scan.csv"))
const SEED_KEY_DIGITS = 6
const ACCEPTABLE_RESIDUAL = 1e-4

const HEADER = join((
    "T_MeV",
    "rho",
    "xi",
    "mu_u_MeV",
    "mu_d_MeV",
    "mu_s_MeV",
    "mu_avg_MeV",
    "pressure_fm4",
    "entropy_fm3",
    "energy_fm4",
    "phi_u",
    "phi_d",
    "phi_s",
    "Phi1",
    "Phi2",
    "M_u_MeV",
    "M_d_MeV",
    "M_s_MeV",
    "iterations",
    "residual_norm",
    "converged",
    "message",
), ",")

# ============================================================================
# 主扫描函数
# ============================================================================

"""
    run_trho_scan(; kwargs...) -> NamedTuple

执行 T-ρ 参数空间扫描。

# 关键字参数
- `T_values`: 温度值列表 (MeV)，默认 50:10:200
- `rho_values`: 归一化密度值列表 (ρ/ρ₀)，默认多分辨率网格
- `xi_values`: 各向异性参数列表，默认 [0.0]
- `output_path`: 输出文件路径
- `overwrite`: 是否覆盖已有文件，默认 false
- `resume`: 是否断点续扫，默认 true
- `reverse_rho`: 是否反向扫描 ρ（从大到小），默认 true
  - 反向扫描可避免 ρ=0 奇异点导致的连续性跟踪失败
- `p_num`, `t_num`: 积分节点数
- `progress_cb`: 进度回调函数 `(point, result) -> nothing`

# 返回
NamedTuple 包含：
- `total`: 总点数
- `success`: 成功点数
- `failure`: 失败点数
- `skipped`: 跳过点数
- `output`: 输出文件路径
"""
function run_trho_scan(;
    T_values=DEFAULT_T_VALUES,
    rho_values=DEFAULT_RHO_VALUES,
    xi_values=[0.0],
    output_path::AbstractString=DEFAULT_OUTPUT_PATH,
    overwrite::Bool=false,
    resume::Bool=true,
    reverse_rho::Bool=true,
    p_num::Int=24,
    t_num::Int=8,
    progress_cb::Union{Nothing, Function}=nothing,
    nlsolve_kwargs...
)
    mkpath(dirname(output_path))
    completed = (resume && !overwrite && isfile(output_path)) ? _load_completed(output_path) : Set{NTuple{3, Float64}}()
    io_mode = (overwrite || !isfile(output_path)) ? "w" : "a"

    stats = Dict(:total => 0, :success => 0, :failure => 0, :skipped => 0)
    
    # 连续性跟踪器（按 T, xi 分组）
    continuation_seeds = Dict{Tuple{Float64, Float64}, Vector{Float64}}()
    
    # 根据 reverse_rho 决定扫描顺序
    rho_scan_order = reverse_rho ? reverse(collect(rho_values)) : collect(rho_values)

    open(output_path, io_mode) do io
        if io_mode == "w"
            println(io, HEADER)
        end

        for xi in xi_values, T in T_values
            # 每个新温度重置连续性种子
            seed_key = _seed_continuation_key(T, xi)
            delete!(continuation_seeds, seed_key)
            
            for rho in rho_scan_order
            stats[:total] += 1
            key = _key(T, rho, xi)
            
            if key in completed
                stats[:skipped] += 1
                continue
            end

            # 构建初值候选
            seed_key = _seed_continuation_key(T, xi)
            candidates = _build_seed_candidates(continuation_seeds, seed_key, T, rho)

            # 转换单位：MeV -> fm⁻¹
            T_fm = T / ħc_MeV_fm

            # 尝试求解
            result, message = _attempt_with_candidates(T_fm, rho, xi, candidates;
                p_num=p_num, t_num=t_num, nlsolve_kwargs...)

            # 更新连续性种子
            if result !== nothing && _is_success(result)
                continuation_seeds[seed_key] = copy(result.solution)
            end

            # 写入结果
            _write_row(io, T, rho, xi, result, message)
            flush(io)
            push!(completed, key)

            # 更新统计
            if _is_success(result)
                stats[:success] += 1
            else
                stats[:failure] += 1
            end

            # 进度回调
            if progress_cb !== nothing
                try
                    progress_cb((T=T, rho=rho, xi=xi), result)
                catch
                    # ignore callback errors
                end
            end
            end  # for rho
        end  # for T, xi
    end  # open

    return (;
        total=stats[:total],
        success=stats[:success],
        failure=stats[:failure],
        skipped=stats[:skipped],
        output=output_path
    )
end

# ============================================================================
# 内部辅助函数
# ============================================================================

_key(T, rho, xi) = (round(Float64(T); digits=6), round(Float64(rho); digits=6), round(Float64(xi); digits=6))
_seed_continuation_key(T, xi) = (round(Float64(T); digits=SEED_KEY_DIGITS), round(Float64(xi); digits=SEED_KEY_DIGITS))

function _load_completed(path::AbstractString)
    completed = Set{NTuple{3, Float64}}()
    open(path, "r") do io
        first_line = true
        for line in eachline(io)
            if first_line
                first_line = false
                continue
            end
            isempty(strip(line)) && continue
            cols = split(line, ',')
            length(cols) < 3 && continue
            try
                T = parse(Float64, strip(cols[1]))
                rho = parse(Float64, strip(cols[2]))
                xi = parse(Float64, strip(cols[3]))
                push!(completed, _key(T, rho, xi))
            catch
                # ignore malformed lines
            end
        end
    end
    return completed
end

"""根据密度选择合适的初值"""
function _select_seed_for_rho(rho::Float64)
    if rho < 0.5
        return HADRON_SEED_8
    elseif rho < 2.0
        return MEDIUM_SEED_8
    else
        return HIGH_DENSITY_SEED_8
    end
end

"""构建初值候选列表"""
function _build_seed_candidates(cache::Dict, seed_key, T, rho)
    candidates = NamedTuple{(:label, :state), Tuple{String, Vector{Float64}}}[]
    
    # 1. 连续性种子（优先）
    if haskey(cache, seed_key)
        cached = cache[seed_key]
        # 确保是 8 维
        if length(cached) == 8
            push!(candidates, (label="continuation", state=copy(cached)))
        elseif length(cached) >= 5
            # 扩展为 8 维
            extended = extend_seed(cached, FixedRho(rho))
            push!(candidates, (label="continuation-ext", state=extended))
        end
    end
    
    # 2. 基于密度的默认种子
    primary = _select_seed_for_rho(rho)
    push!(candidates, (label="density-based", state=copy(primary)))
    
    # 3. 其他候选
    if rho >= 0.5
        push!(candidates, (label="hadron", state=copy(HADRON_SEED_8)))
    end
    if rho < 2.0
        push!(candidates, (label="high-density", state=copy(HIGH_DENSITY_SEED_8)))
    end
    
    return candidates
end

"""尝试多个初值候选"""
function _attempt_with_candidates(T_fm, rho, xi, candidates; p_num, t_num, nlsolve_kwargs...)
    messages = String[]
    
    for candidate in candidates
        result, msg = _solve_point(T_fm, rho, xi, candidate.state; p_num=p_num, t_num=t_num, nlsolve_kwargs...)
        
        if _is_success(result)
            # 尝试精炼
            refined, refine_msg = _refine_result(T_fm, rho, xi, result; p_num=p_num, t_num=t_num, nlsolve_kwargs...)
            result = refined
            
            # 处理近似收敛
            result, promote_msg = _promote_success(result)
            
            if !isempty(msg)
                push!(messages, msg)
            end
            if !isempty(promote_msg)
                push!(messages, promote_msg)
            end
            if !isempty(refine_msg)
                push!(messages, refine_msg)
            end
            
            return result, _join_messages(messages)
        end
        
        push!(messages, _format_candidate_failure(candidate.label, msg, result))
    end
    
    return nothing, _join_messages(messages)
end

"""单点求解"""
function _solve_point(T_fm, rho_target, xi, seed_state; p_num, t_num, nlsolve_kwargs...)
    try
        mode = FixedRho(rho_target)
        
        # 确保种子是 8 维
        if length(seed_state) >= 8
            seed_8 = Float64.(seed_state[1:8])
        else
            # 扩展为 8 维
            seed_5 = seed_state[1:min(5, length(seed_state))]
            seed_8 = extend_seed(seed_5, mode)
        end
        
        # 创建自定义策略，直接返回指定种子
        strategy = _FixedSeedStrategy(seed_8)
        
        result = solve(mode, T_fm;
            xi=xi,
            seed_strategy=strategy,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...
        )
        return result, ""
    catch err
        msg = sprint() do io
            showerror(io, err)
        end
        return nothing, _clean_message(msg)
    end
end

"""固定种子策略（内部使用）"""
struct _FixedSeedStrategy <: SeedStrategy
    seed::Vector{Float64}
end

# 导入 get_seed 以便扩展
import ..SeedStrategies: get_seed

function get_seed(s::_FixedSeedStrategy, θ::AbstractVector, mode::ConstraintMode)
    return copy(s.seed)
end

"""精炼近似收敛的结果"""
function _refine_result(T_fm, rho, xi, result; p_num, t_num, nlsolve_kwargs...)
    result === nothing && return nothing, ""
    if result.converged
        return result, ""
    end
    
    residual = result.residual_norm
    if !isfinite(residual) || residual > ACCEPTABLE_RESIDUAL
        return result, ""
    end
    
    # 用当前解作为初值重新求解
    refined, msg = _solve_point(T_fm, rho, xi, result.solution; p_num=p_num, t_num=t_num, nlsolve_kwargs...)
    if refined !== nothing && refined.converged
        return refined, "refined from near-converged seed"
    end
    return result, msg
end

"""判断是否成功"""
_is_success(::Nothing) = false
function _is_success(result)
    result === nothing && return false
    if result.converged
        return true
    end
    residual = result.residual_norm
    return isfinite(residual) && residual <= ACCEPTABLE_RESIDUAL
end

"""提升近似收敛为成功"""
function _promote_success(result)
    result === nothing && return nothing, ""
    if result.converged
        return result, ""
    end
    
    residual = result.residual_norm
    if !isfinite(residual) || residual > ACCEPTABLE_RESIDUAL
        return result, ""
    end
    
    # 创建标记为收敛的新结果
    promoted = SolverResult(
        result.mode,
        true,  # converged = true
        copy(result.solution),
        result.x_state,
        result.mu_vec,
        result.omega,
        result.pressure,
        result.rho_norm,
        result.entropy,
        result.energy,
        result.masses,
        result.iterations,
        residual,
        result.xi,
    )
    msg = string("force-marked converged (residual ", _fmt(residual), ")")
    return promoted, msg
end

"""写入一行结果"""
function _write_row(io, T, rho, xi, result, message)
    if result === nothing
        values = (
            _fmt(T), _fmt(rho), _fmt(xi),
            "NaN", "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN",
            "-1", "NaN", "false",
            _quote(message),
        )
        println(io, join(values, ','))
        return
    end

    # 提取解
    mu_vec_mev = result.mu_vec .* ħc_MeV_fm
    mu_avg = sum(mu_vec_mev) / 3
    phi = result.x_state[1:3]
    Phi1 = result.x_state[4]
    Phi2 = result.x_state[5]
    masses_mev = result.masses .* ħc_MeV_fm

    values = (
        _fmt(T),
        _fmt(rho),
        _fmt(xi),
        _fmt(mu_vec_mev[1]),
        _fmt(mu_vec_mev[2]),
        _fmt(mu_vec_mev[3]),
        _fmt(mu_avg),
        _fmt(result.pressure),
        _fmt(result.entropy),
        _fmt(result.energy),
        _fmt(phi[1]),
        _fmt(phi[2]),
        _fmt(phi[3]),
        _fmt(Phi1),
        _fmt(Phi2),
        _fmt(masses_mev[1]),
        _fmt(masses_mev[2]),
        _fmt(masses_mev[3]),
        string(result.iterations),
        _fmt(result.residual_norm),
        string(result.converged),
        _quote(message),
    )
    println(io, join(values, ','))
end

# ============================================================================
# 格式化辅助
# ============================================================================

_fmt(x::Float64) = @sprintf("%.6f", x)
_fmt(x::Real) = _fmt(Float64(x))
_fmt(x) = string(x)

function _clean_message(msg::AbstractString)
    stripped = replace(strip(msg), '\n' => ' ')
    return replace(stripped, '"' => '\'')
end

_quote(msg::AbstractString) = isempty(msg) ? "" : string('"', msg, '"')
_quote(::Nothing) = ""

function _join_messages(messages)
    filtered = filter(!isempty, messages)
    return isempty(filtered) ? "" : join(filtered, " | ")
end

function _format_candidate_failure(label, message, result)
    base = "seed[$label] failed"
    if result !== nothing
        base = string(base, " (iterations=", result.iterations, 
                      ", residual=", _fmt(result.residual_norm), 
                      ", converged=", string(result.converged), ")")
    end
    if isempty(message)
        return base
    end
    return string(base, ": ", message)
end

end # module TrhoScan
