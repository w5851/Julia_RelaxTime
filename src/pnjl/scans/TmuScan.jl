"""
    TmuScan

T-μ 参数空间扫描模块（使用新求解器架构）。

## 功能
- 在 (T, μ, ξ) 参数空间进行网格扫描
- 支持断点续扫（resume）
- 连续性跟踪初值策略
- 多初值尝试处理收敛困难点

## 使用示例
```julia
using PNJL.TmuScan

# 基本扫描
result = run_tmu_scan()

# 自定义参数
result = run_tmu_scan(
    T_values = 50.0:10.0:200.0,
    mu_values = 0.0:10.0:400.0,
    xi_values = [0.0],
    output_path = "my_scan.csv"
)
```
"""
module TmuScan

using Printf
using StaticArrays

# 导入新架构模块
using ..Constants_PNJL: ħc_MeV_fm
using ..ConstraintModes: FixedMu, ConstraintMode
using ..SeedStrategies: SeedStrategy, DefaultSeed, ContinuitySeed, MultiSeed
using ..SeedStrategies: get_seed, update!, reset!, HADRON_SEED_5, QUARK_SEED_5
using ..ImplicitSolver: solve, SolverResult

export run_tmu_scan, DEFAULT_T_VALUES, DEFAULT_MU_VALUES, DEFAULT_OUTPUT_PATH

# ============================================================================
# 默认配置
# ============================================================================

const DEFAULT_T_VALUES = collect(50.0:10.0:200.0)
const DEFAULT_MU_VALUES = collect(0.0:10.0:400.0)
const DEFAULT_OUTPUT_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "outputs", "results", "pnjl", "tmu_scan.csv"))
const SEED_KEY_DIGITS = 6
const ACCEPTABLE_RESIDUAL = 1e-4

const HEADER = join((
    "T_MeV",
    "mu_MeV",
    "xi",
    "pressure_fm4",
    "rho",
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
    run_tmu_scan(; kwargs...) -> NamedTuple

执行 T-μ 参数空间扫描。

# 关键字参数
- `T_values`: 温度值列表 (MeV)，默认 50:10:200
- `mu_values`: 化学势值列表 (MeV)，默认 0:10:400
- `xi_values`: 各向异性参数列表，默认 [0.0]
- `output_path`: 输出文件路径
- `overwrite`: 是否覆盖已有文件，默认 false
- `resume`: 是否断点续扫，默认 true
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
function run_tmu_scan(;
    T_values=DEFAULT_T_VALUES,
    mu_values=DEFAULT_MU_VALUES,
    xi_values=[0.0],
    output_path::AbstractString=DEFAULT_OUTPUT_PATH,
    overwrite::Bool=false,
    resume::Bool=true,
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

    open(output_path, io_mode) do io
        if io_mode == "w"
            println(io, HEADER)
        end

        for xi in xi_values, T in T_values, mu in mu_values
            stats[:total] += 1
            key = _key(T, mu, xi)
            
            if key in completed
                stats[:skipped] += 1
                continue
            end

            # 构建初值候选
            seed_key = _seed_continuation_key(T, xi)
            candidates = _build_seed_candidates(continuation_seeds, seed_key, T, mu)

            # 转换单位：MeV -> fm⁻¹
            T_fm = T / ħc_MeV_fm
            μ_fm = mu / ħc_MeV_fm

            # 尝试求解
            result, message = _attempt_with_candidates(T_fm, μ_fm, xi, candidates;
                p_num=p_num, t_num=t_num, nlsolve_kwargs...)

            # 更新连续性种子
            if result !== nothing && _is_success(result)
                continuation_seeds[seed_key] = copy(result.solution)
            end

            # 写入结果
            _write_row(io, T, mu, xi, result, message)
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
                    progress_cb((T=T, mu=mu, xi=xi), result)
                catch
                    # ignore callback errors
                end
            end
        end
    end

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

_key(T, mu, xi) = (round(Float64(T); digits=6), round(Float64(mu); digits=6), round(Float64(xi); digits=6))
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
                mu = parse(Float64, strip(cols[2]))
                xi = parse(Float64, strip(cols[3]))
                push!(completed, _key(T, mu, xi))
            catch
                # ignore malformed lines
            end
        end
    end
    return completed
end

"""构建初值候选列表"""
function _build_seed_candidates(cache::Dict, seed_key, T, mu)
    candidates = NamedTuple{(:label, :state), Tuple{String, Vector{Float64}}}[]
    
    # 1. 连续性种子（优先）
    if haskey(cache, seed_key)
        push!(candidates, (label="continuation", state=copy(cache[seed_key])))
    end
    
    # 2. 基于相位的默认种子
    T_mev = T
    μ_mev = mu
    if T_mev > 150 || μ_mev > 300
        push!(candidates, (label="quark", state=copy(QUARK_SEED_5)))
        push!(candidates, (label="hadron", state=copy(HADRON_SEED_5)))
    else
        push!(candidates, (label="hadron", state=copy(HADRON_SEED_5)))
        push!(candidates, (label="quark", state=copy(QUARK_SEED_5)))
    end
    
    return candidates
end

"""尝试多个初值候选"""
function _attempt_with_candidates(T_fm, μ_fm, xi, candidates; p_num, t_num, nlsolve_kwargs...)
    messages = String[]
    
    for candidate in candidates
        result, msg = _solve_point(T_fm, μ_fm, xi, candidate.state; p_num=p_num, t_num=t_num, nlsolve_kwargs...)
        
        if _is_success(result)
            # 尝试精炼
            refined, refine_msg = _refine_result(T_fm, μ_fm, xi, result; p_num=p_num, t_num=t_num, nlsolve_kwargs...)
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
function _solve_point(T_fm, μ_fm, xi, seed_state; p_num, t_num, nlsolve_kwargs...)
    try
        # 创建固定种子策略
        seed_5 = Float64.(seed_state[1:min(5, length(seed_state))])
        strategy = _FixedSeedStrategy(seed_5)
        
        result = solve(FixedMu(), T_fm, μ_fm;
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
function _refine_result(T_fm, μ_fm, xi, result; p_num, t_num, nlsolve_kwargs...)
    result === nothing && return nothing, ""
    if result.converged
        return result, ""
    end
    
    residual = result.residual_norm
    if !isfinite(residual) || residual > ACCEPTABLE_RESIDUAL
        return result, ""
    end
    
    # 用当前解作为初值重新求解
    refined, msg = _solve_point(T_fm, μ_fm, xi, result.solution; p_num=p_num, t_num=t_num, nlsolve_kwargs...)
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
function _write_row(io, T, mu, xi, result, message)
    if result === nothing
        values = (
            _fmt(T), _fmt(mu), _fmt(xi),
            "NaN", "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN",
            "-1", "NaN", "false",
            _quote(message),
        )
        println(io, join(values, ','))
        return
    end

    # 提取解
    phi = result.x_state[1:3]
    Phi1 = result.x_state[4]
    Phi2 = result.x_state[5]
    masses_mev = result.masses .* ħc_MeV_fm

    values = (
        _fmt(T),
        _fmt(mu),
        _fmt(xi),
        _fmt(result.pressure),
        _fmt(result.rho_norm),
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

end # module TmuScan
