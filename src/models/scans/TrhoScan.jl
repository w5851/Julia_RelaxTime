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
using Main.Constants_PNJL: ħc_MeV_fm
import Main.Models: FixedRho, FixedAsymmetricRho, ConstraintMode
using ..SeedStrategies: SeedStrategy, DefaultSeed, ContinuitySeed, HybridContinuitySeed
using ..SeedStrategies: get_seed, update!, reset!, extend_seed
using ..SeedStrategies: HADRON_SEED_5, QUARK_SEED_5, MEDIUM_SEED_5, HIGH_DENSITY_SEED_5
using ..SeedStrategies: HADRON_SEED_8, MEDIUM_SEED_8, HIGH_DENSITY_SEED_8
using ..ImplicitSolver: solve, SolverResult, solve_weighted_block_fallback
using ..ScanCommon
using ..ScanConfig: TrhoScanConfig, scan_kwargs
using ..ScanResultFinalize: finalize_solver_result, promote_near_converged, is_success, refine_near_converged

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
const DEFAULT_OUTPUT_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "outputs", "results", "pnjl", "scan", "trho", "trho_scan.csv"))
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
    "mu_B_MeV",
    "mu_Q_MeV",
    "mu_S_MeV",
    "pressure_fm4",
    "entropy_fm3",
    "energy_fm4",
    "rho_u_fm3",
    "rho_d_fm3",
    "rho_s_fm3",
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

@inline function _validate_real_vector(name::Symbol, values)
    values isa AbstractVector{<:Real} || throw(ArgumentError("$(name) must be AbstractVector{<:Real}, got $(typeof(values))"))
    isempty(values) && throw(ArgumentError("$(name) must not be empty"))
    for (i, v) in pairs(values)
        isfinite(Float64(v)) || throw(ArgumentError("$(name)[$(i)] must be finite Real, got $(v)"))
    end
    return nothing
end

@inline function _validate_trho_scan_inputs(T_values, rho_values, xi_values,
                                            seed_policy::Symbol, constraint_mode::Symbol,
                                            solver_backend::Symbol)
    _validate_real_vector(:T_values, T_values)
    _validate_real_vector(:rho_values, rho_values)
    _validate_real_vector(:xi_values, xi_values)

    (seed_policy === :hybrid_continuity || seed_policy === :candidates) ||
        throw(ArgumentError("seed_policy must be :hybrid_continuity or :candidates, got $(seed_policy)"))
    (constraint_mode === :fixed_rho || constraint_mode === :fixed_asymmetric_rho) ||
        throw(ArgumentError("constraint_mode must be :fixed_rho or :fixed_asymmetric_rho, got $(constraint_mode)"))
    (solver_backend === :legacy || solver_backend === :models || solver_backend === :auto) ||
        throw(ArgumentError("solver_backend must be :legacy, :models or :auto, got $(solver_backend)"))
    return nothing
end

@inline function _effective_solver_backend(solver_backend::Symbol)::Symbol
    return solver_backend === :auto ? :models : solver_backend
end

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
- `seed_policy`: 初值策略模式，默认 `:hybrid_continuity`
    - `:hybrid_continuity`：连续性优先，失败后回退 MultiSeed（仅 `fixed_asymmetric_rho`）
    - `:candidates`：使用旧的候选初值链路
- `hybrid_weighted_fallback`: 是否在 `hybrid` 失败后启用 weighted-block 兜底（默认 true）
- `hybrid_weighted_max_seed_candidates`: weighted fallback 最多尝试的 seed 数（默认 3）
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
    seed_policy::Symbol=:hybrid_continuity,
    hybrid_weighted_fallback::Bool=true,
    hybrid_weighted_max_seed_candidates::Int=3,
    constraint_mode::Symbol=:fixed_rho,
    asym_ud_ratio_target::Float64=0.876,
    asym_s_target::Float64=0.0,
    solver_backend::Symbol=:legacy,
    model_kind::Symbol=:PNJL,
    p_num::Int=24,
    t_num::Int=8,
    progress_cb::Union{Nothing, Function}=nothing,
    nlsolve_kwargs...
)
    _validate_trho_scan_inputs(T_values, rho_values, xi_values, seed_policy, constraint_mode, solver_backend)

    mkpath(dirname(output_path))
    completed = (resume && !overwrite && isfile(output_path)) ? ScanCommon.load_completed_keys3(output_path; digits=6) : Set{NTuple{3, Float64}}()
    io_mode = (overwrite || !isfile(output_path)) ? "w" : "a"

    stats = Dict(:total => 0, :success => 0, :failure => 0, :skipped => 0)
    
    # 旧候选链路连续性缓存（按 T, xi 分组）
    continuation_seeds = Dict{Tuple{Float64, Float64}, Vector{Float64}}()
    # 新混合策略跟踪器（按 T, xi 分组）
    hybrid_trackers = Dict{Tuple{Float64, Float64}, HybridContinuitySeed}()
    
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
            if seed_policy === :hybrid_continuity
                delete!(hybrid_trackers, seed_key)
                hybrid_trackers[seed_key] = HybridContinuitySeed()
            end
            
            for rho in rho_scan_order
            stats[:total] += 1
            key = ScanCommon.key3(T, rho, xi; digits=6)
            
            if key in completed
                stats[:skipped] += 1
                continue
            end

            # 转换单位：MeV -> fm⁻¹
            T_fm = T / ħc_MeV_fm

            local result, message
            if seed_policy === :hybrid_continuity && constraint_mode === :fixed_asymmetric_rho
                strategy = hybrid_trackers[seed_key]
                result, message = _attempt_with_strategy(T_fm, rho, xi, strategy;
                    constraint_mode=constraint_mode,
                    hybrid_weighted_fallback=hybrid_weighted_fallback,
                    hybrid_weighted_max_seed_candidates=hybrid_weighted_max_seed_candidates,
                    asym_ud_ratio_target=asym_ud_ratio_target,
                    asym_s_target=asym_s_target,
                    solver_backend=solver_backend,
                    model_kind=model_kind,
                    p_num=p_num, t_num=t_num, nlsolve_kwargs...)
            else
                # 兼容旧链路
                candidates = _build_seed_candidates(continuation_seeds, seed_key, T, rho)
                result, message = _attempt_with_candidates(T_fm, rho, xi, candidates;
                    constraint_mode=constraint_mode,
                    asym_ud_ratio_target=asym_ud_ratio_target,
                    asym_s_target=asym_s_target,
                    solver_backend=solver_backend,
                    model_kind=model_kind,
                    p_num=p_num, t_num=t_num, nlsolve_kwargs...)

                if result !== nothing && _is_success(result)
                    continuation_seeds[seed_key] = copy(result.solution)
                end
            end

            # 写入结果
            _write_row(io, T, rho, xi, result, message;
                model_kind=model_kind,
                p_num=p_num,
                t_num=t_num,
            )
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

"""run_trho_scan(config::TrhoScanConfig; kwargs...) -> NamedTuple

结构化配置入口。`kwargs` 会覆盖 config 中同名项，保持与原 kwargs 接口兼容。
"""
function run_trho_scan(config::TrhoScanConfig; kwargs...)
    cfg = scan_kwargs(config)
    return run_trho_scan(; cfg..., config.nlsolve_kwargs..., kwargs...)
end

# ============================================================================
# 内部辅助函数
# ============================================================================

_seed_continuation_key(T, xi) = ScanCommon.key2(T, xi; digits=SEED_KEY_DIGITS)

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
function _attempt_with_candidates(T_fm, rho, xi, candidates;
    constraint_mode::Symbol=:fixed_rho,
    asym_ud_ratio_target::Float64=0.876,
    asym_s_target::Float64=0.0,
    solver_backend::Symbol=:legacy,
    model_kind::Symbol=:PNJL,
    p_num,
    t_num,
    nlsolve_kwargs...)
    return ScanCommon.attempt_with_candidates(candidates;
        solve_point=seed_state -> _solve_point(T_fm, rho, xi, seed_state;
            constraint_mode=constraint_mode,
            asym_ud_ratio_target=asym_ud_ratio_target,
            asym_s_target=asym_s_target,
            solver_backend=solver_backend,
            model_kind=model_kind,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...,
        ),
        refine=result -> _refine_result(T_fm, rho, xi, result;
            constraint_mode=constraint_mode,
            asym_ud_ratio_target=asym_ud_ratio_target,
            asym_s_target=asym_s_target,
            solver_backend=solver_backend,
            model_kind=model_kind,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...,
        ),
        promote=_promote_success,
        is_success=_is_success,
    )
end

"""使用单一策略（如 HybridContinuitySeed）尝试求解"""
function _attempt_with_strategy(T_fm, rho, xi, strategy::SeedStrategy;
    constraint_mode::Symbol=:fixed_rho,
    hybrid_weighted_fallback::Bool=true,
    hybrid_weighted_max_seed_candidates::Int=3,
    asym_ud_ratio_target::Float64=0.876,
    asym_s_target::Float64=0.0,
    solver_backend::Symbol=:legacy,
    model_kind::Symbol=:PNJL,
    p_num,
    t_num,
    nlsolve_kwargs...)

    effective_solver_backend = _effective_solver_backend(solver_backend)

    mode = if constraint_mode === :fixed_rho
        FixedRho(rho)
    elseif constraint_mode === :fixed_asymmetric_rho
        FixedAsymmetricRho(rho, asym_ud_ratio_target, asym_s_target)
    else
        error("unknown constraint_mode=$constraint_mode (expected :fixed_rho or :fixed_asymmetric_rho)")
    end

    primary_err_msg = ""
    result = try
        if effective_solver_backend === :models
            _solve_with_models(mode, T_fm;
                xi=xi,
                model_kind=model_kind,
                seed_strategy=strategy,
                p_num=p_num,
                t_num=t_num,
                nlsolve_kwargs...)
        elseif mode isa FixedAsymmetricRho
            solve(mode, T_fm;
                xi=xi,
                model_kind=model_kind,
                seed_strategy=strategy,
                p_num=p_num,
                t_num=t_num,
                nlsolve_kwargs...
            )
        else
            solve(mode, T_fm;
                xi=xi,
                seed_strategy=strategy,
                p_num=p_num,
                t_num=t_num,
                nlsolve_kwargs...
            )
        end
    catch err
        primary_err_msg = sprint() do io
            showerror(io, err)
        end
        nothing
    end

    if result !== nothing
        result = finalize_solver_result(result, T_fm, xi;
            solver_backend=effective_solver_backend,
            p_num=p_num,
            t_num=t_num,
            model_kind=model_kind,
        )
    end

    if _is_success(result)
        refined, refine_msg = _refine_result(T_fm, rho, xi, result;
            constraint_mode=constraint_mode,
            asym_ud_ratio_target=asym_ud_ratio_target,
            asym_s_target=asym_s_target,
            solver_backend=solver_backend,
            model_kind=model_kind,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...,
        )
        promoted, promote_msg = _promote_success(refined)
        msg = _join_messages(filter(!isempty, [refine_msg, promote_msg]))
        return promoted, msg
    end

    if hybrid_weighted_fallback &&
       constraint_mode === :fixed_asymmetric_rho &&
       strategy isa HybridContinuitySeed

        initial_seed = if result !== nothing && length(result.solution) >= 8
            copy(result.solution)
        elseif strategy.continuity.previous_solution !== nothing
            copy(strategy.continuity.previous_solution)
        else
            get_seed(strategy, [T_fm], mode)
        end

        wb_result = solve_weighted_block_fallback(mode, T_fm;
            initial_seed=initial_seed,
            max_seed_candidates=hybrid_weighted_max_seed_candidates,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            model_kind=model_kind)

        if wb_result !== nothing
            wb_final = finalize_solver_result(wb_result, T_fm, xi;
                solver_backend=effective_solver_backend,
                p_num=p_num,
                t_num=t_num,
                model_kind=model_kind,
            )

            if _is_success(wb_final)
                update!(strategy, wb_final.solution)
                msg = isempty(primary_err_msg) ? "hybrid weighted-block fallback rescued" : "hybrid weighted-block fallback rescued; primary=$( _clean_message(primary_err_msg) )"
                return wb_final, msg
            end
        end
    end

    msg = isempty(primary_err_msg) ? "seed[strategy] failed" : _clean_message(primary_err_msg)
    return result, msg
end

"""单点求解"""
function _solve_point(T_fm, rho_target, xi, seed_state;
    constraint_mode::Symbol=:fixed_rho,
    asym_ud_ratio_target::Float64=0.876,
    asym_s_target::Float64=0.0,
    solver_backend::Symbol=:legacy,
    model_kind::Symbol=:PNJL,
    p_num,
    t_num,
    nlsolve_kwargs...)
    try
        effective_solver_backend = _effective_solver_backend(solver_backend)
        (effective_solver_backend === :legacy || effective_solver_backend === :models) ||
            error("unknown solver_backend=$solver_backend (expected :legacy, :models or :auto)")

        mode = if constraint_mode === :fixed_rho
            FixedRho(rho_target)
        elseif constraint_mode === :fixed_asymmetric_rho
            FixedAsymmetricRho(rho_target, asym_ud_ratio_target, asym_s_target)
        else
            error("unknown constraint_mode=$constraint_mode (expected :fixed_rho or :fixed_asymmetric_rho)")
        end
        
        # 确保种子是 8 维
        if length(seed_state) >= 8
            seed_8 = Float64.(seed_state[1:8])
        else
            # 扩展为 8 维
            seed_5 = seed_state[1:min(5, length(seed_state))]
            seed_8 = extend_seed(seed_5, mode)
        end
        
        # 创建自定义策略，直接返回指定种子
        strategy = ScanCommon.FixedSeedStrategy(seed_8)

        result = if effective_solver_backend === :models
            _solve_with_models(mode, T_fm;
                xi=xi,
                model_kind=model_kind,
                seed_strategy=strategy,
                p_num=p_num,
                t_num=t_num,
                nlsolve_kwargs...)
        elseif mode isa FixedAsymmetricRho
            solve(mode, T_fm;
                xi=xi,
                model_kind=model_kind,
                seed_strategy=strategy,
                p_num=p_num,
                t_num=t_num,
                nlsolve_kwargs...
            )
        else
            solve(mode, T_fm;
                xi=xi,
                seed_strategy=strategy,
                p_num=p_num,
                t_num=t_num,
                nlsolve_kwargs...
            )
        end

        result = finalize_solver_result(result, T_fm, xi;
            solver_backend=effective_solver_backend,
            p_num=p_num,
            t_num=t_num,
            model_kind=model_kind,
        )
        return result, ""
    catch err
        msg = sprint() do io
            showerror(io, err)
        end
        return nothing, _clean_message(msg)
    end
end

"""精炼近似收敛的结果"""
function _refine_result(T_fm, ρ_fm, xi, result;
    constraint_mode::Symbol=:fixed_rho,
    asym_ud_ratio_target::Float64=0.876,
    asym_s_target::Float64=0.0,
    solver_backend::Symbol=:legacy,
    model_kind::Symbol=:PNJL,
    p_num,
    t_num,
    nlsolve_kwargs...)

    return refine_near_converged(result;
        acceptable_residual=ACCEPTABLE_RESIDUAL,
        solve_again=seed -> _solve_point(T_fm, ρ_fm, xi, seed;
            constraint_mode=constraint_mode,
            asym_ud_ratio_target=asym_ud_ratio_target,
            asym_s_target=asym_s_target,
            solver_backend=solver_backend,
            model_kind=model_kind,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...,
        ),
    )
end

"""判断是否成功"""
_is_success(result) = is_success(result; acceptable_residual=ACCEPTABLE_RESIDUAL)

"""提升近似收敛为成功"""
_promote_success(result) = promote_near_converged(result; acceptable_residual=ACCEPTABLE_RESIDUAL)

@inline function _models_mode(mode::FixedRho)
    return Main.Models.FixedRho(mode.rho_target)
end

@inline function _models_mode(mode::FixedAsymmetricRho)
    return Main.Models.FixedAsymmetricRho(mode.rho_target, mode.ud_ratio_target, mode.s_target)
end

function _to_solver_result(mode::ConstraintMode, result, xi::Real)
    return SolverResult(
        mode,
        Bool(result.converged),
        Float64.(result.solution),
        SVector{5, Float64}(Tuple(Float64.(result.x_state))),
        SVector{3, Float64}(Tuple(Float64.(result.mu_vec))),
        Float64(result.omega),
        Float64(result.pressure),
        Float64(result.rho_norm),
        Float64(result.entropy),
        Float64(result.energy),
        SVector{3, Float64}(Tuple(Float64.(result.masses))),
        Int(result.iterations),
        Float64(result.residual_norm),
        Float64(xi),
    )
end

function _solve_with_models(mode::ConstraintMode, T_fm;
    xi::Real,
    model_kind::Symbol,
    seed_strategy::SeedStrategy,
    p_num::Int,
    t_num::Int,
    nlsolve_kwargs...)
    model = Main.Models.create_model(model_kind)
    mapped_mode = _models_mode(mode)
    seed_guess = get_seed(seed_strategy, [T_fm], mode)
    raw = Main.Models.solve_constraint(
        model,
        mapped_mode,
        T_fm;
        seed_guess=seed_guess,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        nlsolve_kwargs...,
    )
    return _to_solver_result(mode, raw, xi)
end

"""写入一行结果"""
function _write_row(io, T, rho, xi, result, message;
    model_kind::Symbol=:PNJL,
    p_num::Int=24,
    t_num::Int=8,
)
    if result === nothing
        values = (
            _fmt(T), _fmt(rho), _fmt(xi),
            "NaN", "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN",
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
    mu_B = mu_vec_mev[1] + 2 * mu_vec_mev[2]
    mu_Q = mu_vec_mev[1] - mu_vec_mev[2]
    mu_S = mu_vec_mev[2] - mu_vec_mev[3]

    T_fm = T / ħc_MeV_fm
    rho_vec = Main.Models.model_rho(
        Main.Models.create_model(model_kind),
        result.x_state,
        result.mu_vec,
        T_fm;
        p_num=p_num,
        t_num=t_num,
        xi=xi,
    )
    rho_u, rho_d, rho_s = rho_vec[1], rho_vec[2], rho_vec[3]

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
        _fmt(mu_B),
        _fmt(mu_Q),
        _fmt(mu_S),
        _fmt(result.pressure),
        _fmt(result.entropy),
        _fmt(result.energy),
        _fmt(rho_u),
        _fmt(rho_d),
        _fmt(rho_s),
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
# 格式化辅助（复用 ScanCommon）
# ============================================================================

const _fmt = ScanCommon.fmt
const _clean_message = ScanCommon.clean_message
const _quote = ScanCommon.quote_csv
const _join_messages = ScanCommon.join_messages
const _format_candidate_failure = ScanCommon.format_candidate_failure

end # module TrhoScan
