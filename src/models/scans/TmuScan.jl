"""
    TmuScan

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

# 使用相变感知策略（推荐用于一阶相变区域）
result = run_tmu_scan(
    T_values = 50.0:10.0:130.0,
    mu_values = 200.0:5.0:400.0,
    xi_values = [0.0],
    use_phase_aware = true
)
```
"""
module TmuScan

using Printf
using StaticArrays

# 导入新架构模块
using Main.Constants_PNJL: ħc_MeV_fm
using ..Models: FixedMu, ConstraintMode, solve, SolverResult, is_physical_solution
using ..Models: build_seed_pool, coerce_solver_result, create_model, solve_constraint
using ..Models: PNJLIntegrals
using ..SeedStrategies: SeedStrategy, DefaultSeed, MultiSeed
using ..SeedStrategies: get_seed, get_all_seeds, HADRON_SEED_5, QUARK_SEED_5, HT_GUESS_0p8_SEED_5, HT_GUESS_0p9_SEED_5, HT_GUESS_0p95_SEED_5, WEAK_CHIRAL_CONF_SEED_5
using ..SeedStrategies: auto_phase_hint
using ..ScanCommon
using ..ScanConfig: TmuScanConfig, scan_kwargs
using ..ScanResultFinalize: finalize_solver_result, promote_near_converged, is_success, refine_near_converged

export run_tmu_scan, DEFAULT_T_VALUES, DEFAULT_MU_VALUES, DEFAULT_OUTPUT_PATH

# ============================================================================
# 默认配置
# ============================================================================

const DEFAULT_T_VALUES = collect(50.0:10.0:200.0)
const DEFAULT_MU_VALUES = collect(0.0:10.0:400.0)
const DEFAULT_OUTPUT_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "outputs", "results", "pnjl", "scan", "tmu", "tmu_scan.csv"))
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

@inline function _validate_real_vector(name::Symbol, values)
    values isa AbstractVector{<:Real} || throw(ArgumentError("$(name) must be AbstractVector{<:Real}, got $(typeof(values))"))
    isempty(values) && throw(ArgumentError("$(name) must not be empty"))
    for (i, v) in pairs(values)
        isfinite(Float64(v)) || throw(ArgumentError("$(name)[$(i)] must be finite Real, got $(v)"))
    end
    return nothing
end

@inline function _validate_tmu_scan_inputs(T_values, mu_values, xi_values, solver_backend::Symbol, model_kind::Symbol)
    _validate_real_vector(:T_values, T_values)
    _validate_real_vector(:mu_values, mu_values)
    _validate_real_vector(:xi_values, xi_values)

    (solver_backend === :models || solver_backend === :auto) ||
        throw(ArgumentError("solver_backend must be :models or :auto (legacy backend removed), got $(solver_backend)"))

    if model_kind in ScanCommon.PARAMETERIZED_MODEL_KIND_ALIASES
        throw(ArgumentError("model_kind=:pnjl_aniso is not supported; use model_kind=:PNJL with profile/xi parameterization"))
    end

    model_kind in ScanCommon.SUPPORTED_SCAN_MODEL_KINDS ||
        throw(ArgumentError("model_kind must be one of $(ScanCommon.SUPPORTED_SCAN_MODEL_KINDS), got $(model_kind)"))

    return nothing
end

@inline function _validate_semantic_mode(semantic_mode::Symbol, selector::Union{Nothing, Function})
    (semantic_mode === :ground_state || semantic_mode === :constrained_manifold) ||
        throw(ArgumentError("semantic_mode must be :ground_state or :constrained_manifold, got $(semantic_mode)"))

    semantic_mode === :ground_state ||
        throw(ArgumentError("run_tmu_scan currently supports semantic_mode=:ground_state only for FixedMu workflow"))

    selector === nothing ||
        throw(ArgumentError("run_tmu_scan does not support custom selector; selector must be nothing"))
    return nothing
end

@inline function _validate_auto_pnjl_backend(auto_pnjl_backend::Symbol)
    (auto_pnjl_backend === :models || auto_pnjl_backend === :legacy) ||
        throw(ArgumentError("auto_pnjl_backend must be :models or :legacy, got $(auto_pnjl_backend)"))
    return nothing
end

@inline function _effective_solver_backend(solver_backend::Symbol, model_kind::Symbol; auto_pnjl_backend::Symbol=:models)::Symbol
    if solver_backend !== :auto
        return solver_backend
    end
    if model_kind === :PNJL
        return auto_pnjl_backend
    end
    return :models
end

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
- `use_phase_aware`: 是否使用相变感知策略，默认 true
- `p_num`, `t_num`: 积分节点数
- `progress_cb`: 进度回调函数 `(point, result) -> nothing`

# 返回
NamedTuple 包含：
- `total`: 总点数
- `success`: 成功点数
- `failure`: 失败点数
- `skipped`: 跳过点数
- `output`: 输出文件路径

# 多初值取优策略
扫描过程默认采用“候选初值全评估 + 治理选优”：
1. 优先尝试连续性种子（若有）
2. 追加多初值候选（`MultiSeed`）
3. 对全部候选求解后按治理规则统一选择

说明：`use_phase_aware` 参数为兼容保留位，不再启用专用相变感知状态机策略。
"""
function run_tmu_scan(;
    T_values=DEFAULT_T_VALUES,
    mu_values=DEFAULT_MU_VALUES,
    xi_values=[0.0],
    output_path::AbstractString=DEFAULT_OUTPUT_PATH,
    overwrite::Bool=false,
    resume::Bool=true,
    use_phase_aware::Bool=true,
    bootstrap_multiseed::Bool=true,
    solver_backend::Symbol=:auto,
    auto_pnjl_backend::Symbol=:models,
    semantic_mode::Symbol=:ground_state,
    selector::Union{Nothing, Function}=nothing,
    model_kind::Symbol=:PNJL,
    p_num::Int=24,
    t_num::Int=8,
    progress_cb::Union{Nothing, Function}=nothing,
    diagnostic_level::Symbol=:none,
    thermo_quadrature_policy::Symbol=:tensor_gauss,
    thermo_quadrature_rtol::Float64=1e-8,
    thermo_quadrature_atol::Float64=1e-10,
    thermo_quadrature_maxevals::Int=10^7,
    nlsolve_kwargs...
)
    _validate_tmu_scan_inputs(T_values, mu_values, xi_values, solver_backend, model_kind)
    _validate_semantic_mode(semantic_mode, selector)
    _validate_auto_pnjl_backend(auto_pnjl_backend)
    model_kind === :PNJL && foreach(PNJLIntegrals.validate_rs_anisotropy, xi_values)
    PNJLIntegrals.validate_thermal_quadrature_policy(thermo_quadrature_policy)
    if thermo_quadrature_policy === :rs_reduced_adaptive && model_kind !== :PNJL
        throw(ArgumentError("rs_reduced_adaptive thermal quadrature is supported only for model_kind=:PNJL"))
    end
    PNJLIntegrals.validate_thermal_quadrature_controls(
        thermo_quadrature_rtol,
        thermo_quadrature_atol,
        thermo_quadrature_maxevals,
    )
    solver_kwargs = (
        thermo_quadrature_policy=thermo_quadrature_policy,
        thermo_quadrature_rtol=thermo_quadrature_rtol,
        thermo_quadrature_atol=thermo_quadrature_atol,
        thermo_quadrature_maxevals=thermo_quadrature_maxevals,
        nlsolve_kwargs...,
    )

    mkpath(dirname(output_path))
    completed = (resume && !overwrite && isfile(output_path)) ? ScanCommon.load_completed_keys3(output_path; digits=6) : Set{NTuple{3, Float64}}()
    io_mode = (overwrite || !isfile(output_path)) ? "w" : "a"

    stats = Dict(:total => 0, :success => 0, :failure => 0, :skipped => 0)
    
    # 连续性种子缓存（按 T, xi 分组）
    continuation_seeds = Dict{Tuple{Float64, Float64}, Vector{Float64}}()

    open(output_path, io_mode) do io
        if io_mode == "w"
            println(io, HEADER)
        end

        for xi in xi_values
            for T in T_values
                # 每个新温度重置连续性缓存
                seed_key = _seed_continuation_key(T, xi)
                delete!(continuation_seeds, seed_key)
                
                for mu in mu_values
                    stats[:total] += 1
                    key = ScanCommon.key3(T, mu, xi; digits=6)
                    
                    if key in completed
                        stats[:skipped] += 1
                        continue
                    end

                    # 转换单位：MeV -> fm⁻¹
                    T_fm = T / ħc_MeV_fm
                    μ_fm = mu / ħc_MeV_fm

                    # 构建候选并统一选优
                    candidates = _build_seed_candidates_v2(
                        continuation_seeds, T, mu, xi, T_fm, μ_fm;
                        bootstrap_multiseed=bootstrap_multiseed,
                    )

                    result, message = _attempt_with_candidates(T_fm, μ_fm, xi, candidates;
                        solver_backend=solver_backend,
                        auto_pnjl_backend=auto_pnjl_backend,
                        model_kind=model_kind,
                        diagnostic_level=diagnostic_level,
                        p_num=p_num,
                        t_num=t_num,
                        solver_kwargs...)

                    # 更新连续性缓存
                    if result !== nothing && _is_success(result)
                        seed_key = _seed_continuation_key(T, xi)
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

"""run_tmu_scan(config::TmuScanConfig; kwargs...) -> NamedTuple

结构化配置入口。`kwargs` 会覆盖 config 中同名项，保持与原 kwargs 接口兼容。
"""
function run_tmu_scan(config::TmuScanConfig; kwargs...)
    cfg = scan_kwargs(config)
    return run_tmu_scan(; cfg..., config.nlsolve_kwargs..., kwargs...)
end

# ============================================================================
# 内部辅助函数
# ============================================================================

_seed_continuation_key(T, xi) = ScanCommon.key2(T, xi; digits=SEED_KEY_DIGITS)

@inline function _push_unique_candidate!(
    candidates::Vector{NamedTuple{(:label, :state), Tuple{String, Vector{Float64}}}},
    seen::Set{String},
    label::String,
    state::AbstractVector{<:Real},
)
    key = join(round.(Float64.(state); digits=8), ",")
    key in seen && return nothing
    push!(seen, key)
    push!(candidates, (label=label, state=Float64.(state)))
    return nothing
end

@inline function _seed_pool_source_to_label(source::Symbol, idx::Int)
    if source == :primary
        return "continuation"
    elseif source == :provided
        return "multiseed_$(idx)"
    elseif source == :extra
        return "fallback_$(idx)"
    end
    return "default_$(idx)"
end

"""构建初值候选列表（连续性 + 多初值全评估）"""
function _build_seed_candidates_v2(cache::Dict, T, mu, xi, T_fm, μ_fm; bootstrap_multiseed::Bool)
    provided_seed_pool = Vector{Vector{Float64}}()
    default_seed_pool = Vector{Vector{Float64}}()
    primary_seed = nothing

    # 1. 连续性种子（最优先）
    seed_key = _seed_continuation_key(T, xi)
    if haskey(cache, seed_key)
        primary_seed = Float64.(cache[seed_key])
    end

    # 2. 多初值候选（可选）
    if bootstrap_multiseed
        multiseeds = get_all_seeds(MultiSeed(), [T_fm, μ_fm], FixedMu())
        append!(provided_seed_pool, [Float64.(seed) for seed in multiseeds])
    end
    
    # 3. 基于相位启发式的基础候选（兜底）
    hint = auto_phase_hint(T_fm, μ_fm)
    if hint === :quark
        push!(default_seed_pool, Float64.(QUARK_SEED_5))
        push!(default_seed_pool, Float64.(HADRON_SEED_5))
    else
        push!(default_seed_pool, Float64.(HADRON_SEED_5))
        push!(default_seed_pool, Float64.(QUARK_SEED_5))
    end

    # 4. 温区兜底候选（高温/弱手征场景）
    append!(default_seed_pool, (
        Float64.(WEAK_CHIRAL_CONF_SEED_5),
        Float64.(HT_GUESS_0p8_SEED_5),
        Float64.(HT_GUESS_0p9_SEED_5),
        Float64.(HT_GUESS_0p95_SEED_5),
    ))

    if primary_seed === nothing
        if !isempty(provided_seed_pool)
            primary_seed = provided_seed_pool[1]
            provided_seed_pool = provided_seed_pool[2:end]
        elseif !isempty(default_seed_pool)
            primary_seed = default_seed_pool[1]
            default_seed_pool = default_seed_pool[2:end]
        else
            primary_seed = Float64.(HADRON_SEED_5)
        end
    end

    seed_pool = build_seed_pool(FixedMu();
        primary_seed=primary_seed,
        provided_seed_pool=provided_seed_pool,
        default_seed_pool=default_seed_pool,
        seed_extend=(seed, _) -> Float64.(seed),
    )

    candidates = NamedTuple{(:label, :state), Tuple{String, Vector{Float64}}}[]
    seen_states = Set{String}()
    for (idx, entry) in enumerate(seed_pool)
        label = _seed_pool_source_to_label(entry.source, idx)
        _push_unique_candidate!(candidates, seen_states, label, entry.seed)
    end

    return candidates
end

"""尝试多个初值候选"""
function _attempt_with_candidates(T_fm, μ_fm, xi, candidates;
    solver_backend::Symbol=:auto,
    auto_pnjl_backend::Symbol=:models,
    model_kind::Symbol=:PNJL,
    diagnostic_level::Symbol=:none,
    p_num,
    t_num,
    nlsolve_kwargs...)
    return ScanCommon.attempt_with_candidates(candidates;
        solve_point=seed_state -> _solve_point(T_fm, μ_fm, xi, seed_state;
            solver_backend=solver_backend,
            auto_pnjl_backend=auto_pnjl_backend,
            model_kind=model_kind,
            diagnostic_level=diagnostic_level,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...,
        ),
        refine=result -> _refine_result(T_fm, μ_fm, xi, result;
            solver_backend=solver_backend,
            auto_pnjl_backend=auto_pnjl_backend,
            model_kind=model_kind,
            diagnostic_level=diagnostic_level,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...,
        ),
        promote=_promote_success,
        is_success=_is_success,
        stop_on_first_success=false,
    )
end

"""单点求解"""
function _solve_point(T_fm, μ_fm, xi, seed_state;
    solver_backend::Symbol=:auto,
    auto_pnjl_backend::Symbol=:models,
    model_kind::Symbol=:PNJL,
    diagnostic_level::Symbol=:none,
    p_num,
    t_num,
    nlsolve_kwargs...)
    try
        effective_solver_backend = _effective_solver_backend(solver_backend, model_kind; auto_pnjl_backend=auto_pnjl_backend)
        effective_solver_backend === :models ||
            throw(ArgumentError("solver_backend=:legacy has been removed from TmuScan models path; use :models or :auto"))

        # 创建固定种子策略
        seed_5 = Float64.(seed_state[1:min(5, length(seed_state))])
        strategy = ScanCommon.FixedSeedStrategy(seed_5)
        
        result = _solve_with_models(FixedMu(), T_fm, μ_fm;
            xi=xi,
            model_kind=model_kind,
            diagnostic_level=diagnostic_level,
            seed_strategy=strategy,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...)

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
function _refine_result(T_fm, μ_fm, xi, result;
    solver_backend::Symbol=:auto,
    auto_pnjl_backend::Symbol=:models,
    model_kind::Symbol=:PNJL,
    diagnostic_level::Symbol=:none,
    p_num,
    t_num,
    nlsolve_kwargs...)

    return refine_near_converged(result;
        acceptable_residual=ACCEPTABLE_RESIDUAL,
        solve_again=seed -> _solve_point(T_fm, μ_fm, xi, seed;
            solver_backend=solver_backend,
            auto_pnjl_backend=auto_pnjl_backend,
            model_kind=model_kind,
            diagnostic_level=diagnostic_level,
            p_num=p_num,
            t_num=t_num,
            nlsolve_kwargs...,
        ),
    )
end

"""判断是否成功"""
function _is_success(result)
    result === nothing && return false
    thermo_finite = isfinite(result.omega) && isfinite(result.pressure) &&
                   isfinite(result.rho_norm) && isfinite(result.entropy) &&
                   isfinite(result.energy)
    phys_ok = thermo_finite && is_physical_solution(result.x_state, result.masses)
    return phys_ok && is_success(result; acceptable_residual=ACCEPTABLE_RESIDUAL)
end

"""提升近似收敛为成功"""
_promote_success(result) = promote_near_converged(result; acceptable_residual=ACCEPTABLE_RESIDUAL)

@inline function _models_mode(mode::FixedMu)
    return FixedMu()
end

function _to_solver_result(mode::ConstraintMode, result, xi::Real)
    return coerce_solver_result(mode, result; xi_override=xi)
end

@inline function _reject_legacy_solver_kwargs(nlsolve_kwargs)
    legacy_switches = (:use_problem_spec, :allow_legacy_path, :warn_on_legacy_path)
    for key in keys(nlsolve_kwargs)
        if key in legacy_switches || key === :problem_spec
            throw(ArgumentError("legacy solver switch '$key' is not allowed from TmuScan models path"))
        end
    end
    return nothing
end

function _solve_with_models(mode::ConstraintMode, T_fm, μ_fm;
    xi::Real,
    model_kind::Symbol,
    diagnostic_level::Symbol=:none,
    seed_strategy::SeedStrategy,
    p_num::Int,
    t_num::Int,
    nlsolve_kwargs...)
    model = create_model(model_kind)
    mapped_mode = _models_mode(mode)
    seed_guess = get_seed(seed_strategy, [T_fm, μ_fm], mode)
    _reject_legacy_solver_kwargs(nlsolve_kwargs)
    models_kwargs = (; (k => v for (k, v) in nlsolve_kwargs if k in (
        :solver,
        :residual_norm_max,
        :physicality_check,
        :thermo_quadrature_policy,
        :thermo_quadrature_rtol,
        :thermo_quadrature_atol,
        :thermo_quadrature_maxevals,
    ))...)
    raw = solve_constraint(
        model,
        mapped_mode,
        T_fm;
        μ_fm=μ_fm,
        seed_guess=seed_guess,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        diagnostic_level=diagnostic_level,
        models_kwargs...,
    )
    return _to_solver_result(mode, raw, xi)
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
# 格式化辅助（复用 ScanCommon）
# ============================================================================

const _fmt = ScanCommon.fmt
const _clean_message = ScanCommon.clean_message
const _quote = ScanCommon.quote_csv
const _join_messages = ScanCommon.join_messages
const _format_candidate_failure = ScanCommon.format_candidate_failure

end # module TmuScan
