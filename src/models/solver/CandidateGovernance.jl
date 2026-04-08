"""
    Candidate governance context and adapters.
"""

const HardRule = Function
const CandidateSelector = Function

@inline function classify_attempt_error(err)::Symbol
    err isa InterruptException && return :interrupt
    err isa ArgumentError && return :constraint_error
    err isa DomainError && return :constraint_error
    return :solver_error
end

@inline function normalize_error_message(err; max_chars::Int=240)::String
    msg = sprint(showerror, err)
    msg = replace(strip(msg), '\n' => ' ', '\r' => ' ')
    if ncodeunits(msg) <= max_chars
        return msg
    end
    return first(msg, max_chars)
end

@inline function _seed_pool_key(seed_vec::AbstractVector{<:Real})
    return join(round.(Float64.(seed_vec); digits=12), ",")
end

@inline function _push_seed_pool_entry!(
    pool::Vector{NamedTuple{(:seed, :source), Tuple{Vector{Float64}, Symbol}}},
    seen::Set{String},
    seed::AbstractVector{<:Real},
    source::Symbol,
    mode::ConstraintMode,
    seed_extend::Function,
)
    normalized = Float64.(seed_extend(Float64.(seed), mode))
    key = _seed_pool_key(normalized)
    key in seen && return nothing
    push!(pool, (seed=normalized, source=source))
    push!(seen, key)
    return nothing
end

function build_seed_pool(
    mode::ConstraintMode;
    primary_seed::Union{Nothing, AbstractVector{<:Real}}=nothing,
    extra_seed_pool=(),
    provided_seed_pool=(),
    default_seed_pool=(),
    seed_extend::Function=(seed, _mode) -> Float64.(seed),
)
    pool = NamedTuple{(:seed, :source), Tuple{Vector{Float64}, Symbol}}[]
    seen = Set{String}()

    if primary_seed !== nothing
        _push_seed_pool_entry!(pool, seen, primary_seed, :primary, mode, seed_extend)
    end
    for seed in extra_seed_pool
        _push_seed_pool_entry!(pool, seen, seed, :extra, mode, seed_extend)
    end
    for seed in provided_seed_pool
        _push_seed_pool_entry!(pool, seen, seed, :provided, mode, seed_extend)
    end
    for seed in default_seed_pool
        _push_seed_pool_entry!(pool, seen, seed, :default, mode, seed_extend)
    end
    return pool
end

@inline function evaluate_candidate_success(candidate; residual_norm_max::Real=1e-6)
    if hasproperty(candidate, :hard_constraint_ok)
        return Bool(getproperty(candidate, :hard_constraint_ok))
    end
    Bool(get(candidate, :converged, false)) && return true
    residual = Float64(get(candidate, :residual_norm, Inf))
    return isfinite(residual) && residual <= Float64(residual_norm_max)
end

@inline function normalize_governance_candidate(candidate;
    seed_index::Int,
    residual_norm_max::Real=1e-6,
    failed_default::Symbol=:residual_too_large,
)
    pressure = Float64(get(candidate, :pressure, -Inf))
    residual_norm = Float64(get(candidate, :residual_norm, Inf))
    converged = Bool(get(candidate, :converged, false))
    hard_ok = if hasproperty(candidate, :hard_constraint_ok)
        Bool(getproperty(candidate, :hard_constraint_ok))
    else
        evaluate_candidate_success(candidate; residual_norm_max=residual_norm_max)
    end
    failed_constraints = if hasproperty(candidate, :failed_constraints)
        Symbol.(getproperty(candidate, :failed_constraints))
    elseif hard_ok
        Symbol[]
    elseif !isfinite(residual_norm) && !converged
        Symbol[:solver_failed]
    elseif failed_default == :solver_failed
        Symbol[:solver_failed]
    else
        Symbol[failed_default]
    end
    return (
        converged=converged,
        pressure=(isfinite(pressure) ? pressure : -Inf),
        residual_norm=(isfinite(residual_norm) ? residual_norm : Inf),
        hard_constraint_ok=hard_ok,
        failed_constraints=failed_constraints,
        seed_index=seed_index,
    )
end

@inline function build_candidate_context(
    mode::ConstraintMode;
    continuity_seed_available::Bool=false,
    phase_hint::Symbol=:unknown,
    residual_norm_max::Real=1e-6,
    prefer_continuity::Bool=false,
)
    return (
        mode=mode,
        continuity_seed_available=Bool(continuity_seed_available),
        phase_hint=phase_hint,
        residual_norm_max=Float64(residual_norm_max),
        prefer_continuity=Bool(prefer_continuity),
    )
end

"""
    execute_attempt_pool(attempts; evaluate_attempt, on_error,
                         stop_on_first_success=true,
                         evaluate_all_attempts=false)

统一 attempt 执行引擎。

# 输入契约
- `attempts`：attempt plan，可包含 seed/method/fallback 等 metadata。
- `evaluate_attempt(attempt_cfg, attempt_index)`：返回 `(candidate, success::Bool)`。
- `on_error(attempt_cfg, attempt_index, err)`：异常降级回调，返回 `(candidate, success::Bool)`。

# 执行语义
- `InterruptException` 直接重抛。
- 其余异常交由 `on_error` 处理。
- 当 `evaluate_all_attempts=false` 且 `stop_on_first_success=true` 时，首个成功候选后提前停止。

# 输出契约
- 返回 `Vector{NamedTuple}`，顺序与已执行 attempt 顺序一致。
- 候选字段由调用方定义，但需保证后续治理选择所需字段可用。
"""
function execute_attempt_pool(
    attempts;
    evaluate_attempt::Function,
    on_error::Function,
    stop_on_first_success::Bool=true,
    evaluate_all_attempts::Bool=false,
)
    candidates = NamedTuple[]
    for (attempt_index, attempt_cfg) in enumerate(attempts)
        candidate, success = try
            evaluate_attempt(attempt_cfg, attempt_index)
        catch err
            err isa InterruptException && rethrow()
            on_error(attempt_cfg, attempt_index, err)
        end
        push!(candidates, candidate)
        if !evaluate_all_attempts && stop_on_first_success && Bool(success)
            break
        end
    end
    return candidates
end
