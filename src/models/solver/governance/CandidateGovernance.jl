"""
    Candidate governance context and adapters.
"""

const HardRule = Function
const CandidateSelector = Function
const SELECTOR_INPUT_REQUIRED_FIELDS = (
    :converged,
    :pressure,
    :residual_norm,
    :hard_constraint_ok,
    :failed_constraints,
    :seed_index,
    :selection_score,
    :quality_tag,
)

@inline function _candidate_seed_index(candidate, fallback::Int)::Int
    return hasproperty(candidate, :seed_index) ? Int(getproperty(candidate, :seed_index)) : fallback
end

@inline function _candidate_attempt_origin(candidate)::Symbol
    if hasproperty(candidate, :governed_attempt_origin)
        return Symbol(getproperty(candidate, :governed_attempt_origin))
    elseif hasproperty(candidate, :fixedrho_joint_attempt_origin)
        return Symbol(getproperty(candidate, :fixedrho_joint_attempt_origin))
    elseif hasproperty(candidate, :attempt_origin)
        return Symbol(getproperty(candidate, :attempt_origin))
    end
    return :unknown
end

@inline function _candidate_selection_score(candidate)::Float64
    if hasproperty(candidate, :selection_score)
        val = Float64(getproperty(candidate, :selection_score))
        return isfinite(val) ? val : NaN
    elseif hasproperty(candidate, :pressure)
        val = Float64(getproperty(candidate, :pressure))
        return isfinite(val) ? val : -Inf
    end
    return -Inf
end

@inline function governance_quality_tag(candidate;
    residual_norm_max::Real=1e-6,
    require_converged::Bool=true,
)::Symbol
    converged = Bool(get(candidate, :converged, false))
    residual = Float64(get(candidate, :residual_norm, Inf))
    hard_ok = Bool(get(candidate, :hard_constraint_ok, false))
    residual_ok = isfinite(residual) && residual <= Float64(residual_norm_max)
    converged_ok = !require_converged || converged
    attempt_origin = _candidate_attempt_origin(candidate)
    is_fallback_attempt = attempt_origin == :method_rescue || attempt_origin == :fallback

    if converged_ok && residual_ok && hard_ok
        return is_fallback_attempt ? :fallback : :good
    end
    if converged || residual_ok || hard_ok
        return :degraded
    end
    return :bad
end

@inline function classify_attempt_error(err)::Symbol
    err isa InterruptException && return :interrupt
    err isa ArgumentError && return :constraint_error
    err isa DomainError && return :constraint_error
    return :solver_error
end

@inline function normalize_error_message(err; max_chars::Int=240)::String
    msg = sprint(showerror, err)
    msg = replace(strip(msg), '\n' => ' ', '\r' => ' ')
    if length(msg) <= max_chars
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
)
    hasproperty(candidate, :converged) || throw(ArgumentError("governance candidate missing field :converged"))
    hasproperty(candidate, :pressure) || throw(ArgumentError("governance candidate missing field :pressure"))
    hasproperty(candidate, :residual_norm) || throw(ArgumentError("governance candidate missing field :residual_norm"))
    hasproperty(candidate, :hard_constraint_ok) || throw(ArgumentError("governance candidate missing field :hard_constraint_ok"))
    hasproperty(candidate, :failed_constraints) || throw(ArgumentError("governance candidate missing field :failed_constraints"))

    pressure = Float64(getproperty(candidate, :pressure))
    residual_norm = Float64(getproperty(candidate, :residual_norm))
    converged = Bool(getproperty(candidate, :converged))
    hard_ok = Bool(getproperty(candidate, :hard_constraint_ok))
    failed_constraints = Symbol.(getproperty(candidate, :failed_constraints))

    if hard_ok && !isempty(failed_constraints)
        throw(ArgumentError("governance candidate hard_constraint_ok=true requires empty failed_constraints"))
    end
    if !hard_ok && isempty(failed_constraints)
        throw(ArgumentError("governance candidate hard_constraint_ok=false requires non-empty failed_constraints"))
    end

    score = _candidate_selection_score(candidate)
    isfinite(score) || (score = NaN)

    return (
        converged=converged,
        pressure=(isfinite(pressure) ? pressure : -Inf),
        residual_norm=(isfinite(residual_norm) ? residual_norm : Inf),
        hard_constraint_ok=hard_ok,
        failed_constraints=failed_constraints,
        seed_index=seed_index,
        selection_score=score,
    )
end

@inline function build_governance_candidate(
    raw_candidate;
    hard_constraint_ok::Bool,
    failed_constraints,
    seed_index::Int,
    residual_norm_max::Real=1e-6,
    error_kind::Symbol=:none,
    error_msg::AbstractString="",
)
    candidate = (
        ; raw_candidate...,
        hard_constraint_ok=hard_constraint_ok,
        failed_constraints=Symbol.(failed_constraints),
    )
    normalized = normalize_governance_candidate(candidate; seed_index=seed_index)
    quality_tag = governance_quality_tag((; candidate...,
        converged=normalized.converged,
        residual_norm=normalized.residual_norm,
        hard_constraint_ok=normalized.hard_constraint_ok,
    ); residual_norm_max=residual_norm_max)
    return (
        ; candidate...,
        converged=normalized.converged,
        pressure=normalized.pressure,
        residual_norm=normalized.residual_norm,
        hard_constraint_ok=normalized.hard_constraint_ok,
        failed_constraints=normalized.failed_constraints,
        seed_index=normalized.seed_index,
        selection_score=normalized.selection_score,
        quality_tag=quality_tag,
        error_kind=error_kind,
        error_msg=String(error_msg),
    )
end

@inline function _ensure_selector_input_contract(candidate, index::Int)
    for field in SELECTOR_INPUT_REQUIRED_FIELDS
        hasproperty(candidate, field) || throw(ArgumentError(
            "selector input candidate #$(index) missing required field :$(field)"
        ))
    end
    return nothing
end

function normalize_selector_candidates(candidates::AbstractVector{<:NamedTuple};
    residual_norm_max::Real=1e-6,
    require_converged::Bool=true,
)
    isempty(candidates) && throw(ArgumentError("candidates must be non-empty"))
    normalized = NamedTuple[]
    for (idx, cand) in enumerate(candidates)
        seed_index = _candidate_seed_index(cand, idx)
        base = normalize_governance_candidate(cand; seed_index=seed_index)
        quality_tag = governance_quality_tag((; cand...,
            converged=base.converged,
            residual_norm=base.residual_norm,
            hard_constraint_ok=base.hard_constraint_ok,
        ); residual_norm_max=residual_norm_max, require_converged=require_converged)
        push!(normalized, (; cand...,
            converged=base.converged,
            pressure=base.pressure,
            residual_norm=base.residual_norm,
            hard_constraint_ok=base.hard_constraint_ok,
            failed_constraints=base.failed_constraints,
            seed_index=base.seed_index,
            selection_score=base.selection_score,
            quality_tag=quality_tag,
        ))
    end
    return normalized
end

function execute_governance_selector(candidates::AbstractVector{<:NamedTuple};
    selector::Function,
    residual_norm_max::Real=1e-6,
    require_converged::Bool=true,
)
    normalized = normalize_selector_candidates(candidates;
        residual_norm_max=residual_norm_max,
        require_converged=require_converged,
    )
    for (idx, candidate) in enumerate(normalized)
        _ensure_selector_input_contract(candidate, idx)
    end
    selected = selector(normalized)
    hasproperty(selected, :selected_index) || throw(ArgumentError("selector must return field :selected_index"))
    hasproperty(selected, :selected_candidate) || throw(ArgumentError("selector must return field :selected_candidate"))
    hasproperty(selected, :selection_reason) || throw(ArgumentError("selector must return field :selection_reason"))
    _ensure_selector_input_contract(getproperty(selected, :selected_candidate), Int(getproperty(selected, :selected_index)))
    return (; selected..., normalized_candidates=normalized)
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
