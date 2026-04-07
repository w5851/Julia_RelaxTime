"""
    Candidate governance context and adapters.
"""

const HardRule = Function
const CandidateSelector = Function

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
