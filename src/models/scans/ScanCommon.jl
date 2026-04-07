module ScanCommon

using Printf

import Main.Models: ConstraintMode
using ..SeedStrategies: SeedStrategy
import ..SeedStrategies: get_seed
import Main.Models: SolverResult

export fmt
export clean_message
export quote_csv
export join_messages
export format_candidate_failure
export attempt_with_candidates
export FixedSeedStrategy
export key3
export key2
export load_completed_keys3
export SUPPORTED_SCAN_MODEL_KINDS
export PARAMETERIZED_MODEL_KIND_ALIASES

const SUPPORTED_SCAN_MODEL_KINDS = (:PNJL, :NJL, :RPNJL, :PNJLMagnetic, :Rotation, :GasLiquid)
const PARAMETERIZED_MODEL_KIND_ALIASES = (:pnjl_aniso, :PNJL_ANISO)

fmt(x::Float64) = @sprintf("%.6f", x)
fmt(x::Real) = fmt(Float64(x))
fmt(x) = string(x)

function clean_message(msg::AbstractString)
    stripped = replace(strip(msg), '\n' => ' ')
    return replace(stripped, '"' => '\'')
end

quote_csv(msg::AbstractString) = isempty(msg) ? "" : string('"', msg, '"')
quote_csv(::Nothing) = ""

@inline function _normalized_fallback_reason(msg::AbstractString)
    if occursin("weighted-block fallback rescued", msg)
        return "weighted_block_rescued"
    elseif occursin("seed[", msg) && occursin(" failed", msg)
        return "seed_failed"
    end
    return ""
end

@inline function _message_rank(msg::AbstractString)
    if startswith(msg, "governance.selection=")
        return 1
    elseif startswith(msg, "fallback.reason=")
        return 2
    end
    return 3
end

function join_messages(messages)
    filtered = filter(!isempty, String.(messages))
    isempty(filtered) && return ""

    normalized = String[]
    for msg in filtered
        if startswith(msg, "governance.selection=")
            push!(normalized, msg)
            continue
        end
        fallback_reason = _normalized_fallback_reason(msg)
        if !isempty(fallback_reason)
            push!(normalized, "fallback.reason=$(fallback_reason)")
            continue
        end
        push!(normalized, msg)
    end

    sort!(normalized; by=m -> (_message_rank(m), m))
    unique_msgs = unique(normalized)
    return join(unique_msgs, " | ")
end

function format_candidate_failure(label, message, result)
    base = "seed[$label] failed"
    if result !== nothing
        base = string(
            base,
            " (iterations=",
            result.iterations,
            ", residual=",
            fmt(result.residual_norm),
            ", converged=",
            string(result.converged),
            ")",
        )
    end
    if isempty(message)
        return base
    end
    return string(base, ": ", message)
end

"""key3(a, b, c; digits=6) -> NTuple{3,Float64}

Scan 用的三元 key（用于 completed set）。
"""
function key3(a, b, c; digits::Integer=6)
    return (
        round(Float64(a); digits=digits),
        round(Float64(b); digits=digits),
        round(Float64(c); digits=digits),
    )
end

"""key2(a, b; digits=6) -> Tuple{Float64,Float64}

Scan 用的二元 key（例如连续性 seed cache）。
"""
function key2(a, b; digits::Integer=6)
    return (
        round(Float64(a); digits=digits),
        round(Float64(b); digits=digits),
    )
end

"""load_completed_keys3(path; digits=6) -> Set{NTuple{3,Float64}}

从 scan CSV（默认 header + 前三列为 key）加载已完成点，供 resume 使用。
约定：前三列为 (T, x, xi)，其中 x 可为 mu 或 rho。
"""
function load_completed_keys3(path::AbstractString; digits::Integer=6)
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
                a = parse(Float64, strip(cols[1]))
                b = parse(Float64, strip(cols[2]))
                c = parse(Float64, strip(cols[3]))
                push!(completed, key3(a, b, c; digits=digits))
            catch
                # ignore malformed lines
            end
        end
    end
    return completed
end

"""attempt_with_candidates(candidates; solve_point, refine, promote, is_success) -> (result, message)

通用 scan glue：依次尝试候选 seed。
- `solve_point(seed_state)` -> `(result, msg)`
- `refine(result)` -> `(result, msg)`
- `promote(result)` -> `(result, msg)`
- `is_success(result)` -> `Bool`

要求 candidates 元素包含字段 `label` 与 `state`（如 NamedTuple）。
"""
function attempt_with_candidates(candidates;
    solve_point::Function,
    refine::Function,
    promote::Function,
    is_success::Function,
    stop_on_first_success::Bool=true,
    evaluate_all_attempts::Bool=false,
)
    messages = String[]
    governance_candidates = NamedTuple{(:pressure, :residual_norm, :hard_constraint_ok, :failed_constraints, :converged), Tuple{Float64, Float64, Bool, Vector{Symbol}, Bool}}[]
    point_results = Union{Nothing, SolverResult}[]
    candidate_labels = String[]

    scanned = Main.Models.execute_attempt_pool(candidates;
        stop_on_first_success=stop_on_first_success,
        evaluate_all_attempts=evaluate_all_attempts,
        evaluate_attempt=(candidate, _) -> begin
            result, msg = solve_point(candidate.state)
            refine_msg = ""
            promote_msg = ""

            if is_success(result)
                result, refine_msg = refine(result)
                result, promote_msg = promote(result)
            end

            success = is_success(result)
            if success
                !isempty(msg) && push!(messages, msg)
                !isempty(refine_msg) && push!(messages, refine_msg)
                !isempty(promote_msg) && push!(messages, promote_msg)
            else
                push!(messages, format_candidate_failure(candidate.label, msg, result))
            end

            scanned_candidate = (
                label=String(candidate.label),
                result=result,
                governance=(
                    pressure=(result === nothing ? -Inf : result.pressure),
                    residual_norm=(result === nothing ? Inf : result.residual_norm),
                    hard_constraint_ok=success,
                    failed_constraints=(success ? Symbol[] : Symbol[:scan_candidate_failed]),
                    converged=success,
                ),
            )
            return scanned_candidate, success
        end,
        on_error=(candidate, _) -> begin
            push!(messages, format_candidate_failure(candidate.label, "", nothing))
            scanned_candidate = (
                label=String(candidate.label),
                result=nothing,
                governance=(
                    pressure=-Inf,
                    residual_norm=Inf,
                    hard_constraint_ok=false,
                    failed_constraints=Symbol[:scan_candidate_failed],
                    converged=false,
                ),
            )
            return scanned_candidate, false
        end,
    )

    for scanned_candidate in scanned
        push!(candidate_labels, scanned_candidate.label)
        push!(point_results, scanned_candidate.result)
        push!(governance_candidates, scanned_candidate.governance)
    end

    isempty(governance_candidates) && return nothing, join_messages(messages)

    has_success = any(c -> c.converged, governance_candidates)
    selected = Main.Models.select_pressure_max_candidate(governance_candidates, nothing, nothing)
    selected_idx = Int(selected.selected_index)
    selection_note = "governance.selection=$(selected.selection_reason);seed=$(candidate_labels[selected_idx])"
    push!(messages, selection_note)

    selected_result = has_success ? point_results[selected_idx] : nothing
    return selected_result, join_messages(messages)
end

"""FixedSeedStrategy(seed)

Scan 内部 glue：固定返回指定 seed（复制），用于把某个向量当作一次求解的初值。
"""
struct FixedSeedStrategy <: SeedStrategy
    seed::Vector{Float64}
end

FixedSeedStrategy(seed::AbstractVector) = FixedSeedStrategy(Float64.(seed))

function get_seed(s::FixedSeedStrategy, θ::AbstractVector, mode::ConstraintMode)
    return copy(s.seed)
end

end # module ScanCommon
