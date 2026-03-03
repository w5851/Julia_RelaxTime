module ScanCommon

using Printf

import Main.Models: ConstraintMode
using ..SeedStrategies: SeedStrategy
import ..SeedStrategies: get_seed
using ..ImplicitSolver: SolverResult

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

fmt(x::Float64) = @sprintf("%.6f", x)
fmt(x::Real) = fmt(Float64(x))
fmt(x) = string(x)

function clean_message(msg::AbstractString)
    stripped = replace(strip(msg), '\n' => ' ')
    return replace(stripped, '"' => '\'')
end

quote_csv(msg::AbstractString) = isempty(msg) ? "" : string('"', msg, '"')
quote_csv(::Nothing) = ""

function join_messages(messages)
    filtered = filter(!isempty, messages)
    return isempty(filtered) ? "" : join(filtered, " | ")
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
)
    messages = String[]

    for candidate in candidates
        result, msg = solve_point(candidate.state)

        if is_success(result)
            refined, refine_msg = refine(result)
            result = refined

            result, promote_msg = promote(result)

            !isempty(msg) && push!(messages, msg)
            !isempty(promote_msg) && push!(messages, promote_msg)
            !isempty(refine_msg) && push!(messages, refine_msg)

            return result, join_messages(messages)
        end

        push!(messages, format_candidate_failure(candidate.label, msg, result))
    end

    return nothing, join_messages(messages)
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
