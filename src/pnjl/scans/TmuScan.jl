module TmuScan

using Printf
using LineSearches
using ..Constants_PNJL: ħc_MeV_fm
using ..SeedCache: find_initial_seed, DEFAULT_SEED_PATH
using ..AnisoGapSolver: solve_fixed_mu, DEFAULT_MU_GUESS, SolverResult

export run_tmu_scan, DEFAULT_T_VALUES, DEFAULT_MU_VALUES, DEFAULT_OUTPUT_PATH

const DEFAULT_T_VALUES = collect(50.0:10.0:200.0)
const DEFAULT_MU_VALUES = collect(0.0:10.0:400.0)
const DEFAULT_OUTPUT_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "outputs", "results", "pnjl", "tmu_scan.csv"))
const DEFAULT_SOLVER_KWARGS = (; linesearch = LineSearches.BackTracking())
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
    "iterations",
    "residual_norm",
    "converged",
    "message",
), ",")

function run_tmu_scan(; T_values=DEFAULT_T_VALUES, mu_values=DEFAULT_MU_VALUES, xi_values=[0.0],
    output_path::AbstractString=DEFAULT_OUTPUT_PATH, seed_path::AbstractString=DEFAULT_SEED_PATH,
    overwrite::Bool=false, resume::Bool=true, p_num::Int=24, t_num::Int=8,
    progress_cb::Union{Nothing, Function}=nothing, solver_kwargs...)

    mkpath(dirname(output_path))
    completed = (resume && !overwrite && isfile(output_path)) ? _load_completed(output_path) : Set{NTuple{3, Float64}}()
    io_mode = (overwrite || !isfile(output_path)) ? "w" : "a"

    stats = Dict(:total => 0, :success => 0, :failure => 0, :skipped => 0)
    solver_opts = merge(DEFAULT_SOLVER_KWARGS, (; solver_kwargs...))
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

            seed_key = _seed_continuation_key(T, xi)
            seed_candidates, seed_msg = _build_seed_candidates(continuation_seeds, seed_key, T, mu, xi, seed_path)

            mu_fm = mu / ħc_MeV_fm
            result, message = _attempt_with_candidates(T, mu_fm, xi, seed_candidates;
                p_num=p_num, t_num=t_num, solver_opts..., initial_msg=seed_msg)
            if _is_success(result)
                continuation_seeds[seed_key] = copy(result.solution)
            end
            _write_row(io, T, mu, xi, result, message)
            flush(io)
            push!(completed, key)

            if _is_success(result)
                stats[:success] += 1
            else
                stats[:failure] += 1
            end

            if progress_cb !== nothing
                try
                    progress_cb((T=T, mu=mu, xi=xi), result)
                catch
                    # ignore callback errors
                end
            end
        end
    end

    return (; total=stats[:total], success=stats[:success], failure=stats[:failure], skipped=stats[:skipped], output=output_path)
end

# ---------------- internal helpers ----------------

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

function _extract_mu_seed(nearest)
    if nearest === nothing
        return copy(DEFAULT_MU_GUESS)
    end
    blended = hasproperty(nearest, :blended_state) ? nearest.blended_state : nothing
    if blended !== nothing && hasproperty(blended, :mu_seed)
        return Float64.(blended.mu_seed)
    end
    state = hasproperty(nearest, :state) ? nearest.state : nothing
    if state !== nothing && hasproperty(state, :x) && length(state.x) >= 5
        return Float64.(state.x[1:5])
    end
    return copy(DEFAULT_MU_GUESS)
end

function _solve_point(T_mev, mu_fm, xi, seed_state; solver_kwargs...)
    try
        res = solve_fixed_mu(T_mev, mu_fm; xi=xi, seed_state=seed_state, solver_kwargs...)
        return res, ""
    catch err
        msg = sprint() do io
            showerror(io, err)
        end
        return nothing, _clean_message(msg)
    end
end

function _build_seed_candidates(cache, key, T, mu, xi, seed_path)
    candidates = NamedTuple{(:label, :state), Tuple{String, Vector{Float64}}}[]
    if haskey(cache, key)
        push!(candidates, (label="continuation", state=copy(cache[key])))
    end
    seed_msg = ""
    try
        nearest = find_initial_seed((T_mev=T, mu_mev=mu, xi=xi); path=seed_path)
        push!(candidates, (label="seed-cache", state=_extract_mu_seed(nearest)))
    catch err
        seed_msg = _clean_message("seed lookup failed: $(err)")
    end
    push!(candidates, (label="default", state=copy(DEFAULT_MU_GUESS)))
    return candidates, seed_msg
end

function _attempt_with_candidates(T, mu_fm, xi, candidates; p_num, t_num, initial_msg="", solver_opts...)
    messages = String[]
    if !isempty(initial_msg)
        push!(messages, initial_msg)
    end
    for candidate in candidates
        result, msg = _solve_point(T, mu_fm, xi, candidate.state; p_num=p_num, t_num=t_num, solver_opts...)
        if _is_success(result)
            refined, refine_msg = _refine_result(T, mu_fm, xi, result; p_num=p_num, t_num=t_num, solver_opts...)
            result = refined
            result, promote_msg = _promote_success(result)
            if !isempty(msg)
                push!(messages, msg)
            end
            if !isempty(promote_msg)
                push!(messages, promote_msg)
            end
            approx_msg = _approx_success_note(result, refine_msg)
            if !isempty(approx_msg)
                push!(messages, approx_msg)
            end
            return result, _join_messages(messages)
        end
        push!(messages, _format_candidate_failure(candidate.label, msg, result))
    end
    return nothing, _join_messages(messages)
end

function _refine_result(T, mu_fm, xi, result; p_num, t_num, solver_opts...)
    result === nothing && return nothing, ""
    if result.converged
        return result, ""
    end
    residual = result.residual_norm
    if !isfinite(residual) || residual > ACCEPTABLE_RESIDUAL
        return result, ""
    end
    refined, msg = _solve_point(T, mu_fm, xi, result.solution; p_num=p_num, t_num=t_num, solver_opts...)
    if refined !== nothing && refined.converged
        return refined, "refined from near-converged seed"
    end
    return result, msg
end

function _join_messages(messages)
    filtered = filter(!isempty, messages)
    return isempty(filtered) ? "" : join(filtered, " | ")
end

function _format_candidate_failure(label, message, result)
    base = "seed[$label] failed"
    if result !== nothing
        base = string(base, " (iterations=", result.iterations, ", residual=", _fmt(result.residual_norm), ", converged=", string(result.converged), ")")
    end
    if isempty(message)
        return base
    end
    return string(base, ": ", message)
end

_is_success(::Nothing) = false
function _is_success(result)
    result === nothing && return false
    if result.converged
        return true
    end
    residual = result.residual_norm
    return isfinite(residual) && residual <= ACCEPTABLE_RESIDUAL
end

function _approx_success_note(result, refine_msg)
    result === nothing && return ""
    notes = String[]
    if !isempty(refine_msg)
        push!(notes, refine_msg)
    end
    if !result.converged
        residual = result.residual_norm
        if isfinite(residual) && residual <= ACCEPTABLE_RESIDUAL
            push!(notes, string("accepted with residual ", _fmt(residual), " (<= ", _fmt(ACCEPTABLE_RESIDUAL), ")"))
        end
    end
    return isempty(notes) ? "" : join(notes, "; ")
end

function _promote_success(result)
    result === nothing && return nothing, ""
    if result.converged
        return result, ""
    end
    residual = result.residual_norm
    if !isfinite(residual) || residual > ACCEPTABLE_RESIDUAL
        return result, ""
    end
    promoted = SolverResult(
        result.mode,
        true,
        copy(result.solution),
        result.mu_vec,
        result.pressure,
        result.rho,
        result.entropy,
        result.energy,
        result.iterations,
        residual,
        result.xi,
    )
    msg = string("force-marked converged (residual ", _fmt(residual), ")")
    return promoted, msg
end

function _write_row(io, T, mu, xi, result, message)
    if result === nothing
        values = (
            _fmt(T),
            _fmt(mu),
            _fmt(xi),
            "NaN", "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN",
            "NaN", "NaN",
            "-1",
            "NaN",
            "false",
            _quote(message),
        )
        println(io, join(values, ','))
        return
    end

    phi = result.solution[1:3]
    Phi1 = result.solution[4]
    Phi2 = result.solution[5]
    values = (
        _fmt(T),
        _fmt(mu),
        _fmt(xi),
        _fmt(result.pressure),
        _fmt(result.rho),
        _fmt(result.entropy),
        _fmt(result.energy),
        _fmt(phi[1]),
        _fmt(phi[2]),
        _fmt(phi[3]),
        _fmt(Phi1),
        _fmt(Phi2),
        string(result.iterations),
        _fmt(result.residual_norm),
        string(result.converged),
        _quote(message),
    )
    println(io, join(values, ','))
end

_fmt(x::Float64) = @sprintf("%.6f", x)
_fmt(x::Real) = _fmt(Float64(x))
_fmt(x) = string(x)

function _clean_message(msg::AbstractString)
    stripped = replace(strip(msg), '\n' => ' ')
    return replace(stripped, '"' => '\'')
end

_quote(msg::AbstractString) = isempty(msg) ? "" : string('"', msg, '"')
_quote(::Nothing) = ""

end # module
