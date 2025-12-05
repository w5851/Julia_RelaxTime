module TrhoScan

using Printf
using LineSearches
using ..Constants_PNJL: ħc_MeV_fm
using ..SeedCache: find_initial_seed, DEFAULT_SEED_PATH
using ..AnisoGapSolver: solve_fixed_rho, DEFAULT_RHO_GUESS, SolverResult

export run_trho_scan, DEFAULT_T_VALUES, DEFAULT_RHO_VALUES, DEFAULT_OUTPUT_PATH

const DEFAULT_T_VALUES = collect(50.0:10.0:200.0)

function build_default_rho_grid(; rho_max::Float64=3.0,
        coarse_step::Float64=0.05,
        medium_switch::Float64=1.0,
        medium_step::Float64=0.02,
        fine_switch::Float64=0.3,
        fine_step::Float64=0.01,
        ultra_fine_switch::Float64=0.15,
        ultra_fine_step::Float64=0.005)
    rho_max > 0 || error("rho_max must be positive")
    values = Float64[]
    function _append_range!(container, step, stop_value)
        step > 0 || error("step must be positive")
        stop_value <= rho_max || error("switch value exceeds rho_max")
        append!(container, collect(0.0:step:stop_value))
    end
    append!(values, collect(0.0:coarse_step:rho_max))
    _append_range!(values, medium_step, medium_switch)
    _append_range!(values, fine_step, fine_switch)
    _append_range!(values, ultra_fine_step, ultra_fine_switch)
    unique_vals = unique(round.(values; digits=6))
    sort!(unique_vals)
    return unique_vals
end

const DEFAULT_RHO_VALUES = build_default_rho_grid()
const DEFAULT_OUTPUT_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "outputs", "results", "pnjl", "trho_scan.csv"))
const DEFAULT_SOLVER_KWARGS = (; linesearch = LineSearches.BackTracking())
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
    "iterations",
    "residual_norm",
    "converged",
    "message",
), ",")

function run_trho_scan(; T_values=DEFAULT_T_VALUES, rho_values=DEFAULT_RHO_VALUES, xi_values=[0.0],
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
        for xi in xi_values, T in T_values, rho in rho_values
            stats[:total] += 1
            key = _key(T, rho, xi)
            if key in completed
                stats[:skipped] += 1
                continue
            end

            seed_key = _seed_continuation_key(T, xi)
            candidates, seed_msg = _build_seed_candidates(continuation_seeds, seed_key, T, rho, xi, seed_path)
            result, message = _attempt_with_candidates(T, rho, xi, candidates; p_num=p_num, t_num=t_num, solver_opts..., initial_msg=seed_msg)
            if _is_success(result)
                continuation_seeds[seed_key] = copy(result.solution)
            end
            _write_row(io, T, rho, xi, result, message)
            flush(io)
            push!(completed, key)

            if _is_success(result)
                stats[:success] += 1
            else
                stats[:failure] += 1
            end

            if progress_cb !== nothing
                try
                    progress_cb((T=T, rho=rho, xi=xi), result)
                catch
                    # ignore callback errors
                end
            end
        end
    end

    return (; total=stats[:total], success=stats[:success], failure=stats[:failure], skipped=stats[:skipped], output=output_path)
end

# ---------------- internal helpers ----------------

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

function _extract_rho_seed(nearest)
    if nearest === nothing
        return copy(DEFAULT_RHO_GUESS)
    end
    blended = hasproperty(nearest, :blended_state) ? nearest.blended_state : nothing
    if blended !== nothing && hasproperty(blended, :rho_seed)
        return Float64.(blended.rho_seed)
    end
    state = hasproperty(nearest, :state) ? nearest.state : nothing
    if state !== nothing && hasproperty(state, :x)
        vec = state.x
        if length(vec) >= 8
            return Float64.(vec[1:8])
        end
    end
    return copy(DEFAULT_RHO_GUESS)
end

function _build_seed_candidates(cache::Dict, seed_key, T, rho, xi, seed_path)
    candidates = NamedTuple{(:label, :state), Tuple{String, Vector{Float64}}}[]
    if haskey(cache, seed_key)
        push!(candidates, (label = "continuation", state = copy(cache[seed_key])))
    end
    seed_msg = ""
    try
        nearest = find_initial_seed((T_mev = T, rho = rho, xi = xi); path = seed_path)
        push!(candidates, (label = "seed-cache", state = _extract_rho_seed(nearest)))
    catch err
        seed_msg = _clean_message("seed lookup failed: $(err)")
    end
    push!(candidates, (label = "default", state = copy(DEFAULT_RHO_GUESS)))
    return candidates, seed_msg
end

function _attempt_with_candidates(T, rho, xi, candidates; p_num, t_num, initial_msg = "", solver_opts...)
    messages = String[]
    if !isempty(initial_msg)
        push!(messages, initial_msg)
    end
    for candidate in candidates
        result, msg = _solve_point(T, rho, xi, candidate.state; p_num = p_num, t_num = t_num, solver_opts...)
        if _is_success(result)
            refined, refine_msg = _refine_result(T, rho, xi, result; p_num=p_num, t_num=t_num, solver_opts...)
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

function _join_messages(messages)
    filtered = filter(!isempty, messages)
    return isempty(filtered) ? "" : join(filtered, " | ")
end

function _format_candidate_failure(label, message, result)
    base = "seed[$label] failed"
    if result !== nothing
        base = string(
            base,
            " (iterations=",
            result.iterations,
            ", residual=",
            _fmt(result.residual_norm),
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

_is_success(::Nothing) = false
function _is_success(result)
    result === nothing && return false
    residual = result.residual_norm
    if result.converged
        return true
    end
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
            push!(notes, string(
                "accepted with residual ",
                _fmt(residual),
                " (<= ",
                _fmt(ACCEPTABLE_RESIDUAL),
                ")",
            ))
        end
    end
    return isempty(notes) ? "" : join(notes, "; ")
end

function _refine_result(T, rho, xi, result; p_num, t_num, solver_opts...)
    result === nothing && return nothing, ""
    if result.converged
        return result, ""
    end
    residual = result.residual_norm
    if !isfinite(residual) || residual > ACCEPTABLE_RESIDUAL
        return result, ""
    end
    refined, msg = _solve_point(T, rho, xi, result.solution; p_num=p_num, t_num=t_num, solver_opts...)
    if refined !== nothing && refined.converged
        return refined, "refined from near-converged seed"
    end
    return result, msg
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

function _solve_point(T_mev, rho_target, xi, seed_state; solver_kwargs...)
    try
        res = solve_fixed_rho(T_mev, rho_target; xi=xi, seed_state=seed_state, solver_kwargs...)
        return res, ""
    catch err
        msg = sprint() do io
            showerror(io, err)
        end
        return nothing, _clean_message(msg)
    end
end

function _write_row(io, T, rho, xi, result, message)
    if result === nothing
        values = (
            _fmt(T),
            _fmt(rho),
            _fmt(xi),
            "NaN", "NaN", "NaN", "NaN",
            "NaN", "NaN", "NaN",
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

    mu_vec = result.mu_vec .* ħc_MeV_fm
    phi = result.solution[1:3]
    Phi1 = result.solution[4]
    Phi2 = result.solution[5]
    mu_avg = sum(mu_vec) / 3
    values = (
        _fmt(T),
        _fmt(rho),
        _fmt(xi),
        _fmt(mu_vec[1]),
        _fmt(mu_vec[2]),
        _fmt(mu_vec[3]),
        _fmt(mu_avg),
        _fmt(result.pressure),
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
