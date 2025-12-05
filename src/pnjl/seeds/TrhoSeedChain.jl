module TrhoSeedChain

using Printf
using LineSearches
using Statistics: mean
using JLD2: jldsave, load
using ..Constants_PNJL: ħc_MeV_fm
using ..AnisoGapSolver: solve_fixed_rho, DEFAULT_RHO_GUESS, SolverResult

export SeedChain, DEFAULT_CHAIN_DIR, DEFAULT_CHAIN_RHO_GRID,
       build_seed_chains, list_chain_files, find_chain_seed, build_dense_rho_grid,
       load_chain

function build_dense_rho_grid(; rho_max::Float64=3.0, coarse_step::Float64=0.05,
        medium_switch::Float64=0.5, medium_step::Float64=0.02,
        fine_switch::Float64=0.3, fine_step::Float64=0.01)
    coarse = collect(range(rho_max, stop=medium_switch, step=-coarse_step))
    medium = collect(range(medium_switch - medium_step, stop=fine_switch, step=-medium_step))
    fine = collect(range(fine_switch - fine_step, stop=0.0, step=-fine_step))
    grid = vcat(coarse, medium, fine)
    unique_vals = unique(round.(grid; digits=4))
    sort!(unique_vals; rev=true)
    return unique_vals
end

const DEFAULT_CHAIN_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "raw", "pnjl", "seeds", "trho_continuation"))
const DEFAULT_CHAIN_RHO_GRID = build_dense_rho_grid()
const DEFAULT_SOLVER_KWARGS = (; linesearch = LineSearches.BackTracking())
const ACCEPTABLE_RESIDUAL = 1e-4
const T_KEY_DIGITS = 3
const XI_KEY_DIGITS = 4

struct SeedChain
    T::Float64
    xi::Float64
    rho_values::Vector{Float64}
    mu_values::Vector{Float64}
    states::Vector{Vector{Float64}}  # length-8 state vectors
end

function build_seed_chains(; T_values::AbstractVector{<:Real}, xi_values::AbstractVector{<:Real},
        rho_grid::AbstractVector{<:Real}=DEFAULT_CHAIN_RHO_GRID,
        output_dir::AbstractString=DEFAULT_CHAIN_DIR, solver_kwargs...)
    mkpath(output_dir)
    solver_opts = merge(DEFAULT_SOLVER_KWARGS, (; solver_kwargs...))
    chains = SeedChain[]
    for xi in xi_values, T in T_values
        chain = _scan_chain(Float64(T), Float64(xi), rho_grid; solver_opts...)
        isempty(chain.rho_values) && continue
        path = joinpath(output_dir, _chain_filename(chain.T, chain.xi))
        _save_chain(path, chain)
        push!(chains, chain)
    end
    return chains
end

function _scan_chain(T::Float64, xi::Float64, rho_grid; solver_opts...)
    rho_values = Float64[]
    mu_values = Float64[]
    states = Vector{Vector{Float64}}()
    last_state = copy(DEFAULT_RHO_GUESS)
    for rho in sort(Float64.(rho_grid); rev=true)
        result, _ = _solve_point(T, rho, xi, last_state; solver_opts...)
        if !_is_success(result)
            continue
        end
        result, _ = _promote_success(result)
        push!(rho_values, rho)
        push!(mu_values, _mean_mu(result))
        push!(states, copy(result.solution))
        last_state = copy(result.solution)
    end
    order = sortperm(rho_values)
    rho_sorted = rho_values[order]
    mu_sorted = mu_values[order]
    states_sorted = states[order]
    return SeedChain(T, xi, rho_sorted, mu_sorted, states_sorted)
end

function _solve_point(T, rho, xi, seed_state; solver_opts...)
    try
        res = solve_fixed_rho(T, rho; xi=xi, seed_state=seed_state, solver_opts...)
        return res, ""
    catch err
        msg = sprint() do io
            showerror(io, err)
        end
        return nothing, msg
    end
end

function _is_success(result)
    result === nothing && return false
    if result.converged
        return true
    end
    residual = result.residual_norm
    return isfinite(residual) && residual <= ACCEPTABLE_RESIDUAL
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
    note = string("force-marked converged (residual ", _fmt(residual), ")")
    return promoted, note
end

_mean_mu(res) = mean(res.mu_vec) * ħc_MeV_fm

function _save_chain(path, chain::SeedChain)
    jldsave(path; T=chain.T, xi=chain.xi, rho=chain.rho_values, mu=chain.mu_values, states=chain.states)
end

const _CHAIN_CACHE = Dict{String, SeedChain}()

function list_chain_files(dir::AbstractString=DEFAULT_CHAIN_DIR)
    if !isdir(dir)
        return String[]
    end
    paths = String[]
    for file in readdir(dir; join=true)
        endswith(file, ".jld2") || continue
        push!(paths, file)
    end
    sort!(paths)
    return paths
end

function _load_chain(path::AbstractString)
    if haskey(_CHAIN_CACHE, path)
        return _CHAIN_CACHE[path]
    end
    data = load(path)
    required = ("T", "xi", "rho", "mu", "states")
    for key in required
        haskey(data, key) || error("Chain file $(path) missing field $(key)")
    end
    chain = SeedChain(Float64(data["T"]), Float64(data["xi"]),
        Float64.(data["rho"]), Float64.(data["mu"]), _convert_states(data["states"]))
    _CHAIN_CACHE[path] = chain
    return chain
end

load_chain(path::AbstractString) = _load_chain(path)

function _convert_states(states)
    return [Float64.(vec) for vec in states]
end

function find_chain_seed(T::Real, xi::Real, mu_target::Real;
        chain_dir::AbstractString=DEFAULT_CHAIN_DIR, T_tol::Real=5.0, xi_tol::Real=1e-3)
    chain = _nearest_chain(Float64(T), Float64(xi); dir=chain_dir, T_tol=T_tol, xi_tol=xi_tol)
    chain === nothing && return nothing
    isempty(chain.mu_values) && return nothing
    idx = _nearest_index(chain.mu_values, Float64(mu_target))
    return copy(chain.states[idx][1:5])
end

function _nearest_chain(T, xi; dir, T_tol, xi_tol)
    best = nothing
    best_dist = Inf
    for path in list_chain_files(dir)
        chain = _load_chain(path)
        ΔT = abs(chain.T - T)
        Δxi = abs(chain.xi - xi)
        if ΔT <= T_tol && Δxi <= xi_tol
            dist = hypot(ΔT, Δxi)
            if dist < best_dist
                best = chain
                best_dist = dist
            end
        end
    end
    return best
end

function _nearest_index(values::Vector{Float64}, target::Float64)
    best = 1
    best_diff = abs(values[1] - target)
    for i in 2:length(values)
        diff = abs(values[i] - target)
        if diff < best_diff
            best = i
            best_diff = diff
        end
    end
    return best
end

function _chain_filename(T, xi)
    return @sprintf("chain_T%05.1f_xi%+.3f.jld2", T, xi)
end

_fmt(x) = @sprintf("%.6f", x)

end # module
