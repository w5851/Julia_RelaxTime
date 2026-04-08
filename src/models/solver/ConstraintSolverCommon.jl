using StaticArrays
using NLsolve
using ForwardDiff
using Statistics: mean
using Main.Constants_PNJL: ρ0_inv_fm3

const rho0 = ρ0_inv_fm3

const HardConstraintRule = Function

@inline function _pack_solution(x_state::AbstractVector, mu_vec::AbstractVector)
    return vcat(Float64.(x_state), Float64.(mu_vec))
end

@inline _to_state_svec(st) = SVector{5}(Tuple(state_vector(st)))
@inline _to_mu_svec(mu_vec) = SVector{3}(Tuple(mu_vec))
@inline _to_chiral_triplet(x_state) = SVector{3}(x_state[1], x_state[2], x_state[3])
@inline _mass_from_state(model::AbstractQCDModel, x_state) = calculate_mass_vec(model, _to_chiral_triplet(x_state))

@inline function _model_kind_for_shared_core(model::AbstractQCDModel)
    if model isa RPNJLModel
        return :RPNJL
    elseif model isa NJL2Model
        return :NJL2
    elseif model isa AbstractNJLModel
        return :NJL
    end
    return :PNJL
end

@inline function _gap_norm_from_state(
    model::AbstractQCDModel,
    x_state::AbstractVector,
    mu_vec::AbstractVector,
    T_fm::Real;
    xi::Real,
    p_num::Int,
    t_num::Int,
)
    params = GapParams(Float64(T_fm), cached_nodes(p_num, t_num), Float64(xi);
        p_num=p_num,
        t_num=t_num,
        model_kind=_model_kind_for_shared_core(model),
    )
    Tout = promote_type(eltype(x_state), eltype(mu_vec), Float64)
    out = Vector{Tout}(undef, length(x_state))
    gap_core_residual!(out, model, x_state, mu_vec, params)
    return sqrt(sum(abs2, out))
end

@inline function _unpack_solution(solution::AbstractVector; state_n::Int=5, mu_n::Int=3)
    return _unpack_solution(solution, Val(state_n), Val(mu_n))
end

@inline function _unpack_solution(solution::AbstractVector, ::Val{STATE_N}, ::Val{MU_N}) where {STATE_N, MU_N}
    expected = STATE_N + MU_N
    length(solution) == expected || throw(ArgumentError("solution length mismatch: expected $expected, got $(length(solution))"))
    x_state = SVector{STATE_N}(Tuple(solution[1:STATE_N]))
    mu_vec = SVector{MU_N}(Tuple(solution[(STATE_N + 1):expected]))
    return x_state, mu_vec
end

@inline function _empty_candidate(; state_n::Int=5, mu_n::Int=3, solution_n::Int=(state_n + mu_n), residual_norm_max::Real)
    return (
        solution=fill(NaN, solution_n),
        x_state=SVector{state_n}(fill(NaN, state_n)),
        mu_vec=SVector{mu_n}(fill(NaN, mu_n)),
        omega=NaN,
        pressure=-Inf,
        rho_norm=NaN,
        entropy=NaN,
        energy=NaN,
        masses=SVector{3}(fill(NaN, 3)),
        iterations=0,
        residual_norm=Inf,
        residual_norm_max=Float64(residual_norm_max),
        hard_constraint_ok=false,
        failed_constraints=Symbol[:solver_failed],
        converged=false,
    )
end

@inline function _build_mode_failure_candidate(; state_n::Int=5, mu_n::Int=3, solution_n::Int=(state_n + mu_n), residual_norm_max::Real, seed_index::Union{Nothing, Int}=nothing)
    raw = _empty_candidate(state_n=state_n, mu_n=mu_n, solution_n=solution_n, residual_norm_max=residual_norm_max)
    return seed_index === nothing ? raw : (; raw..., seed_index=Int(seed_index))
end

@inline function default_mu0_from_seed(seed)::Float64
    if length(seed) >= 8
        return mean(Float64.(seed[6:8]))
    elseif length(seed) >= 6
        return Float64(seed[6])
    end
    return 0.2
end

@inline function default_muvec0_from_seed(seed)::Vector{Float64}
    if length(seed) >= 8
        return Float64.(seed[6:8])
    end
    μ0 = default_mu0_from_seed(seed)
    return [μ0, μ0, μ0]
end

@inline function _constraint_failure!(F)
    @inbounds for i in eachindex(F)
        F[i] = convert(eltype(F), 1e6)
    end
    return nothing
end

@inline function _refresh_constraint_state!(residual_fn!::Function, root)
    tmpF = zeros(length(root))
    residual_fn!(tmpF, root)
    return nothing
end

function _run_outer_nlsolve(
    residual_fn!::Function,
    x0::AbstractVector;
    nlsolve_method::Symbol,
    nlsolve_kwargs...,
)
    res = nlsolve(
        residual_fn!,
        Float64.(x0);
        autodiff=:forward,
        method=nlsolve_method,
        xtol=1e-9,
        ftol=1e-9,
        nlsolve_kwargs...,
    )
    _refresh_constraint_state!(residual_fn!, res.zero)
    return res
end

@inline function _mode_outer_state_ready(
    st_ref,
    x_state_ref,
    mu_vec_ref,
    pressure_ref,
    rho_norm_ref,
    entropy_ref,
    energy_ref,
    masses_ref,
)
    return !(st_ref[] === nothing || x_state_ref[] === nothing || mu_vec_ref[] === nothing ||
             pressure_ref[] === nothing || rho_norm_ref[] === nothing || entropy_ref[] === nothing ||
             energy_ref[] === nothing || masses_ref[] === nothing)
end

@inline function _build_mode_result_from_outer_state(
    x_state,
    mu_vec,
    pressure,
    rho_norm,
    entropy,
    energy,
    masses,
    residual_norm;
    iterations::Integer,
    converged::Bool,
    legacy_fallback_used::Bool=false,
)
    return (
        converged=converged,
        solution=_pack_solution(x_state, mu_vec),
        x_state=x_state,
        mu_vec=mu_vec,
        omega=-pressure,
        pressure=pressure,
        rho_norm=rho_norm,
        entropy=entropy,
        energy=energy,
        masses=masses,
        iterations=Int(iterations),
        residual_norm=residual_norm,
        legacy_fallback_used=legacy_fallback_used,
    )
end

@inline function _is_mass_positive(masses)::Bool
    return all(m -> isfinite(m) && m > 0.0, masses)
end

@inline function _is_phi_in_range(x_state)::Bool
    Φ = x_state[4]
    Φbar = x_state[5]
    return isfinite(Φ) && isfinite(Φbar) && (0.0 <= Φ <= 1.0) && (0.0 <= Φbar <= 1.0)
end

function default_hard_constraint_rules(; physicality_check::Function=((_, _) -> true))
    return HardConstraintRule[
        c -> (isfinite(c.residual_norm), :residual_nonfinite),
        c -> (c.residual_norm <= c.residual_norm_max, :residual_too_large),
        c -> (_is_mass_positive(c.masses), :mass_nonpositive),
        c -> (_is_phi_in_range(c.x_state), :phi_out_of_range),
        c -> (isfinite(c.pressure), :pressure_nonfinite),
        c -> (isfinite(c.omega), :omega_nonfinite),
        c -> (isfinite(c.rho_norm), :rho_nonfinite),
        c -> (isfinite(c.entropy), :entropy_nonfinite),
        c -> (isfinite(c.energy), :energy_nonfinite),
        c -> (physicality_check(c.x_state, c.masses), :physicality_check_failed),
    ]
end

function evaluate_hard_constraints(candidate, rules::AbstractVector{<:HardConstraintRule}, params=nothing, context=nothing)
    prepared_rules = if params === nothing && context === nothing
        rules
    else
        map(rules) do rule
            if applicable(rule, candidate, params, context)
                (c -> rule(c, params, context))
            else
                rule
            end
        end
    end

    failed = Symbol[]
    for rule in prepared_rules
        ok, reason = rule(candidate)
        if !ok
            push!(failed, Symbol(reason))
        end
    end
    return isempty(failed), failed
end

function select_pressure_max_candidate(candidates::AbstractVector, params=nothing, context=nothing)
    isempty(candidates) && throw(ArgumentError("candidates must be non-empty"))

    candidate_seed_index(cand, fallback::Int) = if hasproperty(cand, :seed_index)
        Int(getproperty(cand, :seed_index))
    else
        fallback
    end

    candidate_residual(cand) = begin
        value = Float64(cand.residual_norm)
        isfinite(value) ? value : Inf
    end

    candidate_pressure(cand) = begin
        value = Float64(cand.pressure)
        isfinite(value) ? value : -Inf
    end

    function better_candidate(cand, cand_idx::Int, best, best_idx::Int)
        cand_ok = Bool(cand.hard_constraint_ok)
        best_ok = Bool(best.hard_constraint_ok)
        if cand_ok != best_ok
            return cand_ok
        end

        cand_residual = candidate_residual(cand)
        best_residual = candidate_residual(best)
        if cand_residual != best_residual
            return cand_residual < best_residual
        end

        cand_pressure = candidate_pressure(cand)
        best_pressure = candidate_pressure(best)
        if cand_pressure != best_pressure
            return cand_pressure > best_pressure
        end

        cand_seed_index = candidate_seed_index(cand, cand_idx)
        best_seed_index = candidate_seed_index(best, best_idx)
        return cand_seed_index < best_seed_index
    end

    passed = Int[]
    for i in eachindex(candidates)
        if Bool(candidates[i].hard_constraint_ok)
            push!(passed, i)
        end
    end

    if !isempty(passed)
        selected_idx = passed[1]
        selected = candidates[selected_idx]
        for i in passed[2:end]
            cand = candidates[i]
            if better_candidate(cand, i, selected, selected_idx)
                selected_idx = i
                selected = cand
            end
        end
        return (
            selected_index=selected_idx,
            selected_candidate=selected,
            selection_reason=:pressure_max_under_constraints,
            eligible_indices=passed,
        )
    end

    selected_idx = 1
    selected = candidates[selected_idx]
    for i in eachindex(candidates)
        cand = candidates[i]
        if better_candidate(cand, i, selected, selected_idx)
            selected_idx = i
            selected = cand
        end
    end
    return (
        selected_index=selected_idx,
        selected_candidate=selected,
        selection_reason=:no_candidate_passed_constraints,
        eligible_indices=Int[],
    )
end

function select_residual_min_candidate(candidates::AbstractVector, params=nothing, context=nothing)
    isempty(candidates) && throw(ArgumentError("candidates must be non-empty"))

    candidate_seed_index(cand, fallback::Int) = if hasproperty(cand, :seed_index)
        Int(getproperty(cand, :seed_index))
    else
        fallback
    end

    candidate_residual(cand) = begin
        value = Float64(cand.residual_norm)
        isfinite(value) ? value : Inf
    end

    candidate_pressure(cand) = begin
        value = Float64(cand.pressure)
        isfinite(value) ? value : -Inf
    end

    function better_candidate(cand, cand_idx::Int, best, best_idx::Int)
        cand_ok = Bool(cand.hard_constraint_ok)
        best_ok = Bool(best.hard_constraint_ok)
        if cand_ok != best_ok
            return cand_ok
        end

        cand_residual = candidate_residual(cand)
        best_residual = candidate_residual(best)
        if cand_residual != best_residual
            return cand_residual < best_residual
        end

        cand_seed_index = candidate_seed_index(cand, cand_idx)
        best_seed_index = candidate_seed_index(best, best_idx)
        if cand_seed_index != best_seed_index
            return cand_seed_index < best_seed_index
        end

        cand_pressure = candidate_pressure(cand)
        best_pressure = candidate_pressure(best)
        return cand_pressure > best_pressure
    end

    passed = Int[]
    for i in eachindex(candidates)
        if Bool(candidates[i].hard_constraint_ok)
            push!(passed, i)
        end
    end

    if !isempty(passed)
        selected_idx = passed[1]
        selected = candidates[selected_idx]
        for i in passed[2:end]
            cand = candidates[i]
            if better_candidate(cand, i, selected, selected_idx)
                selected_idx = i
                selected = cand
            end
        end
        return (
            selected_index=selected_idx,
            selected_candidate=selected,
            selection_reason=:residual_min_under_constraints,
            eligible_indices=passed,
        )
    end

    selected_idx = 1
    selected = candidates[selected_idx]
    for i in eachindex(candidates)
        cand = candidates[i]
        if better_candidate(cand, i, selected, selected_idx)
            selected_idx = i
            selected = cand
        end
    end
    return (
        selected_index=selected_idx,
        selected_candidate=selected,
        selection_reason=:no_candidate_passed_constraints,
        eligible_indices=Int[],
    )
end

function _compute_fixedmu_candidate(
    model::AbstractQCDModel,
    T_fm::Real,
    μ_fm::Real,
    st,
    residual_norm_max::Real;
    xi::Real,
    p_num::Int,
    t_num::Int,
)
    x_state = _to_state_svec(st)
    mu_vec = normalize_mu_vec(μ_fm)

    thermo = _compute_mode_thermo_quantities(
        model,
        x_state,
        T_fm,
        mu_vec;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        rho0_scale=nothing,
    )
    residual_norm = _compose_mode_residual_norm(
        model,
        x_state,
        mu_vec,
        T_fm;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )

    return (
        solution=Vector{Float64}(x_state),
        x_state=x_state,
        mu_vec=SVector{3}(Tuple(mu_vec)),
        omega=thermo.omega,
        pressure=thermo.pressure,
        rho_norm=thermo.rho_norm,
        entropy=thermo.entropy,
        energy=thermo.energy,
        masses=thermo.masses,
        iterations=0,
        residual_norm=residual_norm,
        residual_norm_max=Float64(residual_norm_max),
    )
end

function _compute_mode_thermo_quantities(
    model::AbstractQCDModel,
    x_state::AbstractVector,
    T_fm::Real,
    mu_vec::AbstractVector;
    xi::Real,
    p_num::Int,
    t_num::Int,
    rho0_scale::Union{Nothing, Real}=rho0,
)
    pressure = -omega(model, x_state, T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
    pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
    rho_vec = ForwardDiff.gradient(pressure_mu, mu_vec)
    rho_norm = if rho0_scale === nothing
        sum(rho_vec) / 3.0
    else
        sum(rho_vec) / (3.0 * Float64(rho0_scale))
    end
    pressure_T = τ -> -omega(model, x_state, τ, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
    entropy = ForwardDiff.derivative(pressure_T, T_fm)
    energy = -pressure + sum(mu_vec .* rho_vec) + T_fm * entropy
    masses = _mass_from_state(model, x_state)

    return (
        pressure=pressure,
        omega=-pressure,
        rho_vec=rho_vec,
        rho_norm=rho_norm,
        entropy=entropy,
        energy=energy,
        masses=masses,
    )
end

@inline function _residual_component_value(v::Real)
    value = abs(Float64(v))
    return isfinite(value) ? value : Inf
end

@inline function _residual_component_value(v::Tuple{<:Real, <:Real})
    value = abs(Float64(v[1]) - Float64(v[2]))
    return isfinite(value) ? value : Inf
end

function _compose_mode_residual_norm(
    model::AbstractQCDModel,
    x_state::AbstractVector,
    mu_vec::AbstractVector,
    T_fm::Real,
    components...;
    xi::Real,
    p_num::Int,
    t_num::Int,
)
    gap_norm = _gap_norm_from_state(model, x_state, mu_vec, T_fm; xi=xi, p_num=p_num, t_num=t_num)
    residual_norm = gap_norm
    for component in components
        comp = _residual_component_value(component)
        residual_norm = max(residual_norm, comp)
    end
    return residual_norm
end

function _build_default_seed_candidates(seed_guess::AbstractVector)
    base = Float64.(seed_guess)
    seeds = Vector{Vector{Float64}}()
    push!(seeds, copy(base))

    length(base) >= 5 || return seeds

    perturb_specs = (
        (1, 1.01),
        (1, 0.99),
        (2, 1.01),
        (2, 0.99),
        (3, 1.01),
        (3, 0.99),
        (4, 1.05),
        (5, 1.05),
    )

    for (idx, scale) in perturb_specs
        s = copy(base)
        s[idx] *= scale
        push!(seeds, s)
    end

    uniq = Vector{Vector{Float64}}()
    seen = Set{String}()
    for s in seeds
        key = join(round.(s; digits=10), ",")
        key in seen && continue
        push!(uniq, s)
        push!(seen, key)
    end
    return uniq
end

@inline function _mode_seed_key(seed::AbstractVector{<:Real})
    return join(round.(Float64.(seed); digits=10), ",")
end

@inline function _push_unique_seed!(pool::Vector{Vector{Float64}}, seen::Set{String}, seed::AbstractVector{<:Real})
    normalized = Float64.(seed)
    key = _mode_seed_key(normalized)
    key in seen && return nothing
    push!(pool, normalized)
    push!(seen, key)
    return nothing
end

function _build_default_seed_candidates(mode, seed_guess::AbstractVector)
    seeds = _build_default_seed_candidates(seed_guess)
    seen = Set{String}(_mode_seed_key(s) for s in seeds)

    base5 = if length(seed_guess) >= 5
        Float64.(seed_guess[1:5])
    else
        copy(HADRON_SEED_5)
    end

    if mode isa FixedEntropy || mode isa FixedSigma
        mu_anchors = (0.1, 1.0, 1.7)
        base_candidates = (base5, HADRON_SEED_5, QUARK_SEED_5)
        for b in base_candidates
            for μ0 in mu_anchors
                _push_unique_seed!(seeds, seen, Float64[b[1:5]..., μ0, μ0, μ0])
            end
        end
    elseif mode isa FixedAsymmetricRho
        for μ0 in (0.5, 1.0)
            _push_unique_seed!(seeds, seen, Float64[base5..., μ0, μ0, μ0])
        end
    elseif mode isa FixedRho
        for μ0 in (0.2, 1.0)
            _push_unique_seed!(seeds, seen, Float64[base5..., μ0, μ0, μ0])
        end
    end

    return seeds
end

function _solve_gap_with_outer_fallback(
    model::AbstractQCDModel,
    T_fm,
    μ_input;
    st_prev,
    seed_guess,
    solver_primary::AbstractGapSolver,
    solver_secondary::Union{Nothing, AbstractGapSolver},
    residual_norm_max::Real,
    xi,
    p_num::Int,
    t_num::Int,
)
    initial_guess = if st_prev === nothing
        Float64.(seed_guess)
    else
        Vector{Float64}(state_vector(st_prev))
    end

    try
        return solve_gap(
            model,
            T_fm,
            μ_input;
            solver_backend=:models,
            solver=solver_primary,
            initial_guess=initial_guess,
            residual_norm_max=residual_norm_max,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )
    catch err
        err isa InterruptException && rethrow()
        solver_secondary === nothing && return nothing
        try
            return solve_gap(
                model,
                T_fm,
                μ_input;
                solver_backend=:models,
                solver=solver_secondary,
                initial_guess=initial_guess,
                residual_norm_max=residual_norm_max,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
        catch err2
            err2 isa InterruptException && rethrow()
            return nothing
        end
    end
end
