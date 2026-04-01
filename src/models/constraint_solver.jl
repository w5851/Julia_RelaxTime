"""constraint_solver.jl

Models 侧约束求解（阶段 B 内核下沉）：
- 先以 FixedEntropy 为起点，将单变量 μ 外层约束求解从 PNJL 兼容层下沉到 models。
"""

using StaticArrays
using NLsolve
using ForwardDiff
using Statistics: mean

const HardConstraintRule = Function

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
            if cand.pressure > selected.pressure || (cand.pressure == selected.pressure && cand.residual_norm < selected.residual_norm)
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
        if cand.pressure > selected.pressure || (cand.pressure == selected.pressure && cand.residual_norm < selected.residual_norm)
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
    x_state = SVector{5}(Tuple(state_vector(st)))
    mu_vec = normalize_mu_vec(μ_fm)

    pressure = -omega(model, x_state, T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
    pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
    rho_vec = ForwardDiff.gradient(pressure_mu, mu_vec)
    rho_norm = sum(rho_vec) / 3.0
    pressure_T = τ -> -omega(model, x_state, τ, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
    entropy = ForwardDiff.derivative(pressure_T, T_fm)
    energy = -pressure + sum(mu_vec .* rho_vec) + T_fm * entropy
    omega_val = -pressure
    masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

    residual_vec = gap_residual(model, st, T_fm, mu_vec; xi=xi, p_num=p_num, t_num=t_num)
    residual_norm = sqrt(sum(abs2, residual_vec))

    return (
        solution=Vector{Float64}(x_state),
        x_state=x_state,
        mu_vec=SVector{3}(Tuple(mu_vec)),
        omega=omega_val,
        pressure=pressure,
        rho_norm=rho_norm,
        entropy=entropy,
        energy=energy,
        masses=masses,
        iterations=0,
        residual_norm=residual_norm,
        residual_norm_max=Float64(residual_norm_max),
    )
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

    uniq = Dict{String,Vector{Float64}}()
    for s in seeds
        key = join(round.(s; digits=10), ",")
        uniq[key] = s
    end
    return collect(values(uniq))
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
    catch
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
        catch
            return nothing
        end
    end
end

function solve_fixedmu_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    μ_fm::Real;
    seed_guess::AbstractVector,
    solver::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    physicality_check::Function=((_, _) -> true),
    seed_candidates::Union{Nothing, AbstractVector}=nothing,
    hard_constraints::Union{Nothing, AbstractVector{<:HardConstraintRule}}=nothing,
    __allow_hard_deprecated_internal::Bool=false,
)
    __allow_hard_deprecated_internal || throw(ArgumentError("Models.solve_fixedmu_constraint is hard-deprecated in Wave-D; use Models.solve_constraint(model, FixedMu(), T; μ_fm=...)."))
    seed_pool = if seed_candidates === nothing
        _build_default_seed_candidates(seed_guess)
    else
        [Float64.(s) for s in seed_candidates]
    end

    rules = hard_constraints === nothing ? default_hard_constraint_rules(; physicality_check=physicality_check) : hard_constraints

    candidates = NamedTuple[]
    for seed in seed_pool
        local raw
        try
            st = solve_gap(
                model,
                T_fm,
                μ_fm;
                solver_backend=:models,
                solver=solver,
                initial_guess=seed,
                residual_norm_max=residual_norm_max,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )

            raw = _compute_fixedmu_candidate(
                model,
                T_fm,
                μ_fm,
                st,
                residual_norm_max;
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
            ok, failed = evaluate_hard_constraints(raw, rules)
            push!(candidates, (; raw..., hard_constraint_ok=ok, failed_constraints=failed, converged=ok))
        catch
            raw = (
                solution=Float64[],
                x_state=SVector{5}(fill(NaN, 5)),
                mu_vec=SVector{3}(fill(NaN, 3)),
                omega=NaN,
                pressure=-Inf,
                rho_norm=NaN,
                entropy=NaN,
                energy=NaN,
                masses=SVector{3}(fill(NaN, 3)),
                iterations=0,
                residual_norm=Inf,
                residual_norm_max=Float64(residual_norm_max),
            )
            push!(candidates, (; raw..., hard_constraint_ok=false, failed_constraints=Symbol[:solver_failed], converged=false))
        end
    end
    selected = select_pressure_max_candidate(candidates)
    s = selected.selected_candidate

    return (
        converged=s.converged,
        solution=s.solution,
        x_state=s.x_state,
        mu_vec=s.mu_vec,
        omega=s.omega,
        pressure=s.pressure,
        rho_norm=s.rho_norm,
        entropy=s.entropy,
        energy=s.energy,
        masses=s.masses,
        iterations=s.iterations,
        residual_norm=s.residual_norm,
        hard_constraint_ok=s.hard_constraint_ok,
        failed_constraints=s.failed_constraints,
        selection_reason=selected.selection_reason,
        selected_index=selected.selected_index,
        candidate_count=length(candidates),
    )
end

function solve_fixedrho_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    rho_target::Real;
    seed_guess::AbstractVector,
    mu0::Real=0.2,
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    nlsolve_method::Symbol=:newton,
    physicality_check::Function=((_, _) -> true),
    __allow_hard_deprecated_internal::Bool=false,
    nlsolve_kwargs...,
)
    __allow_hard_deprecated_internal || throw(ArgumentError("Models.solve_fixedrho_constraint is hard-deprecated in Wave-D; use Models.solve_constraint(model, FixedRho(...), T)."))
    st_ref = Ref{Any}(nothing)
    x_state_ref = Ref{Any}(nothing)
    mu_vec_ref = Ref{Any}(nothing)
    pressure_ref = Ref{Any}(nothing)
    rho_norm_ref = Ref{Any}(nothing)
    entropy_ref = Ref{Any}(nothing)
    energy_ref = Ref{Any}(nothing)
    masses_ref = Ref{Any}(nothing)

    residual_fn! = (F, x) -> begin
        μ = x[1]
        μ_vec = normalize_mu_vec(μ)

        st = _solve_gap_with_outer_fallback(
            model,
            T_fm,
            μ;
            st_prev=st_ref[],
            seed_guess=seed_guess,
            solver_primary=solver_primary,
            solver_secondary=solver_secondary,
            residual_norm_max=residual_norm_max,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )

        if st === nothing
            _constraint_failure!(F)
            return nothing
        end

        x_state = SVector{5}(Tuple(state_vector(st)))
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)

        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / 3.0

        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)

        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = SVector{3}(Tuple(μ_vec))
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses

        F[1] = convert(eltype(F), rho_norm - rho_target)
        return nothing
    end

    res = try
        nlsolve(
            residual_fn!,
            [Float64(mu0)];
            autodiff=:forward,
            method=nlsolve_method,
            xtol=1e-9,
            ftol=1e-9,
            nlsolve_kwargs...,
        )
    catch
        if nlsolve_method === :trust_region
            rethrow()
        end
        nlsolve(
            residual_fn!,
            [Float64(mu0)];
            autodiff=:forward,
            method=:trust_region,
            xtol=1e-9,
            ftol=1e-9,
            nlsolve_kwargs...,
        )
    end

    _refresh_constraint_state!(residual_fn!, res.zero)

    gap_vec = gap_residual(model, st_ref[], T_fm, mu_vec_ref[]; xi=xi, p_num=p_num, t_num=t_num)
    gap_norm = sqrt(sum(abs2, gap_vec))
    rho_residual = abs(rho_norm_ref[] - rho_target)
    residual_norm = max(gap_norm, rho_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
    converged = res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-4)

    return (
        converged=converged,
        solution=Float64[
            x_state_ref[][1], x_state_ref[][2], x_state_ref[][3], x_state_ref[][4], x_state_ref[][5],
            mu_vec_ref[][1], mu_vec_ref[][2], mu_vec_ref[][3],
        ],
        x_state=x_state_ref[],
        mu_vec=mu_vec_ref[],
        omega=omega_val,
        pressure=pressure_ref[],
        rho_norm=rho_norm_ref[],
        entropy=entropy_ref[],
        energy=energy_ref[],
        masses=masses_ref[],
        iterations=res.iterations,
        residual_norm=residual_norm,
    )
end

function solve_fixedentropy_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    s_target::Real;
    seed_guess::AbstractVector,
    mu0::Real=0.2,
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    rho0::Real,
    physicality_check::Function=((_, _) -> true),
    mass_positive_constraint::Bool=true,
    __allow_hard_deprecated_internal::Bool=false,
    nlsolve_kwargs...,
)
    __allow_hard_deprecated_internal || throw(ArgumentError("Models.solve_fixedentropy_constraint is hard-deprecated in Wave-D; use Models.solve_constraint(model, FixedEntropy(...), T)."))
    st_ref = Ref{Any}(nothing)
    x_state_ref = Ref{Any}(nothing)
    mu_vec_ref = Ref{Any}(nothing)
    pressure_ref = Ref{Any}(nothing)
    rho_norm_ref = Ref{Any}(nothing)
    entropy_ref = Ref{Any}(nothing)
    energy_ref = Ref{Any}(nothing)
    masses_ref = Ref{Any}(nothing)

    residual_fn! = (F, x) -> begin
        μ = x[1]
        μ_vec = normalize_mu_vec(μ)

        st = _solve_gap_with_outer_fallback(
            model,
            T_fm,
            μ;
            st_prev=st_ref[],
            seed_guess=seed_guess,
            solver_primary=solver_primary,
            solver_secondary=solver_secondary,
            residual_norm_max=residual_norm_max,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )

        if st === nothing
            _constraint_failure!(F)
            return nothing
        end

        x_state = SVector{5}(Tuple(state_vector(st)))
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)

        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / (3.0 * rho0)

        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)

        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

        if mass_positive_constraint && any(m -> !isfinite(m) || m <= 0.0, masses)
            _constraint_failure!(F)
            return nothing
        end

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = SVector{3}(Tuple(μ_vec))
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses

        F[1] = convert(eltype(F), entropy - s_target)
        return nothing
    end

    res = nlsolve(
        residual_fn!,
        [Float64(mu0)];
        autodiff=:forward,
        method=:trust_region,
        xtol=1e-9,
        ftol=1e-9,
        nlsolve_kwargs...,
    )

    _refresh_constraint_state!(residual_fn!, res.zero)

    gap_vec = gap_residual(model, st_ref[], T_fm, mu_vec_ref[]; xi=xi, p_num=p_num, t_num=t_num)
    gap_norm = sqrt(sum(abs2, gap_vec))
    entropy_residual = abs(entropy_ref[] - s_target)
    residual_norm = max(gap_norm, entropy_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite

    Φ = x_state_ref[][4]
    Φbar = x_state_ref[][5]
    soft_phys = thermo_finite &&
        isfinite(Φ) && isfinite(Φbar) &&
        (-1e-3 <= Φ <= 1 + 1e-3) && (-1e-3 <= Φbar <= 1 + 1e-3) &&
        all(isfinite, masses_ref[]) && all(m -> m >= -1e-8, masses_ref[])

    accept_tol = max(Float64(residual_norm_max), 1e-3)
    near_accept_tol = min(1e-4, 10 * accept_tol)
    converged = isfinite(residual_norm) && (
        (res.f_converged && phys && residual_norm <= accept_tol) ||
        ((phys || soft_phys) && residual_norm <= near_accept_tol)
    )

    return (
        converged=converged,
        solution=Float64[
            x_state_ref[][1], x_state_ref[][2], x_state_ref[][3], x_state_ref[][4], x_state_ref[][5],
            mu_vec_ref[][1], mu_vec_ref[][2], mu_vec_ref[][3],
        ],
        x_state=x_state_ref[],
        mu_vec=mu_vec_ref[],
        omega=omega_val,
        pressure=pressure_ref[],
        rho_norm=rho_norm_ref[],
        entropy=entropy_ref[],
        energy=energy_ref[],
        masses=masses_ref[],
        iterations=res.iterations,
        residual_norm=residual_norm,
    )
end

function solve_fixedsigma_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    sigma_target::Real;
    seed_guess::AbstractVector,
    mu0::Real=0.2,
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    rho0::Real,
    physicality_check::Function=((_, _) -> true),
    __allow_hard_deprecated_internal::Bool=false,
    nlsolve_kwargs...,
)
    __allow_hard_deprecated_internal || throw(ArgumentError("Models.solve_fixedsigma_constraint is hard-deprecated in Wave-D; use Models.solve_constraint(model, FixedSigma(...), T)."))
    st_ref = Ref{Any}(nothing)
    x_state_ref = Ref{Any}(nothing)
    mu_vec_ref = Ref{Any}(nothing)
    pressure_ref = Ref{Any}(nothing)
    rho_norm_ref = Ref{Any}(nothing)
    entropy_ref = Ref{Any}(nothing)
    energy_ref = Ref{Any}(nothing)
    masses_ref = Ref{Any}(nothing)

    residual_fn! = (F, x) -> begin
        μ = x[1]
        μ_vec = normalize_mu_vec(μ)

        st = _solve_gap_with_outer_fallback(
            model,
            T_fm,
            μ;
            st_prev=st_ref[],
            seed_guess=seed_guess,
            solver_primary=solver_primary,
            solver_secondary=solver_secondary,
            residual_norm_max=residual_norm_max,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )

        if st === nothing
            _constraint_failure!(F)
            return nothing
        end

        x_state = SVector{5}(Tuple(state_vector(st)))
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)

        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / (3.0 * rho0)

        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)

        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = SVector{3}(Tuple(μ_vec))
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses

        n_B = rho_norm * rho0
        if !isfinite(n_B) || abs(n_B) <= 1e-12
            _constraint_failure!(F)
            return nothing
        end

        F[1] = convert(eltype(F), entropy / n_B - sigma_target)
        return nothing
    end

    res = nlsolve(
        residual_fn!,
        [Float64(mu0)];
        autodiff=:forward,
        method=:trust_region,
        xtol=1e-9,
        ftol=1e-9,
        nlsolve_kwargs...,
    )

    _refresh_constraint_state!(residual_fn!, res.zero)

    gap_vec = gap_residual(model, st_ref[], T_fm, mu_vec_ref[]; xi=xi, p_num=p_num, t_num=t_num)
    gap_norm = sqrt(sum(abs2, gap_vec))
    n_B_ref = rho_norm_ref[] * rho0
    sigma_residual = if isfinite(n_B_ref) && abs(n_B_ref) > 1e-12
        abs(entropy_ref[] / n_B_ref - sigma_target)
    else
        Inf
    end
    residual_norm = max(gap_norm, sigma_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
    converged = res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)

    return (
        converged=converged,
        solution=Float64[
            x_state_ref[][1], x_state_ref[][2], x_state_ref[][3], x_state_ref[][4], x_state_ref[][5],
            mu_vec_ref[][1], mu_vec_ref[][2], mu_vec_ref[][3],
        ],
        x_state=x_state_ref[],
        mu_vec=mu_vec_ref[],
        omega=omega_val,
        pressure=pressure_ref[],
        rho_norm=rho_norm_ref[],
        entropy=entropy_ref[],
        energy=energy_ref[],
        masses=masses_ref[],
        iterations=res.iterations,
        residual_norm=residual_norm,
    )
end

function solve_fixedasymrho_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    rho_target::Real,
    ud_ratio_target::Real,
    s_target::Real;
    seed_guess::AbstractVector,
    mu0::AbstractVector=Float64[0.2, 0.2, 0.2],
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    rho0::Real,
    enforce_physicality::Bool=false,
    physicality_check::Function=((_, _) -> true),
    __allow_hard_deprecated_internal::Bool=false,
    nlsolve_kwargs...,
)
    __allow_hard_deprecated_internal || throw(ArgumentError("Models.solve_fixedasymrho_constraint is hard-deprecated in Wave-D; use Models.solve_constraint(model, FixedAsymmetricRho(...), T)."))
    st_ref = Ref{Any}(nothing)
    x_state_ref = Ref{Any}(nothing)
    mu_vec_ref = Ref{Any}(nothing)
    pressure_ref = Ref{Any}(nothing)
    rho_norm_ref = Ref{Any}(nothing)
    entropy_ref = Ref{Any}(nothing)
    energy_ref = Ref{Any}(nothing)
    masses_ref = Ref{Any}(nothing)
    rho_vec_ref = Ref{Any}(nothing)

    residual_fn! = (F, x) -> begin
        μ_vec = SVector{3}(x[1], x[2], x[3])

        st = _solve_gap_with_outer_fallback(
            model,
            T_fm,
            μ_vec;
            st_prev=st_ref[],
            seed_guess=seed_guess,
            solver_primary=solver_primary,
            solver_secondary=solver_secondary,
            residual_norm_max=residual_norm_max,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )

        if st === nothing
            _constraint_failure!(F)
            return nothing
        end

        x_state = SVector{5}(Tuple(state_vector(st)))
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)

        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / (3.0 * rho0)

        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)

        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

        rho_u, rho_d, rho_s = rho_vec[1], rho_vec[2], rho_vec[3]
        nB = sum(rho_vec) / (3.0 * rho0)
        ud_ratio = if abs(rho_d) > 1e-12
            rho_u / rho_d
        else
            rho_u / (rho_d >= 0 ? 1e-12 : -1e-12)
        end

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = μ_vec
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses
        rho_vec_ref[] = SVector{3}(Tuple(rho_vec))

        F[1] = convert(eltype(F), nB - rho_target)
        F[2] = convert(eltype(F), ud_ratio - ud_ratio_target)
        F[3] = convert(eltype(F), rho_s - s_target)
        return nothing
    end

    res = nlsolve(
        residual_fn!,
        Float64.(mu0);
        autodiff=:forward,
        method=:trust_region,
        xtol=1e-9,
        ftol=1e-9,
        nlsolve_kwargs...,
    )

    _refresh_constraint_state!(residual_fn!, res.zero)

    gap_vec = gap_residual(model, st_ref[], T_fm, mu_vec_ref[]; xi=xi, p_num=p_num, t_num=t_num)
    gap_norm = sqrt(sum(abs2, gap_vec))
    rho_u, rho_d, rho_s = rho_vec_ref[][1], rho_vec_ref[][2], rho_vec_ref[][3]
    nB = sum(rho_vec_ref[]) / (3.0 * rho0)
    ud_ratio = if abs(rho_d) > 1e-12
        rho_u / rho_d
    else
        rho_u / (rho_d >= 0 ? 1e-12 : -1e-12)
    end
    nB_residual = abs(nB - rho_target)
    ud_residual = abs(ud_ratio - ud_ratio_target)
    s_residual = abs(rho_s - s_target)
    residual_norm = max(gap_norm, nB_residual, ud_residual, s_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
    converged = if enforce_physicality
        res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)
    else
        isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)
    end

    return (
        converged=converged,
        solution=Float64[
            x_state_ref[][1], x_state_ref[][2], x_state_ref[][3], x_state_ref[][4], x_state_ref[][5],
            mu_vec_ref[][1], mu_vec_ref[][2], mu_vec_ref[][3],
        ],
        x_state=x_state_ref[],
        mu_vec=mu_vec_ref[],
        omega=omega_val,
        pressure=pressure_ref[],
        rho_norm=rho_norm_ref[],
        entropy=entropy_ref[],
        energy=energy_ref[],
        masses=masses_ref[],
        iterations=res.iterations,
        residual_norm=residual_norm,
    )
end
