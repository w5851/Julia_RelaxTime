"""constraint_solver.jl

Models 侧约束求解（阶段 B 内核下沉）：
- 先以 FixedEntropy 为起点，将单变量 μ 外层约束求解从 PNJL 兼容层下沉到 models。
"""

using StaticArrays
using NLsolve
using ForwardDiff
using Statistics: mean
using Main.Constants_PNJL: ρ0_inv_fm3

const rho0 = ρ0_inv_fm3

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

function _solve_constraint_fixedmu(
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
)
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

function _solve_constraint_fixedrho(
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
    nlsolve_kwargs...,
)
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

        F[1] = convert(eltype(F), rho_norm - rho_target)
        return nothing
    end

    function run_outer_attempt(mu_init::Real, method::Symbol)
        local res
        try
            res = nlsolve(
                residual_fn!,
                [Float64(mu_init)];
                autodiff=:forward,
                method=method,
                xtol=1e-9,
                ftol=1e-9,
                nlsolve_kwargs...,
            )
        catch
            return nothing
        end

        _refresh_constraint_state!(residual_fn!, res.zero)

        if st_ref[] === nothing || x_state_ref[] === nothing || mu_vec_ref[] === nothing ||
           pressure_ref[] === nothing || rho_norm_ref[] === nothing || entropy_ref[] === nothing ||
           energy_ref[] === nothing || masses_ref[] === nothing
            return nothing
        end

        gap_vec = gap_residual(model, st_ref[], T_fm, mu_vec_ref[]; xi=xi, p_num=p_num, t_num=t_num)
        gap_norm = sqrt(sum(abs2, gap_vec))
        rho_residual = abs(rho_norm_ref[] - rho_target)
        residual_norm = max(gap_norm, rho_residual)

        omega_val = -pressure_ref[]
        thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
        phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
        converged = res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-4)

        return (
            mode=method,
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

    function evaluate_mu_candidate(mu_eval::Real)
        μ_vec = normalize_mu_vec(mu_eval)
        st = _solve_gap_with_outer_fallback(
            model,
            T_fm,
            mu_eval;
            st_prev=nothing,
            seed_guess=seed_guess,
            solver_primary=solver_primary,
            solver_secondary=solver_secondary,
            residual_norm_max=residual_norm_max,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )
        st === nothing && return nothing

        x_state = SVector{5}(Tuple(state_vector(st)))
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / (3.0 * rho0)
        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)
        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

        gap_vec = gap_residual(model, st, T_fm, SVector{3}(Tuple(μ_vec)); xi=xi, p_num=p_num, t_num=t_num)
        gap_norm = sqrt(sum(abs2, gap_vec))
        rho_residual = abs(rho_norm - rho_target)
        residual_norm = max(gap_norm, rho_residual)

        omega_val = -pressure
        thermo_finite = isfinite(omega_val) && isfinite(pressure) && isfinite(rho_norm) && isfinite(entropy) && isfinite(energy)
        phys = physicality_check(x_state, masses) && thermo_finite
        converged = phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-4)

        return (
            mode=:direct_mu,
            converged=converged,
            solution=Float64[
                x_state[1], x_state[2], x_state[3], x_state[4], x_state[5],
                μ_vec[1], μ_vec[2], μ_vec[3],
            ],
            x_state=x_state,
            mu_vec=SVector{3}(Tuple(μ_vec)),
            omega=omega_val,
            pressure=pressure,
            rho_norm=rho_norm,
            entropy=entropy,
            energy=energy,
            masses=masses,
            iterations=0,
            residual_norm=residual_norm,
        )
    end

    mu_seed = default_mu0_from_seed(seed_guess)
    attempt_specs = Tuple{Float64, Symbol}[]
    push!(attempt_specs, (Float64(mu0), nlsolve_method))
    nlsolve_method !== :trust_region && push!(attempt_specs, (Float64(mu0), :trust_region))
    push!(attempt_specs, (Float64(mu_seed), :trust_region))
    for μ0 in (0.0, 0.2, 0.5, 0.8, 1.2, 1.6, 2.0)
        push!(attempt_specs, (Float64(μ0), :trust_region))
    end

    best = nothing
    for (mu_init, method) in attempt_specs
        cand = run_outer_attempt(mu_init, method)
        cand === nothing && continue
        if best === nothing
            best = cand
            continue
        end
        if cand.converged && !best.converged
            best = cand
        elseif cand.converged == best.converged
            if cand.residual_norm < best.residual_norm || (cand.residual_norm == best.residual_norm && cand.pressure > best.pressure)
                best = cand
            end
        end
    end

    if best === nothing || !best.converged
        mu_min = 0.0
        mu_max = max(Float64(mu_seed), Float64(mu0), 2.0)
        mu_grid = collect(range(mu_min, mu_max; length=25))
        grid_candidates = NamedTuple[]
        for μ in mu_grid
            cand = evaluate_mu_candidate(μ)
            cand === nothing && continue
            push!(grid_candidates, cand)
        end

        if !isempty(grid_candidates)
            for cand in grid_candidates
                if best === nothing
                    best = cand
                elseif cand.converged && !best.converged
                    best = cand
                elseif cand.converged == best.converged
                    if cand.residual_norm < best.residual_norm || (cand.residual_norm == best.residual_norm && cand.pressure > best.pressure)
                        best = cand
                    end
                end
            end

            sort!(grid_candidates; by=c -> c.mu_vec[1])
            for i in 1:(length(grid_candidates)-1)
                c1 = grid_candidates[i]
                c2 = grid_candidates[i + 1]
                f1 = c1.rho_norm - rho_target
                f2 = c2.rho_norm - rho_target
                if !(isfinite(f1) && isfinite(f2))
                    continue
                end
                if f1 == 0.0
                    best = c1
                    break
                end
                if f1 * f2 > 0.0
                    continue
                end

                left = c1.mu_vec[1]
                right = c2.mu_vec[1]
                f_left = f1
                for _ in 1:16
                    mid = 0.5 * (left + right)
                    cmid = evaluate_mu_candidate(mid)
                    cmid === nothing && break

                    if best === nothing
                        best = cmid
                    elseif cmid.converged && !best.converged
                        best = cmid
                    elseif cmid.converged == best.converged
                        if cmid.residual_norm < best.residual_norm || (cmid.residual_norm == best.residual_norm && cmid.pressure > best.pressure)
                            best = cmid
                        end
                    end

                    f_mid = cmid.rho_norm - rho_target
                    if abs(f_mid) <= max(Float64(residual_norm_max), 1e-4)
                        break
                    end
                    if f_left * f_mid <= 0.0
                        right = mid
                    else
                        left = mid
                        f_left = f_mid
                    end
                end
                break
            end
        end
    end

    best === nothing && throw(ArgumentError("FixedRho outer solve failed for all attempts"))

    if !best.converged && model isa AbstractPNJLModel
        seed_legacy = if length(seed_guess) >= 8
            Float64.(seed_guess[1:8])
        else
            Main.Models.extend_seed(Float64.(seed_guess), Main.Models.FixedRho(rho_target))
        end

        legacy_result = try
            Main.Models.solve(
                Main.Models.FixedRho(rho_target),
                T_fm;
                xi=xi,
                seed_strategy=Main.Models.DefaultSeed(seed_legacy, seed_legacy, :hadron),
                p_num=p_num,
                t_num=t_num,
                nlsolve_method=nlsolve_method,
                trust_region_fallback=true,
                residual_norm_max=residual_norm_max,
                nlsolve_kwargs...,
            )
        catch
            nothing
        end

        if legacy_result !== nothing && legacy_result.converged
            return (
                converged=true,
                solution=Float64.(legacy_result.solution),
                x_state=legacy_result.x_state,
                mu_vec=legacy_result.mu_vec,
                omega=legacy_result.omega,
                pressure=legacy_result.pressure,
                rho_norm=legacy_result.rho_norm,
                entropy=legacy_result.entropy,
                energy=legacy_result.energy,
                masses=legacy_result.masses,
                iterations=legacy_result.iterations,
                residual_norm=legacy_result.residual_norm,
            )
        end
    end

    return (
        converged=best.converged,
        solution=best.solution,
        x_state=best.x_state,
        mu_vec=best.mu_vec,
        omega=best.omega,
        pressure=best.pressure,
        rho_norm=best.rho_norm,
        entropy=best.entropy,
        energy=best.energy,
        masses=best.masses,
        iterations=best.iterations,
        residual_norm=best.residual_norm,
    )
end

function _solve_constraint_fixedentropy(
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
    nlsolve_kwargs...,
)
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

function _solve_constraint_fixedsigma(
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
    nlsolve_kwargs...,
)
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

function _solve_constraint_fixedasymrho(
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
    nlsolve_kwargs...,
)
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
