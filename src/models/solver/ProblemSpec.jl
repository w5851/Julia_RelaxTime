"""
    ProblemSpec

约束问题契约：统一描述 mode、维度以及求解/后处理回调。
"""
struct ProblemSpec{M,R,F,C,U,P,H,S}
    mode::M
    x_dim::Int
    theta_dim::Int
    residual!::R
    forward_solve::F
    conditions::C
    unpack_solution::U
    postprocess::P
    hard_rules::H
    selector::S
end

@inline function ProblemSpec(
    mode::ConstraintMode;
    x_dim::Int=state_dim(mode),
    theta_dim::Int=param_dim(mode),
    residual! = _unimplemented_residual,
    forward_solve=_unimplemented_forward_solve,
    conditions=_unimplemented_conditions,
    unpack_solution=identity,
    postprocess=identity,
    hard_rules=default_hard_constraint_rules(),
    selector=select_pressure_max_candidate,
)
    x_dim > 0 || throw(ArgumentError("x_dim must be positive, got $(x_dim)"))
    theta_dim > 0 || throw(ArgumentError("theta_dim must be positive, got $(theta_dim)"))
    return ProblemSpec(mode, x_dim, theta_dim, residual!, forward_solve, conditions, unpack_solution, postprocess, hard_rules, selector)
end

@inline _unimplemented_residual(F, x, theta, cfg, mode) = throw(ArgumentError("residual! is not configured for $(typeof(mode))"))
@inline _unimplemented_forward_solve(theta, cfg, mode) = throw(ArgumentError("forward_solve is not configured for $(typeof(mode))"))
@inline _unimplemented_conditions(theta, x, meta, cfg, mode) = throw(ArgumentError("conditions is not configured for $(typeof(mode))"))

@inline function _resolve_candidate_selector(kwargs::Dict{Symbol,Any})::Function
    selector = get(kwargs, :selector, nothing)
    if selector !== nothing
        selector isa Function || throw(ArgumentError("selector must be Function or nothing, got $(typeof(selector))"))
        return selector
    end

    semantic_mode = get(kwargs, :semantic_mode, :ground_state)
    if semantic_mode === :ground_state
        return select_pressure_max_candidate
    elseif semantic_mode === :constrained_manifold
        return select_residual_min_candidate
    end
    throw(ArgumentError("semantic_mode must be :ground_state or :constrained_manifold, got $(semantic_mode)"))
end

@inline function _strip_problemspec_forwardsolve_kwargs!(kwargs::Dict{Symbol,Any})
    for key in (:seed_candidates, :hard_constraints, :semantic_mode, :selector, :fixedrho_joint_solve)
        delete!(kwargs, key)
    end
    return kwargs
end

@inline function _joint_model_kind(::RPNJLModel)
    return :RPNJL
end

@inline function _joint_model_kind(::AbstractQCDModel)
    return :PNJL
end

@inline function _strip_fixedrho_joint_nlsolve_kwargs!(kwargs::Dict{Symbol,Any})
    for key in (
        :seed_guess,
        :mu0,
        :solver_primary,
        :solver_secondary,
        :physicality_check,
        :fixedrho_joint_solve,
        :nlsolve_method,
        :residual_norm_max,
        :xi,
        :p_num,
        :t_num,
    )
        delete!(kwargs, key)
    end
    return kwargs
end

function _fixedrho_joint_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedRho, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))

    seed_guess = get(kwargs, :seed_guess, nothing)
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedRho joint forward_solve"))

    xi = Float64(get(kwargs, :xi, 0.0))
    p_num = Int(get(kwargs, :p_num, 24))
    t_num = Int(get(kwargs, :t_num, 8))
    residual_norm_max = Float64(get(kwargs, :residual_norm_max, 1e-6))
    nlsolve_method = get(kwargs, :nlsolve_method, :newton)

    x0 = if length(seed_guess) >= 8
        Float64.(seed_guess[1:8])
    else
        extend_seed(Float64.(seed_guess), mode)
    end

    thermal_nodes = cached_nodes(p_num, t_num)
    params = GapParams(Float64(T_fm), thermal_nodes, xi;
        p_num=p_num,
        t_num=t_num,
        model_kind=_joint_model_kind(model),
    )
    residual_fn! = build_residual!(mode, params)

    local_nls_kwargs = Dict{Symbol,Any}(kwargs)
    _strip_problemspec_forwardsolve_kwargs!(local_nls_kwargs)
    _strip_fixedrho_joint_nlsolve_kwargs!(local_nls_kwargs)

    res = nlsolve(
        residual_fn!,
        x0;
        autodiff=:forward,
        method=nlsolve_method,
        xtol=1e-9,
        ftol=1e-9,
        pairs(local_nls_kwargs)...,
    )

    solution = Float64.(res.zero)
    x_state, mu_vec = _unpack_solution(solution; state_n=5, mu_n=3)

    pressure = -omega(model, x_state, T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
    pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
    rho_vec = ForwardDiff.gradient(pressure_mu, mu_vec)
    rho_norm = sum(rho_vec) / (3.0 * rho0)
    pressure_T = τ -> -omega(model, x_state, τ, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
    entropy = ForwardDiff.derivative(pressure_T, T_fm)
    energy = -pressure + sum(mu_vec .* rho_vec) + T_fm * entropy
    omega_val = -pressure
    masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

    residual_vec = zeros(Float64, 8)
    residual_fn!(residual_vec, solution)
    residual_norm = sqrt(sum(abs2, residual_vec))

    converged = res.f_converged && isfinite(residual_norm) && residual_norm <= residual_norm_max

    return (
        converged=converged,
        solution=solution,
        x_state=x_state,
        mu_vec=mu_vec,
        omega=omega_val,
        pressure=pressure,
        rho_norm=rho_norm,
        entropy=entropy,
        energy=energy,
        masses=masses,
        iterations=res.iterations,
        residual_norm=residual_norm,
        fixedrho_joint_solve_active=true,
    )
end

function _fixedrho_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedRho, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    fixedrho_joint_solve = get(kwargs, :fixedrho_joint_solve, false)
    fixedrho_joint_solve isa Bool || throw(ArgumentError("fixedrho_joint_solve must be Bool, got $(typeof(fixedrho_joint_solve))"))

    seed_guess = get(kwargs, :seed_guess, nothing)
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedRho forward_solve"))

    seed_pool = if haskey(kwargs, :seed_candidates)
        [Float64.(s) for s in kwargs[:seed_candidates]]
    else
        _build_default_seed_candidates(seed_guess)
    end

    physicality_check = get(kwargs, :physicality_check, ((_, _) -> true))
    hard_constraints = get(kwargs, :hard_constraints, default_hard_constraint_rules(; physicality_check=physicality_check))

    candidates = NamedTuple[]
    for (seed_index, seed) in enumerate(seed_pool)
        local raw
        try
            local_kwargs = Dict{Symbol,Any}(kwargs)
            local_kwargs[:seed_guess] = seed
            _strip_problemspec_forwardsolve_kwargs!(local_kwargs)

            solved = if fixedrho_joint_solve
                joint = try
                    _fixedrho_joint_problem_spec_forward_solve(model, mode, T_fm; pairs(local_kwargs)...)
                catch
                    nothing
                end

                joint_res = joint === nothing ? Inf : get(joint, :residual_norm, Inf)
                joint_ok = joint !== nothing && isfinite(joint_res) && joint_res <= max(get(local_kwargs, :residual_norm_max, 1e-6), 1e-3)

                if joint_ok
                    (; joint..., fixedrho_joint_solve_active=true, fixedrho_joint_fallback=false)
                else
                    legacy = _solve_constraint_fixedrho(model, T_fm, mode.rho_target; pairs(local_kwargs)...)
                    (; legacy..., fixedrho_joint_solve_active=false, fixedrho_joint_fallback=true, fixedrho_joint_attempt_residual=joint_res)
                end
            else
                legacy = _solve_constraint_fixedrho(model, T_fm, mode.rho_target; pairs(local_kwargs)...)
                (; legacy..., fixedrho_joint_solve_active=false, fixedrho_joint_fallback=false)
            end
            raw = (
                ; solved...,
                residual_norm_max=get(local_kwargs, :residual_norm_max, 1e-6),
                fixedrho_joint_solve_requested=fixedrho_joint_solve,
                fixedrho_joint_solve_active=get(solved, :fixedrho_joint_solve_active, false),
                fixedrho_joint_fallback=get(solved, :fixedrho_joint_fallback, false),
            )
            ok, failed = evaluate_hard_constraints(raw, hard_constraints)
            push!(candidates, (; raw..., hard_constraint_ok=ok, failed_constraints=failed, converged=ok, seed_index=Int(seed_index)))
        catch
            raw = (
                converged=false,
                solution=Float64[],
                x_state=zeros(5),
                mu_vec=zeros(3),
                omega=NaN,
                pressure=-Inf,
                rho_norm=NaN,
                entropy=NaN,
                energy=NaN,
                masses=zeros(3),
                iterations=0,
                residual_norm=Inf,
                residual_norm_max=Float64(get(kwargs, :residual_norm_max, 1e-6)),
                fixedrho_joint_solve_requested=fixedrho_joint_solve,
                fixedrho_joint_solve_active=false,
                fixedrho_joint_fallback=false,
            )
            push!(candidates, (; raw..., hard_constraint_ok=false, failed_constraints=Symbol[:solver_failed], seed_index=Int(seed_index)))
        end
    end

    selector_fn = _resolve_candidate_selector(kwargs)
    selected = selector_fn(candidates)
    s = selected.selected_candidate
    return (
        converged=Bool(s.converged),
        solution=Vector{Float64}(s.solution),
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
        fixedrho_joint_solve_requested=fixedrho_joint_solve,
        fixedrho_joint_solve_active=get(s, :fixedrho_joint_solve_active, false),
        fixedrho_joint_fallback=get(s, :fixedrho_joint_fallback, false),
    )
end

function _governed_mode_forward_solve(
    model::AbstractQCDModel,
    T_fm::Real,
    seed_pool,
    kwargs::Dict{Symbol,Any},
    hard_constraints,
    solve_one::Function,
)
    candidates = NamedTuple[]
    for (seed_index, seed) in enumerate(seed_pool)
        local raw
        try
            local_kwargs = Dict{Symbol,Any}(kwargs)
            local_kwargs[:seed_guess] = seed
            _strip_problemspec_forwardsolve_kwargs!(local_kwargs)

            solved = solve_one(model, T_fm, local_kwargs)
            raw = (; solved..., residual_norm_max=get(local_kwargs, :residual_norm_max, 1e-6))
            ok, failed = evaluate_hard_constraints(raw, hard_constraints)
            push!(candidates, (; raw..., hard_constraint_ok=ok, failed_constraints=failed, converged=ok, seed_index=Int(seed_index)))
        catch
            raw = (
                converged=false,
                solution=Float64[],
                x_state=zeros(5),
                mu_vec=zeros(3),
                omega=NaN,
                pressure=-Inf,
                rho_norm=NaN,
                entropy=NaN,
                energy=NaN,
                masses=zeros(3),
                iterations=0,
                residual_norm=Inf,
                residual_norm_max=Float64(get(kwargs, :residual_norm_max, 1e-6)),
            )
            push!(candidates, (; raw..., hard_constraint_ok=false, failed_constraints=Symbol[:solver_failed], seed_index=Int(seed_index)))
        end
    end

    selector_fn = _resolve_candidate_selector(kwargs)
    selected = selector_fn(candidates)
    s = selected.selected_candidate
    return (
        converged=Bool(s.converged),
        solution=Vector{Float64}(s.solution),
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

function _fixedentropy_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedEntropy, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    seed_guess = get(kwargs, :seed_guess, nothing)
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedEntropy forward_solve"))

    seed_pool = if haskey(kwargs, :seed_candidates)
        [Float64.(s) for s in kwargs[:seed_candidates]]
    else
        _build_default_seed_candidates(seed_guess)
    end
    physicality_check = get(kwargs, :physicality_check, ((_, _) -> true))
    hard_constraints = get(kwargs, :hard_constraints, default_hard_constraint_rules(; physicality_check=physicality_check))

    solve_one = (m, t, local_kwargs) -> begin
        haskey(local_kwargs, :rho0) || throw(ArgumentError("rho0 is required for ProblemSpec FixedEntropy forward_solve"))
        return _solve_constraint_fixedentropy(m, t, mode.s_target; pairs(local_kwargs)...)
    end
    return _governed_mode_forward_solve(model, T_fm, seed_pool, kwargs, hard_constraints, solve_one)
end

function _fixedsigma_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedSigma, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    seed_guess = get(kwargs, :seed_guess, nothing)
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedSigma forward_solve"))

    seed_pool = if haskey(kwargs, :seed_candidates)
        [Float64.(s) for s in kwargs[:seed_candidates]]
    else
        _build_default_seed_candidates(seed_guess)
    end
    physicality_check = get(kwargs, :physicality_check, ((_, _) -> true))
    hard_constraints = get(kwargs, :hard_constraints, default_hard_constraint_rules(; physicality_check=physicality_check))

    solve_one = (m, t, local_kwargs) -> begin
        haskey(local_kwargs, :rho0) || throw(ArgumentError("rho0 is required for ProblemSpec FixedSigma forward_solve"))
        return _solve_constraint_fixedsigma(m, t, mode.sigma_target; pairs(local_kwargs)...)
    end
    return _governed_mode_forward_solve(model, T_fm, seed_pool, kwargs, hard_constraints, solve_one)
end

function _fixedasymrho_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedAsymmetricRho, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    seed_guess = get(kwargs, :seed_guess, nothing)
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedAsymmetricRho forward_solve"))

    seed_pool = if haskey(kwargs, :seed_candidates)
        [Float64.(s) for s in kwargs[:seed_candidates]]
    else
        _build_default_seed_candidates(seed_guess)
    end
    physicality_check = get(kwargs, :physicality_check, ((_, _) -> true))
    hard_constraints = get(kwargs, :hard_constraints, default_hard_constraint_rules(; physicality_check=physicality_check))

    solve_one = (m, t, local_kwargs) -> begin
        haskey(local_kwargs, :rho0) || throw(ArgumentError("rho0 is required for ProblemSpec FixedAsymmetricRho forward_solve"))
        return _solve_constraint_fixedasymrho(m, t, mode.rho_target, mode.ud_ratio_target, mode.s_target; pairs(local_kwargs)...)
    end
    return _governed_mode_forward_solve(model, T_fm, seed_pool, kwargs, hard_constraints, solve_one)
end

@inline function build_problem_spec(mode::ConstraintMode; kwargs...)
    if isempty(kwargs)
        components = build_constraint_components(mode)
        forward_solve = if mode isa FixedRho
            (model, T_fm; fwd_kwargs...) -> _fixedrho_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        elseif mode isa FixedEntropy
            (model, T_fm; fwd_kwargs...) -> _fixedentropy_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        elseif mode isa FixedSigma
            (model, T_fm; fwd_kwargs...) -> _fixedsigma_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        elseif mode isa FixedAsymmetricRho
            (model, T_fm; fwd_kwargs...) -> _fixedasymrho_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        else
            (model, T_fm; fwd_kwargs...) -> solve_constraint(model, mode, T_fm; fwd_kwargs...)
        end
        return ProblemSpec(
            mode;
            x_dim=constraint_total_dim(components),
            conditions=params -> build_conditions(mode, params),
            forward_solve=forward_solve,
        )
    end
    return ProblemSpec(mode; kwargs...)
end
