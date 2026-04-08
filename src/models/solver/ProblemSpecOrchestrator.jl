@inline function _build_governed_attempt_plan(
    mode::ConstraintMode,
    primary_seed::AbstractVector{<:Real},
    provided_seed_pool,
    default_seed_pool,
    extra_constraints::ExtraConstraints;
    primary_method,
    primary_use_fallback::Bool,
    fallback_method,
    extra_fallback_seeds=Vector{Vector{Float64}}(),
)
    seed_extend = (seed, local_mode) -> _extend_seed_with_extra(seed, local_mode, extra_constraints)
    seed_pool = build_seed_pool(mode;
        primary_seed=primary_seed,
        extra_seed_pool=extra_fallback_seeds,
        provided_seed_pool=provided_seed_pool,
        default_seed_pool=default_seed_pool,
        seed_extend=seed_extend,
    )

    isempty(seed_pool) && throw(ArgumentError("build_seed_pool produced empty pool for $(typeof(mode))"))
    primary_seed_vec = seed_pool[1].seed
    fallback_seeds = [entry.seed for entry in seed_pool[2:end]]

    attempt_plan = NamedTuple[]
    push!(attempt_plan, (
        seed=primary_seed_vec,
        method=primary_method,
        use_fallback=primary_use_fallback,
        fallback_method=fallback_method,
        attempt_origin=:primary,
    ))

    if primary_method != :trust_region
        push!(attempt_plan, (
            seed=primary_seed_vec,
            method=:trust_region,
            use_fallback=false,
            fallback_method=:trust_region,
            attempt_origin=:method_rescue,
        ))
    end

    for seed in fallback_seeds
        push!(attempt_plan, (
            seed=seed,
            method=:trust_region,
            use_fallback=false,
            fallback_method=:trust_region,
            attempt_origin=:seed_rescue,
        ))
    end

    return attempt_plan
end

function _execute_governed_attempt_plan(
    attempt_plan,
    kwargs::Dict{Symbol,Any},
    hard_constraints,
    selector_fn::Function,
    solve_attempt::Function,
    failure_attempt::Function,
)
    evaluate_all_attempts = Bool(get(kwargs, :evaluate_all_attempts, false))
    candidates = execute_attempt_pool(attempt_plan;
        stop_on_first_success=true,
        evaluate_all_attempts=evaluate_all_attempts,
        evaluate_attempt=(attempt_cfg, attempt_index) -> begin
            local_kwargs = Dict{Symbol,Any}(kwargs)
            local_kwargs[:seed_guess] = attempt_cfg.seed
            local_kwargs[:nlsolve_method] = attempt_cfg.method
            local_kwargs[:trust_region_fallback] = attempt_cfg.use_fallback
            local_kwargs[:fallback_method] = attempt_cfg.fallback_method

            raw = solve_attempt(local_kwargs, attempt_cfg, attempt_index)
            ok, failed = evaluate_hard_constraints(raw, hard_constraints)
            candidate = (; raw..., hard_constraint_ok=ok, failed_constraints=failed, converged=ok, seed_index=Int(attempt_index))
            normalized = normalize_governance_candidate(candidate;
                seed_index=Int(attempt_index),
            )
            merged = (; candidate...,
                converged=normalized.converged,
                pressure=normalized.pressure,
                residual_norm=normalized.residual_norm,
                hard_constraint_ok=normalized.hard_constraint_ok,
                failed_constraints=normalized.failed_constraints,
                seed_index=normalized.seed_index,
            )
            success = evaluate_candidate_success(merged;
                residual_norm_max=Float64(get(local_kwargs, :residual_norm_max, 1e-6)),
            )
            return merged, success
        end,
        on_error=(attempt_cfg, attempt_index, err) -> begin
            raw = failure_attempt(kwargs, attempt_cfg, attempt_index)
            err_kind = classify_attempt_error(err)
            err_msg = normalize_error_message(err)
            candidate = (; raw..., hard_constraint_ok=false, failed_constraints=Symbol[:solver_failed], seed_index=Int(attempt_index))
            normalized = normalize_governance_candidate(candidate;
                seed_index=Int(attempt_index),
            )
            merged = (; candidate...,
                converged=normalized.converged,
                pressure=normalized.pressure,
                residual_norm=normalized.residual_norm,
                hard_constraint_ok=normalized.hard_constraint_ok,
                failed_constraints=normalized.failed_constraints,
                seed_index=normalized.seed_index,
                error_kind=err_kind,
                error_msg=err_msg,
            )
            success = evaluate_candidate_success(merged;
                residual_norm_max=Float64(get(kwargs, :residual_norm_max, 1e-6)),
            )
            return merged, success
        end,
    )

    selected = selector_fn(candidates)
    return selected, candidates
end

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
    for key in (:seed_candidates, :hard_constraints, :semantic_mode, :selector, :fixedrho_joint_solve, :continuity_seed, :extra_constraints, :fixedmu_use_problem_spec, :legacy_fallback_plugin, :diagnostic_level)
        delete!(kwargs, key)
    end
    return kwargs
end

@inline function _resolve_optional_seed_candidates(kwargs::Dict{Symbol,Any})::Vector{Vector{Float64}}
    if !haskey(kwargs, :seed_candidates)
        return Vector{Vector{Float64}}()
    end
    raw = kwargs[:seed_candidates]
    raw === nothing && return Vector{Vector{Float64}}()
    return [Float64.(s) for s in raw]
end

function _fixedmu_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedMu, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    diagnostic_level = _resolve_diagnostic_level(kwargs)

    μ_fm = get(kwargs, :μ_fm, nothing)
    μ_fm === nothing && throw(ArgumentError("μ_fm is required for ProblemSpec FixedMu forward_solve"))

    seed_guess = get(kwargs, :seed_guess, nothing)
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedMu forward_solve"))

    local_kwargs = Dict{Symbol,Any}(kwargs)
    _strip_problemspec_forwardsolve_kwargs!(local_kwargs)

    solved = _solve_constraint_fixedmu(model, T_fm, μ_fm; pairs(local_kwargs)...)
    result = (; solved..., fixedmu_problem_spec_active=true)
    if diagnostic_level === :none
        return result
    end

    seed_source = Bool(get(kwargs, :continuity_seed, false)) ? :warm_start : :seed
    summary_candidate = (
        ; solved...,
        attempt_origin=:seed,
        hard_constraint_ok=Bool(get(solved, :hard_constraint_ok, false)),
        failed_constraints=Symbol.(get(solved, :failed_constraints, Symbol[])),
        error_kind=:none,
        error_msg="",
        selection_reason=Symbol(get(solved, :selection_reason, :none)),
        endpoint_cause=Bool(get(solved, :converged, false)) ? :converged : :nonconvergence,
        continuity_distance=nothing,
    )
    summary = _solver_diagnostic_from_candidate(summary_candidate; seed_source=seed_source)
    if diagnostic_level === :summary
        return merge(result, (diagnostic=summary,))
    end
    return merge(result, (diagnostic=(; summary..., candidates=[summary]),))
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
        :continuity_seed,
        :trust_region_fallback,
        :fallback_method,
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
    nlsolve_method = get(kwargs, :nlsolve_method, :trust_region)
    trust_region_fallback = Bool(get(kwargs, :trust_region_fallback, true))
    fallback_method = get(kwargs, :fallback_method, :trust_region)
    physicality_check = get(kwargs, :physicality_check, ((_, _) -> true))

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

    postprocess_solution = function (solution)
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
        thermo_finite = isfinite(omega_val) && isfinite(pressure) && isfinite(rho_norm) && isfinite(entropy) && isfinite(energy)
        phys = physicality_check(x_state, masses) && thermo_finite

        residual_vec = zeros(Float64, 8)
        residual_fn!(residual_vec, solution)
        residual_norm = sqrt(sum(abs2, residual_vec))

        return (
            x_state=x_state,
            mu_vec=mu_vec,
            omega=omega_val,
            pressure=pressure,
            rho_norm=rho_norm,
            entropy=entropy,
            energy=energy,
            masses=masses,
            residual_norm=residual_norm,
            phys=phys,
        )
    end

    cache = Dict{Symbol,NamedTuple}()
    solve_once = function (method::Symbol, seed::Vector{Float64})
        res = nlsolve(
            residual_fn!,
            seed;
            autodiff=:forward,
            method=method,
            xtol=1e-9,
            ftol=1e-9,
            pairs(local_nls_kwargs)...,
        )

        solution = Float64.(res.zero)
        pp = postprocess_solution(solution)
        converged = Bool(res.f_converged) && pp.phys && isfinite(pp.residual_norm) && pp.residual_norm <= residual_norm_max

        cache[method] = (
            res=res,
            solution=solution,
            x_state=pp.x_state,
            mu_vec=pp.mu_vec,
            omega=pp.omega,
            pressure=pp.pressure,
            rho_norm=pp.rho_norm,
            entropy=pp.entropy,
            energy=pp.energy,
            masses=pp.masses,
            residual_norm=pp.residual_norm,
            converged=converged,
        )

        return (
            x=solution,
            converged=converged,
            residual_norm=pp.residual_norm,
            score=Float64(pp.omega),
        )
    end

    policy = RootPolicy(
        primary_method=nlsolve_method,
        fallback_method=fallback_method,
        use_fallback=trust_region_fallback,
        use_multiseed=false,
        residual_norm_max=residual_norm_max,
        require_converged=true,
        diagnostics_level=:basic,
    )

    solved = solve_root_with_policy(solve_once, x0; policy=policy)
    selected_attempt = solved.diagnostics.attempts[solved.diagnostics.selected_attempt]
    selected_method = selected_attempt.method

    if !haskey(cache, selected_method)
        _ = solve_once(selected_method, copy(x0))
    end
    picked = cache[selected_method]

    return (
        converged=picked.converged,
        solution=picked.solution,
        x_state=picked.x_state,
        mu_vec=picked.mu_vec,
        omega=picked.omega,
        pressure=picked.pressure,
        rho_norm=picked.rho_norm,
        entropy=picked.entropy,
        energy=picked.energy,
        masses=picked.masses,
        iterations=picked.res.iterations,
        residual_norm=picked.residual_norm,
        fixedrho_joint_solve_active=true,
        fixedrho_joint_selected_method=selected_method,
        fixedrho_joint_selected_quality=selected_attempt.quality_tag,
        fixedrho_joint_fallback_used=(selected_attempt.quality_tag == :fallback),
    )
end

function _fixedrho_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedRho, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    diagnostic_level = _resolve_diagnostic_level(kwargs)
    extra_constraints = _resolve_extra_constraints(kwargs)
    fixedrho_joint_solve = get(kwargs, :fixedrho_joint_solve, true)
    fixedrho_joint_solve isa Bool || throw(ArgumentError("fixedrho_joint_solve must be Bool, got $(typeof(fixedrho_joint_solve))"))
    fixedrho_joint_solve || throw(ArgumentError("fixedrho_joint_solve=false is no longer supported; FixedRho uses joint solve only"))

    seed_guess = get(kwargs, :seed_guess, nothing)
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedRho forward_solve"))

    provided_seed_pool = _resolve_optional_seed_candidates(kwargs)
    default_seed_pool = _build_default_seed_candidates(mode, seed_guess)

    primary_seed = Float64.(seed_guess)
    primary_method = if haskey(kwargs, :nlsolve_method)
        kwargs[:nlsolve_method]
    else
        Bool(get(kwargs, :continuity_seed, false)) ? :newton : :trust_region
    end
    primary_use_fallback = Bool(get(kwargs, :trust_region_fallback, true))
    fallback_method = get(kwargs, :fallback_method, :trust_region)

    legacy_seed = if length(primary_seed) >= 5
        _extend_seed_with_extra(Float64.(extend_seed(Float64.(primary_seed[1:5]), mode)), mode, extra_constraints)
    else
        _extend_seed_with_extra(Float64.(extend_seed(primary_seed, mode)), mode, extra_constraints)
    end
    attempt_plan = _build_governed_attempt_plan(
        mode,
        primary_seed,
        provided_seed_pool,
        default_seed_pool,
        extra_constraints;
        primary_method=primary_method,
        primary_use_fallback=primary_use_fallback,
        fallback_method=fallback_method,
        extra_fallback_seeds=[legacy_seed],
    )

    physicality_check = get(kwargs, :physicality_check, ((_, _) -> true))
    hard_constraints = get(kwargs, :hard_constraints, default_hard_constraint_rules(; physicality_check=physicality_check))
    hard_constraints = _append_extra_feasible_rule(hard_constraints, extra_constraints, mode)

    selector_fn = _resolve_candidate_selector(kwargs)
    selected, candidates = _execute_governed_attempt_plan(
        attempt_plan,
        kwargs,
        hard_constraints,
        selector_fn,
        function (local_kwargs, attempt_cfg, _)
            _strip_problemspec_forwardsolve_kwargs!(local_kwargs)
            solved = _fixedrho_joint_problem_spec_forward_solve(model, mode, T_fm; pairs(local_kwargs)...)
            return (
                converged=Bool(solved.converged),
                solution=Float64.(solved.solution),
                x_state=solved.x_state,
                mu_vec=solved.mu_vec,
                omega=Float64(solved.omega),
                pressure=Float64(solved.pressure),
                rho_norm=Float64(solved.rho_norm),
                entropy=Float64(solved.entropy),
                energy=Float64(solved.energy),
                masses=solved.masses,
                iterations=Int(solved.iterations),
                residual_norm=Float64(solved.residual_norm),
                residual_norm_max=get(local_kwargs, :residual_norm_max, 1e-6),
                fixedrho_joint_solve_requested=fixedrho_joint_solve,
                fixedrho_joint_solve_active=get(solved, :fixedrho_joint_solve_active, false),
                fixedrho_joint_fallback=false,
                fixedrho_joint_selected_method=get(solved, :fixedrho_joint_selected_method, :none),
                fixedrho_joint_selected_quality=get(solved, :fixedrho_joint_selected_quality, :bad),
                fixedrho_joint_fallback_used=get(solved, :fixedrho_joint_fallback_used, false),
                fixedrho_joint_attempt_origin=attempt_cfg.attempt_origin,
            )
        end,
        function (base_kwargs, attempt_cfg, _)
            return (
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
                residual_norm_max=Float64(get(base_kwargs, :residual_norm_max, 1e-6)),
                fixedrho_joint_solve_requested=fixedrho_joint_solve,
                fixedrho_joint_solve_active=false,
                fixedrho_joint_fallback=false,
                fixedrho_joint_selected_method=:none,
                fixedrho_joint_selected_quality=:bad,
                fixedrho_joint_fallback_used=false,
                fixedrho_joint_attempt_origin=attempt_cfg.attempt_origin,
            )
        end,
    )
    s = selected.selected_candidate
    result = (
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
        error_kind=(hasproperty(s, :error_kind) ? Symbol(getproperty(s, :error_kind)) : :none),
        error_msg=(hasproperty(s, :error_msg) ? String(getproperty(s, :error_msg)) : ""),
        fixedrho_joint_solve_requested=fixedrho_joint_solve,
        fixedrho_joint_solve_active=get(s, :fixedrho_joint_solve_active, false),
        fixedrho_joint_fallback=get(s, :fixedrho_joint_fallback, false),
        fixedrho_joint_selected_method=get(s, :fixedrho_joint_selected_method, :none),
        fixedrho_joint_selected_quality=get(s, :fixedrho_joint_selected_quality, :bad),
        fixedrho_joint_fallback_used=get(s, :fixedrho_joint_fallback_used, false),
    )
    seed_source = Bool(get(kwargs, :continuity_seed, false)) ? :warm_start : :seed
    return _attach_solver_diagnostic(result, s, candidates, diagnostic_level;
        seed_source=seed_source,
        selection_reason=selected.selection_reason,
    )
end

function _governed_nonrho_problem_spec_forward_solve(
    model::AbstractQCDModel,
    mode::ConstraintMode,
    T_fm::Real,
    kwargs::Dict{Symbol,Any},
    mode_label::AbstractString,
    attempt_origin_key::Symbol,
    solve_mode_constraint::Function,
)
    extra_constraints = _resolve_extra_constraints(kwargs)
    diagnostic_level = _resolve_diagnostic_level(kwargs)
    seed_guess = get(kwargs, :seed_guess, nothing)
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec $(mode_label) forward_solve"))

    provided_seed_pool = _resolve_optional_seed_candidates(kwargs)
    default_seed_pool = _build_default_seed_candidates(mode, seed_guess)

    primary_seed = Float64.(seed_guess)
    primary_method = if haskey(kwargs, :nlsolve_method)
        kwargs[:nlsolve_method]
    else
        Bool(get(kwargs, :continuity_seed, false)) ? :newton : :trust_region
    end
    primary_use_fallback = Bool(get(kwargs, :trust_region_fallback, true))
    fallback_method = get(kwargs, :fallback_method, :trust_region)

    attempt_plan = _build_governed_attempt_plan(
        mode,
        primary_seed,
        provided_seed_pool,
        default_seed_pool,
        extra_constraints;
        primary_method=primary_method,
        primary_use_fallback=primary_use_fallback,
        fallback_method=fallback_method,
    )

    physicality_check = get(kwargs, :physicality_check, ((_, _) -> true))
    hard_constraints = get(kwargs, :hard_constraints, default_hard_constraint_rules(; physicality_check=physicality_check))
    hard_constraints = _append_extra_feasible_rule(hard_constraints, extra_constraints, mode)

    selector_fn = _resolve_candidate_selector(kwargs)
    selected, candidates = _execute_governed_attempt_plan(
        attempt_plan,
        kwargs,
        hard_constraints,
        selector_fn,
        function (local_kwargs, attempt_cfg, _)
            _strip_problemspec_forwardsolve_kwargs!(local_kwargs)
            delete!(local_kwargs, :trust_region_fallback)
            delete!(local_kwargs, :fallback_method)

            haskey(local_kwargs, :rho0) || throw(ArgumentError("rho0 is required for ProblemSpec $(mode_label) forward_solve"))
            allow_legacy_fallback = Bool(get(kwargs, :allow_legacy_fallback, false))
            if haskey(kwargs, :legacy_fallback_plugin)
                allow_legacy_fallback = Bool(kwargs[:legacy_fallback_plugin])
            end
            local_kwargs[:allow_legacy_fallback] = allow_legacy_fallback

            solved = solve_mode_constraint(local_kwargs)
            solver_converged = Bool(solved.converged)
            attempt_quality = if solver_converged
                :degraded
            else
                :bad
            end
            raw = (
                ; solved...,
                residual_norm_max=get(local_kwargs, :residual_norm_max, 1e-6),
                solver_converged=solver_converged,
                governed_selected_method=attempt_cfg.method,
                governed_selected_quality=attempt_quality,
                governed_fallback_used=(attempt_cfg.attempt_origin == :method_rescue),
                governed_attempt_origin=attempt_cfg.attempt_origin,
            )
            return (; raw..., attempt_origin_key => attempt_cfg.attempt_origin)
        end,
        function (base_kwargs, attempt_cfg, _)
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
                residual_norm_max=Float64(get(base_kwargs, :residual_norm_max, 1e-6)),
                solver_converged=false,
                governed_selected_method=attempt_cfg.method,
                governed_selected_quality=:bad,
                governed_fallback_used=(attempt_cfg.attempt_origin == :method_rescue),
                governed_attempt_origin=attempt_cfg.attempt_origin,
            )
            return (; raw..., attempt_origin_key => attempt_cfg.attempt_origin)
        end,
    )
    s = selected.selected_candidate
    result = (
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
        error_kind=(hasproperty(s, :error_kind) ? Symbol(getproperty(s, :error_kind)) : :none),
        error_msg=(hasproperty(s, :error_msg) ? String(getproperty(s, :error_msg)) : ""),
        governed_selected_method=(hasproperty(s, :governed_selected_method) ? getproperty(s, :governed_selected_method) : :none),
        governed_selected_quality=(hasproperty(s, :governed_selected_quality) ? getproperty(s, :governed_selected_quality) : :bad),
        governed_fallback_used=(hasproperty(s, :governed_fallback_used) ? Bool(getproperty(s, :governed_fallback_used)) : false),
        legacy_fallback_used=(hasproperty(s, :legacy_fallback_used) ? Bool(getproperty(s, :legacy_fallback_used)) : false),
    )
    seed_source = Bool(get(kwargs, :continuity_seed, false)) ? :warm_start : :seed
    return _attach_solver_diagnostic(result, s, candidates, diagnostic_level;
        seed_source=seed_source,
        selection_reason=selected.selection_reason,
    )
end

function _fixedentropy_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedEntropy, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    return _governed_nonrho_problem_spec_forward_solve(
        model,
        mode,
        T_fm,
        kwargs,
        "FixedEntropy",
        :entropy_attempt_origin,
        local_kwargs -> _solve_constraint_fixedentropy(model, T_fm, mode.s_target; pairs(local_kwargs)...),
    )
end

function _fixedsigma_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedSigma, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    return _governed_nonrho_problem_spec_forward_solve(
        model,
        mode,
        T_fm,
        kwargs,
        "FixedSigma",
        :sigma_attempt_origin,
        local_kwargs -> _solve_constraint_fixedsigma(model, T_fm, mode.sigma_target; pairs(local_kwargs)...),
    )
end

function _fixedasymrho_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedAsymmetricRho, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    return _governed_nonrho_problem_spec_forward_solve(
        model,
        mode,
        T_fm,
        kwargs,
        "FixedAsymmetricRho",
        :asym_attempt_origin,
        local_kwargs -> _solve_constraint_fixedasymrho(model, T_fm, mode.rho_target, mode.ud_ratio_target, mode.s_target; pairs(local_kwargs)...),
    )
end
