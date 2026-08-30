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
    raw_candidates = execute_attempt_pool(attempt_plan;
        stop_on_first_success=true,
        evaluate_all_attempts=evaluate_all_attempts,
        evaluate_attempt=(attempt_cfg, attempt_index) -> begin
            local_kwargs = Dict{Symbol,Any}(kwargs)
            local_kwargs[:seed_guess] = attempt_cfg.seed
            local_kwargs[:nlsolve_method] = attempt_cfg.method
            local_kwargs[:trust_region_fallback] = attempt_cfg.use_fallback
            local_kwargs[:fallback_method] = attempt_cfg.fallback_method
            record_governed_attempt!(get(local_kwargs, :work_telemetry, nothing); origin=attempt_cfg.attempt_origin)

            raw = solve_attempt(local_kwargs, attempt_cfg, attempt_index)
            ok, failed = evaluate_hard_constraints(raw, hard_constraints)
            # A governed attempt may deliberately return a structured
            # incomplete-output candidate (for example when a backend
            # produces `nothing` for one thermodynamic field).  Preserve
            # that diagnostic reason alongside the ordinary hard-rule
            # failures instead of allowing it to be overwritten here.
            raw_failed_value = hasproperty(raw, :failed_constraints) ?
                getproperty(raw, :failed_constraints) : Symbol[]
            raw_failed = raw_failed_value === nothing ? Symbol[] : Symbol.(raw_failed_value)
            failed = unique(vcat(raw_failed, failed))
            ok = ok && isempty(raw_failed)
            raw_error_kind_value = hasproperty(raw, :error_kind) ? getproperty(raw, :error_kind) : :none
            raw_error_kind = raw_error_kind_value === nothing ? :none : Symbol(raw_error_kind_value)
            raw_error_msg_value = hasproperty(raw, :error_msg) ? getproperty(raw, :error_msg) : ""
            raw_error_msg = raw_error_msg_value === nothing ? "" : String(raw_error_msg_value)
            merged = build_governance_candidate(raw;
                hard_constraint_ok=ok,
                failed_constraints=failed,
                seed_index=Int(attempt_index),
                residual_norm_max=Float64(get(local_kwargs, :residual_norm_max, 1e-6)),
                error_kind=raw_error_kind,
                error_msg=raw_error_msg,
            )
            success = evaluate_candidate_success(merged;
                residual_norm_max=Float64(get(local_kwargs, :residual_norm_max, 1e-6)),
            )
            return merged, success
        end,
        on_error=(attempt_cfg, attempt_index, err) -> begin
            record_solver_exception!(get(kwargs, :work_telemetry, nothing))
            raw = failure_attempt(kwargs, attempt_cfg, attempt_index)
            err_kind = classify_attempt_error(err)
            err_msg = normalize_error_message(err)
            merged = build_governance_candidate(raw;
                hard_constraint_ok=false,
                failed_constraints=Symbol[:solver_failed],
                seed_index=Int(attempt_index),
                residual_norm_max=Float64(get(kwargs, :residual_norm_max, 1e-6)),
                error_kind=err_kind,
                error_msg=err_msg,
            )
            success = evaluate_candidate_success(merged;
                residual_norm_max=Float64(get(kwargs, :residual_norm_max, 1e-6)),
            )
            return merged, success
        end,
    )

    residual_norm_max = Float64(get(kwargs, :residual_norm_max, 1e-6))
    selected = execute_governance_selector(raw_candidates;
        selector=selector_fn,
        residual_norm_max=residual_norm_max,
        require_converged=true,
    )
    return selected, selected.normalized_candidates
end

const _FIXEDRHO_JOINT_REQUIRED_SCALARS = (
    :omega,
    :pressure,
    :rho_norm,
    :entropy,
    :energy,
    :residual_norm,
)
const _FIXEDRHO_JOINT_REQUIRED_VECTORS = (
    (:solution, 8),
    (:x_state, 5),
    (:mu_vec, 3),
    (:masses, 3),
)

"""Return explicit reasons for a malformed FixedRho joint-solve payload.

The solver callback is allowed to report a failed attempt, but it must not
leak a partial payload into the outer orchestration layer.  In particular,
`Float64(nothing)` is an orchestration error, not a physical solver result.
"""
function _fixedrho_joint_output_issues(solved)
    solved === nothing && return Symbol[:solver_returned_nothing]
    issues = Symbol[]
    if !hasproperty(solved, :converged) || !(getproperty(solved, :converged) isa Bool)
        push!(issues, :invalid_converged)
    end
    if !hasproperty(solved, :iterations) || !(getproperty(solved, :iterations) isa Integer)
        push!(issues, :invalid_iterations)
    end
    for field in _FIXEDRHO_JOINT_REQUIRED_SCALARS
        if !hasproperty(solved, field)
            push!(issues, Symbol("missing_$(field)"))
        elseif getproperty(solved, field) === nothing || !(getproperty(solved, field) isa Real)
            push!(issues, Symbol("invalid_$(field)"))
        end
    end
    for (field, expected_length) in _FIXEDRHO_JOINT_REQUIRED_VECTORS
        if !hasproperty(solved, field)
            push!(issues, Symbol("missing_$(field)"))
            continue
        end
        value = getproperty(solved, field)
        if value === nothing || !(value isa AbstractVector) || length(value) != expected_length ||
                any(entry -> !(entry isa Real), value)
            push!(issues, Symbol("invalid_$(field)"))
        end
    end
    return issues
end

function _fixedrho_incomplete_output_candidate(
    attempt_cfg,
    residual_norm_max::Real,
    issues::AbstractVector{<:Symbol},
)
    issue_text = isempty(issues) ? "unknown incomplete solver output" : join(string.(issues), ", ")
    return (
        converged=false,
        solution=Float64[],
        x_state=zeros(5),
        mu_vec=zeros(3),
        omega=NaN,
        pressure=NaN,
        rho_norm=NaN,
        entropy=NaN,
        energy=NaN,
        masses=zeros(3),
        iterations=0,
        residual_norm=Inf,
        residual_norm_max=Float64(residual_norm_max),
        fixedrho_joint_solve_requested=true,
        fixedrho_joint_solve_active=false,
        fixedrho_joint_fallback=false,
        fixedrho_joint_selected_method=:none,
        fixedrho_joint_selected_quality=:bad,
        fixedrho_joint_fallback_used=false,
        fixedrho_joint_attempt_origin=attempt_cfg.attempt_origin,
        failed_constraints=vcat(Symbol[:solver_incomplete_output], collect(issues)),
        error_kind=:solver_incomplete_output,
        error_msg="FixedRho joint solver returned incomplete output ($(issue_text))",
    )
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
    for key in (:seed_candidates, :hard_constraints, :semantic_mode, :selector, :fixedrho_joint_solve, :continuity_seed, :evaluate_all_attempts, :extra_constraints, :diagnostic_level)
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
    delete!(local_kwargs, :work_telemetry)

    solved = _solve_constraint_fixedmu(model, T_fm, μ_fm; pairs(local_kwargs)...)
    result = solved
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
        :thermo_quadrature_policy,
        :thermo_quadrature_rtol,
        :thermo_quadrature_atol,
        :thermo_quadrature_maxevals,
        :work_telemetry,
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
    thermo_quadrature_policy = Symbol(get(kwargs, :thermo_quadrature_policy, :tensor_gauss))
    thermo_quadrature_rtol = Float64(get(kwargs, :thermo_quadrature_rtol, 1e-8))
    thermo_quadrature_atol = Float64(get(kwargs, :thermo_quadrature_atol, 1e-10))
    thermo_quadrature_maxevals = Int(get(kwargs, :thermo_quadrature_maxevals, 10^7))
    residual_norm_max = Float64(get(kwargs, :residual_norm_max, 1e-6))
    nlsolve_method = get(kwargs, :nlsolve_method, :trust_region)
    trust_region_fallback = Bool(get(kwargs, :trust_region_fallback, true))
    fallback_method = get(kwargs, :fallback_method, :trust_region)
    physicality_check = get(kwargs, :physicality_check, ((_, _) -> true))
    work_telemetry = get(kwargs, :work_telemetry, nothing)

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
        thermo_quadrature_policy=thermo_quadrature_policy,
        thermo_quadrature_rtol=thermo_quadrature_rtol,
        thermo_quadrature_atol=thermo_quadrature_atol,
        thermo_quadrature_maxevals=thermo_quadrature_maxevals,
    )
    residual_fn! = build_residual!(mode, params)

    local_nls_kwargs = Dict{Symbol,Any}(kwargs)
    _strip_problemspec_forwardsolve_kwargs!(local_nls_kwargs)
    _strip_fixedrho_joint_nlsolve_kwargs!(local_nls_kwargs)

    postprocess_solution = function (solution)
        thermo = compute_thermo_from_solution(
            model,
            solution,
            T_fm;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            rho0_scale=rho0,
            state_n=5,
            mu_n=3,
            thermo_quadrature_policy=thermo_quadrature_policy,
            thermo_quadrature_rtol=thermo_quadrature_rtol,
            thermo_quadrature_atol=thermo_quadrature_atol,
            thermo_quadrature_maxevals=thermo_quadrature_maxevals,
        )
        residual_norm = compute_residual_norm_from_solution(
            model,
            solution,
            T_fm;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            state_n=5,
            mu_n=3,
            residual_fn=residual_fn!,
            thermo_quadrature_policy=thermo_quadrature_policy,
            thermo_quadrature_rtol=thermo_quadrature_rtol,
            thermo_quadrature_atol=thermo_quadrature_atol,
            thermo_quadrature_maxevals=thermo_quadrature_maxevals,
        )
        record_postprocess_residual!(work_telemetry)
        phys = physicality_check(thermo.x_state, thermo.masses) && _thermo_quantities_finite(thermo)

        return (
            x_state=thermo.x_state,
            mu_vec=thermo.mu_vec,
            omega=thermo.omega,
            pressure=thermo.pressure,
            rho_norm=thermo.rho_norm,
            entropy=thermo.entropy,
            energy=thermo.energy,
            masses=thermo.masses,
            residual_norm=residual_norm,
            phys=phys,
        )
    end

    cache = Dict{Symbol,NamedTuple}()
    solve_once = function (method::Symbol, seed::Vector{Float64})
        local res
        try
            res = nlsolve(
                residual_fn!,
                seed;
                autodiff=:forward,
                method=method,
                xtol=1e-9,
                ftol=1e-9,
                pairs(local_nls_kwargs)...,
            )
        catch
            record_solver_exception!(work_telemetry)
            rethrow()
        end
        record_nlsolve_work!(work_telemetry, res, method)

        solution = Float64.(res.zero)
        pp = postprocess_solution(solution)
        converged = Bool(res.f_converged) && pp.phys && isfinite(pp.residual_norm) && pp.residual_norm <= residual_norm_max
        record_attempt_outcome!(work_telemetry; converged=converged)

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
    record_solver_request!(get(kwargs, :work_telemetry, nothing); fixedrho=true)
    cfg = _fixedrho_runtime_config_from_kwargs(mode, kwargs)
    diagnostic_level = _resolve_diagnostic_level(kwargs)
    extra_constraints = _resolve_extra_constraints(kwargs)
    fixedrho_joint_solve = cfg.fixedrho_joint_solve
    fixedrho_joint_solve || throw(ArgumentError("fixedrho_joint_solve=false is no longer supported; FixedRho uses joint solve only"))

    provided_seed_pool = cfg.seed_candidates
    default_seed_pool = _build_default_seed_candidates(mode, cfg.seed_guess)

    primary_seed = cfg.seed_guess
    primary_method = cfg.primary_method
    primary_use_fallback = cfg.trust_region_fallback
    fallback_method = cfg.fallback_method

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
    kwargs[:evaluate_all_attempts] = cfg.evaluate_all_attempts
    kwargs[:residual_norm_max] = cfg.residual_norm_max
    kwargs[:xi] = cfg.xi
    kwargs[:p_num] = cfg.p_num
    kwargs[:t_num] = cfg.t_num
    kwargs[:continuity_seed] = cfg.continuity_seed
    selected, candidates = _execute_governed_attempt_plan(
        attempt_plan,
        kwargs,
        hard_constraints,
        selector_fn,
        function (local_kwargs, attempt_cfg, _)
            _strip_problemspec_forwardsolve_kwargs!(local_kwargs)
            solved = _fixedrho_joint_problem_spec_forward_solve(model, mode, T_fm; pairs(local_kwargs)...)
            output_issues = _fixedrho_joint_output_issues(solved)
            if !isempty(output_issues)
                return _fixedrho_incomplete_output_candidate(
                    attempt_cfg,
                    get(local_kwargs, :residual_norm_max, 1e-6),
                    output_issues,
                )
            end
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
        selection_score=(hasproperty(s, :selection_score) ? Float64(getproperty(s, :selection_score)) : NaN),
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
    runtime_config=nothing,
)
    extra_constraints = _resolve_extra_constraints(kwargs)
    diagnostic_level = _resolve_diagnostic_level(kwargs)
    seed_guess = runtime_config === nothing ? get(kwargs, :seed_guess, nothing) : runtime_config.seed_guess
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec $(mode_label) forward_solve"))

    provided_seed_pool = runtime_config === nothing ? _resolve_optional_seed_candidates(kwargs) : runtime_config.seed_candidates
    default_seed_pool = _build_default_seed_candidates(mode, seed_guess)

    primary_seed = Float64.(seed_guess)
    primary_method = if runtime_config !== nothing
        runtime_config.primary_method
    elseif haskey(kwargs, :nlsolve_method)
        kwargs[:nlsolve_method]
    else
        Bool(get(kwargs, :continuity_seed, false)) ? :newton : :trust_region
    end
    primary_use_fallback = runtime_config === nothing ? Bool(get(kwargs, :trust_region_fallback, true)) : runtime_config.trust_region_fallback
    fallback_method = runtime_config === nothing ? get(kwargs, :fallback_method, :trust_region) : runtime_config.fallback_method

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
    if runtime_config !== nothing
        kwargs[:evaluate_all_attempts] = runtime_config.evaluate_all_attempts
        kwargs[:residual_norm_max] = runtime_config.residual_norm_max
        kwargs[:xi] = runtime_config.xi
        kwargs[:p_num] = runtime_config.p_num
        kwargs[:t_num] = runtime_config.t_num
        kwargs[:continuity_seed] = runtime_config.continuity_seed
    end
    selected, candidates = _execute_governed_attempt_plan(
        attempt_plan,
        kwargs,
        hard_constraints,
        selector_fn,
        function (local_kwargs, attempt_cfg, _)
            _strip_problemspec_forwardsolve_kwargs!(local_kwargs)
            delete!(local_kwargs, :work_telemetry)
            delete!(local_kwargs, :trust_region_fallback)
            delete!(local_kwargs, :fallback_method)

            haskey(local_kwargs, :rho0) || throw(ArgumentError("rho0 is required for ProblemSpec $(mode_label) forward_solve"))
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
        selection_score=(hasproperty(s, :selection_score) ? Float64(getproperty(s, :selection_score)) : NaN),
        hard_constraint_ok=s.hard_constraint_ok,
        failed_constraints=s.failed_constraints,
        selection_reason=selected.selection_reason,
        selected_index=selected.selected_index,
        candidate_count=length(candidates),
        error_kind=(hasproperty(s, :error_kind) ? Symbol(getproperty(s, :error_kind)) : :none),
        error_msg=(hasproperty(s, :error_msg) ? String(getproperty(s, :error_msg)) : ""),
        governed_selected_method=(hasproperty(s, :governed_selected_method) ? getproperty(s, :governed_selected_method) : :none),
        governed_selected_quality=(hasproperty(s, :quality_tag) ? Symbol(getproperty(s, :quality_tag)) :
            (hasproperty(s, :governed_selected_quality) ? getproperty(s, :governed_selected_quality) : :bad)),
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
    cfg = _fixedentropy_runtime_config_from_kwargs(mode, kwargs)
    kwargs[:rho0] = cfg.rho0
    kwargs[:mu0] = cfg.mu0
    return _governed_nonrho_problem_spec_forward_solve(
        model,
        mode,
        T_fm,
        kwargs,
        "FixedEntropy",
        :entropy_attempt_origin,
        local_kwargs -> _solve_constraint_fixedentropy(model, T_fm, mode.s_target; pairs(local_kwargs)...),
        cfg,
    )
end

function _fixedsigma_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedSigma, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    cfg = _fixedsigma_runtime_config_from_kwargs(mode, kwargs)
    kwargs[:rho0] = cfg.rho0
    kwargs[:mu0] = cfg.mu0
    return _governed_nonrho_problem_spec_forward_solve(
        model,
        mode,
        T_fm,
        kwargs,
        "FixedSigma",
        :sigma_attempt_origin,
        local_kwargs -> _solve_constraint_fixedsigma(model, T_fm, mode.sigma_target; pairs(local_kwargs)...),
        cfg,
    )
end

function _fixedasymrho_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedAsymmetricRho, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    cfg = _fixedasymrho_runtime_config_from_kwargs(mode, kwargs)
    kwargs[:rho0] = cfg.rho0
    kwargs[:mu0] = cfg.mu0
    return _governed_nonrho_problem_spec_forward_solve(
        model,
        mode,
        T_fm,
        kwargs,
        "FixedAsymmetricRho",
        :asym_attempt_origin,
        local_kwargs -> _solve_constraint_fixedasymrho(model, T_fm, mode.rho_target, mode.ud_ratio_target, mode.s_target; pairs(local_kwargs)...),
        cfg,
    )
end

function _fixedmub_conserved_problem_spec_forward_solve(
    model::AbstractQCDModel,
    mode::FixedMuBConservedCharges,
    T_fm::Real;
    fwd_kwargs...,
)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))
    kwargs[:rho0] = Float64(get(kwargs, :rho0, rho0))
    isfinite(kwargs[:rho0]) && kwargs[:rho0] > 0 || throw(ArgumentError(
        "rho0 must be finite and positive, got $(kwargs[:rho0])",
    ))
    if haskey(kwargs, :nlsolve_method)
        method = kwargs[:nlsolve_method]
        method in (:newton, :trust_region) || throw(ArgumentError(
            "nlsolve_method must be :newton or :trust_region, got $(method)",
        ))
    end
    if !haskey(kwargs, :mu0)
        flavor = flavor_mu_from_bqs(mode.muB_fm, 0.0, mode.muB_fm / 3)
        kwargs[:mu0] = Float64[flavor.mu_u, flavor.mu_d, flavor.mu_s]
    end
    return _governed_nonrho_problem_spec_forward_solve(
        model,
        mode,
        T_fm,
        kwargs,
        "FixedMuBConservedCharges",
        :conserved_charge_attempt_origin,
        local_kwargs -> _solve_constraint_fixedmub_conserved(
            model,
            T_fm,
            mode;
            pairs(local_kwargs)...,
        ),
    )
end
