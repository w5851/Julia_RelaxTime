"""
    solve_constraint(model, mode, T_fm; kwargs...)

统一约束求解入口。根据 `mode` 的类型分发到 models 域约束求解内核实现。
"""
struct SolverResult{VX<:AbstractVector{<:Real}, VM<:AbstractVector{<:Real}, MM<:AbstractVector{<:Real}}
    mode::ConstraintMode
    converged::Bool
    solution::Vector{Float64}
    x_state::VX
    mu_vec::VM
    omega::Float64
    pressure::Float64
    rho_norm::Float64
    entropy::Float64
    energy::Float64
    masses::MM
    iterations::Int
    residual_norm::Float64
    xi::Float64
end

@inline function _coerce_solver_result(mode::ConstraintMode, raw_result; xi_override=nothing)
    xi_val = xi_override === nothing ? Float64(getproperty(raw_result, :xi)) : Float64(xi_override)
    return SolverResult(
        mode,
        Bool(getproperty(raw_result, :converged)),
        Float64.(getproperty(raw_result, :solution)),
        getproperty(raw_result, :x_state),
        getproperty(raw_result, :mu_vec),
        Float64(getproperty(raw_result, :omega)),
        Float64(getproperty(raw_result, :pressure)),
        Float64(getproperty(raw_result, :rho_norm)),
        Float64(getproperty(raw_result, :entropy)),
        Float64(getproperty(raw_result, :energy)),
        getproperty(raw_result, :masses),
        Int(getproperty(raw_result, :iterations)),
        Float64(getproperty(raw_result, :residual_norm)),
        xi_val,
    )
end

@inline function _strip_forward_kwargs(kwargs, blocked::Tuple)
    return (; (k => v for (k, v) in pairs(kwargs) if !(k in blocked))...)
end

@inline function _resolve_solver_model(model::AbstractPNJLModel, kwargs)
    if haskey(kwargs, :model_kind)
        model_kind = kwargs[:model_kind]
        model_kind isa Symbol || throw(ArgumentError("model_kind must be Symbol, got $(typeof(model_kind))"))
        return create_model(model_kind)
    end
    return model
end

@inline function _resolve_nonfixedmu_bridge(mode::ConstraintMode, T_fm::Real, kwargs)
    xi = get(kwargs, :xi, 0.0)
    p_num = get(kwargs, :p_num, default_momentum_count())
    t_num = get(kwargs, :t_num, default_theta_count())
    residual_norm_max = get(kwargs, :residual_norm_max, 1e-6)
    rho0 = get(kwargs, :rho0, Main.Constants_PNJL.ρ0_inv_fm3)

    semantic_mode = get(kwargs, :semantic_mode, :ground_state)
    (semantic_mode === :ground_state || semantic_mode === :constrained_manifold) ||
        throw(ArgumentError("semantic_mode must be :ground_state or :constrained_manifold, got $(semantic_mode)"))

    selector = get(kwargs, :selector, nothing)
    selector === nothing || selector isa Function ||
        throw(ArgumentError("selector must be Function or nothing, got $(typeof(selector))"))

    seed_strategy = get(kwargs, :seed_strategy, DefaultSeed())
    seed_guess = haskey(kwargs, :seed_guess) ? kwargs[:seed_guess] : get_seed(seed_strategy, [T_fm], mode)
    seed_candidates = if haskey(kwargs, :seed_candidates)
        kwargs[:seed_candidates]
    else
        (seed_guess,)
    end

    use_problem_spec_chain =
        haskey(kwargs, :problem_spec) ||
        haskey(kwargs, :seed_guess) ||
        haskey(kwargs, :seed_candidates) ||
        semantic_mode !== :ground_state ||
        selector !== nothing

    return (
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        residual_norm_max=residual_norm_max,
        rho0=rho0,
        semantic_mode=semantic_mode,
        selector=selector,
        seed_guess=seed_guess,
        seed_candidates=seed_candidates,
        use_problem_spec_chain=use_problem_spec_chain,
    )
end

@inline function _solve_with_problem_spec_default(
    model::AbstractQCDModel,
    mode::ConstraintMode,
    T_fm::Real,
    kwargs,
)
    haskey(kwargs, :use_problem_spec) && throw(ArgumentError("use_problem_spec has been removed; solve_constraint always uses ProblemSpec chain"))
    haskey(kwargs, :allow_legacy_path) && throw(ArgumentError("allow_legacy_path has been removed together with legacy fallback path"))
    haskey(kwargs, :warn_on_legacy_path) && throw(ArgumentError("warn_on_legacy_path has been removed together with legacy fallback path"))

    spec = get(kwargs, :problem_spec, nothing)
    forwarded = (; (k => v for (k, v) in pairs(kwargs) if k != :problem_spec)...)
    if spec === nothing
        spec = build_problem_spec(mode)
    else
        spec isa ProblemSpec || throw(ArgumentError("problem_spec must be ProblemSpec or nothing, got $(typeof(spec))"))
    end
    return spec.forward_solve(model, T_fm; forwarded...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedMu, T_fm::Real; μ_fm::Real, kwargs...)
    return _solve_constraint_fixedmu(model, T_fm, μ_fm; kwargs...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedRho, T_fm::Real;
    problem_spec::Union{Nothing, ProblemSpec}=nothing,
    kwargs...)
    merged = (; kwargs..., problem_spec=problem_spec)
    return _solve_with_problem_spec_default(model, mode, T_fm, merged)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedEntropy, T_fm::Real; kwargs...)
    return _solve_with_problem_spec_default(model, mode, T_fm, kwargs)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedSigma, T_fm::Real; kwargs...)
    return _solve_with_problem_spec_default(model, mode, T_fm, kwargs)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedAsymmetricRho, T_fm::Real; kwargs...)
    return _solve_with_problem_spec_default(model, mode, T_fm, kwargs)
end

function solve(model::AbstractPNJLModel, mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    xi = get(kwargs, :xi, 0.0)
    p_num = get(kwargs, :p_num, default_momentum_count())
    t_num = get(kwargs, :t_num, default_theta_count())
    residual_norm_max = get(kwargs, :residual_norm_max, 1e-6)
    auto_multiseed_fallback = get(kwargs, :auto_multiseed_fallback, true)
    seed_strategy = get(kwargs, :seed_strategy, DefaultSeed())

    if seed_strategy isa MultiSeed
        return solve_multi(model, mode, T_fm, μ_fm; kwargs...)
    end

    seed_guess = haskey(kwargs, :seed_guess) ? kwargs[:seed_guess] : get_seed(seed_strategy, [T_fm, μ_fm], mode)
    forwarded = _strip_forward_kwargs(kwargs, (
        :seed_strategy,
        :seed_guess,
        :xi,
        :p_num,
        :t_num,
        :residual_norm_max,
        :auto_multiseed_fallback,
    ))
    raw = solve_constraint(
        model,
        mode,
        T_fm;
        μ_fm=μ_fm,
        seed_guess=seed_guess,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        residual_norm_max=residual_norm_max,
        forwarded...,
    )

    single = SolverResult(
        mode,
        Bool(raw.converged),
        Float64.(raw.solution),
        raw.x_state,
        raw.mu_vec,
        Float64(raw.omega),
        Float64(raw.pressure),
        Float64(raw.rho_norm) / Float64(rho0),
        Float64(raw.entropy),
        Float64(raw.energy),
        raw.masses,
        Int(raw.iterations),
        Float64(raw.residual_norm),
        Float64(xi),
    )

    if single.converged || !Bool(auto_multiseed_fallback)
        return single
    end

    try
        return solve_multi(
            model,
            mode,
            T_fm,
            μ_fm;
            seed_strategy=MultiSeed(),
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            residual_norm_max=residual_norm_max,
            forwarded...,
        )
    catch
        return single
    end
end

function solve(model::AbstractPNJLModel, mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, T_fm::Real; kwargs...)
    effective_model = _resolve_solver_model(model, kwargs)
    bridge = _resolve_nonfixedmu_bridge(mode, T_fm, kwargs)

    if bridge.use_problem_spec_chain || mode isa FixedRho || mode isa FixedEntropy || mode isa FixedSigma || mode isa FixedAsymmetricRho
        problem_spec = get(kwargs, :problem_spec, nothing)
        rho0_kwargs = (mode isa FixedRho) ? NamedTuple() : (; rho0=bridge.rho0)
        forwarded = _strip_forward_kwargs(kwargs, (
            :problem_spec,
            :seed_strategy,
            :seed_guess,
            :seed_candidates,
            :semantic_mode,
            :selector,
            :rho0,
            :xi,
            :p_num,
            :t_num,
            :residual_norm_max,
            :model_kind,
        ))
        raw = solve_constraint(
            effective_model,
            mode,
            T_fm;
            problem_spec=problem_spec,
            seed_guess=bridge.seed_guess,
            seed_candidates=bridge.seed_candidates,
            semantic_mode=bridge.semantic_mode,
            selector=bridge.selector,
            rho0_kwargs...,
            xi=bridge.xi,
            p_num=bridge.p_num,
            t_num=bridge.t_num,
            residual_norm_max=bridge.residual_norm_max,
            forwarded...,
        )
        return SolverResult(
            mode,
            Bool(raw.converged),
            Float64.(raw.solution),
            raw.x_state,
            raw.mu_vec,
            Float64(raw.omega),
            Float64(raw.pressure),
            Float64(raw.rho_norm),
            Float64(raw.entropy),
            Float64(raw.energy),
            raw.masses,
            Int(raw.iterations),
            Float64(raw.residual_norm),
            Float64(bridge.xi),
        )
    end

    throw(ArgumentError("unsupported mode for solve bridge: $(typeof(mode))"))
end

function solve_multi(model::AbstractPNJLModel, mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    xi = get(kwargs, :xi, 0.0)
    p_num = get(kwargs, :p_num, default_momentum_count())
    t_num = get(kwargs, :t_num, default_theta_count())
    residual_norm_max = get(kwargs, :residual_norm_max, 1e-6)
    seed_strategy = get(kwargs, :seed_strategy, MultiSeed())

    seeds = if haskey(kwargs, :seeds)
        [Float64.(s) for s in kwargs[:seeds]]
    else
        get_all_seeds(seed_strategy, [T_fm, μ_fm], mode)
    end
    isempty(seeds) && (seeds = [get_seed(DefaultSeed(), [T_fm, μ_fm], mode)])

    forwarded = _strip_forward_kwargs(kwargs, (
        :seed_strategy,
        :seeds,
        :xi,
        :p_num,
        :t_num,
        :residual_norm_max,
    ))
    candidates = NamedTuple[]
    for (seed_index, seed) in enumerate(seeds)
        local raw
        try
            raw = solve_constraint(
                model,
                mode,
                T_fm;
                μ_fm=μ_fm,
                seed_guess=seed,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
                residual_norm_max=residual_norm_max,
                forwarded...,
            )
            ok = Bool(raw.converged) && isfinite(raw.residual_norm) && raw.residual_norm <= residual_norm_max
            candidate = (
                converged=Bool(raw.converged),
                solution=Float64.(raw.solution),
                x_state=raw.x_state,
                mu_vec=raw.mu_vec,
                omega=Float64(raw.omega),
                pressure=Float64(raw.pressure),
                rho_norm=Float64(raw.rho_norm) / Float64(rho0),
                entropy=Float64(raw.entropy),
                energy=Float64(raw.energy),
                masses=raw.masses,
                iterations=Int(raw.iterations),
                residual_norm=Float64(raw.residual_norm),
                hard_constraint_ok=ok,
                failed_constraints=(ok ? Symbol[] : Symbol[:residual_too_large]),
                seed_index=Int(seed_index),
            )
            push!(candidates, candidate)
        catch
            push!(candidates, (
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
                hard_constraint_ok=false,
                failed_constraints=Symbol[:solver_failed],
                seed_index=Int(seed_index),
            ))
        end
    end

    selected = select_pressure_max_candidate(candidates)
    s = selected.selected_candidate
    return SolverResult(
        mode,
        Bool(s.converged),
        Float64.(s.solution),
        s.x_state,
        s.mu_vec,
        Float64(s.omega),
        Float64(s.pressure),
        Float64(s.rho_norm),
        Float64(s.entropy),
        Float64(s.energy),
        s.masses,
        Int(s.iterations),
        Float64(s.residual_norm),
        Float64(xi),
    )
end

function solve_multi(model::AbstractPNJLModel, mode::Union{FixedRho, FixedAsymmetricRho}, T_fm::Real; kwargs...)
    effective_model = _resolve_solver_model(model, kwargs)
    bridge = _resolve_nonfixedmu_bridge(mode, T_fm, kwargs)

    if bridge.use_problem_spec_chain
        seed_strategy = get(kwargs, :seed_strategy, MultiSeed())
        seeds = if haskey(kwargs, :seeds)
            [Float64.(s) for s in kwargs[:seeds]]
        elseif haskey(kwargs, :seed_candidates)
            [Float64.(s) for s in kwargs[:seed_candidates]]
        else
            get_all_seeds(seed_strategy, [T_fm], mode)
        end
        isempty(seeds) && (seeds = [Float64.(bridge.seed_guess)])

        forwarded = _strip_forward_kwargs(kwargs, (
            :seed_strategy,
            :seeds,
            :seed_candidates,
            :seed_guess,
            :semantic_mode,
            :selector,
            :rho0,
            :xi,
            :p_num,
            :t_num,
            :residual_norm_max,
            :model_kind,
            :problem_spec,
        ))

        candidates = NamedTuple[]
        for (seed_index, seed) in enumerate(seeds)
            local raw
            try
                rho0_kwargs = (mode isa FixedRho) ? NamedTuple() : (; rho0=bridge.rho0)
                raw = solve_constraint(
                    effective_model,
                    mode,
                    T_fm;
                    problem_spec=get(kwargs, :problem_spec, nothing),
                    seed_guess=seed,
                    seed_candidates=(seed,),
                    semantic_mode=bridge.semantic_mode,
                    selector=bridge.selector,
                    rho0_kwargs...,
                    xi=bridge.xi,
                    p_num=bridge.p_num,
                    t_num=bridge.t_num,
                    residual_norm_max=bridge.residual_norm_max,
                    forwarded...,
                )
                ok = Bool(raw.converged) && isfinite(raw.residual_norm) && raw.residual_norm <= max(bridge.residual_norm_max, 1e-3)
                push!(candidates, (
                    converged=Bool(raw.converged),
                    solution=Float64.(raw.solution),
                    x_state=raw.x_state,
                    mu_vec=raw.mu_vec,
                    omega=Float64(raw.omega),
                    pressure=Float64(raw.pressure),
                    rho_norm=Float64(raw.rho_norm),
                    entropy=Float64(raw.entropy),
                    energy=Float64(raw.energy),
                    masses=raw.masses,
                    iterations=Int(raw.iterations),
                    residual_norm=Float64(raw.residual_norm),
                    hard_constraint_ok=ok,
                    failed_constraints=(ok ? Symbol[] : Symbol[:residual_too_large]),
                    seed_index=Int(seed_index),
                ))
            catch
                push!(candidates, (
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
                    hard_constraint_ok=false,
                    failed_constraints=Symbol[:solver_failed],
                    seed_index=Int(seed_index),
                ))
            end
        end

        selected = select_pressure_max_candidate(candidates)
        s = selected.selected_candidate
        return SolverResult(
            mode,
            Bool(s.converged),
            Float64.(s.solution),
            s.x_state,
            s.mu_vec,
            Float64(s.omega),
            Float64(s.pressure),
            Float64(s.rho_norm),
            Float64(s.entropy),
            Float64(s.energy),
            s.masses,
            Int(s.iterations),
            Float64(s.residual_norm),
            Float64(bridge.xi),
        )
    end

    if mode isa FixedRho || mode isa FixedAsymmetricRho
        seed_strategy = get(kwargs, :seed_strategy, MultiSeed())
        seeds = if haskey(kwargs, :seeds)
            [Float64.(s) for s in kwargs[:seeds]]
        else
            get_all_seeds(seed_strategy, [T_fm], mode)
        end
        isempty(seeds) && (seeds = [Float64.(bridge.seed_guess)])

        selector_fn = haskey(kwargs, :selector) ? kwargs[:selector] : SeedStrategies.default_omega_selector
        forwarded = _strip_forward_kwargs(kwargs, (
            :seed_strategy,
            :seeds,
            :seed_guess,
            :seed_candidates,
            :semantic_mode,
            :selector,
        ))

        results = SolverResult[]
        for seed in seeds
            try
                push!(results, solve(effective_model, mode, T_fm;
                    seed_strategy=DefaultSeed(Float64.(seed), Float64.(seed), :hadron),
                    forwarded...,
                ))
            catch
            end
        end

        isempty(results) && error("All seeds failed (exceptions) in solve_multi")
        converged = filter(r -> r.converged, results)
        isempty(converged) && error("All seeds failed to converge to a physical solution")
        return selector_fn(converged)
    end

    throw(ArgumentError("unsupported mode for solve_multi bridge: $(typeof(mode))"))
end

function solve(mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    model = create_model(:PNJL)
    return solve(model, mode, T_fm, μ_fm; kwargs...)
end

function solve(mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, T_fm::Real; kwargs...)
    model = create_model(:PNJL)
    return solve(model, mode, T_fm; kwargs...)
end

function solve_multi(mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    model = create_model(:PNJL)
    if haskey(kwargs, :seeds)
        seeds = kwargs[:seeds]
        candidates = [DefaultSeed(Float64.(seed), Float64.(seed), :hadron) for seed in seeds]
        selector = haskey(kwargs, :selector) ? kwargs[:selector] : SeedStrategies.default_omega_selector
        forward_kwargs = (; (k => v for (k, v) in pairs(kwargs) if k != :seeds && k != :selector)...)
        return solve_multi(model, mode, T_fm, μ_fm;
            seed_strategy=MultiSeed(candidates, selector),
            forward_kwargs...)
    end
    return solve_multi(model, mode, T_fm, μ_fm; kwargs...)
end

function solve_multi(mode::Union{FixedRho, FixedAsymmetricRho}, T_fm::Real; kwargs...)
    model = create_model(:PNJL)
    return solve_multi(model, mode, T_fm; kwargs...)
end

@inline function is_physical_solution(x_state::AbstractVector{<:Real}, masses::AbstractVector{<:Real}; phi_tol::Float64=1e-8)
    if length(x_state) < 5 || length(masses) < 3
        return false
    end
    Φ = x_state[4]
    Φbar = x_state[5]
    if !(isfinite(Φ) && isfinite(Φbar) && (-phi_tol <= Φ <= 1 + phi_tol) && (-phi_tol <= Φbar <= 1 + phi_tol))
        return false
    end
    if any(!isfinite, masses) || any(m -> m <= 0.0, masses)
        return false
    end
    return true
end
@inline solve_with_derivatives(T_fm::Real, μ_fm::Real; kwargs...) =
    solve_pnjl_with_derivatives(T_fm, μ_fm; kwargs...)

@inline function solve_with_derivatives(model::AbstractPNJLModel, mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    _ = model, mode
    return solve_pnjl_with_derivatives(T_fm, μ_fm; kwargs...)
end
