"""
    solve_constraint(model, mode, T_fm; kwargs...)

统一约束求解入口，固定单主链：
`solve* -> solve_constraint -> ProblemSpec.forward_solve -> mode executor`。
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
    contract_version::Symbol
end

const SOLVER_CONTRACT_VERSION_V1 = :v1
const SOLVER_RESULT_REQUIRED_FIELDS = (
    :contract_version,
    :mode,
    :converged,
    :solution,
    :x_state,
    :mu_vec,
    :omega,
    :pressure,
    :rho_norm,
    :entropy,
    :energy,
    :masses,
    :iterations,
    :residual_norm,
    :xi,
)

@inline function _normalize_solver_contract_version(version)::Symbol
    if version isa Symbol
        version == SOLVER_CONTRACT_VERSION_V1 || throw(ArgumentError("unsupported solver contract_version: $(version), expected $(SOLVER_CONTRACT_VERSION_V1)"))
        return version
    end
    version isa AbstractString || throw(ArgumentError("contract_version must be Symbol or AbstractString, got $(typeof(version))"))
    normalized = Symbol(version)
    normalized == SOLVER_CONTRACT_VERSION_V1 || throw(ArgumentError("unsupported solver contract_version: $(normalized), expected $(SOLVER_CONTRACT_VERSION_V1)"))
    return normalized
end

function SolverResult(
    mode::ConstraintMode,
    converged::Bool,
    solution::Vector{Float64},
    x_state::VX,
    mu_vec::VM,
    omega::Float64,
    pressure::Float64,
    rho_norm::Float64,
    entropy::Float64,
    energy::Float64,
    masses::MM,
    iterations::Int,
    residual_norm::Float64,
    xi::Float64,
) where {VX<:AbstractVector{<:Real}, VM<:AbstractVector{<:Real}, MM<:AbstractVector{<:Real}}
    return SolverResult(
        mode,
        converged,
        solution,
        x_state,
        mu_vec,
        omega,
        pressure,
        rho_norm,
        entropy,
        energy,
        masses,
        iterations,
        residual_norm,
        xi,
        SOLVER_CONTRACT_VERSION_V1,
    )
end

function SolverResult(
    mode::ConstraintMode,
    converged::Bool,
    solution::Vector{Float64},
    x_state::VX,
    mu_vec::VM,
    omega::Float64,
    pressure::Float64,
    rho_norm::Float64,
    entropy::Float64,
    energy::Float64,
    masses::MM,
    iterations::Int,
    residual_norm::Float64,
    xi::Float64,
    contract_version,
) where {VX<:AbstractVector{<:Real}, VM<:AbstractVector{<:Real}, MM<:AbstractVector{<:Real}}
    return SolverResult{VX, VM, MM}(
        mode,
        converged,
        solution,
        x_state,
        mu_vec,
        omega,
        pressure,
        rho_norm,
        entropy,
        energy,
        masses,
        iterations,
        residual_norm,
        xi,
        _normalize_solver_contract_version(contract_version),
    )
end

@inline solver_contract_version(result::SolverResult) = result.contract_version

@inline function to_namedtuple(result::SolverResult)
    return (
        contract_version=result.contract_version,
        mode=result.mode,
        converged=result.converged,
        solution=result.solution,
        x_state=result.x_state,
        mu_vec=result.mu_vec,
        omega=result.omega,
        pressure=result.pressure,
        rho_norm=result.rho_norm,
        entropy=result.entropy,
        energy=result.energy,
        masses=result.masses,
        iterations=result.iterations,
        residual_norm=result.residual_norm,
        xi=result.xi,
    )
end

@inline function _coerce_solver_result(mode::ConstraintMode, raw_result; xi_override=nothing)
    xi_val = if xi_override === nothing
        hasproperty(raw_result, :xi) ? Float64(getproperty(raw_result, :xi)) : 0.0
    else
        Float64(xi_override)
    end
    contract_version = hasproperty(raw_result, :contract_version) ?
        _normalize_solver_contract_version(getproperty(raw_result, :contract_version)) :
        SOLVER_CONTRACT_VERSION_V1
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
        contract_version,
    )
end

@inline coerce_solver_result(mode::ConstraintMode, raw_result; xi_override=nothing) =
    _coerce_solver_result(mode, raw_result; xi_override=xi_override)

@inline solver_result_view(result::SolverResult) = to_namedtuple(result)

@inline function solver_result_is_success(result::SolverResult; residual_norm_max::Real=1e-6, require_converged::Bool=true)
    require_converged && !result.converged && return false
    return isfinite(result.residual_norm) && result.residual_norm <= Float64(residual_norm_max)
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

@inline function _normalize_seed_vector(seed, key::Symbol)::Vector{Float64}
    seed isa AbstractVector || throw(ArgumentError("$(key) must be an AbstractVector, got $(typeof(seed))"))
    return Float64.(seed)
end

@inline function _normalize_seed_pool(seed_input, key::Symbol)::Vector{Vector{Float64}}
    raw = try
        collect(seed_input)
    catch
        throw(ArgumentError("$(key) must be an iterable of AbstractVector, got $(typeof(seed_input))"))
    end
    seeds = Vector{Vector{Float64}}(undef, length(raw))
    for i in eachindex(raw)
        seeds[i] = _normalize_seed_vector(raw[i], key)
    end
    return seeds
end

@inline function _resolve_common_runtime_options(kwargs)
    xi_raw = get(kwargs, :xi, 0.0)
    xi_raw isa Real || throw(ArgumentError("xi must be Real, got $(typeof(xi_raw))"))

    p_num_raw = get(kwargs, :p_num, default_momentum_count())
    p_num_raw isa Integer || throw(ArgumentError("p_num must be Integer, got $(typeof(p_num_raw))"))
    p_num_raw > 0 || throw(ArgumentError("p_num must be positive, got $(p_num_raw)"))

    t_num_raw = get(kwargs, :t_num, default_theta_count())
    t_num_raw isa Integer || throw(ArgumentError("t_num must be Integer, got $(typeof(t_num_raw))"))
    t_num_raw > 0 || throw(ArgumentError("t_num must be positive, got $(t_num_raw)"))

    residual_norm_max_raw = get(kwargs, :residual_norm_max, 1e-6)
    residual_norm_max_raw isa Real || throw(ArgumentError("residual_norm_max must be Real, got $(typeof(residual_norm_max_raw))"))
    residual_norm_max_raw > 0 || throw(ArgumentError("residual_norm_max must be positive, got $(residual_norm_max_raw)"))

    return (
        xi=Float64(xi_raw),
        p_num=Int(p_num_raw),
        t_num=Int(t_num_raw),
        residual_norm_max=Float64(residual_norm_max_raw),
    )
end

@inline function _resolve_bool_option(kwargs, key::Symbol, default::Bool)::Bool
    value_raw = get(kwargs, key, default)
    value_raw isa Bool || throw(ArgumentError("$(key) must be Bool, got $(typeof(value_raw))"))
    return value_raw
end

@inline function _candidate_is_physical_for_selection(cand)::Bool
    if !Bool(get(cand, :converged, false))
        return false
    end
    if !isfinite(get(cand, :omega, NaN)) || !isfinite(get(cand, :pressure, NaN)) ||
       !isfinite(get(cand, :rho_norm, NaN)) || !isfinite(get(cand, :entropy, NaN)) ||
       !isfinite(get(cand, :energy, NaN))
        return false
    end
    if !haskey(cand, :x_state) || !haskey(cand, :masses)
        return false
    end
    return is_physical_solution(cand.x_state, cand.masses)
end

@inline function _select_pressure_max_local(candidates::AbstractVector)
    isempty(candidates) && throw(ArgumentError("candidates must be non-empty"))
    candidate_required(cand, field::Symbol) = begin
        if hasproperty(cand, field)
            return getproperty(cand, field)
        end
        if hasmethod(haskey, Tuple{typeof(cand), Symbol}) && haskey(cand, field)
            return get(cand, field, nothing)
        end
        throw(ArgumentError("candidate missing required field :$(field) for _select_pressure_max_local"))
    end
    candidate_seed_index(cand, fallback::Int) = haskey(cand, :seed_index) ? Int(get(cand, :seed_index, fallback)) : fallback
    candidate_residual(cand) = begin
        value = Float64(candidate_required(cand, :residual_norm))
        isfinite(value) ? value : Inf
    end
    candidate_pressure(cand) = begin
        value = Float64(candidate_required(cand, :pressure))
        isfinite(value) ? value : -Inf
    end

    function better_candidate(cand, cand_idx::Int, best, best_idx::Int)
        cand_ok = Bool(candidate_required(cand, :hard_constraint_ok))
        best_ok = Bool(candidate_required(best, :hard_constraint_ok))
        if cand_ok != best_ok
            return cand_ok
        end

        cand_pressure = candidate_pressure(cand)
        best_pressure = candidate_pressure(best)
        if cand_pressure != best_pressure
            return cand_pressure > best_pressure
        end

        cand_residual = candidate_residual(cand)
        best_residual = candidate_residual(best)
        if cand_residual != best_residual
            return cand_residual < best_residual
        end

        cand_seed_index = candidate_seed_index(cand, cand_idx)
        best_seed_index = candidate_seed_index(best, best_idx)
        return cand_seed_index < best_seed_index
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
    return selected
end

@inline function _fixedmu_multiseed_selector_adapter(candidates::AbstractVector)
    converged_physical = [c for c in candidates if _candidate_is_physical_for_selection(c)]
    if !isempty(converged_physical)
        return _select_pressure_max_local(converged_physical)
    end
    return _select_pressure_max_local(candidates)
end

@inline function _resolve_nonfixedmu_bridge(mode::ConstraintMode, T_fm::Real, kwargs)
    common = _resolve_common_runtime_options(kwargs)
    rho0_raw = get(kwargs, :rho0, Main.Constants_PNJL.ρ0_inv_fm3)
    rho0_raw isa Real || throw(ArgumentError("rho0 must be Real, got $(typeof(rho0_raw))"))
    rho0 = Float64(rho0_raw)

    semantic_mode = get(kwargs, :semantic_mode, :ground_state)
    (semantic_mode === :ground_state || semantic_mode === :constrained_manifold) ||
        throw(ArgumentError("semantic_mode must be :ground_state or :constrained_manifold, got $(semantic_mode)"))

    selector = get(kwargs, :selector, nothing)
    selector === nothing || selector isa Function ||
        throw(ArgumentError("selector must be Function or nothing, got $(typeof(selector))"))

    seed_strategy = get(kwargs, :seed_strategy, DefaultSeed())
    seed_guess = haskey(kwargs, :seed_guess) ? kwargs[:seed_guess] : get_seed(seed_strategy, [T_fm], mode)
    seed_guess = _normalize_seed_vector(seed_guess, :seed_guess)
    seed_candidates = if haskey(kwargs, :seed_candidates)
        _normalize_seed_pool(kwargs[:seed_candidates], :seed_candidates)
    else
        [seed_guess]
    end

    return (
        xi=common.xi,
        p_num=common.p_num,
        t_num=common.t_num,
        residual_norm_max=common.residual_norm_max,
        rho0=rho0,
        semantic_mode=semantic_mode,
        selector=selector,
        seed_strategy=seed_strategy,
        seed_guess=seed_guess,
        seed_candidates=seed_candidates,
    )
end

@inline function _resolve_fixedmu_runtime_options(mode::FixedMu, T_fm::Real, μ_fm::Real, kwargs)
    common = _resolve_common_runtime_options(kwargs)
    auto_multiseed_fallback = _resolve_bool_option(kwargs, :auto_multiseed_fallback, true)
    seed_strategy = get(kwargs, :seed_strategy, DefaultSeed())
    seed_guess_raw = haskey(kwargs, :seed_guess) ? kwargs[:seed_guess] : get_seed(seed_strategy, [T_fm, μ_fm], mode)
    seed_guess = _normalize_seed_vector(seed_guess_raw, :seed_guess)
    return (
        xi=common.xi,
        p_num=common.p_num,
        t_num=common.t_num,
        residual_norm_max=common.residual_norm_max,
        auto_multiseed_fallback=auto_multiseed_fallback,
        seed_strategy=seed_strategy,
        seed_guess=seed_guess,
    )
end

@inline function _resolve_fixedmu_multi_runtime_options(mode::FixedMu, T_fm::Real, μ_fm::Real, kwargs)
    common = _resolve_common_runtime_options(kwargs)
    evaluate_all_attempts = _resolve_bool_option(kwargs, :evaluate_all_attempts, true)
    seed_strategy = get(kwargs, :seed_strategy, MultiSeed())
    explicit_seeds = if haskey(kwargs, :seeds)
        _normalize_seed_pool(kwargs[:seeds], :seeds)
    else
        Vector{Vector{Float64}}()
    end
    return (
        xi=common.xi,
        p_num=common.p_num,
        t_num=common.t_num,
        residual_norm_max=common.residual_norm_max,
        evaluate_all_attempts=evaluate_all_attempts,
        seed_strategy=seed_strategy,
        explicit_seeds=explicit_seeds,
    )
end

@inline function _solve_with_problem_spec_default(
    model::AbstractQCDModel,
    mode::ConstraintMode,
    T_fm::Real,
    kwargs,
)
    haskey(kwargs, :use_problem_spec) && throw(ArgumentError("use_problem_spec has been removed; use solve_constraint(...; problem_spec=...) or default ProblemSpec chain"))
    haskey(kwargs, :allow_legacy_path) && throw(ArgumentError("allow_legacy_path has been removed; solve_constraint now only uses ProblemSpec chain"))
    haskey(kwargs, :warn_on_legacy_path) && throw(ArgumentError("warn_on_legacy_path has been removed; solve_constraint now only uses ProblemSpec chain"))
    haskey(kwargs, :fixedmu_use_problem_spec) && throw(ArgumentError("fixedmu_use_problem_spec has been removed; solve_constraint now always uses ProblemSpec chain"))
    haskey(kwargs, :legacy_fallback_plugin) && throw(ArgumentError("legacy_fallback_plugin has been removed; governed solver no longer accepts legacy fallback plugin path"))

    spec = get(kwargs, :problem_spec, nothing)
    forwarded = (; (k => v for (k, v) in pairs(kwargs) if k != :problem_spec)...)
    if spec === nothing
        spec = build_problem_spec(mode)
    else
        spec isa ProblemSpec || throw(ArgumentError("problem_spec must be ProblemSpec or nothing, got $(typeof(spec))"))
    end
    return spec.forward_solve(model, T_fm; forwarded...)
end

function solve_constraint(model::AbstractQCDModel, mode::FixedMu, T_fm::Real;
    problem_spec::Union{Nothing, ProblemSpec}=nothing,
    μ_fm::Real,
    kwargs...)
    # FixedMu 与其他 mode 一致，统一走 ProblemSpec 主链。
    merged = (; kwargs..., μ_fm=μ_fm, problem_spec=problem_spec)
    return _solve_with_problem_spec_default(model, mode, T_fm, merged)
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
    kwargs = _resolve_primary_strategy_kwargs(kwargs)
    opts = _resolve_fixedmu_runtime_options(mode, T_fm, μ_fm, kwargs)

    if opts.seed_strategy isa MultiSeed
        return solve_multi(model, mode, T_fm, μ_fm; kwargs...)
    end

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
        seed_guess=opts.seed_guess,
        xi=opts.xi,
        p_num=opts.p_num,
        t_num=opts.t_num,
        residual_norm_max=opts.residual_norm_max,
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
        Float64(raw.rho_norm) / Float64(Main.Constants_PNJL.ρ0_inv_fm3),
        Float64(raw.entropy),
        Float64(raw.energy),
        raw.masses,
        Int(raw.iterations),
        Float64(raw.residual_norm),
        Float64(opts.xi),
    )

    if single.converged || !opts.auto_multiseed_fallback
        return single
    end

    try
        return solve_multi(
            model,
            mode,
            T_fm,
            μ_fm;
            seed_strategy=MultiSeed(),
            xi=opts.xi,
            p_num=opts.p_num,
            t_num=opts.t_num,
            residual_norm_max=opts.residual_norm_max,
            forwarded...,
        )
    catch
        return single
    end
end

function solve(model::AbstractPNJLModel, mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, T_fm::Real; kwargs...)
    kwargs = _resolve_primary_strategy_kwargs(kwargs)
    effective_model = _resolve_solver_model(model, kwargs)
    bridge = _resolve_nonfixedmu_bridge(mode, T_fm, kwargs)

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

function solve_multi(model::AbstractPNJLModel, mode::FixedMu, T_fm::Real, μ_fm::Real; kwargs...)
    kwargs = _resolve_primary_strategy_kwargs(kwargs)
    opts = _resolve_fixedmu_multi_runtime_options(mode, T_fm, μ_fm, kwargs)

    strategy_seeds = if isempty(opts.explicit_seeds)
        if opts.seed_strategy isa MultiSeed
            get_all_seeds(opts.seed_strategy, [T_fm, μ_fm], mode)
        else
            [get_seed(opts.seed_strategy, [T_fm, μ_fm], mode)]
        end
    else
        Vector{Vector{Float64}}()
    end
    fallback_default = [get_seed(DefaultSeed(), [T_fm, μ_fm], mode)]
    seed_pool = if !isempty(opts.explicit_seeds)
        build_seed_pool(mode;
            primary_seed=opts.explicit_seeds[1],
            extra_seed_pool=opts.explicit_seeds[2:end],
            seed_extend=(seed, _) -> Float64.(seed),
        )
    elseif !isempty(strategy_seeds)
        build_seed_pool(mode;
            primary_seed=strategy_seeds[1],
            extra_seed_pool=strategy_seeds[2:end],
            default_seed_pool=fallback_default,
            seed_extend=(seed, _) -> Float64.(seed),
        )
    else
        build_seed_pool(mode;
            primary_seed=fallback_default[1],
            seed_extend=(seed, _) -> Float64.(seed),
        )
    end
    seeds = [entry.seed for entry in seed_pool]

    forwarded = _strip_forward_kwargs(kwargs, (
        :seed_strategy,
        :seeds,
        :xi,
        :p_num,
        :t_num,
        :residual_norm_max,
        :evaluate_all_attempts,
    ))
    candidates = execute_attempt_pool(seeds;
        stop_on_first_success=true,
        evaluate_all_attempts=opts.evaluate_all_attempts,
        evaluate_attempt=(seed, seed_index) -> begin
            raw = solve_constraint(
                model,
                mode,
                T_fm;
                μ_fm=μ_fm,
                seed_guess=seed,
                xi=opts.xi,
                p_num=opts.p_num,
                t_num=opts.t_num,
                residual_norm_max=opts.residual_norm_max,
                forwarded...,
            )
            thermo_finite = isfinite(raw.omega) && isfinite(raw.pressure) && isfinite(raw.rho_norm) && isfinite(raw.entropy) && isfinite(raw.energy)
            phys_ok = thermo_finite && is_physical_solution(raw.x_state, raw.masses)
            residual_ok = isfinite(raw.residual_norm) && raw.residual_norm <= opts.residual_norm_max
            ok = (Bool(raw.converged) || (residual_ok && phys_ok))
            raw_candidate = (
                converged=ok,
                solution=Float64.(raw.solution),
                x_state=raw.x_state,
                mu_vec=raw.mu_vec,
                omega=Float64(raw.omega),
                pressure=Float64(raw.pressure),
                rho_norm=Float64(raw.rho_norm) / Float64(Main.Constants_PNJL.ρ0_inv_fm3),
                entropy=Float64(raw.entropy),
                energy=Float64(raw.energy),
                masses=raw.masses,
                iterations=Int(raw.iterations),
                residual_norm=Float64(raw.residual_norm),
                hard_constraint_ok=ok,
                failed_constraints=(ok ? Symbol[] : Symbol[:residual_too_large]),
            )
            merged = build_governance_candidate(raw_candidate;
                hard_constraint_ok=ok,
                failed_constraints=(ok ? Symbol[] : Symbol[:residual_too_large]),
                seed_index=Int(seed_index),
                residual_norm_max=opts.residual_norm_max,
            )
            return merged, evaluate_candidate_success(merged; residual_norm_max=opts.residual_norm_max)
        end,
        on_error=(_, seed_index, err) -> begin
            err_kind = classify_attempt_error(err)
            err_msg = normalize_error_message(err)
            raw_candidate = (
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
            )
            merged = build_governance_candidate(raw_candidate;
                hard_constraint_ok=false,
                failed_constraints=Symbol[:solver_failed],
                seed_index=Int(seed_index),
                residual_norm_max=opts.residual_norm_max,
                error_kind=err_kind,
                error_msg=err_msg,
            )
            return merged, evaluate_candidate_success(merged; residual_norm_max=opts.residual_norm_max)
        end,
    )

    s = _fixedmu_multiseed_selector_adapter(candidates)
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
        Float64(opts.xi),
    )
end

function solve_multi(model::AbstractPNJLModel, mode::Union{FixedRho, FixedAsymmetricRho}, T_fm::Real; kwargs...)
    kwargs = _resolve_primary_strategy_kwargs(kwargs)
    effective_model = _resolve_solver_model(model, kwargs)
    bridge = _resolve_nonfixedmu_bridge(mode, T_fm, kwargs)
    evaluate_all_attempts = _resolve_bool_option(kwargs, :evaluate_all_attempts, true)

    seed_strategy = haskey(kwargs, :seed_strategy) ? bridge.seed_strategy : MultiSeed()
    explicit_seeds = if haskey(kwargs, :seeds)
        _normalize_seed_pool(kwargs[:seeds], :seeds)
    else
        Vector{Vector{Float64}}()
    end
    provided_seed_pool = if haskey(kwargs, :seed_candidates)
        _normalize_seed_pool(kwargs[:seed_candidates], :seed_candidates)
    else
        Vector{Vector{Float64}}()
    end
    strategy_seeds = if isempty(explicit_seeds) && isempty(provided_seed_pool)
        if seed_strategy isa MultiSeed
            get_all_seeds(seed_strategy, [T_fm], mode)
        else
            [get_seed(seed_strategy, [T_fm], mode)]
        end
    else
        Vector{Vector{Float64}}()
    end
    seed_pool = if !isempty(explicit_seeds)
        build_seed_pool(mode;
            primary_seed=explicit_seeds[1],
            extra_seed_pool=explicit_seeds[2:end],
            provided_seed_pool=provided_seed_pool,
            default_seed_pool=[Float64.(bridge.seed_guess)],
            seed_extend=(seed, _) -> Float64.(seed),
        )
    elseif !isempty(provided_seed_pool)
        build_seed_pool(mode;
            primary_seed=provided_seed_pool[1],
            extra_seed_pool=provided_seed_pool[2:end],
            default_seed_pool=[Float64.(bridge.seed_guess)],
            seed_extend=(seed, _) -> Float64.(seed),
        )
    elseif !isempty(strategy_seeds)
        build_seed_pool(mode;
            primary_seed=strategy_seeds[1],
            extra_seed_pool=strategy_seeds[2:end],
            default_seed_pool=[Float64.(bridge.seed_guess)],
            seed_extend=(seed, _) -> Float64.(seed),
        )
    else
        build_seed_pool(mode;
            primary_seed=Float64.(bridge.seed_guess),
            seed_extend=(seed, _) -> Float64.(seed),
        )
    end
    seeds = [entry.seed for entry in seed_pool]

    forwarded = _strip_forward_kwargs(kwargs, (
        :seed_strategy,
        :seeds,
        :seed_guess,
        :seed_candidates,
        :semantic_mode,
        :selector,
        :rho0,
        :xi,
        :p_num,
        :t_num,
        :residual_norm_max,
        :evaluate_all_attempts,
        :model_kind,
        :problem_spec,
    ))

    rho0_kwargs = (mode isa FixedRho) ? NamedTuple() : (; rho0=bridge.rho0)
    raw = solve_constraint(
        effective_model,
        mode,
        T_fm;
        problem_spec=get(kwargs, :problem_spec, nothing),
        seed_guess=seeds[1],
        seed_candidates=seeds,
        semantic_mode=bridge.semantic_mode,
        selector=bridge.selector,
        rho0_kwargs...,
        xi=bridge.xi,
        p_num=bridge.p_num,
        t_num=bridge.t_num,
        residual_norm_max=bridge.residual_norm_max,
        evaluate_all_attempts=evaluate_all_attempts,
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
        forward_kwargs = (; (k => v for (k, v) in pairs(kwargs) if k != :seeds && k != :selector)...)
        return solve_multi(model, mode, T_fm, μ_fm;
            seed_strategy=MultiSeed(candidates),
            forward_kwargs...)
    end
    return solve_multi(model, mode, T_fm, μ_fm; kwargs...)
end

function solve_multi(mode::Union{FixedRho, FixedAsymmetricRho}, T_fm::Real; kwargs...)
    model = create_model(:PNJL)
    return solve_multi(model, mode, T_fm; kwargs...)
end

@inline function _extract_theta_fixedmu(theta_vec::AbstractVector{<:Real})
    length(theta_vec) == 2 || throw(ArgumentError("FixedMu expects theta_vec length 2 ([T_fm, μ_fm]), got $(length(theta_vec))"))
    return Float64(theta_vec[1]), Float64(theta_vec[2])
end

@inline function _extract_theta_nonfixedmu(theta_vec::AbstractVector{<:Real}, mode::ConstraintMode)
    _ = mode
    length(theta_vec) == 1 || throw(ArgumentError("$(typeof(mode)) expects theta_vec length 1 ([T_fm]), got $(length(theta_vec))"))
    return Float64(theta_vec[1])
end

function solve_vec(model::AbstractPNJLModel, mode::FixedMu, theta_vec::AbstractVector{<:Real}; kwargs...)
    T_fm, μ_fm = _extract_theta_fixedmu(theta_vec)
    return solve(model, mode, T_fm, μ_fm; kwargs...)
end

function solve_vec(model::AbstractPNJLModel, mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, theta_vec::AbstractVector{<:Real}; kwargs...)
    T_fm = _extract_theta_nonfixedmu(theta_vec, mode)
    return solve(model, mode, T_fm; kwargs...)
end

function solve_named(model::AbstractPNJLModel, mode::FixedMu, theta_named::NamedTuple; kwargs...)
    haskey(theta_named, :T_fm) || throw(ArgumentError("FixedMu solve_named requires :T_fm"))
    haskey(theta_named, :μ_fm) || throw(ArgumentError("FixedMu solve_named requires :μ_fm"))
    theta_vec = [theta_named[:T_fm], theta_named[:μ_fm]]
    return solve_vec(model, mode, theta_vec; kwargs...)
end

function solve_named(model::AbstractPNJLModel, mode::Union{FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma}, theta_named::NamedTuple; kwargs...)
    haskey(theta_named, :T_fm) || throw(ArgumentError("$(typeof(mode)) solve_named requires :T_fm"))
    theta_vec = [theta_named[:T_fm]]
    return solve_vec(model, mode, theta_vec; kwargs...)
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
