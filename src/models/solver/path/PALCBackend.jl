import BifurcationKit
import BifurcationKit: @optic
import LinearAlgebra

"""Per-path counters, including anchor discovery and target-point correction."""
Base.@kwdef mutable struct PALCWorkTelemetry
    root_solves::Int = 0
    root_residual_calls::Int = 0
    postprocess_residual_calls::Int = 0
    continuation_residual_calls::Int = 0
    continuation_jacobians::Int = 0
    continuation_points::Int = 0
    anchor_solves::Int = 0
    correction_solves::Int = 0
    fallback_solves::Int = 0
    rejected_roots::Int = 0
    exceptions::Int = 0
    root_iterations::Int = 0
    root_jacobian_calls::Int = 0
    workspace_builds::Int = 0
    residual_first_accepts::Int = 0
    reused_anchor_solves::Int = 0
    newton_corrections::Int = 0
    corrector_fallbacks::Int = 0
    root_time_s::Float64 = 0.0
    root_candidate_time_s::Float64 = 0.0
    continuation_time_s::Float64 = 0.0
    selected_rebuild_time_s::Float64 = 0.0
end

palc_work_snapshot(work::PALCWorkTelemetry) =
    (; (key => getfield(work, key) for key in fieldnames(PALCWorkTelemetry))...)

const PATH_SOLVER_KEYWORDS = (:residual_norm_max, :iterations,
    :thermo_quadrature_policy, :thermo_quadrature_rtol, :thermo_quadrature_atol,
    :thermo_quadrature_maxevals, :physicality_check, :hard_constraints,
    :work_telemetry, :path_work, :evaluate_all_attempts,
    :reuse_corrector, :corrector_residual_first, :corrector_method)

function _scan_path_controls_supported(kwargs)
    return all(key -> key in PATH_SOLVER_KEYWORDS, keys(kwargs)) &&
        get(kwargs, :evaluate_all_attempts, true) === true
end

function _path_options(model, path; residual_norm_max::Real=1e-6,
    iterations::Integer=300, thermo_quadrature_policy::Symbol=:tensor_gauss,
    thermo_quadrature_rtol::Real=1e-8, thermo_quadrature_atol::Real=1e-10,
    thermo_quadrature_maxevals::Integer=10^7, physicality_check::Function=((_, _) -> true),
    hard_constraints=nothing, reuse_corrector::Bool=false,
    corrector_residual_first::Bool=false, corrector_method::Symbol=:trust_region, kwargs...)
    unknown = setdiff(collect(keys(kwargs)), PATH_SOLVER_KEYWORDS)
    isempty(unknown) || throw(ArgumentError("unsupported path solver keywords: $(join(unknown, ", "))"))
    get(kwargs, :evaluate_all_attempts, true) === true || throw(ArgumentError(
        "path candidate discovery requires evaluate_all_attempts=true"))
    model isa Union{PNJLModel,RPNJLModel} ||
        throw(UnsupportedCapabilityError(model_kind_symbol(model), :palc_continuation))
    if !(path isa FixedMuPath)
        cached_model = get_cached_model(_joint_model_kind(model))
        same_model = model isa PNJLModel ? model.consts == cached_model.consts :
            model.base.consts == cached_model.base.consts && model.ext == cached_model.ext &&
            model.use_extensions == cached_model.use_extensions
        same_model || throw(ArgumentError(
            "joint constraint paths require the configured model profile used by the joint residual; " *
            "select the profile before constructing the path"))
    end
    PNJLIntegrals.validate_rs_anisotropy(path.xi)
    PNJLIntegrals.validate_thermal_quadrature_policy(thermo_quadrature_policy)
    thermo_quadrature_policy === :tensor_gauss || throw(ArgumentError(
        "PALC paths currently support thermo_quadrature_policy=:tensor_gauss"))
    residual_tol = _path_positive_real(residual_norm_max, :residual_norm_max)
    iteration_limit = _path_positive_int(iterations, :iterations)
    PNJLIntegrals.validate_thermal_quadrature_controls(
        Float64(thermo_quadrature_rtol), Float64(thermo_quadrature_atol),
        Int(thermo_quadrature_maxevals))
    rules = hard_constraints === nothing ? default_hard_constraint_rules(; physicality_check) : hard_constraints
    corrector_method in (:trust_region, :newton) ||
        throw(ArgumentError("corrector_method must be :trust_region or :newton"))
    return (; residual_tol, iteration_limit, rules, reuse_corrector, corrector_residual_first, corrector_method,
        thermo=(; xi=path.xi, p_num=path.p_num, t_num=path.t_num,
            thermo_quadrature_policy, thermo_quadrature_rtol=Float64(thermo_quadrature_rtol),
            thermo_quadrature_atol=Float64(thermo_quadrature_atol),
            thermo_quadrature_maxevals=Int(thermo_quadrature_maxevals)))
end

struct FixedMuPathResidual{M,P,O}
    model::M
    path::P
    options::O
end

function (r::FixedMuPathResidual)(x, value)
    return collect(gap_residual(r.model, x, r.path.T_fm,
        normalize_mu_vec(value); r.options.thermo...))
end

function _path_residual!(out, r::FixedMuPathResidual, x, value)
    copyto!(out, gap_residual(r.model, x, r.path.T_fm,
        normalize_mu_vec(value); r.options.thermo...))
    return out
end

struct JointPathResidual{P,G}
    path::P
    params::G
end

function _path_residual!(out, r::JointPathResidual, x, value)
    explicit_residual!(out, x, SVector(_path_temperature(r.path, value)),
        r.params, _path_mode(r.path, value))
    return out
end

(r::JointPathResidual)(x, value) = _path_residual!(similar(x), r, x, value)
_path_residual!(out, residual, x, value) = copyto!(out, residual(x, value))

function _path_residual(model, path, options)
    if path isa FixedMuPath
        # Use the same model object and stationarity equation as FixedMu production.
        return FixedMuPathResidual(model, path, options)
    end
    nodes = cached_nodes(path.p_num, path.t_num;
        p_max_inv_fm=thermal_p_max_inv_fm(model))
    params = GapParams(_path_temperature(path, first(_path_values(path))), nodes, path.xi;
        p_num=path.p_num, t_num=path.t_num, model_kind=_joint_model_kind(model))
    return JointPathResidual(path, params)
end

mutable struct PathRootFunction{R}
    residual::R
    value::Float64
    work::PALCWorkTelemetry
end

function (f::PathRootFunction)(out, x)
    f.work.root_residual_calls += 1
    _path_residual!(out, f.residual, x, f.value)
    return nothing
end

struct PathRootWorkspace{F,D,C}
    f!::F
    objective::D
    cache::C
end

function PathRootWorkspace(residual, value, seed, work)
    f! = PathRootFunction(residual, Float64(value), work)
    x = Float64.(seed)
    objective = NLsolve.OnceDifferentiable(f!, x, similar(x); autodiff=:forward)
    work.workspace_builds += 1
    return PathRootWorkspace(f!, objective, NLsolve.NewtonTrustRegionCache(objective))
end

function _path_workspace(residual, path, options, work)
    options === nothing && return nothing
    get(options, :reuse_corrector, false) || get(options, :corrector_method, :trust_region) === :newton || return nothing
    attempt = first(_endpoint_anchor_attempts(path))
    return PathRootWorkspace(residual, attempt.value, attempt.seed, work)
end

"""Local Newton correction with residual acceptance before building the next Jacobian."""
function _newton_corrector!(workspace::PathRootWorkspace, value, seed, iteration_limit)
    workspace.f!.value = Float64(value)
    df, cache = workspace.objective, workspace.cache
    NLsolve.NLSolversBase.clear!(df)
    copyto!(cache.x, seed)
    NLsolve.NLSolversBase.value!!(df, cache.x)
    iterations = 0
    accepted = false
    for iteration in 0:min(iteration_limit, 8)
        residual_norm = maximum(abs, df.F)
        if isfinite(residual_norm) && residual_norm <= 1e-9
            accepted = true
            break
        end
        isfinite(residual_norm) && iteration < min(iteration_limit, 8) || break
        NLsolve.NLSolversBase.jacobian!!(df, cache.x)
        all(isfinite, df.DF) || break
        try
            copyto!(cache.p, df.DF \ df.F)
        catch err
            err isa Union{LinearAlgebra.SingularException,LinearAlgebra.LAPACKException} || rethrow()
            break
        end
        all(isfinite, cache.p) || break
        copyto!(cache.xold, cache.x)
        improved = false
        for reduction in 0:5
            @. cache.x = cache.xold - (0.5^reduction) * cache.p
            NLsolve.NLSolversBase.value!!(df, cache.x)
            trial_norm = maximum(abs, df.F)
            if isfinite(trial_norm) && trial_norm < residual_norm
                improved = true
                break
            end
        end
        iterations += 1
        improved || break
    end
    return (; zero=copy(cache.x), iterations, f_calls=first(df.f_calls),
        g_calls=first(df.df_calls), accepted)
end

function _workspace_root!(workspace::PathRootWorkspace, value, seed, iterations)
    workspace.f!.value = Float64(value)
    # The external parameter may change at exactly the same x. Invalidate both
    # value and Jacobian caches as well as the per-solve counters.
    NLsolve.NLSolversBase.clear!(workspace.objective)
    return NLsolve.trust_region(workspace.objective, Float64.(seed), 1e-9, 1e-9,
        iterations, false, false, false, 1.0, true, workspace.cache)
end

function _path_candidate(model, path, value, solution, residual, iterations, options)
    mu_vec = path isa FixedMuPath ? normalize_mu_vec(value) : solution[6:8]
    full_solution = path isa FixedMuPath ? vcat(solution, mu_vec) : solution
    thermo = compute_thermo_from_solution(model, full_solution, _path_temperature(path, value);
        options.thermo..., rho0_scale=rho0)
    residual_norm = LinearAlgebra.norm(residual(solution, value))
    raw = build_solver_candidate(solution, thermo, residual_norm;
        converged=isfinite(residual_norm) && residual_norm <= options.residual_tol,
        iterations, residual_norm_max=options.residual_tol)
    hard_ok, failed = evaluate_hard_constraints(raw, options.rules)
    return (; raw..., converged=raw.converged && hard_ok,
        hard_constraint_ok=hard_ok, failed_constraints=failed)
end

function _path_root(model, path, value, seed, residual, options, work;
    source::Symbol=:correction, work_telemetry=nothing, workspace=nothing)
    started = time_ns()
    work.root_solves += 1
    source === :anchor && (work.anchor_solves += 1)
    source === :correction && (work.correction_solves += 1)
    source === :fallback && (work.fallback_solves += 1)
    record_governed_attempt!(work_telemetry)
    path isa FixedAsymmetricRhoPath && record_fixedasymrho_inner_solve!(work_telemetry)
    f! = function (F, x)
        work.root_residual_calls += 1
        copyto!(F, residual(x, value))
        return nothing
    end
    try
        local_root = nothing
        if source === :correction && options.corrector_method === :newton
            workspace = workspace === nothing ? PathRootWorkspace(residual, value, seed, work) : workspace
            local_root = try
                _newton_corrector!(workspace, value, seed, options.iteration_limit)
            catch err
                err isa Union{InterruptException,MethodError,UndefVarError} && rethrow()
                nothing
            end
            work.newton_corrections += 1
            if local_root !== nothing
                work.root_iterations += local_root.iterations
                work.root_jacobian_calls += local_root.g_calls
                record_nlsolve_work!(work_telemetry, local_root, :newton)
            end
            if local_root === nothing || !local_root.accepted
                work.corrector_fallbacks += 1
                local_root = nothing
            end
        end
        accepted_seed = local_root !== nothing
        if !accepted_seed && source === :correction && options.corrector_residual_first
            work.root_residual_calls += 1
            seed_residual = residual(seed, value)
            accepted_seed = all(isfinite, seed_residual) && maximum(abs, seed_residual) <= 1e-9
        end
        root = if local_root !== nothing
            local_root
        elseif accepted_seed
            work.residual_first_accepts += 1
            (; zero=Float64.(seed), iterations=0)
        elseif workspace === nothing
            nlsolve(f!, Float64.(seed); autodiff=:forward, method=:trust_region,
                ftol=1e-9, xtol=1e-9, iterations=options.iteration_limit)
        else
            _workspace_root!(workspace, value, seed, options.iteration_limit)
        end
        if !accepted_seed
            record_nlsolve_work!(work_telemetry, root, :trust_region)
            work.root_iterations += root.iterations
            work.root_jacobian_calls += root.g_calls
        end
        work.postprocess_residual_calls += 1
        candidate_started = time_ns()
        candidate = _path_candidate(model, path, value, Float64.(root.zero),
            residual, root.iterations, options)
        work.root_candidate_time_s += (time_ns() - candidate_started) / 1e9
        if !candidate.converged && local_root !== nothing
            work.corrector_fallbacks += 1
            root = _workspace_root!(workspace, value, seed, options.iteration_limit)
            record_nlsolve_work!(work_telemetry, root, :trust_region)
            work.root_iterations += root.iterations
            work.root_jacobian_calls += root.g_calls
            work.postprocess_residual_calls += 1
            candidate_started = time_ns()
            candidate = _path_candidate(model, path, value, Float64.(root.zero),
                residual, root.iterations, options)
            work.root_candidate_time_s += (time_ns() - candidate_started) / 1e9
        end
        !candidate.converged && (work.rejected_roots += 1)
        return candidate
    catch err
        err isa InterruptException && rethrow()
        err isa Union{MethodError,UndefVarError} && rethrow()
        work.exceptions += 1
        return nothing
    finally
        work.root_time_s += (time_ns() - started) / 1e9
    end
end

function _path_catalog(path, value; previous=nothing)
    mode = _path_mode(path, value)
    temperature = _path_temperature(path, value)
    theta = path isa FixedMuPath ? [temperature, value] : [temperature]
    seeds = Vector{Vector{Float64}}()
    previous === nothing || push!(seeds, Float64.(previous))
    append!(seeds, get_all_seeds(MultiSeed(), theta, mode))
    if path isa FixedAsymmetricRhoPath
        expanded = Vector{Vector{Float64}}()
        for seed in seeds
            append!(expanded, _fixedasymrho_inner_seeds(mode, path.T_fm, seed, [0.2, 0.2, 0.2]))
        end
        seeds = expanded
    elseif path isa FixedMuPath
        append!(seeds, [SeedStrategies.HADRON_SEED_5, SeedStrategies.QUARK_SEED_5,
            SeedStrategies.HT_GUESS_0p8_SEED_5, SeedStrategies.HT_GUESS_0p9_SEED_5,
            SeedStrategies.HT_GUESS_0p95_SEED_5, SeedStrategies.WEAK_CHIRAL_CONF_SEED_5])
        # FixedMu production adds small perturbations inside each outer request.
        seeds = reduce(vcat, [_build_default_seed_candidates(seed) for seed in seeds])
    elseif path isa Union{ConstraintTPath{FixedEntropy},ConstraintTPath{FixedSigma}}
        # Include the same alternate chemical-potential start used by the
        # existing joint solvers, once per catalog state instead of nesting it.
        append!(seeds, [vcat(seed[1:5], fill(0.2, 3)) for seed in copy(seeds)])
    end
    return sort!(_unique_path_seeds(seeds); by=Tuple)
end

function _single_pass_path_result(model, path, value, previous; kwargs...)
    options = _path_options(model, path; kwargs...)
    residual = _path_residual(model, path, options)
    work = get(kwargs, :path_work, PALCWorkTelemetry())
    roots = NamedTuple[]
    for seed in _path_catalog(path, value; previous)
        candidate = _path_root(model, path, value, seed, residual, options, work;
            work_telemetry=get(kwargs, :work_telemetry, nothing))
        _push_distinct_path_root!(roots, candidate)
    end
    isempty(roots) && return nothing
    selected = select_pressure_max_candidate(roots).selected_candidate
    return SolverResult(_path_mode(path, value), selected.converged, selected.solution,
        selected.x_state, selected.mu_vec, selected.omega, selected.pressure,
        selected.rho_norm, selected.entropy, selected.energy, selected.masses,
        selected.iterations, selected.residual_norm, path.xi)
end

"""Resolve the scan default; explicit strategies never silently change route."""
function resolve_scan_continuation(strategy::Symbol, model_kind::Symbol,
    thermo_policy::Symbol, semantic_mode::Symbol, selector, point_count::Integer;
    constraint_mode::Symbol=:fixed_rho, controls_supported::Bool=true)
    strategy in (:auto, :legacy, :single_pass, :palc, :natural) || throw(ArgumentError(
        "continuation must be :auto, :legacy, :single_pass, :palc or :natural, got $(strategy)"))
    supported = model_kind in (:PNJL, :RPNJL) && thermo_policy === :tensor_gauss &&
        semantic_mode === :ground_state && selector === nothing && controls_supported &&
        constraint_mode in (:fixed_mu, :fixed_rho, :fixed_asymmetric_rho,
            :fixed_entropy, :fixed_sigma, :fixed_mub_conserved_charges)
    if strategy === :auto
        # The matched workloads cover 21-point density and 31-point chemical-
        # potential paths. Keep shorter paths and unqualified models on legacy.
        minimum_points = constraint_mode === :fixed_mu ? 31 : 21
        qualified = supported && model_kind === :PNJL && point_count >= minimum_points &&
            constraint_mode in (:fixed_mu, :fixed_rho, :fixed_asymmetric_rho,
                :fixed_sigma, :fixed_mub_conserved_charges)
        qualified || return :legacy
        return constraint_mode in (:fixed_rho,:fixed_asymmetric_rho) ? :natural : :palc
    end
    strategy === :natural && !(constraint_mode in (:fixed_rho,:fixed_asymmetric_rho)) &&
        throw(ArgumentError("natural continuation requires a density path"))
    strategy in (:palc, :single_pass, :natural) && !supported && throw(ArgumentError(
        "$(strategy) scan requires PNJL/RPNJL, tensor_gauss, the standard ground_state selector " *
        "and supported path solver controls"))
    return strategy
end

function _push_distinct_path_root!(roots, candidate; tol=1e-5)
    candidate === nothing && return false
    candidate.converged || return false
    any(root -> _solution_distance(root.solution, candidate.solution) <= tol, roots) && return false
    push!(roots, candidate)
    return true
end

function _constraint_path_route(strategy, model, path::ConstraintTPath; kwargs...)
    kind = model isa PNJLModel ? :PNJL : model isa RPNJLModel ? :RPNJL : model_kind_symbol(model)
    return resolve_scan_continuation(strategy, kind,
        get(kwargs, :thermo_quadrature_policy, :tensor_gauss), :ground_state, nothing,
        length(path.T_values); constraint_mode=_path_constraint_symbol(path),
        controls_supported=_scan_path_controls_supported(kwargs))
end

"""Return every crossing of a target, including folds; interpolation is only a seed."""
function _palc_interpolate(points, i, target, predictor)
    a, b = points[i-1], points[i]
    weight = (target-a.param)/(b.param-a.param)
    linear = (1-weight) .* a.solution .+ weight .* b.solution
    predictor === :linear && return linear
    # Interpolate within one monotone piece only. Stencils spanning a turning
    # point would mix different roots of the same target parameter.
    3 <= i <= length(points)-1 || return linear
    stencil = points[i-2:i+1]
    steps = diff([p.param for p in stencil])
    all(d -> d * (b.param-a.param) > 0, steps) || return linear
    maximum(abs, steps) <= 20 * minimum(abs, steps) || return linear
    cubic = zero(linear)
    for j in 1:4
        w = prod((target-stencil[k].param)/(stencil[j].param-stencil[k].param)
            for k in 1:4 if k != j)
        cubic .+= w .* stencil[j].solution
    end
    # Retain a local seed even when the inverse curve has large curvature.
    all(isfinite, cubic) && LinearAlgebra.norm(cubic-linear) <=
        LinearAlgebra.norm(b.solution-a.solution) || return linear
    return cubic
end

"""A previous slice supplies an initial guess only; every new slice is solved again."""
function _reusable_endpoint_anchor(previous, path)
    previous === nothing && return EndpointAnchor()
    old = previous.path
    typeof(old) === typeof(path) || return EndpointAnchor()
    old.xi == path.xi && old.p_num == path.p_num && old.t_num == path.t_num || return EndpointAnchor()
    if path isa FixedAsymmetricRhoPath
        old.ud_ratio_target == path.ud_ratio_target && old.s_target == path.s_target || return EndpointAnchor()
    end
    path isa Union{FixedMuPath,FixedRhoPath,FixedAsymmetricRhoPath} || return EndpointAnchor()
    high = maximum(_path_values(path))
    points = [point for branch in previous.branches for point in branch.points
        if point.param_value == high && point.converged && all(isfinite, point.solution)]
    isempty(points) && return EndpointAnchor()
    point = points[argmax(p.pressure for p in points)]
    return EndpointAnchor(point.solution)
end

function _palc_target_seeds(points, target; nearest::Bool=false, predictor::Symbol=:linear)
    predictor in (:linear, :cubic) || throw(ArgumentError("unknown target predictor $(predictor)"))
    seeds = Vector{Vector{Float64}}()
    for point in points
        abs(point.param - target) <= 1e-12 && push!(seeds, copy(point.solution))
    end
    for i in 2:length(points)
        a, b = points[i-1], points[i]
        a.param == b.param && continue
        if min(a.param, b.param) < target < max(a.param, b.param)
            push!(seeds, _palc_interpolate(points, i, target, predictor))
        end
    end
    if nearest && isempty(seeds) && !isempty(points)
        push!(seeds, copy(points[argmin(abs(p.param - target) for p in points)].solution))
    end
    return _unique_path_seeds(seeds)
end

mutable struct PALCJacobianFunction{R}
    residual::R
    value::Float64
    work::PALCWorkTelemetry
end

function (f::PALCJacobianFunction)(out, x)
    f.work.continuation_residual_calls += 1
    _path_residual!(out, f.residual, x, f.value)
    return nothing
end

function _palc_jacobian(residual, seed, value, work, backend)
    f! = PALCJacobianFunction(residual, Float64(value), work)
    x, out = Float64.(seed), similar(seed, Float64)
    matrix = Matrix{Float64}(undef, length(x), length(x))
    if backend === :forwarddiff
        config = ForwardDiff.JacobianConfig(f!, out, x)
        return function (x, p)
            work.continuation_jacobians += 1
            f!.value = p.value
            ForwardDiff.jacobian!(matrix, f!, out, x, config)
            return matrix
        end
    end
    plus, minus, fplus, fminus = copy(x), copy(x), similar(x), similar(x)
    return function (x, p)
        work.continuation_jacobians += 1
        f!.value = p.value
        copyto!(plus, x)
        copyto!(minus, x)
        for j in eachindex(x)
            delta = cbrt(eps(Float64)) * max(1.0, abs(x[j]))
            plus[j], minus[j] = x[j]+delta, x[j]-delta
            f!(fplus, plus)
            f!(fminus, minus)
            for i in eachindex(x)
                matrix[i,j] = (fplus[i]-fminus[i])/(2delta)
            end
            plus[j] = minus[j] = x[j]
        end
        return matrix
    end
end

function _palc_curve(residual, root, value, bounds, strategy, work; positive_parameter::Bool=false)
    f = function (x, p)
        work.continuation_residual_calls += 1
        residual(x, p.value)
    end
    jac = _palc_jacobian(residual, root.solution, value, work, strategy.jacobian)
    record = (x, p; kwargs...) -> (; solution=Float64.(x))
    problem = BifurcationKit.BifurcationProblem(f, copy(root.solution), (; value),
        (@optic _.value); J=jac, record_from_solution=record)
    span = bounds[2] - bounds[1]
    pad = max(strategy.dsmax, 0.05 * span)
    lower = positive_parameter ? max(bounds[1]-pad, eps(Float64)) : bounds[1]-pad
    options = BifurcationKit.ContinuationPar(p_min=lower, p_max=bounds[2]+pad,
        ds=strategy.ds, dsmax=strategy.dsmax, max_steps=strategy.max_steps,
        detect_bifurcation=0, newton_options=BifurcationKit.NewtonPar(
            tol=strategy.newton_tol, max_iterations=strategy.newton_max_iterations))
    started = time_ns()
    curve = try
        BifurcationKit.continuation(problem, BifurcationKit.PALC(), options;
            bothside=true, verbosity=0)
    finally
        work.continuation_time_s += (time_ns() - started) / 1e9
    end
    points = [(; param=Float64(point.param), solution=Float64.(point.solution)) for point in curve.branch]
    work.continuation_points += length(points)
    return points
end

function _default_path_anchor(model, path)
    qualified = model isa PNJLModel && length(unique(_path_values(path))) > 1 &&
        !(path isa ConstraintTPath{FixedEntropy})
    return qualified ? EndpointAnchor() : SeedPoolAnchor()
end

"""Cheap starting guesses at the high endpoint, then the opposite endpoint."""
function _endpoint_anchor_attempts(path, anchor::EndpointAnchor=EndpointAnchor())
    low, high = extrema(_path_values(path))
    temperature = _path_temperature(path, high)
    hot = temperature * 197.327 > 150.0
    hint = if path isa FixedMuPath
        SeedStrategies.auto_phase_hint(temperature, high)
    elseif path isa Union{FixedRhoPath,FixedAsymmetricRhoPath}
        hot || high >= 1.0 ? :quark : :hadron
    else
        hot ? :quark : :hadron
    end
    attempts = NamedTuple[]
    if anchor.seed !== nothing
        length(anchor.seed) == state_dim(_path_mode(path, high)) ||
            throw(ArgumentError("EndpointAnchor seed length must match the path state dimension"))
        push!(attempts, (; value=high, seed=copy(anchor.seed)))
    end
    for (value, phase) in ((high, hint), (high, hint === :quark ? :hadron : :quark),
                           (low, :auto), (low,:quark))
        theta = path isa FixedMuPath ? [_path_temperature(path, value), value] :
            [_path_temperature(path, value)]
        seed = get_seed(DefaultSeed(; phase_hint=phase), theta, _path_mode(path, value))
        if path isa ConstraintTPath{FixedSigma} && value == high && phase === :hadron
            # A finite entropy per baryon benefits from a seed with nonzero
            # thermal occupation, rather than the nearly vacuum Polyakov loop.
            seed = SeedStrategies.extend_seed(SeedStrategies.MEDIUM_SEED_5, _path_mode(path,value))
        end
        any(a -> a.value == value && a.seed == seed, attempts) && continue
        push!(attempts, (; value, seed))
    end
    return attempts
end

function _solve_endpoint_palc_path(model, path, strategy;
    branch_policy, diagnostic_level, branch_jump_tol, path_work, work_telemetry,
    anchor_strategy::EndpointAnchor=EndpointAnchor(), kwargs...)
    started = time()
    options = _path_options(model, path; kwargs...)
    residual = _path_residual(model, path, options)
    workspace = _path_workspace(residual, path, options, path_work)
    messages = String[]
    attempted_roots = NamedTuple[]
    for (attempt_index, attempt) in enumerate(_endpoint_anchor_attempts(path, anchor_strategy))
        attempt_index == 1 && anchor_strategy.seed !== nothing && (path_work.reused_anchor_solves += 1)
        root = _path_root(model, path, attempt.value, attempt.seed, residual, options, path_work;
            source=:anchor, work_telemetry, workspace)
        root === nothing && continue
        root.converged || continue
        any(r -> r.value == attempt.value && _solution_distance(r.solution,root.solution) <= 1e-5,
            attempted_roots) && continue
        push!(attempted_roots,(;value=attempt.value,solution=root.solution))
        exceptions_before = path_work.exceptions
        result = _solve_palc_path(model, path, strategy;
            anchor_strategy=SingleSeedAnchor(attempt.value, root.solution),
            initial_anchor=root, initial_residual=residual, initial_workspace=workspace,
            boundary_correction=true,
            branch_policy, diagnostic_level, branch_jump_tol,
            path_work, work_telemetry, kwargs...)
        if result.diagnostics.failure_count > 0 && path_work.exceptions > exceptions_before
            append!(messages,result.diagnostics.messages)
            push!(messages,"Retrying endpoint PALC from the same root with smaller steps.")
            reduced = PALCContinuation(ds=strategy.ds/10,dsmax=strategy.dsmax/5,
                max_steps=strategy.max_steps,newton_tol=strategy.newton_tol,
                newton_max_iterations=strategy.newton_max_iterations,
                jacobian=strategy.jacobian, target_predictor=strategy.target_predictor)
            result = _solve_palc_path(model,path,reduced;
                anchor_strategy=SingleSeedAnchor(attempt.value,root.solution),
                initial_anchor=root,initial_residual=residual,initial_workspace=workspace,
                boundary_correction=true,
                branch_policy,diagnostic_level,branch_jump_tol,path_work,work_telemetry,kwargs...)
        end
        if result.diagnostics.failure_count == 0
            diagnostics = PathDiagnostics(:bifurcationkit, :endpoint_anchor,
                _policy_symbol(branch_policy), time()-started, path_work.anchor_solves,
                path_work.continuation_points, 0, result.diagnostics.branch_jump_count,
                vcat(messages,result.diagnostics.messages))
            return PathSolveResult(path, result.branches, result.selections, result.anchors, diagnostics)
        end
        append!(messages, result.diagnostics.messages)
        push!(messages, "Endpoint continuation left $(result.diagnostics.failure_count) targets unresolved; trying another endpoint root.")
    end
    isempty(messages) && push!(messages, "Endpoint guesses did not establish a complete path; using catalog recovery.")
    result = _solve_palc_path(model, path, strategy; anchor_strategy=SeedPoolAnchor(),
        branch_policy, diagnostic_level, branch_jump_tol, path_work, work_telemetry, kwargs...)
    append!(messages, result.diagnostics.messages)
    diagnostics = PathDiagnostics(:bifurcationkit, :endpoint_catalog_recovery,
        _policy_symbol(branch_policy), time()-started, path_work.anchor_solves,
        path_work.continuation_points, result.diagnostics.failure_count,
        result.diagnostics.branch_jump_count, messages)
    return PathSolveResult(path, result.branches, result.selections, result.anchors, diagnostics)
end

function _candidate_path_point(value, candidate, source, step)
    return BranchPoint(Float64(value), copy(candidate.solution), collect(candidate.x_state),
        collect(candidate.mu_vec), Float64(candidate.pressure), Float64(candidate.rho_norm),
        Float64(candidate.residual_norm), Bool(candidate.converged), source, step)
end

function _solve_palc_path(model, path::EquilibriumPath, strategy::PALCContinuation;
    anchor_strategy::Union{Nothing,AbstractAnchorStrategy}=nothing,
    branch_policy::AbstractBranchPolicy=GroundStateBranchPolicy(),
    diagnostic_level::Symbol=:summary, branch_jump_tol::Real=DEFAULT_PATH_BRANCH_JUMP_TOL,
    path_work::PALCWorkTelemetry=PALCWorkTelemetry(), work_telemetry=nothing,
    anchor_values=nothing, initial_anchor=nothing, initial_residual=nothing,
    initial_workspace=nothing, boundary_correction::Bool=false, kwargs...)
    diagnostic_level in (:summary, :full) || throw(ArgumentError("diagnostic_level must be :summary or :full"))
    if anchor_strategy === nothing
        anchor_strategy = anchor_values === nothing ? _default_path_anchor(model, path) : SeedPoolAnchor()
    end
    if anchor_strategy isa EndpointAnchor
        anchor_values === nothing || throw(ArgumentError("EndpointAnchor selects its own endpoint; omit anchor_values"))
        return _solve_endpoint_palc_path(model, path, strategy;
            branch_policy, diagnostic_level, branch_jump_tol, path_work, work_telemetry,
            anchor_strategy, kwargs...)
    end
    anchor_strategy isa Union{SeedPoolAnchor,MultiSeedAnchor,SingleSeedAnchor} || throw(ArgumentError(
        "select SeedPoolAnchor, MultiSeedAnchor or SingleSeedAnchor"))
    anchor_strategy isa MultiSeedAnchor && !anchor_strategy.evaluate_all_attempts &&
        throw(ArgumentError("PALC anchor discovery requires evaluate_all_attempts=true"))
    jump_tol = _path_positive_real(branch_jump_tol, :branch_jump_tol)
    started = time()
    options = _path_options(model, path; kwargs...)
    residual = initial_residual === nothing ? _path_residual(model, path, options) : initial_residual
    workspace = initial_workspace === nothing ? _path_workspace(residual, path, options, path_work) : initial_workspace
    values = sort!(unique(_path_values(path)))
    single_seed = anchor_strategy isa SingleSeedAnchor
    single_seed && anchor_values !== nothing && throw(ArgumentError(
        "SingleSeedAnchor already specifies its parameter; omit anchor_values"))
    if single_seed
        length(anchor_strategy.seed) == state_dim(_path_mode(path, anchor_strategy.param_value)) ||
            throw(ArgumentError("SingleSeedAnchor seed length must match the path state dimension"))
    end
    positions = single_seed ? [anchor_strategy.param_value] :
        anchor_values === nothing ? unique([values[cld(length(values), 2)], first(values), last(values)]) :
            _path_float_vector(anchor_values, :anchor_values)
    all(value -> first(values) <= value <= last(values), positions) ||
        throw(ArgumentError("anchor_values must lie within the requested path"))
    curves = Vector{NamedTuple{(:param,:solution),Tuple{Float64,Vector{Float64}}}}[]
    anchors = AnchorRoot[]
    messages = String[]
    for value in positions
        roots = NamedTuple[]
        seeds = single_seed ? [anchor_strategy.seed] : _path_catalog(path, value)
        for seed in seeds
            candidate = single_seed && initial_anchor !== nothing ? initial_anchor :
                _path_root(model, path, value, seed, residual, options, path_work;
                    source=:anchor, work_telemetry, workspace)
            _push_distinct_path_root!(roots, candidate)
        end
        covered = NamedTuple[]
        for curve in curves, seed in _palc_target_seeds(curve, value; predictor=strategy.target_predictor)
            candidate = _path_root(model, path, value, seed, residual, options, path_work; work_telemetry, workspace)
            _push_distinct_path_root!(covered, candidate)
        end
        sort!(roots; by=root -> -root.pressure)
        for root in roots
            any(other -> _solution_distance(root.solution, other.solution) <= 1e-4, covered) && continue
            anchor_id = Symbol("anchor_", length(anchors) + 1)
            point = _candidate_path_point(value, root, single_seed ? :single_seed_anchor : :catalog_anchor, 0)
            push!(anchors, _anchor_from_branch_point(anchor_id, point, length(roots)))
            try
                curve = length(values) == 1 ? [(; param=value, solution=copy(root.solution))] :
                    _palc_curve(residual, root, value, extrema(values), strategy, path_work;
                        positive_parameter=path isa ConstraintTPath)
                push!(curves, curve)
            catch err
                err isa InterruptException && rethrow()
                err isa Union{MethodError,UndefVarError} && rethrow()
                path_work.exceptions += 1
                push!(messages, "continuation failed: " * normalize_error_message(err))
                push!(curves, [(; param=value, solution=copy(root.solution))])
            end
        end
    end
    branches = ContinuationBranch[]
    branch_points = [BranchPoint[] for _ in curves]
    fallback_points = BranchPoint[]
    failure_count = 0
    for (step, value) in enumerate(values)
        accepted = NamedTuple[]
        for (i, curve) in enumerate(curves)
            nearest = boundary_correction && value in extrema(values)
            for seed in _palc_target_seeds(curve, value; nearest, predictor=strategy.target_predictor)
                candidate = _path_root(model, path, value, seed, residual, options, path_work; work_telemetry, workspace)
                if _push_distinct_path_root!(accepted, candidate)
                    push!(branch_points[i], _candidate_path_point(value, candidate, :palc_corrected, step))
                end
            end
        end
        if isempty(accepted) && !single_seed
            for seed in _path_catalog(path, value)
                candidate = _path_root(model, path, value, seed, residual, options, path_work;
                    source=:fallback, work_telemetry, workspace)
                if _push_distinct_path_root!(accepted, candidate)
                    push!(fallback_points, _candidate_path_point(value, candidate, :catalog_fallback, step))
                end
            end
        end
        if isempty(accepted)
            failure_count += 1
            push!(fallback_points, BranchPoint(value, Float64[], zeros(5), zeros(3),
                NaN, NaN, Inf, false, :failed, step))
        end
    end
    snapshot = palc_work_snapshot(path_work)
    for (i, points) in enumerate(branch_points)
        push!(branches, ContinuationBranch(Symbol("palc_", i), anchors[i].anchor_id,
            :bifurcationkit, points, :complete, (; curve=curves[i], work=snapshot)))
    end
    if !isempty(fallback_points)
        push!(branches, ContinuationBranch(:fallback, :fallback_anchor,
            single_seed ? :uncovered_single_seed : :catalog_fallback,
            fallback_points, failure_count == 0 ? :complete : :partial, (; work=snapshot)))
    end
    selected = apply_branch_policy(branches, branch_policy)
    jumps = _path_branch_jump_metrics(branches; tol=jump_tol)
    diagnostics = PathDiagnostics(:bifurcationkit,
        single_seed ? :single_seed_anchor : :catalog_anchors, _policy_symbol(branch_policy),
        time()-started, path_work.anchor_solves, path_work.continuation_points,
        failure_count, jumps.count, messages)
    return PathSolveResult(path, selected.branches, selected.selections, anchors, diagnostics)
end

function _solve_natural_path(model, path; anchor_strategy=nothing,
    branch_policy=GroundStateBranchPolicy(),diagnostic_level=:summary,
    branch_jump_tol=DEFAULT_PATH_BRANCH_JUMP_TOL,path_work=PALCWorkTelemetry(),
    work_telemetry=nothing,recovery_strategy::PALCContinuation=PALCContinuation(),kwargs...)
    path isa Union{FixedRhoPath,FixedAsymmetricRhoPath} ||
        throw(ArgumentError("NaturalContinuation requires FixedRhoPath or FixedAsymmetricRhoPath"))
    anchor_strategy === nothing || anchor_strategy isa EndpointAnchor ||
        throw(ArgumentError("NaturalContinuation selects an endpoint; use EndpointAnchor or omit anchor_strategy"))
    diagnostic_level in (:summary,:full) || throw(ArgumentError("invalid diagnostic_level"))
    jump_tol = _path_positive_real(branch_jump_tol,:branch_jump_tol)
    started = time()
    options = _path_options(model,path;kwargs...)
    residual = _path_residual(model,path,options)
    workspace = _path_workspace(residual, path, options, path_work)
    values = sort!(unique(_path_values(path));rev=true)
    root = nothing
    for (attempt_index, attempt) in enumerate(_endpoint_anchor_attempts(path, something(anchor_strategy, EndpointAnchor())))
        attempt.value == first(values) || continue
        attempt_index == 1 && anchor_strategy !== nothing && anchor_strategy.seed !== nothing &&
            (path_work.reused_anchor_solves += 1)
        candidate = _path_root(model,path,attempt.value,attempt.seed,residual,options,path_work;
            source=:anchor,work_telemetry,workspace)
        candidate === nothing && continue
        candidate.converged || continue
        root = candidate
        break
    end
    points = BranchPoint[]
    anchors = AnchorRoot[]
    for (step,value) in enumerate(values)
        root === nothing && break
        if step > 1
            root = _path_root(model,path,value,root.solution,residual,options,path_work;work_telemetry,workspace)
            root === nothing && break
            root.converged || break
        end
        point = _candidate_path_point(value,root,:natural_corrected,step)
        push!(points,point)
        step == 1 && push!(anchors,_anchor_from_branch_point(:natural_anchor,point,1))
    end
    if length(points) != length(values)
        result = _solve_endpoint_palc_path(model,path,recovery_strategy;
            branch_policy,diagnostic_level,branch_jump_tol,path_work,work_telemetry,kwargs...)
        diagnostics = PathDiagnostics(:natural_palc_recovery,result.diagnostics.anchor_strategy,
            _policy_symbol(branch_policy),time()-started,path_work.anchor_solves,
            path_work.continuation_points,result.diagnostics.failure_count,
            result.diagnostics.branch_jump_count,
            vcat(["Descending density continuation stopped; recovered with PALC."],result.diagnostics.messages))
        return PathSolveResult(path,result.branches,result.selections,result.anchors,diagnostics)
    end
    branches = [ContinuationBranch(:natural,:natural_anchor,:natural_parameter,points,:complete,
        (;work=palc_work_snapshot(path_work)))]
    selected = apply_branch_policy(branches,branch_policy)
    diagnostics = PathDiagnostics(:natural_parameter,:endpoint_anchor,_policy_symbol(branch_policy),
        time()-started,path_work.anchor_solves,length(points),0,
        _path_branch_jump_metrics(branches;tol=jump_tol).count,String[])
    return PathSolveResult(path,selected.branches,selected.selections,anchors,diagnostics)
end

"""Rebuild SolverResult v1 at exactly the requested points, preserving input order."""
function selected_path_results(model, result::PathSolveResult; kwargs...)
    started = time_ns()
    path = result.path
    options = _path_options(model, path; kwargs...)
    path_work = get(kwargs, :path_work, nothing)
    residual = _path_residual(model, path, options)
    selected = Dict(selection.param_value => selection for selection in result.selections)
    out = Union{Nothing,SolverResult}[]
    for value in _path_values(path)
        choice = get(selected, value, nothing)
        if choice === nothing || choice.selected_branch_id === nothing
            push!(out, nothing)
            continue
        end
        branch = only(filter(branch -> branch.branch_id === choice.selected_branch_id, result.branches))
        matches = filter(point -> point.param_value == value && point.converged, branch.points)
        isempty(matches) && (push!(out, nothing); continue)
        point = matches[argmax(point.pressure for point in matches)]
        path_work === nothing || (path_work.postprocess_residual_calls += 1)
        candidate = _path_candidate(model, path, value, point.solution, residual, 0, options)
        push!(out, SolverResult(_path_mode(path, value), candidate.converged, candidate.solution,
            candidate.x_state, candidate.mu_vec, candidate.omega, candidate.pressure,
            candidate.rho_norm, candidate.entropy, candidate.energy, candidate.masses,
            candidate.iterations, candidate.residual_norm, path.xi))
    end
    path_work === nothing || (path_work.selected_rebuild_time_s += (time_ns()-started)/1e9)
    return out
end

function _solve_seed_equilibrium_path(model, path, strategy;
    anchor_strategy=SeedPoolAnchor(), branch_policy=GroundStateBranchPolicy(),
    diagnostic_level=:summary, branch_jump_tol=DEFAULT_PATH_BRANCH_JUMP_TOL, kwargs...)
    messages = _seed_continuation_anchor_messages(anchor_strategy)
    if path isa ConstraintTPath
        messages = ["Temperature continuation uses the previous solution and the existing joint point solver's seed policy."]
    end
    diagnostic_level in (:summary, :full) || throw(ArgumentError("invalid diagnostic_level"))
    started = time()
    points = BranchPoint[]
    previous = nothing
    for (step, value) in enumerate(_path_values(path))
        mode = _path_mode(path, value)
        result = if path isa ConstraintTPath
            # Preserve the existing point solver's joint solve and seed policy.
            seed = previous === nothing ? get(kwargs, :seed_guess,
                get_seed(get(kwargs, :seed_strategy, DefaultSeed()), [value], mode)) : previous
            common = merge((; xi=path.xi,
                p_num=path.p_num, t_num=path.t_num,
                evaluate_all_attempts=strategy.evaluate_all_attempts), (;kwargs...), (;seed_guess=seed))
            solve(model, mode, value; common...)
        else
            seeds = _path_catalog(path, value; previous)
            common = (; seeds, xi=path.xi, p_num=path.p_num, t_num=path.t_num,
                evaluate_all_attempts=strategy.evaluate_all_attempts, kwargs...)
            path isa FixedMuPath ? solve_multi(model, mode, path.T_fm, value; common...) :
                solve_multi(model, mode, path.T_fm; common...)
        end
        push!(points, _branch_point_from_solver_result(value, result, :seed_continuation, step))
        result.converged && (previous = copy(result.solution))
    end
    failures = count(point -> !point.converged, points)
    branches = [ContinuationBranch(:seed_ground_state, :seed_anchor_1, :seed_continuation,
        points, failures == 0 ? :complete : :partial, (;))]
    anchors = [_anchor_from_branch_point(:seed_anchor_1, first(points), 1)]
    selected = apply_branch_policy(branches, branch_policy)
    diagnostics = PathDiagnostics(:seed_continuation, :seed_pool_anchor, _policy_symbol(branch_policy),
        time()-started, 0, length(points), failures,
        _path_branch_jump_metrics(branches; tol=branch_jump_tol).count, messages)
    return PathSolveResult(path, selected.branches, selected.selections, anchors, diagnostics)
end
