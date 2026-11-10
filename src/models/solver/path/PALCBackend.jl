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
end

palc_work_snapshot(work::PALCWorkTelemetry) =
    (; (key => getfield(work, key) for key in fieldnames(PALCWorkTelemetry))...)

const PATH_SOLVER_KEYWORDS = (:residual_norm_max, :iterations,
    :thermo_quadrature_policy, :thermo_quadrature_rtol, :thermo_quadrature_atol,
    :thermo_quadrature_maxevals, :physicality_check, :hard_constraints,
    :work_telemetry, :path_work, :evaluate_all_attempts)

function _scan_path_controls_supported(kwargs)
    return all(key -> key in PATH_SOLVER_KEYWORDS, keys(kwargs)) &&
        get(kwargs, :evaluate_all_attempts, true) === true
end

function _path_options(model, path; residual_norm_max::Real=1e-6,
    iterations::Integer=300, thermo_quadrature_policy::Symbol=:tensor_gauss,
    thermo_quadrature_rtol::Real=1e-8, thermo_quadrature_atol::Real=1e-10,
    thermo_quadrature_maxevals::Integer=10^7, physicality_check::Function=((_, _) -> true),
    hard_constraints=nothing, kwargs...)
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
            "density paths require the configured model profile used by the joint residual; " *
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
    return (; residual_tol, iteration_limit, rules,
        thermo=(; xi=path.xi, p_num=path.p_num, t_num=path.t_num,
            thermo_quadrature_policy, thermo_quadrature_rtol=Float64(thermo_quadrature_rtol),
            thermo_quadrature_atol=Float64(thermo_quadrature_atol),
            thermo_quadrature_maxevals=Int(thermo_quadrature_maxevals)))
end

function _path_residual(model, path, options)
    if path isa FixedMuPath
        # Use the same model object and stationarity equation as FixedMu production.
        return (x, value) -> collect(gap_residual(model, x, path.T_fm,
            normalize_mu_vec(value); options.thermo...))
    end
    nodes = cached_nodes(path.p_num, path.t_num;
        p_max_inv_fm=thermal_p_max_inv_fm(model))
    params = GapParams(path.T_fm, nodes, path.xi;
        p_num=path.p_num, t_num=path.t_num, model_kind=_joint_model_kind(model))
    return function (x, value)
        out = similar(x, length(x))
        explicit_residual!(out, x, [path.T_fm], params, _path_mode(path, value))
        return out
    end
end

function _path_candidate(model, path, value, solution, residual, iterations, options)
    mu_vec = path isa FixedMuPath ? normalize_mu_vec(value) : solution[6:8]
    full_solution = path isa FixedMuPath ? vcat(solution, mu_vec) : solution
    thermo = compute_thermo_from_solution(model, full_solution, path.T_fm;
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
    source::Symbol=:correction, work_telemetry=nothing)
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
        root = nlsolve(f!, Float64.(seed); autodiff=:forward, method=:trust_region,
            ftol=1e-9, xtol=1e-9, iterations=options.iteration_limit)
        record_nlsolve_work!(work_telemetry, root, :trust_region)
        work.postprocess_residual_calls += 1
        candidate = _path_candidate(model, path, value, Float64.(root.zero),
            residual, root.iterations, options)
        !candidate.converged && (work.rejected_roots += 1)
        return candidate
    catch err
        err isa InterruptException && rethrow()
        err isa Union{MethodError,UndefVarError} && rethrow()
        work.exceptions += 1
        return nothing
    end
end

function _path_catalog(path, value; previous=nothing)
    mode = _path_mode(path, value)
    theta = path isa FixedMuPath ? [path.T_fm, value] : [path.T_fm]
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
    strategy in (:auto, :legacy, :single_pass, :palc) || throw(ArgumentError(
        "continuation must be :auto, :legacy, :single_pass or :palc, got $(strategy)"))
    supported = model_kind in (:PNJL, :RPNJL) && thermo_policy === :tensor_gauss &&
        semantic_mode === :ground_state && selector === nothing && controls_supported &&
        constraint_mode in (:fixed_mu, :fixed_rho, :fixed_asymmetric_rho)
    if strategy === :auto
        # The matched workloads cover 21-point density and 31-point chemical-
        # potential paths. Keep shorter paths and unqualified models on legacy.
        minimum_points = constraint_mode === :fixed_mu ? 31 : 21
        qualified = supported && model_kind === :PNJL && point_count >= minimum_points
        return qualified ? :palc : :legacy
    end
    strategy in (:palc, :single_pass) && !supported && throw(ArgumentError(
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

"""Return every crossing of a target, including folds; interpolation is only a seed."""
function _palc_target_seeds(points, target; nearest::Bool=false)
    seeds = Vector{Vector{Float64}}()
    for point in points
        abs(point.param - target) <= 1e-12 && push!(seeds, copy(point.solution))
    end
    for i in 2:length(points)
        a, b = points[i-1], points[i]
        a.param == b.param && continue
        if min(a.param, b.param) < target < max(a.param, b.param)
            weight = (target - a.param) / (b.param - a.param)
            push!(seeds, (1 - weight) .* a.solution .+ weight .* b.solution)
        end
    end
    if nearest && isempty(seeds) && !isempty(points)
        push!(seeds, copy(points[argmin(abs(p.param - target) for p in points)].solution))
    end
    return _unique_path_seeds(seeds)
end

function _palc_curve(residual, root, value, bounds, strategy, work)
    f = function (x, p)
        work.continuation_residual_calls += 1
        residual(x, p.value)
    end
    jac = function (x, p)
        work.continuation_jacobians += 1
        # Same centered finite difference as the qualified BifurcationKit pilot.
        matrix = Matrix{Float64}(undef, length(x), length(x))
        for j in eachindex(x)
            delta = cbrt(eps(Float64)) * max(1.0, abs(x[j]))
            plus, minus = copy(x), copy(x)
            plus[j] += delta
            minus[j] -= delta
            matrix[:, j] = (f(plus, p) - f(minus, p)) / (2 * delta)
        end
        return matrix
    end
    record = (x, p; kwargs...) -> (; solution=Float64.(x))
    problem = BifurcationKit.BifurcationProblem(f, copy(root.solution), (; value),
        (@optic _.value); J=jac, record_from_solution=record)
    span = bounds[2] - bounds[1]
    pad = max(strategy.dsmax, 0.05 * span)
    options = BifurcationKit.ContinuationPar(p_min=bounds[1]-pad, p_max=bounds[2]+pad,
        ds=strategy.ds, dsmax=strategy.dsmax, max_steps=strategy.max_steps,
        detect_bifurcation=0, newton_options=BifurcationKit.NewtonPar(
            tol=strategy.newton_tol, max_iterations=strategy.newton_max_iterations))
    curve = BifurcationKit.continuation(problem, BifurcationKit.PALC(), options;
        bothside=true, verbosity=0)
    points = [(; param=Float64(point.param), solution=Float64.(point.solution)) for point in curve.branch]
    work.continuation_points += length(points)
    return points
end

function _candidate_path_point(value, candidate, source, step)
    return BranchPoint(Float64(value), copy(candidate.solution), collect(candidate.x_state),
        collect(candidate.mu_vec), Float64(candidate.pressure), Float64(candidate.rho_norm),
        Float64(candidate.residual_norm), Bool(candidate.converged), source, step)
end

function _solve_palc_path(model, path::EquilibriumPath, strategy::PALCContinuation;
    anchor_strategy::AbstractAnchorStrategy=SeedPoolAnchor(),
    branch_policy::AbstractBranchPolicy=GroundStateBranchPolicy(),
    diagnostic_level::Symbol=:summary, branch_jump_tol::Real=DEFAULT_PATH_BRANCH_JUMP_TOL,
    path_work::PALCWorkTelemetry=PALCWorkTelemetry(), work_telemetry=nothing,
    anchor_values=nothing, kwargs...)
    diagnostic_level in (:summary, :full) || throw(ArgumentError("diagnostic_level must be :summary or :full"))
    anchor_strategy isa Union{SeedPoolAnchor,MultiSeedAnchor} || throw(ArgumentError(
        "PALC production uses catalog anchors; select SeedPoolAnchor or MultiSeedAnchor and optional anchor_values"))
    anchor_strategy isa MultiSeedAnchor && !anchor_strategy.evaluate_all_attempts &&
        throw(ArgumentError("PALC anchor discovery requires evaluate_all_attempts=true"))
    jump_tol = _path_positive_real(branch_jump_tol, :branch_jump_tol)
    started = time()
    options = _path_options(model, path; kwargs...)
    residual = _path_residual(model, path, options)
    values = sort!(unique(_path_values(path)))
    positions = anchor_values === nothing ? unique([values[cld(length(values), 2)], first(values), last(values)]) :
        _path_float_vector(anchor_values, :anchor_values)
    all(value -> first(values) <= value <= last(values), positions) ||
        throw(ArgumentError("anchor_values must lie within the requested path"))
    curves = Vector{NamedTuple{(:param,:solution),Tuple{Float64,Vector{Float64}}}}[]
    anchors = AnchorRoot[]
    messages = String[]
    for value in positions
        roots = NamedTuple[]
        for seed in _path_catalog(path, value)
            candidate = _path_root(model, path, value, seed, residual, options, path_work;
                source=:anchor, work_telemetry)
            _push_distinct_path_root!(roots, candidate)
        end
        covered = NamedTuple[]
        for curve in curves, seed in _palc_target_seeds(curve, value)
            candidate = _path_root(model, path, value, seed, residual, options, path_work; work_telemetry)
            _push_distinct_path_root!(covered, candidate)
        end
        sort!(roots; by=root -> -root.pressure)
        for root in roots
            any(other -> _solution_distance(root.solution, other.solution) <= 1e-4, covered) && continue
            anchor_id = Symbol("anchor_", length(anchors) + 1)
            point = _candidate_path_point(value, root, :catalog_anchor, 0)
            push!(anchors, _anchor_from_branch_point(anchor_id, point, length(roots)))
            try
                curve = length(values) == 1 ? [(; param=value, solution=copy(root.solution))] :
                    _palc_curve(residual, root, value, extrema(values), strategy, path_work)
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
            for seed in _palc_target_seeds(curve, value)
                candidate = _path_root(model, path, value, seed, residual, options, path_work; work_telemetry)
                if _push_distinct_path_root!(accepted, candidate)
                    push!(branch_points[i], _candidate_path_point(value, candidate, :palc_corrected, step))
                end
            end
        end
        if isempty(accepted)
            for seed in _path_catalog(path, value)
                candidate = _path_root(model, path, value, seed, residual, options, path_work;
                    source=:fallback, work_telemetry)
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
        push!(branches, ContinuationBranch(:fallback, :fallback_anchor, :catalog_fallback,
            fallback_points, failure_count == 0 ? :complete : :partial, (; work=snapshot)))
    end
    selected = apply_branch_policy(branches, branch_policy)
    jumps = _path_branch_jump_metrics(branches; tol=jump_tol)
    diagnostics = PathDiagnostics(:bifurcationkit, :catalog_anchors, _policy_symbol(branch_policy),
        time()-started, path_work.anchor_solves, path_work.continuation_points,
        failure_count, jumps.count, messages)
    return PathSolveResult(path, selected.branches, selected.selections, anchors, diagnostics)
end

"""Rebuild SolverResult v1 at exactly the requested points, preserving input order."""
function selected_path_results(model, result::PathSolveResult; kwargs...)
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
    return out
end

function _solve_seed_equilibrium_path(model, path, strategy;
    anchor_strategy=SeedPoolAnchor(), branch_policy=GroundStateBranchPolicy(),
    diagnostic_level=:summary, branch_jump_tol=DEFAULT_PATH_BRANCH_JUMP_TOL, kwargs...)
    messages = _seed_continuation_anchor_messages(anchor_strategy)
    diagnostic_level in (:summary, :full) || throw(ArgumentError("invalid diagnostic_level"))
    started = time()
    points = BranchPoint[]
    previous = nothing
    for (step, value) in enumerate(_path_values(path))
        mode = _path_mode(path, value)
        seeds = _path_catalog(path, value; previous)
        common = (; seeds, xi=path.xi, p_num=path.p_num, t_num=path.t_num,
            evaluate_all_attempts=strategy.evaluate_all_attempts, kwargs...)
        result = path isa FixedMuPath ? solve_multi(model, mode, path.T_fm, value; common...) :
            solve_multi(model, mode, path.T_fm; common...)
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
