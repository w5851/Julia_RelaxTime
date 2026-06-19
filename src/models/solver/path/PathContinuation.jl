"""
Experimental path-level continuation scaffolding.

This file intentionally has no external continuation dependency. PALC can be
configured here, but the BifurcationKit adapter remains outside the root
environment until the experimental gate is passed.
"""

abstract type AbstractPathProblem end
abstract type AbstractAnchorStrategy end
abstract type AbstractContinuationStrategy end
abstract type AbstractBranchPolicy end

const DEFAULT_PATH_BRANCH_JUMP_TOL = 0.25

@inline function _path_finite_real(value, name::Symbol)::Float64
    value isa Real || throw(ArgumentError("$(name) must be Real, got $(typeof(value))"))
    out = Float64(value)
    isfinite(out) || throw(ArgumentError("$(name) must be finite, got $(value)"))
    return out
end

@inline function _path_positive_real(value, name::Symbol)::Float64
    out = _path_finite_real(value, name)
    out > 0.0 || throw(ArgumentError("$(name) must be positive, got $(value)"))
    return out
end

@inline function _path_nonnegative_real(value, name::Symbol)::Float64
    out = _path_finite_real(value, name)
    out >= 0.0 || throw(ArgumentError("$(name) must be non-negative, got $(value)"))
    return out
end

@inline function _path_positive_int(value, name::Symbol)::Int
    value isa Integer || throw(ArgumentError("$(name) must be Integer, got $(typeof(value))"))
    out = Int(value)
    out > 0 || throw(ArgumentError("$(name) must be positive, got $(value)"))
    return out
end

function _path_float_vector(values, name::Symbol; nonnegative::Bool=false)::Vector{Float64}
    values isa AbstractVector || throw(ArgumentError("$(name) must be an AbstractVector, got $(typeof(values))"))
    isempty(values) && throw(ArgumentError("$(name) must be non-empty"))
    out = Float64[]
    sizehint!(out, length(values))
    for value in values
        v = nonnegative ? _path_nonnegative_real(value, name) : _path_finite_real(value, name)
        push!(out, v)
    end
    return out
end

struct FixedAsymmetricRhoPath <: AbstractPathProblem
    T_fm::Float64
    rho_values::Vector{Float64}
    ud_ratio_target::Float64
    s_target::Float64
    xi::Float64
    p_num::Int
    t_num::Int

    function FixedAsymmetricRhoPath(
        T_fm::Real,
        rho_values::AbstractVector{<:Real};
        ud_ratio_target::Real=0.876,
        s_target::Real=0.0,
        xi::Real=0.0,
        p_num::Integer=default_momentum_count(),
        t_num::Integer=default_theta_count(),
    )
        return new(
            _path_positive_real(T_fm, :T_fm),
            _path_float_vector(rho_values, :rho_values; nonnegative=true),
            _path_positive_real(ud_ratio_target, :ud_ratio_target),
            _path_finite_real(s_target, :s_target),
            _path_finite_real(xi, :xi),
            _path_positive_int(p_num, :p_num),
            _path_positive_int(t_num, :t_num),
        )
    end
end

FixedAsymmetricRhoPath(
    T_fm::Real,
    rho_values::AbstractVector{<:Real},
    ud_ratio_target::Real,
    s_target::Real,
    xi::Real,
    p_num::Integer,
    t_num::Integer,
) = FixedAsymmetricRhoPath(
    T_fm,
    rho_values;
    ud_ratio_target=ud_ratio_target,
    s_target=s_target,
    xi=xi,
    p_num=p_num,
    t_num=t_num,
)

struct MultiSeedAnchor <: AbstractAnchorStrategy
    evaluate_all_attempts::Bool
    root_distance_tol::Float64
end

function MultiSeedAnchor(evaluate_all_attempts::Bool=true, root_distance_tol::Real=0.25)
    return MultiSeedAnchor(evaluate_all_attempts, _path_positive_real(root_distance_tol, :root_distance_tol))
end

struct ProductionLikeAnchor <: AbstractAnchorStrategy
    source_param::Float64
    step::Float64
    root_distance_tol::Float64
end

function ProductionLikeAnchor(source_param::Real, step::Real, root_distance_tol::Real=0.25)
    return ProductionLikeAnchor(
        _path_finite_real(source_param, :source_param),
        _path_positive_real(step, :step),
        _path_positive_real(root_distance_tol, :root_distance_tol),
    )
end

struct CompositeAnchor{S<:AbstractVector{<:AbstractAnchorStrategy}} <: AbstractAnchorStrategy
    strategies::S
    root_distance_tol::Float64
end

function CompositeAnchor(strategies::AbstractVector{<:AbstractAnchorStrategy}, root_distance_tol::Real=0.25)
    isempty(strategies) && throw(ArgumentError("strategies must be non-empty"))
    return CompositeAnchor(collect(AbstractAnchorStrategy, strategies), _path_positive_real(root_distance_tol, :root_distance_tol))
end

function default_fixedasymrho_anchor_strategy()
    return CompositeAnchor(AbstractAnchorStrategy[
        MultiSeedAnchor(true, 0.25),
        ProductionLikeAnchor(1.0, 0.05, 0.25),
    ], 0.25)
end

struct SeedPoolAnchor <: AbstractAnchorStrategy end

struct SeedContinuation <: AbstractContinuationStrategy
    evaluate_all_attempts::Bool
end

SeedContinuation() = SeedContinuation(true)

struct PALCContinuation <: AbstractContinuationStrategy
    ds::Float64
    dsmax::Float64
    max_steps::Int
    newton_tol::Float64
    newton_max_iterations::Int
    backend::Symbol
end

function PALCContinuation(;
    ds::Real=0.01,
    dsmax::Real=0.05,
    max_steps::Integer=200,
    newton_tol::Real=1e-8,
    newton_max_iterations::Integer=12,
    backend::Symbol=:bifurcationkit,
)
    ds_val = _path_positive_real(ds, :ds)
    dsmax_val = _path_positive_real(dsmax, :dsmax)
    dsmax_val >= ds_val || throw(ArgumentError("dsmax must be greater than or equal to ds"))
    max_steps_val = _path_positive_int(max_steps, :max_steps)
    newton_max_iterations_val = _path_positive_int(newton_max_iterations, :newton_max_iterations)
    newton_tol_val = _path_positive_real(newton_tol, :newton_tol)
    backend === :bifurcationkit ||
        throw(ArgumentError("PALCContinuation backend must be :bifurcationkit, got $(backend)"))
    return PALCContinuation(ds_val, dsmax_val, max_steps_val, newton_tol_val, newton_max_iterations_val, backend)
end

struct GroundStateBranchPolicy <: AbstractBranchPolicy
    pressure_gap_tol::Float64
    residual_norm_max::Float64
end

function GroundStateBranchPolicy(pressure_gap_tol::Real=1e-3, residual_norm_max::Real=1e-6)
    return GroundStateBranchPolicy(
        _path_nonnegative_real(pressure_gap_tol, :pressure_gap_tol),
        _path_positive_real(residual_norm_max, :residual_norm_max),
    )
end

struct ContinuationBranchPolicy <: AbstractBranchPolicy
    branch_id::Symbol
end

function ContinuationBranchPolicy(branch_id)
    branch_id isa Symbol || throw(ArgumentError("branch_id must be Symbol, got $(typeof(branch_id))"))
    return ContinuationBranchPolicy(branch_id)
end

struct AllBranchesPolicy <: AbstractBranchPolicy end

struct AnchorRoot
    anchor_id::Symbol
    source::Symbol
    param_value::Float64
    solution::Vector{Float64}
    x_state::Vector{Float64}
    mu_vec::Vector{Float64}
    pressure::Float64
    residual_norm::Float64
    converged::Bool
    candidate_count::Int
end

struct BranchPoint
    param_value::Float64
    solution::Vector{Float64}
    x_state::Vector{Float64}
    mu_vec::Vector{Float64}
    pressure::Float64
    rho_norm::Float64
    residual_norm::Float64
    converged::Bool
    source::Symbol
    step_index::Int
end

struct ContinuationBranch
    branch_id::Symbol
    anchor_id::Symbol
    source::Symbol
    points::Vector{BranchPoint}
    status::Symbol
    diagnostics::NamedTuple
end

struct BranchSelection
    param_value::Float64
    selected_branch_id::Union{Nothing, Symbol}
    selection_reason::Symbol
    candidate_branch_count::Int
    selected_pressure::Float64
    runner_up_branch_id::Union{Nothing, Symbol}
    pressure_gap::Float64
    selected_residual_norm::Float64
end

struct PathDiagnostics
    continuation_backend::Symbol
    anchor_strategy::Symbol
    branch_policy::Symbol
    wall_time_s::Float64
    anchor_solve_count::Int
    continuation_step_count::Int
    failure_count::Int
    branch_jump_count::Int
    messages::Vector{String}
end

struct PathSolveResult
    path::AbstractPathProblem
    branches::Vector{ContinuationBranch}
    selections::Vector{BranchSelection}
    anchors::Vector{AnchorRoot}
    diagnostics::PathDiagnostics
end

@inline _strategy_symbol(::MultiSeedAnchor) = :multi_seed_anchor
@inline _strategy_symbol(::ProductionLikeAnchor) = :production_like_anchor
@inline _strategy_symbol(::CompositeAnchor) = :composite_anchor
@inline _strategy_symbol(::SeedPoolAnchor) = :seed_pool_anchor
@inline _strategy_symbol(::SeedContinuation) = :seed_continuation
@inline _strategy_symbol(c::PALCContinuation) = c.backend
@inline _policy_symbol(::GroundStateBranchPolicy) = :ground_state
@inline _policy_symbol(::ContinuationBranchPolicy) = :continuation_branch
@inline _policy_symbol(::AllBranchesPolicy) = :all_branches

@inline function _branch_point_success(point::BranchPoint, residual_norm_max::Float64)::Bool
    return point.converged &&
        isfinite(point.pressure) &&
        isfinite(point.residual_norm) &&
        point.residual_norm <= residual_norm_max
end

function _points_by_param(branches::AbstractVector{<:ContinuationBranch})
    grouped = Dict{Float64, Vector{Tuple{Symbol, BranchPoint}}}()
    for branch in branches
        for point in branch.points
            bucket = get!(grouped, point.param_value, Tuple{Symbol, BranchPoint}[])
            push!(bucket, (branch.branch_id, point))
        end
    end
    return grouped
end

@inline function _better_ground_state_candidate(
    candidate::Tuple{Symbol, BranchPoint},
    incumbent::Tuple{Symbol, BranchPoint},
)::Bool
    cand_id, cand_point = candidate
    best_id, best_point = incumbent
    cand_point.pressure != best_point.pressure && return cand_point.pressure > best_point.pressure
    cand_point.residual_norm != best_point.residual_norm && return cand_point.residual_norm < best_point.residual_norm
    return String(cand_id) < String(best_id)
end

function _runner_up(
    candidates::AbstractVector{Tuple{Symbol, BranchPoint}},
    selected_branch_id::Union{Nothing, Symbol},
)
    runner = nothing
    for candidate in candidates
        candidate[1] === selected_branch_id && continue
        if runner === nothing || _better_ground_state_candidate(candidate, runner)
            runner = candidate
        end
    end
    return runner
end

function _selection_from_point(
    param_value::Float64,
    selected_id::Union{Nothing, Symbol},
    reason::Symbol,
    all_candidates::AbstractVector{Tuple{Symbol, BranchPoint}},
    eligible_candidates::AbstractVector{Tuple{Symbol, BranchPoint}},
    selected_point::Union{Nothing, BranchPoint},
)
    runner = _runner_up(eligible_candidates, selected_id)
    runner_id = runner === nothing ? nothing : runner[1]
    runner_pressure = runner === nothing ? NaN : runner[2].pressure
    selected_pressure = selected_point === nothing ? -Inf : selected_point.pressure
    selected_residual = selected_point === nothing ? Inf : selected_point.residual_norm
    pressure_gap = runner === nothing ? NaN : selected_pressure - runner_pressure
    return BranchSelection(
        param_value,
        selected_id,
        reason,
        length(all_candidates),
        Float64(selected_pressure),
        runner_id,
        Float64(pressure_gap),
        Float64(selected_residual),
    )
end

function _ground_state_selections(
    branches::AbstractVector{<:ContinuationBranch},
    policy::GroundStateBranchPolicy,
)
    grouped = _points_by_param(branches)
    selections = BranchSelection[]
    for param_value in sort(collect(keys(grouped)))
        all_candidates = grouped[param_value]
        eligible = Tuple{Symbol, BranchPoint}[]
        for candidate in all_candidates
            if _branch_point_success(candidate[2], policy.residual_norm_max)
                push!(eligible, candidate)
            end
        end

        if isempty(eligible)
            push!(selections, _selection_from_point(
                param_value,
                nothing,
                :no_candidate_passed_constraints,
                all_candidates,
                eligible,
                nothing,
            ))
            continue
        end

        selected = eligible[1]
        for candidate in eligible[2:end]
            if _better_ground_state_candidate(candidate, selected)
                selected = candidate
            end
        end
        runner = _runner_up(eligible, selected[1])
        reason = :pressure_max_under_constraints
        if runner !== nothing
            pressure_gap = selected[2].pressure - runner[2].pressure
            if isfinite(pressure_gap) && pressure_gap <= policy.pressure_gap_tol
                reason = :pressure_degenerate_under_constraints
            end
        end
        push!(selections, _selection_from_point(
            param_value,
            selected[1],
            reason,
            all_candidates,
            eligible,
            selected[2],
        ))
    end
    return selections
end

function _continuation_branch_selections(
    branches::AbstractVector{<:ContinuationBranch},
    policy::ContinuationBranchPolicy,
)
    target_idx = findfirst(branch -> branch.branch_id === policy.branch_id, branches)
    target_idx === nothing && throw(ArgumentError("branch_id $(policy.branch_id) was not found in continuation branches"))
    target = branches[target_idx]
    grouped = _points_by_param(branches)
    selections = BranchSelection[]
    for point in target.points
        all_candidates = get(grouped, point.param_value, Tuple{Symbol, BranchPoint}[])
        push!(selections, _selection_from_point(
            point.param_value,
            policy.branch_id,
            :continuation_branch,
            all_candidates,
            all_candidates,
            point,
        ))
    end
    return selections
end

function apply_branch_policy(branches::AbstractVector{<:ContinuationBranch}, policy::GroundStateBranchPolicy)
    return (
        branches=ContinuationBranch[branches...],
        selections=_ground_state_selections(branches, policy),
    )
end

function apply_branch_policy(branches::AbstractVector{<:ContinuationBranch}, policy::ContinuationBranchPolicy)
    return (
        branches=ContinuationBranch[branches...],
        selections=_continuation_branch_selections(branches, policy),
    )
end

function apply_branch_policy(branches::AbstractVector{<:ContinuationBranch}, ::AllBranchesPolicy)
    return (
        branches=ContinuationBranch[branches...],
        selections=BranchSelection[],
    )
end

@inline function _branch_point_from_solver_result(param_value::Float64, result::SolverResult, source::Symbol, step_index::Int)
    return BranchPoint(
        param_value,
        Float64.(result.solution),
        Float64.(result.x_state),
        Float64.(result.mu_vec),
        Float64(result.pressure),
        Float64(result.rho_norm),
        Float64(result.residual_norm),
        Bool(result.converged),
        source,
        step_index,
    )
end

@inline function _anchor_from_branch_point(anchor_id::Symbol, point::BranchPoint, candidate_count::Int)
    return AnchorRoot(
        anchor_id,
        point.source,
        point.param_value,
        copy(point.solution),
        copy(point.x_state),
        copy(point.mu_vec),
        point.pressure,
        point.residual_norm,
        point.converged,
        candidate_count,
    )
end

function _solution_distance(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})::Float64
    length(a) == length(b) || return Inf
    acc = 0.0
    for idx in eachindex(a, b)
        delta = Float64(a[idx]) - Float64(b[idx])
        acc += delta * delta
    end
    return sqrt(acc)
end

function _branch_jump_metrics(points::AbstractVector{<:BranchPoint}; tol::Float64)
    length(points) < 2 && return (count=0, max_jump=0.0)
    count = 0
    max_jump = 0.0
    for idx in 2:length(points)
        jump = _solution_distance(points[idx].solution, points[idx - 1].solution)
        isfinite(jump) || continue
        max_jump = max(max_jump, jump)
        jump > tol && (count += 1)
    end
    return (count=count, max_jump=max_jump)
end

function _path_branch_jump_metrics(branches::AbstractVector{<:ContinuationBranch}; tol::Real=DEFAULT_PATH_BRANCH_JUMP_TOL)
    tol_value = _path_positive_real(tol, :branch_jump_tol)
    count = 0
    max_jump = 0.0
    for branch in branches
        jumps = _branch_jump_metrics(branch.points; tol=tol_value)
        count += jumps.count
        max_jump = max(max_jump, jumps.max_jump)
    end
    return (count=count, max_jump=max_jump)
end

function _seed_continuation_anchor_messages(anchor_strategy::AbstractAnchorStrategy)
    if anchor_strategy isa SeedPoolAnchor
        return String["SeedContinuation uses previous-solution plus MultiSeed seed pools; no distinct anchor-root discovery is performed."]
    end
    throw(ArgumentError(
        "anchor_strategy $(typeof(anchor_strategy)) is not supported with SeedContinuation; " *
        "use SeedPoolAnchor() for the built-in seed-pool path, or keep multi-anchor discovery in the isolated PALC adapter.",
    ))
end

@inline function _path_seed_key(seed_vec::AbstractVector{<:Real})
    return join(round.(Float64.(seed_vec); digits=12), ",")
end

function _unique_path_seeds(seeds::AbstractVector{<:AbstractVector{<:Real}})
    out = Vector{Vector{Float64}}()
    seen = Set{String}()
    for seed in seeds
        normalized = Float64.(seed)
        key = _path_seed_key(normalized)
        key in seen && continue
        push!(out, normalized)
        push!(seen, key)
    end
    return out
end

function _seed_continuation_seeds(path::FixedAsymmetricRhoPath, mode::FixedAsymmetricRho, previous_solution)
    seeds = Vector{Vector{Float64}}()
    if previous_solution !== nothing
        push!(seeds, Float64.(previous_solution))
    end
    append!(seeds, get_all_seeds(MultiSeed(), [path.T_fm], mode))
    return _unique_path_seeds(seeds)
end

function _run_seed_continuation_branch(
    model::AbstractQCDModel,
    path::FixedAsymmetricRhoPath,
    continuation::SeedContinuation;
    kwargs...,
)
    points = BranchPoint[]
    failures = 0
    previous_solution = nothing

    for (idx, rho_value) in enumerate(path.rho_values)
        mode = FixedAsymmetricRho(rho_value, path.ud_ratio_target, path.s_target)
        seeds = _seed_continuation_seeds(path, mode, previous_solution)
        result = solve_multi(
            model,
            mode,
            path.T_fm;
            seeds=seeds,
            evaluate_all_attempts=continuation.evaluate_all_attempts,
            semantic_mode=:ground_state,
            xi=path.xi,
            p_num=path.p_num,
            t_num=path.t_num,
            kwargs...,
        )
        point = _branch_point_from_solver_result(rho_value, result, :seed_continuation, idx)
        push!(points, point)
        point.converged || (failures += 1)
        previous_solution = copy(result.solution)
    end

    status = failures == 0 ? :complete : :partial
    branch = ContinuationBranch(
        :seed_ground_state,
        :seed_anchor_1,
        :seed_continuation,
        points,
        status,
        (failure_count=failures, point_count=length(points)),
    )
    anchors = isempty(points) ? AnchorRoot[] : AnchorRoot[_anchor_from_branch_point(:seed_anchor_1, first(points), length(points))]
    return (branches=ContinuationBranch[branch], anchors=anchors, failure_count=failures)
end

function _unsupported_palc_error(model::AbstractQCDModel, continuation::PALCContinuation)
    _ = continuation
    return UnsupportedCapabilityError(model_kind_symbol(model), :bifurcationkit_palc_continuation)
end

function solve_path(
    model::AbstractQCDModel,
    path::FixedAsymmetricRhoPath;
    anchor_strategy::AbstractAnchorStrategy=SeedPoolAnchor(),
    continuation_strategy::AbstractContinuationStrategy=SeedContinuation(),
    branch_policy::AbstractBranchPolicy=GroundStateBranchPolicy(),
    diagnostic_level::Symbol=:summary,
    branch_jump_tol::Real=DEFAULT_PATH_BRANCH_JUMP_TOL,
    kwargs...,
)
    diagnostic_level in (:summary, :full) ||
        throw(ArgumentError("diagnostic_level must be :summary or :full, got $(diagnostic_level)"))

    if continuation_strategy isa PALCContinuation
        throw(_unsupported_palc_error(model, continuation_strategy))
    elseif !(continuation_strategy isa SeedContinuation)
        throw(ArgumentError("unsupported continuation_strategy type $(typeof(continuation_strategy))"))
    end
    messages = _seed_continuation_anchor_messages(anchor_strategy)
    jump_tol = _path_positive_real(branch_jump_tol, :branch_jump_tol)

    started = time()
    run = _run_seed_continuation_branch(model, path, continuation_strategy; kwargs...)
    selected = apply_branch_policy(run.branches, branch_policy)
    wall_time_s = time() - started
    branch_jumps = _path_branch_jump_metrics(run.branches; tol=jump_tol)
    diagnostics = PathDiagnostics(
        _strategy_symbol(continuation_strategy),
        _strategy_symbol(anchor_strategy),
        _policy_symbol(branch_policy),
        Float64(wall_time_s),
        length(run.anchors),
        sum(length(branch.points) for branch in run.branches),
        Int(run.failure_count),
        branch_jumps.count,
        messages,
    )
    return PathSolveResult(path, selected.branches, selected.selections, run.anchors, diagnostics)
end

@inline function to_namedtuple(root::AnchorRoot)
    return (
        anchor_id=root.anchor_id,
        source=root.source,
        param_value=root.param_value,
        solution=root.solution,
        x_state=root.x_state,
        mu_vec=root.mu_vec,
        pressure=root.pressure,
        residual_norm=root.residual_norm,
        converged=root.converged,
        candidate_count=root.candidate_count,
    )
end

@inline function to_namedtuple(point::BranchPoint)
    return (
        param_value=point.param_value,
        solution=point.solution,
        x_state=point.x_state,
        mu_vec=point.mu_vec,
        pressure=point.pressure,
        rho_norm=point.rho_norm,
        residual_norm=point.residual_norm,
        converged=point.converged,
        source=point.source,
        step_index=point.step_index,
    )
end

@inline function to_namedtuple(branch::ContinuationBranch)
    return (
        branch_id=branch.branch_id,
        anchor_id=branch.anchor_id,
        source=branch.source,
        points=[to_namedtuple(point) for point in branch.points],
        status=branch.status,
        diagnostics=branch.diagnostics,
    )
end

@inline function to_namedtuple(selection::BranchSelection)
    return (
        param_value=selection.param_value,
        selected_branch_id=selection.selected_branch_id,
        selection_reason=selection.selection_reason,
        candidate_branch_count=selection.candidate_branch_count,
        selected_pressure=selection.selected_pressure,
        runner_up_branch_id=selection.runner_up_branch_id,
        pressure_gap=selection.pressure_gap,
        selected_residual_norm=selection.selected_residual_norm,
    )
end

@inline function to_namedtuple(diagnostics::PathDiagnostics)
    return (
        continuation_backend=diagnostics.continuation_backend,
        anchor_strategy=diagnostics.anchor_strategy,
        branch_policy=diagnostics.branch_policy,
        wall_time_s=diagnostics.wall_time_s,
        anchor_solve_count=diagnostics.anchor_solve_count,
        continuation_step_count=diagnostics.continuation_step_count,
        failure_count=diagnostics.failure_count,
        branch_jump_count=diagnostics.branch_jump_count,
        messages=diagnostics.messages,
    )
end

@inline function to_namedtuple(path::FixedAsymmetricRhoPath)
    return (
        path_kind=:fixed_asymmetric_rho,
        T_fm=path.T_fm,
        rho_values=path.rho_values,
        ud_ratio_target=path.ud_ratio_target,
        s_target=path.s_target,
        xi=path.xi,
        p_num=path.p_num,
        t_num=path.t_num,
    )
end

@inline function to_namedtuple(result::PathSolveResult)
    return (
        path=to_namedtuple(result.path),
        branches=[to_namedtuple(branch) for branch in result.branches],
        selections=[to_namedtuple(selection) for selection in result.selections],
        anchors=[to_namedtuple(anchor) for anchor in result.anchors],
        diagnostics=to_namedtuple(result.diagnostics),
    )
end
