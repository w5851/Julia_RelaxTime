#!/usr/bin/env julia

"""Solver-free feasibility replay for the Stage-C extrema outer-sample guard.

The replay consumes the immutable numerical curves from the 2026-08-02
hybrid shadow.  It never calls an equilibrium solver.  The only route to a
Stage-C point is a feature ranked from the complete Stage-B curve; oracle
labels are consulted only after the replay for gates and reporting.
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using CSV
using JSON3
using SHA
using Statistics
using Printf

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
end

module StageCExtremaGuardFeasibility

using CSV
using JSON3
using SHA
using Statistics
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS = Main.Models

const SOURCE_RUN_ID = "30737739707"
const SOURCE_CALCULATION_SHA = "467be1fce847a9c991ec362c3335be07fccbe604"
const REVALIDATION_RUN_ID = "30730990835"
const METHODS = ("independent_oracle", "production_hybrid")
const XIS = (-0.5, 0.0, 0.5)
const CAPS = (12, 24, 48, 96, 160)
const STAGE_A_COARSE = 0.05
const STAGE_A_FINE = 0.025
const STAGE_B_COARSE = 0.0125
const STAGE_B_FINE = 0.00625
const STAGE_C_FINE = 0.003125
const POSITION_TOL = 0.025
const DENSITY_TOL = 0.0025
const AREA_TOL = 5e-5
const AREA_TOL_GOOD = 1e-4
const AREA_TOL_BAD = 5e-4
const MAXWELL_SOLVER_TOL = 5e-6
const COMPARISON_EPS = 32eps(Float64)
const AUTHOR_FIRST_ORDER = Set([
    (-0.5, 147.0947265625),
    (0.5, 106.9599609375),
    (-0.5, 20.0),
    (0.0, 5.0),
])
const CONSENSUS_MONOTONE = Set([
    (-0.5, 147.2197265625),
    (0.5, 107.0849609375),
])

const GEOMETRY_TOL = MODELS.PhaseGeometryTolerances(
    position_MeV=POSITION_TOL,
    density=DENSITY_TOL,
    maxwell_area=AREA_TOL,
)

struct Point
    xi::Float64
    T::Float64
    rho::Float64
    mu::Float64
    residual::Float64
end

struct Candidate
    rho_low::Float64
    rho_high::Float64
    drop_mu::Float64
    width::Float64
    negative_secants::Int
    level::Symbol
end

@inline _finite(value) = value isa Real && isfinite(Float64(value))

function _field(row, name::Symbol, default=nothing)
    hasproperty(row, name) || return default
    value = getproperty(row, name)
    value === missing && return default
    return value
end

function _float(value, default=NaN)
    value === nothing || value === missing ? default : try
        result = Float64(value)
        isfinite(result) ? result : default
    catch
        default
    end
end

function _bool(value, default=false)
    value === nothing || value === missing ? default :
        value isa Bool ? value : lowercase(strip(String(value))) in ("true", "1", "yes")
end

function _sha256_file(path::String)
    open(path, "r") do io
        bytes2hex(sha256(io))
    end
end

@inline function _same(a::Float64, b::Float64; atol=2e-4)
    isfinite(a) && isfinite(b) && isapprox(a, b; atol=atol, rtol=0.0)
end

@inline function _on_grid(rho::Float64, step::Float64; atol=3e-6)
    isfinite(rho) && abs(rho / step - round(rho / step)) <= atol / step
end

function _load_input(input_dir::String)
    required = ("manifest.json", "curve_points.csv", "slice_metrics.csv", "method_costs.csv")
    missing_files = [name for name in required if !isfile(joinpath(input_dir, name))]
    isempty(missing_files) || error("missing replay files: $(join(missing_files, ", "))")
    manifest = JSON3.read(read(joinpath(input_dir, "manifest.json"), String))
    curves = collect(CSV.File(joinpath(input_dir, "curve_points.csv")))
    slices = collect(CSV.File(joinpath(input_dir, "slice_metrics.csv")))
    costs = collect(CSV.File(joinpath(input_dir, "method_costs.csv")))
    isempty(curves) && error("curve_points.csv is empty")
    isempty(slices) && error("slice_metrics.csv is empty")
    isempty(costs) && error("method_costs.csv is empty")
    return (; manifest, curves, slices, costs)
end

function _validate_input(input_dir::String, data)
    manifest = data.manifest
    expected = string(_field(manifest, :expected_calculation_sha, ""))
    expected == SOURCE_CALCULATION_SHA || error("input calculation SHA mismatch: $expected")
    string(_field(manifest, :evidence_state, "")) == "final" ||
        error("input aggregate must have evidence_state=final")
    actions = _field(manifest, :actions, nothing)
    actions === nothing && error("input manifest lacks actions provenance")
    string(_field(actions, :run_id, "")) == SOURCE_RUN_ID || error("source run mismatch")
    string(_field(actions, :headSha, "")) == SOURCE_CALCULATION_SHA || error("source head SHA mismatch")
    string(_field(manifest, :source_run_id, SOURCE_RUN_ID)) == SOURCE_RUN_ID ||
        error("manifest source_run_id mismatch")

    seen = Set{Tuple{String, String, String, String}}()
    for row in data.curves
        key = (
            string(_field(row, :xi, "")), string(_field(row, :method, "")),
            string(_field(row, :T_MeV, "")), string(_field(row, :rho, "")),
        )
        key in seen && error("duplicate curve point key: $key")
        push!(seen, key)
        _bool(_field(row, :converged, false)) || error("non-converged curve point: $key")
        _bool(_field(row, :finite, false)) || error("non-finite curve point: $key")
        _finite(_float(_field(row, :rho))) || error("invalid rho: $key")
        _finite(_float(_field(row, :muq_MeV))) || error("invalid mu: $key")
    end
    return true
end

function _point_rows(rows, method::String, xi::Float64, T::Float64, predicate=(rho -> true))
    selected = Dict{Float64, Point}()
    for row in rows
        string(_field(row, :method, "")) == method || continue
        x = _float(_field(row, :xi))
        temp = _float(_field(row, :T_MeV))
        rho = _float(_field(row, :rho))
        mu = _float(_field(row, :muq_MeV))
        residual = _float(_field(row, :residual_norm), Inf)
        (_same(x, xi; atol=1e-8) && _same(temp, T) && predicate(rho)) || continue
        (_finite(rho) && _finite(mu)) || continue
        point = Point(xi, T, rho, mu, residual)
        if !haskey(selected, rho) || point.residual < selected[rho].residual
            selected[rho] = point
        end
    end
    sort!(collect(values(selected)); by=point -> point.rho)
end

function _curve(points::Vector{Point})
    length(points) >= 6 || return nothing
    ordered = sort(copy(points); by=point -> point.rho)
    length(unique(point.rho for point in ordered)) == length(ordered) || return nothing
    all(point -> _finite(point.rho) && _finite(point.mu), ordered) || return nothing
    return (rho=[point.rho for point in ordered], mu=[point.mu for point in ordered], points=ordered)
end

function _evaluate(curve)
    curve === nothing && return (
        status=:invalid, reason="solver_or_curve_failure", mu_transition=nothing,
        rho_hadron=nothing, rho_quark=nothing, area_residual=Inf,
        sres=MODELS.SShapeResult(), maxwell=MODELS.MaxwellResult(), curve=nothing,
    )
    result = MODELS._classify_s_curve(
        curve.mu, curve.rho;
        maxwell_options=(; tol_area=MAXWELL_SOLVER_TOL),
        area_tol_good=AREA_TOL_GOOD, area_tol_bad=AREA_TOL_BAD,
    )
    maxwell = result.maxwell
    area = result.area_residual === nothing ?
        (maxwell.converged ? maxwell.area_residual : Inf) : result.area_residual
    return (
        status=Symbol(result.status), reason=String(result.reason),
        mu_transition=result.mu_transition,
        rho_hadron=maxwell.converged ? maxwell.rho_hadron : nothing,
        rho_quark=maxwell.converged ? maxwell.rho_quark : nothing,
        area_residual=area, sres=result.sres, maxwell=maxwell, curve=curve,
    )
end

function _geometry(left, right)
    MODELS._compare_phase_geometry(left, right, GEOMETRY_TOL)
end

function _semantic(left, right, geometry)
    if left.status == :invalid && right.status == :invalid &&
            left.reason == "no_s_shape" && right.reason == "no_s_shape" && geometry.converged
        return :confirmed_monotone
    elseif left.status == :valid && right.status == :valid && geometry.converged
        return :confirmed_first_order
    end
    return :ambiguous_near_critical
end

function _segments(curve)
    curve === nothing && return (Int[], Int[])
    signs = Int[]
    indices = Int[]
    for index in 1:(length(curve.points) - 1)
        drho = curve.points[index + 1].rho - curve.points[index].rho
        abs(drho) <= eps(Float64) && continue
        slope = (curve.points[index + 1].mu - curve.points[index].mu) / drho
        sign = slope > 0 ? 1 : slope < 0 ? -1 : 0
        sign == 0 && continue
        if isempty(signs) || sign != last(signs)
            push!(signs, sign)
            push!(indices, index)
        end
    end
    return signs, indices
end

function _topology(curve)
    signs, _ = _segments(curve)
    signs
end

function _candidates(curve; level::Symbol=:unknown)
    curve === nothing && return Candidate[]
    signs, indices = _segments(curve)
    length(signs) < 3 && return Candidate[]
    candidates = Candidate[]
    for run in 2:(length(signs) - 1)
        signs[run] == -1 && signs[run - 1] == 1 && signs[run + 1] == 1 || continue
        first_index = indices[run]
        last_index = run < length(indices) ? indices[run + 1] - 1 : length(curve.points) - 1
        first_index < last_index + 1 || continue
        low = curve.points[first_index].rho
        high = curve.points[last_index + 1].rho
        drop = curve.points[first_index].mu - curve.points[last_index + 1].mu
        width = high - low
        negative_secants = max(1, last_index - first_index + 1)
        drop > 0 && width > 0 && push!(candidates, Candidate(low, high, drop, width, negative_secants, level))
    end
    candidates
end

function _unique_topology(curve)
    signs = _topology(curve)
    signs == [1, -1, 1] && length(_candidates(curve)) == 1
end

function _stable_candidates(coarse, fine)
    cands_coarse = _candidates(coarse; level=:coarse)
    cands_fine = _candidates(fine; level=:fine)
    if length(cands_coarse) != 1 || length(cands_fine) != 1
        return false, vcat(cands_coarse, cands_fine)
    end
    c = cands_coarse[1]
    f = cands_fine[1]
    overlap = c.rho_low <= f.rho_high + STAGE_B_COARSE &&
        f.rho_low <= c.rho_high + STAGE_B_COARSE
    overlap, vcat(cands_coarse, cands_fine)
end

"""Return the first strict outer Stage-B samples around the extrema.

The comparison is deliberately strict.  Samples whose chemical potential is
equal to an extremum are skipped; no interpolation, binary search, or fixed
padding is used.
"""
function _extrema_outer_guard(curve, evaluation)
    curve === nothing && return (status=:ambiguous_near_critical, reason="missing_stage_b_curve")
    evaluation.sres.has_s_shape || return (status=:ambiguous_near_critical, reason="no_unique_s_shape")
    _unique_topology(curve) || return (status=:ambiguous_near_critical, reason="unstable_or_multiple_s_topology")
    rho_spinodal = sort(Float64[
        something(evaluation.sres.rho_spinodal_hadron, NaN),
        something(evaluation.sres.rho_spinodal_quark, NaN),
    ])
    mu_spinodal = Float64[
        something(evaluation.sres.mu_spinodal_hadron, NaN),
        something(evaluation.sres.mu_spinodal_quark, NaN),
    ]
    all(isfinite, rho_spinodal) && all(isfinite, mu_spinodal) ||
        return (status=:ambiguous_near_critical, reason="missing_spinodal_extrema")
    mu_low, mu_high = extrema(mu_spinodal)
    left = sort([point for point in curve.points if
        point.rho < rho_spinodal[1] && point.mu < mu_low - COMPARISON_EPS]; by=point -> point.rho)
    right = sort([point for point in curve.points if
        point.rho > rho_spinodal[2] && point.mu > mu_high + COMPARISON_EPS]; by=point -> point.rho)
    isempty(left) && return (status=:ambiguous_near_critical, reason="missing_left_strict_outer_sample",
        mu_low=mu_low, mu_high=mu_high, rho_spinodal_low=rho_spinodal[1], rho_spinodal_high=rho_spinodal[2])
    isempty(right) && return (status=:ambiguous_near_critical, reason="missing_right_strict_outer_sample",
        mu_low=mu_low, mu_high=mu_high, rho_spinodal_low=rho_spinodal[1], rho_spinodal_high=rho_spinodal[2])
    low_point = last(left)
    high_point = first(right)
    low_point.rho < rho_spinodal[1] && high_point.rho > rho_spinodal[2] ||
        return (status=:ambiguous_near_critical, reason="invalid_guard_order")
    (
        status=:ok, reason="first_strict_outer_samples", guard_low=low_point.rho,
        guard_high=high_point.rho, guard_low_mu=low_point.mu, guard_high_mu=high_point.mu,
        mu_low=mu_low, mu_high=mu_high,
        rho_spinodal_low=rho_spinodal[1], rho_spinodal_high=rho_spinodal[2],
        low_strict_margin=mu_low - low_point.mu, high_strict_margin=high_point.mu - mu_high,
    )
end

function _crossing_targets(curve, target_mu)
    targets = Float64[]
    curve === nothing && return targets
    for index in 1:(length(curve.points) - 1)
        a, b = curve.points[index], curve.points[index + 1]
        (a.mu - target_mu) * (b.mu - target_mu) <= 0 &&
            push!(targets, 0.5 * (a.rho + b.rho))
    end
    targets
end

function _feature_targets(curve, evaluation, guard)
    guard.status == :ok || return Float64[]
    targets = Float64[guard.rho_spinodal_low, guard.rho_spinodal_high]
    for value in (evaluation.rho_hadron, evaluation.rho_quark)
        value === nothing || push!(targets, Float64(value))
    end
    evaluation.mu_transition === nothing || append!(targets, _crossing_targets(curve, evaluation.mu_transition))

    # Rank high curvature and area-contribution segments from Stage-B only.
    scored = NamedTuple[]
    for index in 2:(length(curve.points) - 1)
        left, middle, right = curve.points[index - 1], curve.points[index], curve.points[index + 1]
        hleft = middle.rho - left.rho
        hright = right.rho - middle.rho
        hleft > 0 && hright > 0 || continue
        slope_left = (middle.mu - left.mu) / hleft
        slope_right = (right.mu - middle.mu) / hright
        curvature = abs(slope_right - slope_left) / max(0.5 * (hleft + hright), eps(Float64))
        area_weight = evaluation.mu_transition === nothing ? 0.0 :
            abs((middle.mu - Float64(evaluation.mu_transition)) * 0.5 * (hleft + hright))
        push!(scored, (rho=middle.rho, curvature=curvature, area=area_weight))
    end
    sort!(scored; by=item -> (-item.curvature, -item.area, item.rho))
    append!(targets, item.rho for item in Iterators.take(scored, min(12, length(scored))))
    sort!(unique(filter(rho -> guard.guard_low <= rho <= guard.guard_high, targets)))
end

function _select_stage_c_points(pool, curve, evaluation, guard, cap::Int)
    cap > 0 || return Point[]
    guard.status == :ok || return Point[]
    targets = _feature_targets(curve, evaluation, guard)
    isempty(targets) && return Point[]
    stage_b_rho = Set(point.rho for point in curve.points)
    eligible = [point for point in pool if
        guard.guard_low <= point.rho <= guard.guard_high && !(point.rho in stage_b_rho)]
    isempty(eligible) && return Point[]
    scored = NamedTuple[]
    for point in eligible
        distance = minimum(abs(point.rho - target) for target in targets)
        # A deterministic distance-first ranking keeps the route explainable;
        # rho is the final tie-break and no oracle label is read here.
        push!(scored, (point=point, distance=distance))
    end
    sort!(scored; by=item -> (item.distance, item.point.rho))
    [item.point for item in Iterators.take(scored, min(cap, length(scored)))]
end

function _row_points(points)
    Dict(point.rho => point for point in points)
end

function _union_count(groups...)
    keys = Set{Float64}()
    for group in groups
        for point in group
            push!(keys, point.rho)
        end
    end
    length(keys)
end

function _anchor_rows(slices, xi::Float64, T::Float64)
    Dict(string(_field(row, :method, "")) => row for row in slices if
        _same(_float(_field(row, :xi)), xi; atol=1e-8) && _same(_float(_field(row, :T_MeV)), T))
end

function _oracle_status(slices, xi, T)
    row = get(_anchor_rows(slices, xi, T), "independent_oracle", nothing)
    row === nothing ? "missing" : string(_field(row, :result_status, "missing"))
end

function _stage_a(curves, xi, T)
    coarse_points = _point_rows(curves, "production_hybrid", xi, T,
        rho -> _on_grid(rho, STAGE_A_COARSE))
    fine_points = _point_rows(curves, "production_hybrid", xi, T,
        rho -> _on_grid(rho, STAGE_A_FINE))
    coarse = _evaluate(_curve(coarse_points))
    fine = _evaluate(_curve(fine_points))
    geometry = _geometry(coarse, fine)
    (coarse_points=coarse_points, fine_points=fine_points, coarse=coarse, fine=fine,
        geometry=geometry, status=_semantic(coarse, fine, geometry))
end

function _stage_b(curves, xi, T)
    coarse_points = _point_rows(curves, "independent_oracle", xi, T,
        rho -> _on_grid(rho, STAGE_B_COARSE))
    fine_points = _point_rows(curves, "independent_oracle", xi, T,
        rho -> _on_grid(rho, STAGE_B_FINE))
    coarse = _evaluate(_curve(coarse_points))
    fine = _evaluate(_curve(fine_points))
    geometry = _geometry(coarse, fine)
    (coarse_points=coarse_points, fine_points=fine_points, coarse=coarse, fine=fine,
        geometry=geometry, status=_semantic(coarse, fine, geometry))
end

function _replay_anchor(curves, slices, xi::Float64, T::Float64, cap::Int)
    oracle_status = _oracle_status(slices, xi, T)
    stage_a = _stage_a(curves, xi, T)
    if stage_a.status == :confirmed_monotone
        return (
            xi=xi, T_MeV=T, cap=cap, oracle_status=oracle_status,
            result_status=:confirmed_monotone, reason="stage_a_monotone_certificate",
            stage_used=:stage_a, stage_a=stage_a, stage_b=nothing, stage_c=nothing,
            guard=(status=:not_run, reason="stage_a_certificate"), selected_points=Point[],
            candidates=Candidate[], stable_candidate=true,
            unique_solves=_union_count(stage_a.coarse_points, stage_a.fine_points),
        )
    end

    stage_b = _stage_b(curves, xi, T)
    if stage_b.status in (:confirmed_monotone, :confirmed_first_order)
        return (
            xi=xi, T_MeV=T, cap=cap, oracle_status=oracle_status,
            result_status=stage_b.status, reason="stage_b_certificate",
            stage_used=:stage_b, stage_a=stage_a, stage_b=stage_b, stage_c=nothing,
            guard=(status=:not_run, reason="stage_b_certificate"), selected_points=Point[],
            candidates=_candidates(stage_b.fine.curve; level=:fine), stable_candidate=true,
            unique_solves=_union_count(stage_a.coarse_points, stage_a.fine_points, stage_b.fine_points),
        )
    end

    stable_candidate, candidates = _stable_candidates(stage_b.coarse.curve, stage_b.fine.curve)
    guard = _extrema_outer_guard(stage_b.fine.curve, stage_b.fine)
    pool = _point_rows(curves, "independent_oracle", xi, T,
        rho -> _on_grid(rho, STAGE_C_FINE) && !_on_grid(rho, STAGE_B_FINE))
    selected = _select_stage_c_points(pool, stage_b.fine.curve, stage_b.fine, guard, cap)
    merged_points = vcat(stage_b.fine_points, selected)
    merged = _evaluate(_curve(merged_points))
    geometry = _geometry(stage_b.fine, merged)
    valid = stable_candidate && stage_b.fine.status == :valid && merged.status == :valid && geometry.converged
    status = valid ? :confirmed_first_order : :ambiguous_near_critical
    reason = valid ? "extrema_guard_stage_c_certificate" :
        !stable_candidate ? "unstable_or_multiple_s_topology" :
        guard.status != :ok ? guard.reason :
        stage_b.fine.status != :valid ? "stage_b_unresolved" :
        "stage_c_geometry_unresolved"
    (
        xi=xi, T_MeV=T, cap=cap, oracle_status=oracle_status,
        result_status=status, reason=reason, stage_used=:stage_c,
        stage_a=stage_a, stage_b=stage_b,
        stage_c=(result=merged, geometry=geometry), guard=guard,
        selected_points=selected, candidates=candidates,
        stable_candidate=stable_candidate,
        unique_solves=_union_count(stage_a.coarse_points, stage_a.fine_points,
            stage_b.fine_points, selected),
    )
end

function _author_expected(xi, T)
    (xi, T) in AUTHOR_FIRST_ORDER ? "confirmed_first_order" : nothing
end

function _matches_gate(result)
    expected = _author_expected(result.xi, result.T_MeV)
    expected !== nothing && return String(result.result_status) == expected
    oracle = result.oracle_status
    oracle == "ambiguous_near_critical" && return result.result_status == :ambiguous_near_critical
    oracle in ("confirmed_first_order", "confirmed_monotone") && return String(result.result_status) == oracle
    result.result_status == Symbol(oracle)
end

function _unsupported_confirmation(result)
    (result.xi, result.T_MeV) in AUTHOR_FIRST_ORDER && return false
    result.oracle_status == "ambiguous_near_critical" && result.result_status != :ambiguous_near_critical
end

function _frontier(curves, slices, costs, anchors)
    dense = sum(_float(_field(row, :unique_solves), 0.0) for row in costs if
        string(_field(row, :method, "")) == "memoized_dense")
    frontier = NamedTuple[]
    replays = Dict{Int, Vector{Any}}()
    for cap in CAPS
        results = [_replay_anchor(curves, slices, xi, T, cap) for (xi, T) in anchors]
        replays[cap] = results
        mismatches = count(result -> !_matches_gate(result), results)
        unsupported = count(_unsupported_confirmation, results)
        multiple = count(result -> !result.stable_candidate, results)
        guard_failures = count(result -> result.stage_used == :stage_c && result.guard.status != :ok, results)
        unique_solves = sum(result.unique_solves for result in results)
        state_gate = mismatches == 0 && unsupported == 0
        candidate_gate = multiple == 0
        cost_gate = unique_solves <= dense
        push!(frontier, (
            cap=cap, anchors=length(anchors), classification_mismatches=mismatches,
            unsupported_confirmations=unsupported, multiple_candidate_anchors=multiple,
            guard_failures=guard_failures, unique_solves=unique_solves,
            dense_unique_solves=dense, state_gate=state_gate,
            candidate_gate=candidate_gate, cost_gate=cost_gate,
            feasible=state_gate && candidate_gate && cost_gate,
        ))
    end
    selected = findfirst(row -> row.feasible, frontier)
    (; dense, frontier, replays, selected)
end

function _anchors(slices)
    values = Set{Tuple{Float64, Float64}}()
    for row in slices
        string(_field(row, :method, "")) == "independent_oracle" || continue
        xi = _float(_field(row, :xi)); T = _float(_field(row, :T_MeV))
        _finite(xi) && _finite(T) && push!(values, (xi, T))
    end
    sort!(collect(values); by=item -> (item[1], item[2]))
end

function _write_csv(path, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _json_pretty(path, value)
    open(path, "w") do io
        JSON3.pretty(io, value)
        write(io, '\n')
    end
end

function _result_rows(results)
    rows = NamedTuple[]
    for result in results
        guard = result.guard
        push!(rows, (
            xi=result.xi, T_MeV=result.T_MeV, cap=result.cap,
            oracle_status=result.oracle_status, result_status=String(result.result_status),
            reason=result.reason, stage_used=String(result.stage_used),
            stage_a_status=String(result.stage_a.status),
            stage_b_status=result.stage_b === nothing ? "not_run" : String(result.stage_b.status),
            stage_c_status=result.stage_c === nothing ? "not_run" : String(result.result_status),
            guard_status=String(guard.status), guard_reason=String(guard.reason),
            guard_low=get(guard, :guard_low, missing), guard_high=get(guard, :guard_high, missing),
            support_mu_low=get(guard, :mu_low, missing), support_mu_high=get(guard, :mu_high, missing),
            stage_c_point_count=length(result.selected_points), unique_solves=result.unique_solves,
            stable_candidate=result.stable_candidate,
        ))
    end
    rows
end

function _guard_rows(results)
    rows = NamedTuple[]
    for result in results
        guard = result.guard
        push!(rows, (
            xi=result.xi, T_MeV=result.T_MeV, cap=result.cap,
            status=String(guard.status), reason=String(guard.reason),
            guard_low=get(guard, :guard_low, missing), guard_high=get(guard, :guard_high, missing),
            guard_low_mu=get(guard, :guard_low_mu, missing), guard_high_mu=get(guard, :guard_high_mu, missing),
            mu_low=get(guard, :mu_low, missing), mu_high=get(guard, :mu_high, missing),
            rho_spinodal_low=get(guard, :rho_spinodal_low, missing),
            rho_spinodal_high=get(guard, :rho_spinodal_high, missing),
            low_strict_margin=get(guard, :low_strict_margin, missing),
            high_strict_margin=get(guard, :high_strict_margin, missing),
            strict_order=(get(guard, :guard_low, -Inf) < get(guard, :guard_high, Inf)),
        ))
    end
    rows
end

function _candidate_rows(results)
    rows = NamedTuple[]
    for result in results, candidate in result.candidates
        push!(rows, (
            xi=result.xi, T_MeV=result.T_MeV, cap=result.cap,
            level=String(candidate.level), rho_low=candidate.rho_low,
            rho_high=candidate.rho_high, width=candidate.width,
            drop_mu=candidate.drop_mu, negative_secants=candidate.negative_secants,
        ))
    end
    rows
end

function _selected_rows(results)
    rows = NamedTuple[]
    for result in results, point in result.selected_points
        push!(rows, (xi=result.xi, T_MeV=result.T_MeV, cap=result.cap,
            rho=point.rho, muq_MeV=point.mu, residual_norm=point.residual,
            source="stage_c_guard_pool"))
    end
    rows
end

function _write_docs(output_dir, policy, frontier, input_dir, revalidation_dir)
    write(joinpath(output_dir, "README.md"), """# PNJL Stage-C discrete-extrema guard feasibility v1

verdict: `$(policy.verdict)`。

本目录是对 source run `$(SOURCE_RUN_ID)`（calculation SHA
`$(SOURCE_CALCULATION_SHA)`）的 solver-free replay。Stage-B 使用完整
`0.00625` 全域曲线；Stage-C 只从 Stage-B 派生的特征排序中选择 guard 内的
`0.003125` 已有曲线点。没有调用 equilibrium solver，也没有修改 production、
reference 或历史 v1/v2 evidence。

Guard 定义为两个 μ 极值外侧的首个严格采样点：相等点跳过，不插值、不二分、
不使用固定 padding。缺少任一外侧点、多 S 形或 topology 不稳定时保持
`ambiguous_near_critical`。

- tested caps: `$(join(string.(CAPS), ", "))`
- dense unique-solve reference: `$(frontier.dense)`
- selected cap: `$(something(policy.selected_cap, "none"))`
- solver called: `false`
- five-point revalidation run: `$(REVALIDATION_RUN_ID)` (`$(revalidation_dir === nothing ? "not supplied" : "recorded")`)

`author_adjudication` 只作为显式 provenance 使用，不把其它视觉判断升级为自动
证书。即使 feasibility gate 通过，它也只授权下一步 production focused CI，
不等于 shadow、reference 或正式 production 已通过。
""")
    write(joinpath(output_dir, "AUDIT.md"), """# Stage-C extrema-guard feasibility audit

输入 manifest、curve、source run、calculation SHA 和 replay producer SHA 均写入
`manifest.json`。所有 S-shape、Maxwell 和 geometry 结果由当前 Julia PhaseCore
重新计算；历史 oracle 标签只在 replay 完成后用于 gate。Stage-C 分类曲线始终是
完整 Stage-B 全域曲线与选定局部点的并集。

成本按 Stage-A、完整 Stage-B 和选定 Stage-C rho key 的并集计算，Stage-B 网格不
免费消费。离线 replay 不能证明 residual/Jacobian、fallback/retry 或 runner-minutes；
这些必须由后续 Actions targeted/full shadow 验证。source Actions 的 aggregate
physical gate 失败被保留为诊断事实，但所有输入曲线均需 finite/converged。
""")
end

function _make_manifest(output_dir, policy, frontier, data, input_dir, revalidation_dir)
    files = Dict{String, String}()
    for (root, _, names) in walkdir(output_dir)
        for name in names
            name == "manifest.json" && continue
            path = joinpath(root, name)
            files[replace(relpath(path, output_dir), '\\' => '/')] = _sha256_file(path)
        end
    end
    source_manifest_sha = _sha256_file(joinpath(input_dir, "manifest.json"))
    producer_script_sha = _sha256_file(joinpath(PROJECT_ROOT, "scripts", "analysis", "pnjl_cep_hybrid_extrema_guard_feasibility.jl"))
    manifest = Dict(
        "schema_version" => "cep_hybrid_stagec_extrema_guard_feasibility_v1",
        "verdict" => policy.verdict,
        "source_run_id" => SOURCE_RUN_ID,
        "source_calculation_sha" => SOURCE_CALCULATION_SHA,
        "source_manifest_sha256" => source_manifest_sha,
        "revalidation_run_id" => REVALIDATION_RUN_ID,
        "revalidation_artifact_supplied" => revalidation_dir !== nothing,
        "producer_head_sha" => try readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`) catch; "unknown" end,
        "producer_script_sha256" => producer_script_sha,
        "solver_called" => false,
        "source_actions_conclusion" => string(_field(_field(data.manifest, :actions, nothing), :conclusion, "unknown")),
        "selected_policy" => policy,
        "cost_frontier" => frontier.frontier,
        "input_files" => Dict(
            "manifest.json" => source_manifest_sha,
            "curve_points.csv" => _sha256_file(joinpath(input_dir, "curve_points.csv")),
            "slice_metrics.csv" => _sha256_file(joinpath(input_dir, "slice_metrics.csv")),
            "method_costs.csv" => _sha256_file(joinpath(input_dir, "method_costs.csv")),
        ),
        "files" => files,
        "gate" => Dict(
            "automatic_gate_is_not_promotion" => true,
            "author_first_order" => [Dict("xi" => xi, "T_MeV" => T) for (xi, T) in AUTHOR_FIRST_ORDER],
            "consensus_monotone" => [Dict("xi" => xi, "T_MeV" => T) for (xi, T) in CONSENSUS_MONOTONE],
        ),
    )
    _json_pretty(joinpath(output_dir, "manifest.json"), manifest)
    manifest
end

function run(input_dir::String, output_dir::String; revalidation_dir::Union{Nothing, String}=nothing)
    data = _load_input(input_dir)
    _validate_input(input_dir, data)
    anchors = _anchors(data.slices)
    length(anchors) == 24 || error("expected 24 independent-oracle anchors, got $(length(anchors))")
    frontier = _frontier(data.curves, data.slices, data.costs, anchors)
    selected = frontier.selected
    baseline_results = frontier.replays[first(CAPS)]
    all_results = reduce(vcat, (frontier.replays[cap] for cap in CAPS))
    unknown_oracle = any(result -> result.oracle_status == "ambiguous_near_critical" &&
        !_unsupported_confirmation(result) && result.result_status != :ambiguous_near_critical,
        baseline_results)
    author_missing = any(result -> (result.xi, result.T_MeV) in AUTHOR_FIRST_ORDER &&
        result.result_status != :confirmed_first_order, baseline_results)
    verdict = if selected !== nothing
        "feasible_candidate"
    elseif unknown_oracle || author_missing
        "oracle_inconclusive"
    elseif any(row -> row.multiple_candidate_anchors > 0, frontier.frontier)
        "maxwell_candidate_inconclusive"
    elseif any(row -> !row.cost_gate, frontier.frontier)
        "performance_inconclusive"
    else
        "integration_failed"
    end
    policy = (
        schema_version="cep_hybrid_stagec_extrema_guard_feasibility_v1",
        verdict=verdict,
        selected_cap=selected === nothing ? nothing : frontier.frontier[selected].cap,
        guard_rule="extrema_outer_samples_v1",
        comparison_epsilon=COMPARISON_EPS,
        stage_b_grid_step=STAGE_B_FINE,
        stage_c_local_step=STAGE_C_FINE,
        caps=CAPS,
        dense_unique_solves=frontier.dense,
        reason=selected === nothing ? "no cap satisfies all feasibility gates" : "minimum feasible cap",
    )
    mkpath(output_dir)
    _write_csv(joinpath(output_dir, "anchor_replay.csv"), _result_rows(all_results))
    _write_csv(joinpath(output_dir, "guard_table.csv"), _guard_rows(all_results))
    _write_csv(joinpath(output_dir, "candidate_runs.csv"), _candidate_rows(all_results))
    _write_csv(joinpath(output_dir, "selected_points.csv"), _selected_rows(baseline_results))
    _write_csv(joinpath(output_dir, "cost_frontier.csv"), frontier.frontier)
    deep = [(xi=result.xi, T_MeV=result.T_MeV, oracle_status=result.oracle_status,
        hybrid_status=String(result.result_status), reason="oracle_ambiguous_or_guard_inconclusive")
        for result in baseline_results if result.oracle_status == "ambiguous_near_critical" &&
        result.result_status == :ambiguous_near_critical]
    _write_csv(joinpath(output_dir, "deep_oracle_required.csv"), deep)
    _write_csv(joinpath(output_dir, "author_adjudications.csv"), [
        (xi=xi, T_MeV=T, expected_status="confirmed_first_order", label_source="author_adjudication")
        for (xi, T) in AUTHOR_FIRST_ORDER
    ])
    _write_csv(joinpath(output_dir, "claim_ledger.csv"), [
        (claim_id="guard", claim="guard uses first strict Stage-B samples outside both μ extrema", status="observed", boundary="solver-free replay"),
        (claim_id="classification", claim="hybrid classification matches oracle non-ambiguous anchors and explicit author adjudications", status=verdict == "feasible_candidate" ? "pass" : "inconclusive", boundary="requires targeted/full shadow"),
        (claim_id="geometry", claim="Maxwell and rho geometry are recomputed with Julia PhaseCore", status="observed", boundary="discrete curve evidence only"),
        (claim_id="cost", claim="Stage-A/B/C union unique solves do not exceed dense", status=all(row -> row.cost_gate, frontier.frontier) ? "pass" : "inconclusive", boundary="runner telemetry requires Actions"),
        (claim_id="promotion", claim="feasibility authorizes reference or production promotion", status="not_claimed", boundary="not part of this replay"),
    ])
    _write_csv(joinpath(output_dir, "curve_index.csv"), [
        (source="aggregate_replay", path="curve_points.csv", sha256=_sha256_file(joinpath(input_dir, "curve_points.csv")), raw_curve_copy_in_repository=false),
    ])
    _json_pretty(joinpath(output_dir, "selected_policy.json"), policy)
    _json_pretty(joinpath(output_dir, "provenance.json"), Dict(
        "source_run_id" => SOURCE_RUN_ID,
        "source_calculation_sha" => SOURCE_CALCULATION_SHA,
        "revalidation_run_id" => REVALIDATION_RUN_ID,
        "revalidation_artifact" => revalidation_dir,
        "source_manifest_sha256" => _sha256_file(joinpath(input_dir, "manifest.json")),
        "producer_script_sha256" => _sha256_file(joinpath(PROJECT_ROOT, "scripts", "analysis", "pnjl_cep_hybrid_extrema_guard_feasibility.jl")),
        "solver_called" => false,
        "oracle_label_leakage" => false,
    ))
    _json_pretty(joinpath(output_dir, "plot_manifest.json"), Dict(
        "schema_version" => "cep_hybrid_stagec_extrema_guard_feasibility_v1",
        "figures" => Any[], "raw_curve_artifact" => "external_actions_or_local_artifact",
    ))
    _write_docs(output_dir, policy, frontier, input_dir, revalidation_dir)
    _make_manifest(output_dir, policy, frontier, data, input_dir, revalidation_dir)
    policy
end

function _parse_args(args)
    input_dir = nothing
    output_dir = joinpath(PROJECT_ROOT, "docs", "analysis", "pnjl", "cep_maxwell", "stagec", "cep_hybrid_stagec_extrema_guard_feasibility_v1")
    revalidation_dir = nothing
    for arg in args
        startswith(arg, "--input-dir=") && (input_dir = abspath(split(arg, "="; limit=2)[2]))
        startswith(arg, "--output-dir=") && (output_dir = abspath(split(arg, "="; limit=2)[2]))
        startswith(arg, "--revalidation-dir=") && (revalidation_dir = abspath(split(arg, "="; limit=2)[2]))
        arg in ("-h", "--help") && return nothing
    end
    input_dir === nothing && throw(ArgumentError("--input-dir=PATH is required"))
    (; input_dir, output_dir, revalidation_dir)
end

function main(args=ARGS)
    options = _parse_args(args)
    options === nothing && begin
        println("Usage: julia --project=. scripts/analysis/pnjl_cep_hybrid_extrema_guard_feasibility.jl --input-dir=PATH [--output-dir=PATH] [--revalidation-dir=PATH]")
        return 0
    end
    policy = run(options.input_dir, options.output_dir; revalidation_dir=options.revalidation_dir)
    println(JSON3.write(policy))
    policy.verdict == "feasible_candidate" ? 0 : 2
end

end # module StageCExtremaGuardFeasibility

if abspath(PROGRAM_FILE) == @__FILE__
    exit(StageCExtremaGuardFeasibility.main())
end
