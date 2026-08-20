#!/usr/bin/env julia

"""Solver-free Stage-C certificate feasibility replay (v2).

This replay consumes an immutable hybrid-shadow aggregate artifact.  It does
not call an equilibrium solver.  Stage-B and Stage-C semantic curves are
reconstructed from one independent-oracle session, while the existing Julia
phase classifier, Maxwell construction, and geometry comparison are used for
the actual certificate.  Oracle status/geometry are used only for comparison
and gates; they never choose verification points.
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
const MODELS = Main.Models

module StageCCertificateFeasibilityV2

using CSV
using JSON3
using SHA
using Statistics
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS = Main.Models
const SOURCE_RUN_ID = "30601857594"
const SOURCE_CALCULATION_SHA = "fde2b929b60575f1daacb84a1b9b8ff6e3b0a0cc"
const METHODS = ("independent_oracle", "production_hybrid")
const XIS = (-0.5, 0.0, 0.5)
const CAPS = (12, 24, 48, 96, 160)
const POSITION_TOL = 0.025
const DENSITY_TOL = 0.0025
const AREA_TOL = 5e-5
const CEP_AREA_TOL_GOOD = 1e-4
const CEP_AREA_TOL_BAD = 5e-4
const MAXWELL_SOLVER_TOL_FACTOR = 0.1
const MAXWELL_SOLVER_TOL = MAXWELL_SOLVER_TOL_FACTOR * min(CEP_AREA_TOL_GOOD, AREA_TOL)
const STAGE_A_COARSE = 0.05
const STAGE_A_FINE = 0.025
const STAGE_B_COARSE = 0.0125
const STAGE_B_FINE = 0.00625
const STAGE_C_FINE = 0.003125
const FEATURE_RADIUS = 0.025
const AUTHOR_FIRST_ORDER = Set([(-0.5, 147.0947265625), (0.5, 106.9599609375)])
const CONSENSUS_MONOTONE = Set([(-0.5, 147.2197265625), (0.5, 107.0849609375)])

const TOL = MODELS.PhaseGeometryTolerances(
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

@inline _finite(x) = x isa Real && isfinite(Float64(x))

function _field(row, name::Symbol, default=nothing)
    hasproperty(row, name) || return default
    value = getproperty(row, name)
    value === missing && return default
    return value
end

function _float(value, default=NaN)
    value === nothing && return default
    value === missing && return default
    try
        result = Float64(value)
        return isfinite(result) ? result : default
    catch
        return default
    end
end

function _bool(value, default=false)
    value === nothing && return default
    value === missing && return default
    value isa Bool && return value
    lowercase(strip(String(value))) in ("true", "1", "yes")
end

function _sha256_file(path::String)
    bytes2hex(open(sha256, path))
end

function _on_grid(rho::Float64, step::Float64)
    isapprox(rho / step, round(rho / step); atol=3e-7, rtol=0.0)
end

function _same(value::Float64, target::Float64; atol::Float64=2e-4)
    isapprox(value, target; atol=atol, rtol=0.0)
end

function _manifest_value(manifest, name::Symbol, default=nothing)
    try
        value = getproperty(manifest, name)
        value === nothing && return default
        return value
    catch
        return default
    end
end

function _load_rows(input_dir::String)
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
    expected = string(_manifest_value(manifest, :expected_calculation_sha, ""))
    expected == SOURCE_CALCULATION_SHA || error("input calculation SHA mismatch: $expected")
    string(_manifest_value(manifest, :evidence_state, "")) == "final" ||
        error("input aggregate must be final")
    actions = _manifest_value(manifest, :actions, nothing)
    actions === nothing || Bool(_manifest_value(actions, :source_run_completed_success, false)) ||
        error("source Actions run is not completed successfully")
    run_id = string(_manifest_value(actions, :run_id, SOURCE_RUN_ID))
    run_id == SOURCE_RUN_ID || error("source run mismatch: $run_id")

    seen = Set{Tuple{String, String, String, String}}()
    for row in data.curves
        key = (
            string(_field(row, :xi, "")),
            string(_field(row, :method, "")),
            string(_field(row, :T_MeV, "")),
            string(_field(row, :rho, "")),
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

function _point_rows(rows, method::String, xi::Float64, T::Float64, predicate)
    selected = Dict{Float64, Point}()
    for row in rows
        string(_field(row, :method, "")) == method || continue
        x = _float(_field(row, :xi))
        temp = _float(_field(row, :T_MeV))
        rho = _float(_field(row, :rho))
        mu = _float(_field(row, :muq_MeV))
        residual = _float(_field(row, :residual_norm), Inf)
        _same(x, xi; atol=1e-8) && _same(temp, T) || continue
        predicate(rho) || continue
        (_finite(rho) && _finite(mu)) || continue
        candidate = Point(xi, T, rho, mu, residual)
        if !haskey(selected, rho) || residual < selected[rho].residual
            selected[rho] = candidate
        end
    end
    sort!(collect(values(selected)); by=point -> point.rho)
end

function _curve(points::Vector{Point})
    length(points) >= 6 || return nothing
    rho = [point.rho for point in points]
    mu = [point.mu for point in points]
    all(isfinite, rho) && all(isfinite, mu) || return nothing
    length(unique(rho)) == length(rho) || return nothing
    return (rho=rho, mu=mu, points=points)
end

function _evaluate(curve)
    curve === nothing && return (
        status=:invalid, reason="solver_or_curve_failure", mu_transition=nothing,
        rho_hadron=nothing, rho_quark=nothing, area_residual=Inf,
        sres=MODELS.SShapeResult(), curve=nothing, maxwell=nothing,
    )
    classify = MODELS._classify_s_curve(
        curve.mu, curve.rho;
        maxwell_options=(; tol_area=MAXWELL_SOLVER_TOL),
        area_tol_good=CEP_AREA_TOL_GOOD,
        area_tol_bad=CEP_AREA_TOL_BAD,
    )
    maxwell = if classify.sres.has_s_shape
        MODELS.maxwell_construction(
            curve.mu, curve.rho; spinodal_hint=classify.sres, tol_area=MAXWELL_SOLVER_TOL,
        )
    else
        MODELS.MaxwellResult()
    end
    status = Symbol(classify.status)
    if status == :valid && !maxwell.converged
        status = :invalid
    end
    return (
        status=status,
        reason=String(classify.reason),
        mu_transition=maxwell.converged ? maxwell.mu_transition : classify.mu_transition,
        rho_hadron=maxwell.converged ? maxwell.rho_hadron : nothing,
        rho_quark=maxwell.converged ? maxwell.rho_quark : nothing,
        area_residual=maxwell.converged ? maxwell.area_residual : something(classify.area_residual, Inf),
        sres=classify.sres,
        curve=curve,
        maxwell=maxwell,
    )
end

function _geometry(left, right)
    MODELS._compare_phase_geometry(left, right, TOL)
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

function _candidates(curve::Union{Nothing, NamedTuple}; level::Symbol=:unknown)
    curve === nothing && return Candidate[]
    points = curve.points
    length(points) >= 5 || return Candidate[]
    slopes = [(points[i + 1].mu - points[i].mu) / (points[i + 1].rho - points[i].rho) for i in 1:(length(points) - 1)]
    signed = [(index=i, sign=slopes[i] > 0 ? 1 : slopes[i] < 0 ? -1 : 0) for i in eachindex(slopes) if slopes[i] != 0.0]
    isempty(signed) && return Candidate[]
    candidates = Candidate[]
    cursor = 1
    while cursor <= length(signed)
        signed[cursor].sign == -1 && (cursor += 1; continue)
        start = cursor + 1
        start <= length(signed) && signed[start].sign == -1 || (cursor += 1; continue)
        finish = start
        while finish < length(signed) && signed[finish + 1].sign == -1
            finish += 1
        end
        before_positive = signed[cursor].sign == 1
        after_positive = finish < length(signed) && signed[finish + 1].sign == 1
        if before_positive && after_positive
            first_index = signed[start].index
            last_index = signed[finish].index
            low = points[first_index].rho
            high = points[last_index + 1].rho
            width = high - low
            drop = points[first_index].mu - points[last_index + 1].mu
            step = median(diff([point.rho for point in points]))
            if width >= 2 * step && drop > 0.0
                push!(candidates, Candidate(low, high, drop, width, finish - start + 1, level))
            end
        end
        cursor = finish + 1
    end
    return candidates
end

function _stable_candidate_count(coarse, fine)
    coarse_candidates = _candidates(coarse; level=:coarse)
    fine_candidates = _candidates(fine; level=:fine)
    tolerance = 2 * STAGE_B_COARSE
    clusters = Vector{NamedTuple}()
    for candidate in vcat(coarse_candidates, fine_candidates)
        index = findfirst(cluster -> candidate.rho_low <= cluster.high + tolerance &&
            candidate.rho_high >= cluster.low - tolerance, clusters)
        if index === nothing
            push!(clusters, (low=candidate.rho_low, high=candidate.rho_high, count=1))
        else
            cluster = clusters[index]
            clusters[index] = (
                low=min(cluster.low, candidate.rho_low),
                high=max(cluster.high, candidate.rho_high),
                count=cluster.count + 1,
            )
        end
    end
    stable = [cluster for cluster in clusters if cluster.count >= 2]
    return stable, coarse_candidates, fine_candidates
end

function _drawdown(curve)
    curve === nothing && return (drop=0.0, width=0.0, rho_low=NaN, rho_high=NaN)
    points = curve.points
    best_mu = points[1].mu
    best_rho = points[1].rho
    best = (drop=0.0, width=0.0, rho_low=NaN, rho_high=NaN)
    for point in points[2:end]
        drop = best_mu - point.mu
        if drop > best.drop
            best = (drop=drop, width=point.rho - best_rho, rho_low=best_rho, rho_high=point.rho)
        end
        if point.mu > best_mu
            best_mu = point.mu
            best_rho = point.rho
        end
    end
    return best
end

function _target_values(eval, candidates)
    values = Float64[]
    for value in (
        eval.rho_hadron,
        eval.rho_quark,
        eval.sres.rho_spinodal_hadron,
        eval.sres.rho_spinodal_quark,
    )
        value === nothing || push!(values, Float64(value))
    end
    for candidate in candidates
        push!(values, candidate.rho_low)
        push!(values, (candidate.rho_low + candidate.rho_high) / 2)
        push!(values, candidate.rho_high)
    end
    draw = _drawdown(eval.curve)
    isfinite(draw.rho_low) && push!(values, draw.rho_low)
    isfinite(draw.rho_high) && push!(values, draw.rho_high)
    sort!(unique(filter(value -> 0.0 <= value <= 4.0, values)))
end

function _select_points(pool::Vector{Point}, targets::Vector{Float64}, cap::Int)
    cap <= 0 && return Point[]
    isempty(pool) && return Point[]
    eligible = Point[]
    for point in pool
        if isempty(targets) || minimum(abs(point.rho - target) for target in targets) <= FEATURE_RADIUS
            push!(eligible, point)
        end
    end
    isempty(eligible) && (eligible = copy(pool))
    distances = Dict(point.rho => (isempty(targets) ? 0.0 : minimum(abs(point.rho - target) for target in targets)) for point in eligible)
    sort!(eligible; by=point -> (distances[point.rho], point.rho))
    eligible[1:min(cap, length(eligible))]
end

function _anchor_rows(slices, xi::Float64, T::Float64)
    rows = [row for row in slices if _same(_float(_field(row, :xi)), xi; atol=1e-8) &&
        _same(_float(_field(row, :T_MeV)), T)]
    Dict(string(_field(row, :method, "")) => row for row in rows)
end

function _oracle_status(rows, xi, T)
    row = get(rows, "independent_oracle", nothing)
    row === nothing && return "missing"
    string(_field(row, :result_status, "missing"))
end

function _status_string(status)
    String(status)
end

function _build_anchor(curves, slices, xi::Float64, T::Float64, cap::Int)
    oracle = _point_rows(curves, "independent_oracle", xi, T, rho -> _on_grid(rho, STAGE_B_FINE))
    b_coarse = _evaluate(_curve(filter(point -> _on_grid(point.rho, STAGE_B_COARSE), oracle)))
    b_fine = _evaluate(_curve(oracle))
    a_coarse_points = _point_rows(curves, "production_hybrid", xi, T, rho -> _on_grid(rho, STAGE_A_COARSE))
    a_fine_points = _point_rows(curves, "production_hybrid", xi, T, rho -> _on_grid(rho, STAGE_A_FINE))
    a_coarse = _evaluate(_curve(a_coarse_points))
    a_fine = _evaluate(_curve(a_fine_points))
    a_geometry = _geometry(a_coarse, a_fine)
    a_status = _semantic(a_coarse, a_fine, a_geometry)
    if a_status == :confirmed_monotone
        return (
            xi=xi, T_MeV=T, cap=cap, result_status=:confirmed_monotone,
            reason="stage_a_certificate", stage_used=:stage_a,
            stage_a_status=a_status, stage_b_status=:not_run, stage_c_status=:not_run,
            stage_a_geometry=a_geometry, stage_b_geometry=nothing, stage_c_geometry=nothing,
            selected_points=Point[], candidates=Candidate[], stable_candidate_count=0,
            drawdown=_drawdown(a_fine.curve), oracle_status=_oracle_status(_anchor_rows(slices, xi, T), xi, T),
            cost_unique=length(_point_rows(curves, "production_hybrid", xi, T, _ -> true)),
            stage_c_point_count=0, finite_and_converged=true,
        )
    end
    stable, coarse_candidates, fine_candidates = _stable_candidate_count(b_coarse.curve, b_fine.curve)
    targets = _target_values(b_fine, fine_candidates)
    pool = _point_rows(curves, "independent_oracle", xi, T,
        rho -> _on_grid(rho, STAGE_C_FINE) && !_on_grid(rho, STAGE_B_FINE))
    selected = _select_points(pool, targets, cap)
    selected_rho = Set(point.rho for point in selected)
    merged_points = vcat(b_fine.curve === nothing ? Point[] : b_fine.curve.points, selected)
    merged_map = Dict(point.rho => point for point in merged_points)
    merged = _evaluate(_curve(sort!(collect(values(merged_map)); by=point -> point.rho)))
    b_geometry = _geometry(b_coarse, b_fine)
    c_geometry = _geometry(b_fine, merged)
    c_valid = b_fine.status == :valid && merged.status == :valid && c_geometry.converged && length(stable) == 1
    status = c_valid ? :confirmed_first_order : :ambiguous_near_critical
    reason = c_valid ? "stage_c_certificate" :
        length(stable) > 1 ? "multiple_cross_resolution_candidates" :
        b_fine.status != :valid ? "stage_b_unresolved" : "stage_c_geometry_unresolved"
    all_points = Set(point.rho for point in vcat(a_coarse_points, a_fine_points, b_coarse.curve === nothing ? Point[] : b_coarse.curve.points, selected))
    targeted = get(_anchor_rows(slices, xi, T), "production_hybrid", nothing)
    stage_a_targeted = targeted === nothing ? 0 : Int(round(_float(_field(targeted, :targeted_additions), 0.0)))
    return (
        xi=xi, T_MeV=T, cap=cap, result_status=status, reason=reason,
        stage_used=:stage_c, stage_a_status=a_status, stage_b_status=_semantic(b_coarse, b_fine, b_geometry),
        stage_c_status=status, stage_a_geometry=a_geometry, stage_b_geometry=b_geometry,
        stage_c_geometry=c_geometry, selected_points=selected,
        candidates=vcat(coarse_candidates, fine_candidates), stable_candidate_count=length(stable),
        drawdown=_drawdown(b_fine.curve), oracle_status=_oracle_status(_anchor_rows(slices, xi, T), xi, T),
        cost_unique=length(all_points) + stage_a_targeted, stage_c_point_count=length(selected),
        finite_and_converged=true,
    )
end

function _anchor_set(slices)
    sort!(unique(((_float(_field(row, :xi)), _float(_field(row, :T_MeV))) for row in slices if
        string(_field(row, :method, "")) == "independent_oracle")); by=first)
end

function _csv_rows(path::String)
    isfile(path) ? collect(CSV.File(path)) : NamedTuple[]
end

function _write_csv(path::String, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _curve_rows_for_output(result)
    rows = NamedTuple[]
    for point in result.selected_points
        push!(rows, (xi=result.xi, T_MeV=result.T_MeV, rho=point.rho, muq_MeV=point.mu,
            source="stage_c_selected", cap=result.cap, result_status=String(result.result_status)))
    end
    return rows
end

function _cost_frontier(curves, slices, costs, anchors)
    dense = sum(_float(_field(row, :unique_solves), 0.0) for row in costs if string(_field(row, :method, "")) == "memoized_dense")
    frontier = NamedTuple[]
    replays = Dict{Int, Vector{Any}}()
    for cap in CAPS
        results = [_build_anchor(curves, slices, xi, T, cap) for (xi, T) in anchors]
        replays[cap] = results
        mismatches = count(result -> begin
            expected = result.oracle_status
            actual = String(result.result_status)
            if (result.xi, result.T_MeV) in AUTHOR_FIRST_ORDER
                actual != "confirmed_first_order"
            elseif (result.xi, result.T_MeV) in CONSENSUS_MONOTONE
                actual != "confirmed_monotone"
            elseif expected == "ambiguous_near_critical"
                actual != "ambiguous_near_critical"
            else
                actual != expected
            end
        end, results)
        unsupported = count(result -> result.oracle_status == "ambiguous_near_critical" &&
            result.result_status != :ambiguous_near_critical, results)
        geometry_failures = count(result -> begin
            result.result_status == :confirmed_first_order &&
                (result.stage_used == :stage_a ? !result.stage_a_geometry.converged :
                 result.stage_c_geometry === nothing ? true : !result.stage_c_geometry.converged)
        end, results)
        stable_candidate_failures = count(result -> result.stable_candidate_count > 1, results)
        unique_solves = sum(result.cost_unique for result in results)
        state_gate = mismatches == 0 && unsupported == 0
        geometry_gate = geometry_failures == 0
        candidate_gate = stable_candidate_failures == 0
        cost_gate = unique_solves <= dense
        push!(frontier, (
            cap=cap, anchors=length(anchors), classification_mismatches=mismatches,
            unsupported_confirmations=unsupported, geometry_failures=geometry_failures,
            multiple_candidate_anchors=stable_candidate_failures, unique_solves=unique_solves,
            dense_unique_solves=dense, cost_gate=cost_gate, state_gate=state_gate,
            geometry_gate=geometry_gate, candidate_gate=candidate_gate,
            feasible=state_gate && geometry_gate && candidate_gate && cost_gate,
        ))
    end
    selected = findfirst(row -> row.feasible, frontier)
    return (; dense, frontier, replays, selected)
end

function _author_rows()
    [
        (xi=-0.5, T_MeV=147.0947265625, expected_status="confirmed_first_order", label_source="author_adjudication"),
        (xi=0.5, T_MeV=106.9599609375, expected_status="confirmed_first_order", label_source="author_adjudication"),
        (xi=-0.5, T_MeV=147.2197265625, expected_status="confirmed_monotone", label_source="three_method_consensus"),
        (xi=0.5, T_MeV=107.0849609375, expected_status="confirmed_monotone", label_source="three_method_consensus"),
    ]
end

function _write_docs(output_dir::String, policy, frontier, source_manifest)
    verdict = policy.verdict
    write(joinpath(output_dir, "README.md"), """# PNJL Hybrid Stage-C certificate feasibility v2

verdict: `$(verdict)`。本目录使用 source run `$(SOURCE_RUN_ID)` 的 final aggregate
artifact 做 solver-free replay；没有调用 equilibrium solver，不修改 production、reference 或历史 v1 evidence。

v2 直接复用 Julia `PhaseCore`、Maxwell construction 和 phase geometry comparison。Stage-C
路由使用 Stage-B 曲线的 drawdown、density support、turning-point 和 Maxwell 特征；oracle
状态和 geometry 只用于比较与 gate，不能选择补点。固定 slope margin 不再作为弱 S 的必要条件。

- tested caps: `$(join(string.(CAPS), ", "))`
- dense unique-solve reference: `$(frontier.dense)`
- selected cap: `$(policy.selected_cap)`
- solver called: `false`

`author_adjudication` 只适用于两个作者确认的一阶点；高温 monotone 点标记为
`three_method_consensus`。如果 verdict 不是 `feasible_candidate`，本结果不授权 production
集成或 reference 重放。
""")
    write(joinpath(output_dir, "AUDIT.md"), """# Stage-C certificate feasibility v2 audit

输入 manifest SHA、source run ID、calculation SHA 和本地 replay producer SHA 均写入
`manifest.json`。所有最终 curve 证书均由 Julia 核心重算；本审计不把历史 oracle row 的
geometry 直接当作 hybrid 证书。完整 raw `curve_points.csv` 继续留在 Actions/local artifact，
仓库只保留 selected-point index 与代表性曲线。

成本按 Stage-A、完整 Stage-B 和 selected Stage-C rho key 的并集计算，Stage-B semantic
grid 不免费计入。离线 replay 不能证明 residual/Jacobian/runner-minutes，须由后续 Actions
targeted/full shadow 验证。
""")
    _write_csv(joinpath(output_dir, "author_adjudications.csv"), _author_rows())
    _write_csv(joinpath(output_dir, "claim_ledger.csv"), [
        (claim_id="classification", claim="v2 replay matches author anchors and oracle non-ambiguous anchors", status=verdict == "feasible_candidate" ? "pass" : "inconclusive", boundary="offline semantic replay only"),
        (claim_id="geometry", claim="Maxwell and rho geometry are recomputed with Julia core", status=verdict == "feasible_candidate" ? "pass" : "inconclusive", boundary="no new equilibrium solves"),
        (claim_id="cost", claim="Stage-A/B/C union unique solves do not exceed dense", status=all(row -> row.cost_gate, frontier.frontier) ? "pass" : "inconclusive", boundary="residual/Jacobian and runner costs require Actions"),
        (claim_id="promotion", claim="v2 authorizes production/reference promotion", status="not_claimed", boundary="requires targeted/full shadow and author review"),
    ])
    hashes = Dict{String, String}()
    for path in readdir(output_dir; join=true),
            name in ("README.md", "AUDIT.md", "author_adjudications.csv", "claim_ledger.csv")
        basename(path) == name && (hashes[name] = _sha256_file(path))
    end
    return hashes
end

function _write_manifest(output_dir::String, data, frontier, policy, input_dir::String)
    files = Dict{String, String}()
    for path in walkdir(output_dir)
        root, dirs, names = path
        for name in names
            name == "manifest.json" && continue
            full = joinpath(root, name)
            files[replace(relpath(full, output_dir), '\\' => '/')] = _sha256_file(full)
        end
    end
    manifest = Dict(
        "schema_version" => "cep_hybrid_stagec_certificate_feasibility_v2",
        "verdict" => policy.verdict,
        "source_run_id" => SOURCE_RUN_ID,
        "source_calculation_sha" => SOURCE_CALCULATION_SHA,
        "source_manifest_sha256" => _sha256_file(joinpath(input_dir, "manifest.json")),
        "producer_head_sha" => try readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`) catch; "unknown" end,
        "solver_called" => false,
        "selected_policy" => policy,
        "cost_frontier" => frontier.frontier,
        "files" => files,
        "gate" => Dict(
            "automatic_gate_is_not_promotion" => true,
            "author_first_order" => [Dict("xi" => xi, "T_MeV" => T) for (xi, T) in AUTHOR_FIRST_ORDER],
            "consensus_monotone" => [Dict("xi" => xi, "T_MeV" => T) for (xi, T) in CONSENSUS_MONOTONE],
        ),
    )
    open(joinpath(output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, manifest)
        write(io, '\n')
    end
    return manifest
end

function run(input_dir::String, output_dir::String)
    data = _load_rows(input_dir)
    _validate_input(input_dir, data)
    anchors = _anchor_set(data.slices)
    length(anchors) == 24 || error("expected 24 independent-oracle anchors, got $(length(anchors))")
    frontier = _cost_frontier(data.curves, data.slices, data.costs, anchors)
    selected = frontier.selected
    ambiguous_confirmation = any(result -> result.oracle_status == "ambiguous_near_critical" &&
        result.result_status != :ambiguous_near_critical, frontier.replays[first(CAPS)])
    verdict = if selected === nothing
        ambiguous_confirmation ? "oracle_inconclusive" :
        any(row -> !row.state_gate, frontier.frontier) ? "integration_failed" :
        any(row -> !row.candidate_gate || row.multiple_candidate_anchors > 0, frontier.frontier) ? "maxwell_candidate_inconclusive" :
        any(row -> !row.cost_gate, frontier.frontier) ? "performance_inconclusive" : "integration_failed"
    else
        "feasible_candidate"
    end
    policy = (
        schema_version="cep_hybrid_stagec_certificate_feasibility_v2",
        verdict=verdict,
        selected_cap=selected === nothing ? nothing : frontier.frontier[selected].cap,
        stage_b_grid_step=STAGE_B_FINE,
        stage_c_local_step=STAGE_C_FINE,
        feature_radius=FEATURE_RADIUS,
        caps=CAPS,
        dense_unique_solves=frontier.dense,
        reason=selected === nothing ? "no cap satisfies all v2 gates" : "minimum feasible cap",
    )
    mkpath(output_dir)
    replay_rows = NamedTuple[]
    candidate_rows = NamedTuple[]
    selected_rows = NamedTuple[]
    representative_rows = NamedTuple[]
    for cap in CAPS
        for result in frontier.replays[cap]
            push!(replay_rows, (
                xi=result.xi, T_MeV=result.T_MeV, cap=cap,
                oracle_status=result.oracle_status, result_status=String(result.result_status),
                reason=result.reason, stage_used=String(result.stage_used),
                stage_a_status=String(result.stage_a_status), stage_b_status=String(result.stage_b_status),
                stage_c_status=String(result.stage_c_status), stable_candidate_count=result.stable_candidate_count,
                stage_c_point_count=result.stage_c_point_count, unique_solves=result.cost_unique,
                drawdown_mu=result.drawdown.drop, drawdown_width_rho=result.drawdown.width,
                stage_b_geometry=result.stage_b_geometry === nothing ? "not_run" : result.stage_b_geometry.reason,
                stage_c_geometry=result.stage_c_geometry === nothing ? "not_run" : result.stage_c_geometry.reason,
            ))
            for candidate in result.candidates
                push!(candidate_rows, (
                    xi=result.xi, T_MeV=result.T_MeV, cap=cap, level=String(candidate.level),
                    rho_low=candidate.rho_low, rho_high=candidate.rho_high,
                    width=candidate.width, drop_mu=candidate.drop_mu,
                    negative_secants=candidate.negative_secants,
                ))
            end
            if cap == first(CAPS)
                append!(selected_rows, _curve_rows_for_output(result))
            end
        end
    end
    _write_csv(joinpath(output_dir, "anchor_replay.csv"), replay_rows)
    _write_csv(joinpath(output_dir, "candidate_runs.csv"), candidate_rows)
    _write_csv(joinpath(output_dir, "selected_points.csv"), selected_rows)
    _write_csv(joinpath(output_dir, "cost_frontier.csv"), frontier.frontier)
    deep_rows = NamedTuple[]
    for result in frontier.replays[first(CAPS)]
        if result.oracle_status == "ambiguous_near_critical" && result.result_status != :ambiguous_near_critical
            push!(deep_rows, (xi=result.xi, T_MeV=result.T_MeV,
                oracle_status=result.oracle_status, hybrid_status=String(result.result_status),
                reason="deep_oracle_required"))
        end
    end
    _write_csv(joinpath(output_dir, "deep_oracle_required.csv"), deep_rows)
    _write_docs(output_dir, policy, frontier, data.manifest)
    _write_csv(joinpath(output_dir, "curve_index.csv"), [
        (source="aggregate_replay", path="curve_points.csv", sha256=_sha256_file(joinpath(input_dir, "curve_points.csv")), raw_curve_copy_in_repository=false),
    ])
    open(joinpath(output_dir, "plot_manifest.json"), "w") do io
        JSON3.pretty(io, Dict("schema_version" => "cep_hybrid_stagec_certificate_feasibility_v2", "figures" => Any[], "source_curve_artifact" => "external_actions_or_local_artifact"))
        write(io, '\n')
    end
    _write_manifest(output_dir, data, frontier, policy, input_dir)
    return policy
end

function _parse_args(args)
    input_dir = nothing
    output_dir = joinpath(PROJECT_ROOT, "docs", "analysis", "pnjl", "cep_maxwell", "stagec", "cep_hybrid_stagec_certificate_feasibility_v2")
    for arg in args
        startswith(arg, "--input-dir=") && (input_dir = abspath(split(arg, "="; limit=2)[2]))
        startswith(arg, "--output-dir=") && (output_dir = abspath(split(arg, "="; limit=2)[2]))
        arg in ("-h", "--help") && return nothing
    end
    input_dir === nothing && throw(ArgumentError("--input-dir=PATH is required"))
    return (; input_dir, output_dir)
end

function main(args=ARGS)
    options = _parse_args(args)
    options === nothing && (println("Usage: julia --project=. scripts/analysis/pnjl_cep_hybrid_stagec_certificate_feasibility_v2.jl --input-dir=PATH [--output-dir=PATH]"); return 0)
    policy = run(options.input_dir, options.output_dir)
    println(JSON3.write(policy))
    policy.verdict == "feasible_candidate" ? 0 : 2
end

end # module StageCCertificateFeasibilityV2

if abspath(PROGRAM_FILE) == @__FILE__
    exit(StageCCertificateFeasibilityV2.main())
end
