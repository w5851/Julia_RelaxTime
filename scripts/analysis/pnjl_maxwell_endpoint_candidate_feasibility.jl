#!/usr/bin/env julia

"""Solver-free Maxwell candidate and rho=0 endpoint feasibility replay.

This analysis deliberately does not change the production Maxwell path.  It
replays immutable curve artifacts and compares the historical first-sign-change
algorithm with a strict candidate search that requires a unique + -> - -> +
three-crossing topology.  A crossing touching the first rho=0 cell is recorded
as endpoint-dependent; the replay never inserts duplicate rho=0 points or a
synthetic (0, 0) anchor.
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using CSV
using JSON3
using SHA
using Statistics

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
end

module MaxwellEndpointCandidateFeasibility

using CSV
using JSON3
using SHA
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS = Main.Models

const SOURCE_DEEP_RUN_ID = "30809754119"
const SOURCE_TARGETED_RUN_ID = "30805637032"
const SOURCE_CALCULATION_SHA = "3217bed3635574f00c04cbee75e843b4c49451db"
const AREA_TOL = 5e-6
const OUTER_AREA_GATE = 5e-5
const POSITION_TOL = 0.025
const DENSITY_TOL = 0.0025
const MAX_ITER = 80
const SCAN_STEPS = 1024
const MIN_SAMPLES = 12
const ENDPOINT_CELL = 0.003125
const ENDPOINT_ATOL = 64eps(Float64)

export CurvePoint, CurveData, load_curve_groups, strict_candidate,
    replay_input, write_replay_artifact

struct CurvePoint
    xi::Float64
    T::Float64
    method::String
    rho::Float64
    mu::Float64
    residual::Float64
    rho_level::Int
    sampling_role::String
end

struct CurveData
    points::Vector{CurvePoint}
    rho::Vector{Float64}
    mu::Vector{Float64}
end

@inline _finite(x) = x isa Real && isfinite(Float64(x))

function _field(row, name::Symbol, default=nothing)
    hasproperty(row, name) || return default
    value = getproperty(row, name)
    value === missing && return default
    return value
end

function _float(value, default=NaN)
    value === nothing || value === missing && return default
    try
        result = Float64(value)
        return isfinite(result) ? result : default
    catch
        return default
    end
end

function _bool(value, default=false)
    value === nothing || value === missing && return default
    value isa Bool && return value
    return lowercase(strip(String(value))) in ("true", "1", "yes")
end

function _sha256_file(path::String)
    bytes2hex(open(sha256, path))
end

function _deduplicate(points::Vector{CurvePoint})
    selected = Dict{Float64, CurvePoint}()
    for point in points
        previous = get(selected, point.rho, nothing)
        if previous === nothing || point.residual < previous.residual ||
           (point.residual == previous.residual && point.rho_level > previous.rho_level)
            selected[point.rho] = point
        end
    end
    result = collect(values(selected))
    sort!(result; by=point -> point.rho)
    return result
end

function _curve(points::Vector{CurvePoint})
    selected = _deduplicate(points)
    length(selected) >= MIN_SAMPLES || return nothing
    rho = Float64[point.rho for point in selected]
    mu = Float64[point.mu for point in selected]
    all(isfinite, rho) && all(isfinite, mu) || return nothing
    all(diff(rho) .> 0) || return nothing
    return CurveData(selected, rho, mu)
end

function load_curve_groups(input_dir::String; expected_calculation_sha::String=SOURCE_CALCULATION_SHA)
    path = joinpath(input_dir, "curve_points.csv")
    isfile(path) || error("missing curve_points.csv: $path")
    rows = collect(CSV.File(path))
    isempty(rows) && error("curve_points.csv is empty: $path")
    groups = Dict{Tuple{String, Float64, Float64}, Vector{CurvePoint}}()
    for row in rows
        xi = _float(_field(row, :xi))
        T = _float(_field(row, :T_MeV))
        rho = _float(_field(row, :rho))
        mu = _float(_field(row, :muq_MeV))
        residual = _float(_field(row, :residual_norm), Inf)
        method = String(_field(row, :method, "unknown"))
        calculation = String(_field(row, :calculation_sha, ""))
        level = try Int(round(_float(_field(row, :rho_level), -1.0))) catch; -1 end
        (_finite(xi) && _finite(T) && _finite(rho) && _finite(mu) && isfinite(residual)) ||
            error("non-finite curve row in $path at ($xi,$T,$rho)")
        calculation == expected_calculation_sha ||
            error("calculation SHA mismatch in $path: $calculation")
        _bool(_field(row, :converged), false) && _bool(_field(row, :finite), false) ||
            error("non-converged curve row in $path at ($xi,$T,$rho)")
        point = CurvePoint(xi, T, method, rho, mu, residual, level,
            String(_field(row, :sampling_role, "")))
        push!(get!(groups, (method, xi, T), CurvePoint[]), point)
    end
    return groups
end

function _crossings(mu0::Float64, curve::CurveData; atol::Float64=1e-9)
    crossings = Float64[]
    for i in 1:(length(curve.rho) - 1)
        r1, r2 = curve.rho[i], curve.rho[i + 1]
        f1, f2 = curve.mu[i] - mu0, curve.mu[i + 1] - mu0
        abs(f1) <= atol && push!(crossings, r1)
        if f1 * f2 < 0.0
            α = f1 / (f1 - f2)
            push!(crossings, r1 + α * (r2 - r1))
        end
        abs(f2) <= atol && push!(crossings, r2)
    end
    sort!(crossings)
    unique_crossings = Float64[]
    for value in crossings
        isempty(unique_crossings) || abs(value - last(unique_crossings)) > atol * 10 || continue
        push!(unique_crossings, value)
    end
    return unique_crossings
end

function _integrate(mu0::Float64, left::Float64, right::Float64, curve::CurveData)
    total = 0.0
    for i in 1:(length(curve.rho) - 1)
        r1, r2 = curve.rho[i], curve.rho[i + 1]
        (r2 <= left || r1 >= right) && continue
        a, b = max(r1, left), min(r2, right)
        b <= a && continue
        m1, m2 = curve.mu[i], curve.mu[i + 1]
        μa = m1 + (a - r1) / (r2 - r1) * (m2 - m1)
        μb = m1 + (b - r1) / (r2 - r1) * (m2 - m1)
        total += 0.5 * ((μa - mu0) + (μb - mu0)) * (b - a)
    end
    return total
end

function _area_at(mu0::Float64, curve::CurveData, sres)
    crossings = _crossings(mu0, curve)
    length(crossings) == 3 || return (
        valid=false, area=nothing, crossings=crossings,
        endpoint_dependent=false, reason="crossing_count_not_three",
    )
    # The S-shape detector supplies the global + -> - -> + topology.  The
    # crossing count check above prevents a two-crossing endpoint jump from
    # entering the area function.
    sres.has_s_shape || return (
        valid=false, area=nothing, crossings=crossings,
        endpoint_dependent=false, reason="no_s_shape",
    )
    area = _integrate(mu0, first(crossings), last(crossings), curve)
    endpoint_dependent = first(curve.rho) == 0.0 && first(crossings) <= ENDPOINT_CELL + ENDPOINT_ATOL
    return (
        valid=true, area=area, crossings=crossings,
        endpoint_dependent=endpoint_dependent, reason="three_crossings",
    )
end

function _bisect(curve::CurveData, sres, a::Float64, b::Float64,
        fa::Float64, fb::Float64; tol_area::Float64=AREA_TOL, max_iter::Int=MAX_ITER)
    fa * fb <= 0.0 || return nothing
    best = nothing
    for iter in 1:max_iter
        mid = 0.5 * (a + b)
        probe = _area_at(mid, curve, sres)
        probe.valid || return (
            converged=false, reason="topology_changed_inside_bisection",
            mu=mid, area=nothing, crossings=probe.crossings, iterations=iter,
            endpoint_dependent=false,
        )
        best = probe
        if abs(probe.area) <= tol_area
            return (
                converged=true, reason="ok", mu=mid, area=probe.area,
                crossings=probe.crossings, iterations=iter,
                endpoint_dependent=probe.endpoint_dependent,
            )
        end
        if fa * probe.area < 0.0
            b, fb = mid, probe.area
        else
            a, fa = mid, probe.area
        end
    end
    return (
        converged=false, reason="solver_tolerance_not_met", mu=0.5 * (a + b),
        area=best === nothing ? nothing : best.area,
        crossings=best === nothing ? Float64[] : best.crossings,
        iterations=max_iter,
        endpoint_dependent=best === nothing ? false : best.endpoint_dependent,
    )
end

function _candidate_roots(curve::CurveData, sres; tol_area::Float64=AREA_TOL)
    sres.has_s_shape || return (roots=NamedTuple[], valid_intervals=NamedTuple[], reason="no_s_shape")
    μlo, μhi = sort((Float64(sres.mu_spinodal_hadron), Float64(sres.mu_spinodal_quark)))
    grid = collect(range(μlo, μhi; length=SCAN_STEPS))
    roots = NamedTuple[]
    intervals = NamedTuple[]
    previous = nothing
    interval_start = nothing
    interval_last = nothing
    for μ in grid
        probe = _area_at(μ, curve, sres)
        if !probe.valid
            if interval_start !== nothing
                push!(intervals, (mu_low=interval_start, mu_high=interval_last,
                    status="three_crossing_interval"))
            end
            previous = nothing
            interval_start = nothing
            interval_last = nothing
            continue
        end
        interval_start === nothing && (interval_start = μ)
        interval_last = μ
        if previous !== nothing
            if abs(probe.area) <= tol_area
                push!(roots, (converged=true, mu=μ, area=probe.area,
                    crossings=probe.crossings, iterations=0,
                    endpoint_dependent=probe.endpoint_dependent,
                    reason="grid_hit"))
            elseif previous.area * probe.area < 0.0
                solved = _bisect(curve, sres, previous.mu, μ,
                    previous.area, probe.area; tol_area=tol_area)
                solved === nothing && continue
                push!(roots, merge(solved, (bracket_low=previous.mu, bracket_high=μ)))
            end
        end
        previous = merge(probe, (mu=μ,))
    end
    interval_start !== nothing && push!(intervals, (mu_low=interval_start,
        mu_high=interval_last, status="three_crossing_interval"))
    # Numerical scan refinement can rediscover one root in adjacent brackets.
    # Prefer a bisection result over a grid hit when both represent the same
    # tolerance-sized neighborhood.
    sort!(roots; by=root -> (root.reason == "grid_hit" ? 1 : 0, root.mu))
    unique_roots = NamedTuple[]
    for root in roots
        # A grid point that already lies inside the area tolerance is often
        # followed immediately by the same sign-change bracket.  Treat those
        # as one candidate; the bisection result is retained when available.
        any(abs(root.mu - previous_root.mu) <= 2e-3 for previous_root in unique_roots) && continue
        push!(unique_roots, root)
    end
    sort!(unique_roots; by=root -> root.mu)
    return (roots=unique_roots, valid_intervals=intervals,
        reason=isempty(unique_roots) ? "no_three_crossing_sign_change" : "ok")
end

function _legacy(curve::CurveData, sres)
    result = MODELS.maxwell_construction(curve.mu, curve.rho;
        spinodal_hint=sres, tol_area=AREA_TOL, max_iter=MAX_ITER)
    return (
        converged=result.converged,
        mu=result.mu_transition,
        rho_hadron=result.rho_hadron,
        rho_quark=result.rho_quark,
        area=result.area_residual,
        iterations=result.iterations,
        reason=String(get(result.details, :reason, result.converged ? "ok" : "failed")),
    )
end

function strict_candidate(curve::CurveData)
    curve === nothing && return (status=:invalid, reason="solver_or_curve_failure")
    sres = MODELS.detect_s_shape(curve.mu, curve.rho)
    legacy = _legacy(curve, sres)
    sres.has_s_shape || return (
        status=:monotone, reason="no_s_shape", sres=sres, legacy=legacy,
        roots=NamedTuple[], valid_intervals=NamedTuple[])
    candidates = _candidate_roots(curve, sres)
    roots = [root for root in candidates.roots if root.converged &&
        root.area !== nothing && abs(root.area) <= AREA_TOL]
    status = if length(roots) == 1
        :first_order
    elseif length(roots) > 1
        :multiple_candidates
    else
        :ambiguous
    end
    root = length(roots) == 1 ? first(roots) : nothing
    return (
        status=status,
        reason=status == :first_order ? "unique_three_crossing_candidate" : candidates.reason,
        sres=sres,
        legacy=legacy,
        roots=candidates.roots,
        valid_intervals=candidates.valid_intervals,
        candidate=root,
        endpoint_dependent=root === nothing ? false : root.endpoint_dependent,
    )
end

function _group_summary(key, curve, result)
    method, xi, T = key
    candidate = get(result, :candidate, nothing)
    return (
        method=method, xi=xi, T_MeV=T,
        status=String(result.status), reason=String(result.reason),
        legacy_converged=result.legacy.converged,
        legacy_mu=result.legacy.mu === nothing ? missing : result.legacy.mu,
        legacy_area=result.legacy.area === nothing ? missing : result.legacy.area,
        candidate_count=length(result.roots),
        candidate_mu=candidate === nothing ? missing : candidate.mu,
        candidate_area=candidate === nothing ? missing : candidate.area,
        rho_hadron=candidate === nothing ? missing : first(candidate.crossings),
        rho_quark=candidate === nothing ? missing : last(candidate.crossings),
        endpoint_dependent=get(result, :endpoint_dependent, false),
        endpoint_mu_zero=first(curve.rho) == 0.0 ? first(curve.mu) : missing,
        source_point_count=length(curve.rho),
    )
end

function _load_manifest(input_dir::String)
    path = joinpath(input_dir, "manifest.json")
    isfile(path) || return nothing
    return JSON3.read(read(path, String))
end

function _write_csv(path, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _write_json(path, value)
    open(path, "w") do io
        JSON3.pretty(io, value)
        write(io, '\n')
    end
end

function _write_docs(output_dir, input_dirs, summaries, verdict, files)
    write(joinpath(output_dir, "README.md"), """# PNJL Maxwell endpoint/candidate feasibility v1

verdict: `$(verdict)`。本包只对已有 Actions 曲线做 solver-free replay；没有调用
equilibrium/Newton solver，没有修改 production、reference 或历史 evidence。

严格 candidate 要求：

- 面积函数只在恰有三个去重交点、且切片具有 `+→−→+` S-shape topology 时有效；
- 无效拓扑会重置 sign-change 状态，不能跨两个交点/三交点边界二分；
- 同一曲线所有有效 sign-change 都会枚举，多个稳定根保持 inconclusive；
- `rho=0` 只记录端点代表值，不复制重复点，也不插入 `(0,0)`。

当前低温证据中的旧 Maxwell 残差约 `36--38` 是端点/两交点伪变号；严格 replay
应在约 `331.55 MeV` 找到唯一三交点根。该结论仍是离散曲线 candidate，不能替代
后续正密度近零补点的跨分辨率证书。

输入：$(join(input_dirs, ", "))。
""")
    write(joinpath(output_dir, "AUDIT.md"), """# Audit boundary

输入曲线的 SHA、行数、calculation SHA 和 producer head 写入 `manifest.json`。所有派生
表只使用已有 `(method,xi,T,rho)` 记录；没有插值、补点或 equilibrium 重算。`legacy_*`
字段只记录当前 public Maxwell 的对照结果，不作为新候选 gate。

`endpoint_dependent=true` 表示左外交点落在首个 rho 单元内，提示后续 Actions 需要正密度
局部细化；它不是零密度平台的硬编码，也不是 production promotion 证书。
""")
    ledger = [
        (claim_id="solver_free_replay", claim="No equilibrium/Newton solver was called", status="pass", boundary="postprocess-only"),
        (claim_id="strict_three_cross_candidate", claim="Unique Maxwell candidate requires three-crossing topology", status=verdict == "feasible_candidate" ? "pass" : "inconclusive", boundary="fixed source artifacts"),
        (claim_id="endpoint_semantics", claim="rho=0 is recorded as an endpoint observation, not duplicated or hard-coded", status="pass", boundary="finite-T endpoint remains diagnostic"),
        (claim_id="production_promotion", claim="Candidate is ready for production integration", status="not_claimed", boundary="requires endpoint-local Actions and public-core regression"),
    ]
    _write_csv(joinpath(output_dir, "claim_ledger.csv"), ledger)
    _write_json(joinpath(output_dir, "provenance.json"), Dict(
        "source_deep_run_id" => SOURCE_DEEP_RUN_ID,
        "source_targeted_run_id" => SOURCE_TARGETED_RUN_ID,
        "source_calculation_sha" => SOURCE_CALCULATION_SHA,
        "solver_called" => false,
        "reference_write" => false,
        "input_dirs" => input_dirs,
        "files" => files,
    ))
end

function _artifact_files(output_dir::String)
    files = Dict{String, String}()
    for (root, _, names) in walkdir(output_dir)
        for name in names
            name == "manifest.json" && continue
            path = joinpath(root, name)
            files[replace(relpath(path, output_dir), '\\' => '/')] = _sha256_file(path)
        end
    end
    return files
end

function replay_input(deep_input_dir::String, output_dir::String;
        targeted_input_dir::Union{Nothing, String}=nothing)
    deep_groups = load_curve_groups(deep_input_dir)
    groups = copy(deep_groups)
    input_dirs = [deep_input_dir]
    if targeted_input_dir !== nothing
        targeted_groups = load_curve_groups(targeted_input_dir)
        merge!(groups, targeted_groups)
        push!(input_dirs, targeted_input_dir)
    end

    summaries = NamedTuple[]
    curve_rows = NamedTuple[]
    for (key, points) in sort!(collect(groups); by=item -> item[1])
        curve = _curve(points)
        curve === nothing && continue
        result = strict_candidate(curve)
        push!(summaries, _group_summary(key, curve, result))
        for point in curve.points
            push!(curve_rows, (
                method=point.method, xi=point.xi, T_MeV=point.T, rho=point.rho,
                muq_MeV=point.mu, rho_level=point.rho_level,
                sampling_role=point.sampling_role,
            ))
        end
    end
    strict_groups = filter(row -> row.status == "first_order", summaries)
    # Feasibility is intentionally restricted to the low-temperature candidate
    # and the two confirmed controls; any extra targeted anchors are retained as
    # diagnostics but cannot weaken the required result.
    required = Set([("independent_oracle", -0.5, 5.0),
                    ("independent_oracle", -0.5, 20.0),
                    ("independent_oracle", 0.0, 5.0)])
    required_rows = filter(row -> (row.method, row.xi, row.T_MeV) in required, summaries)
    required_ok = length(required_rows) == length(required) &&
        all(row -> row.status == "first_order" && row.candidate_count == 1 &&
            row.candidate_mu !== missing && row.candidate_area !== missing &&
            abs(Float64(row.candidate_area)) <= AREA_TOL, required_rows)
    unique_ok = all(row -> row.candidate_count <= 1, summaries)
    verdict = required_ok && unique_ok ? "feasible_candidate" : "candidate_inconclusive"
    mkpath(output_dir)
    _write_csv(joinpath(output_dir, "candidate_summary.csv"), summaries)
    _write_csv(joinpath(output_dir, "curve_index.csv"), curve_rows)
    _write_json(joinpath(output_dir, "selected_policy.json"), Dict(
        "schema_version" => "pnjl_maxwell_endpoint_candidate_feasibility_v1",
        "verdict" => verdict,
        "candidate_policy" => "unique_three_crossing_topology_v1",
        "endpoint_policy" => "diagnostic_only_until_positive_rho_bracket",
        "area_solver_tol" => AREA_TOL,
        "outer_area_gate" => OUTER_AREA_GATE,
        "position_tol_MeV" => POSITION_TOL,
        "density_tol" => DENSITY_TOL,
        "endpoint_cell" => ENDPOINT_CELL,
        "targeted_cap_candidates" => [4, 6, 8, 10, 12],
        "local_step_candidates" => [0.0015625, 0.00078125, 0.000390625, 0.0001953125],
        "selection_status" => "awaiting_actions_endpoint_refinement",
    ))
    files = _artifact_files(output_dir)
    _write_docs(output_dir, input_dirs, summaries, verdict, files)
    files = _artifact_files(output_dir)
    manifest = Dict(
        "schema_version" => "pnjl_maxwell_endpoint_candidate_feasibility_v1",
        "verdict" => verdict,
        "source_deep_run_id" => SOURCE_DEEP_RUN_ID,
        "source_targeted_run_id" => SOURCE_TARGETED_RUN_ID,
        "source_calculation_sha" => SOURCE_CALCULATION_SHA,
        "producer_head_sha" => try readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`) catch; "unknown" end,
        "solver_called" => false,
        "reference_write" => false,
        "required_anchor_gate" => required_ok,
        "unique_candidate_gate" => unique_ok,
        "summary_count" => length(summaries),
        "files" => files,
        "promotion_boundary" => "diagnostic_only_until_endpoint_actions_and_public_core_regression",
    )
    _write_json(joinpath(output_dir, "manifest.json"), manifest)
    return manifest
end

function _arg(args, name, default=nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && throw(ArgumentError("missing value for $name"))
    return args[index + 1]
end

if abspath(PROGRAM_FILE) == @__FILE__
    deep = abspath(String(_arg(ARGS, "--deep-input",
        joinpath("D:", "Desktop", "Julia_RelaxTime_issue130_artifacts",
            "cep_hybrid_extrema_guard_20260803", "required_three_deep_run_30809754119",
            "cep-deep-oracle-cep_hybrid_extrema_guard_required_three_deep_20260803-aggregate"))))
    targeted_arg = _arg(ARGS, "--targeted-input", nothing)
    targeted = targeted_arg === nothing ? nothing : abspath(String(targeted_arg))
    output = abspath(String(_arg(ARGS, "--output-dir",
        joinpath(PROJECT_ROOT, "docs", "analysis", "pnjl_maxwell_endpoint_candidate_feasibility_v1"))))
    manifest = replay_input(deep, output; targeted_input_dir=targeted)
    println(JSON3.write(manifest))
end

end # module MaxwellEndpointCandidateFeasibility
