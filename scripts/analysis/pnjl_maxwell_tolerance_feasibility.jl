#!/usr/bin/env julia

"""Solver-free Maxwell tolerance-contract feasibility replay.

The runner consumes an immutable deep-oracle aggregate artifact and replays
only the phase post-processing layer.  It never calls an equilibrium/Newton
solver and it does not change the production tolerance contract.  The sweep
separates the tolerance used by the bisection from the existing outer
coarse/fine geometry gate so that the two contracts can be audited together.
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
const MODELS = Main.Models

module MaxwellToleranceFeasibility

using CSV
using JSON3
using SHA
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODELS = Main.Models

const SOURCE_RUN_ID = "30676440627"
const SOURCE_CALCULATION_SHA = "fde2b929b60575f1daacb84a1b9b8ff6e3b0a0cc"
const SOURCE_WORKFLOW_HEAD_SHA = "592374ff23c101575587d1164e6aeca9a7231fc1"
const EXPECTED_POINTS = Set([
    (-0.5, 20.0),
    (0.0, 5.0),
    (0.0, 20.0),
    (0.5, 5.0),
    (0.5, 20.0),
])
const AREA_TOLS = (1e-4, 5e-5, 1e-5, 5e-6)
const STRICT_SOLVER_TOL = last(AREA_TOLS)
const OUTER_AREA_GATE = 5e-5
const POSITION_TOL = 0.025
const DENSITY_TOL = 0.0025
const MAX_ITER = 60
const MIN_SAMPLES = 12

struct Point
    xi::Float64
    T::Float64
    rho::Float64
    mu::Float64
    residual::Float64
    rho_level::Int
    sampling_role::String
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

function _same(a::Float64, b::Float64; atol::Float64=1e-8)
    isfinite(a) && isfinite(b) && isapprox(a, b; atol=atol, rtol=0.0)
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

function _canonical_point_key(rho::Float64)
    # CSV rows are printed to six decimal places.  A 2e-6 bucket removes only
    # formatting duplicates and is much smaller than either rho layer step.
    round(Int, rho / 2e-6)
end

function _deduplicate_points(points::Vector{Point})
    selected = Dict{Int, Point}()
    for point in points
        key = _canonical_point_key(point.rho)
        previous = get(selected, key, nothing)
        if previous === nothing || point.residual < previous.residual
            selected[key] = point
        end
    end
    result = collect(values(selected))
    sort!(result; by=point -> point.rho)
    return result
end

function _load_input(input_dir::String)
    required = ("manifest.json", "curve_points.csv")
    missing_files = [name for name in required if !isfile(joinpath(input_dir, name))]
    isempty(missing_files) || error("missing deep-oracle files: $(join(missing_files, ", "))")

    manifest = JSON3.read(read(joinpath(input_dir, "manifest.json"), String))
    schema = string(_manifest_value(manifest, :schema_version, ""))
    schema == "cep_deep_oracle_v1" || error("unsupported source schema: $schema")
    calculation = string(_manifest_value(manifest, :calculation_sha, ""))
    workflow = string(_manifest_value(manifest, :workflow_head_sha, ""))
    calculation == SOURCE_CALCULATION_SHA || error("calculation SHA mismatch: $calculation")
    workflow == SOURCE_WORKFLOW_HEAD_SHA || error("workflow head SHA mismatch: $workflow")
    string(_manifest_value(manifest, :status, "complete")) == "complete" ||
        error("source deep-oracle artifact is not complete")

    rows = collect(CSV.File(joinpath(input_dir, "curve_points.csv")))
    isempty(rows) && error("curve_points.csv is empty")
    groups = Dict{Tuple{Float64, Float64}, Vector{Point}}()
    for row in rows
        xi = _float(_field(row, :xi))
        T = _float(_field(row, :T_MeV))
        rho = _float(_field(row, :rho))
        mu = _float(_field(row, :muq_MeV))
        residual = _float(_field(row, :residual_norm), Inf)
        level = try Int(round(_float(_field(row, :rho_level), -1.0))) catch; -1 end
        key = (xi, T)
        (_finite(xi) && _finite(T) && _finite(rho) && _finite(mu) && isfinite(residual)) ||
            error("non-finite source curve row for $key at rho=$rho")
        _bool(_field(row, :converged), false) && _bool(_field(row, :finite), false) ||
            error("source curve row is not finite/converged for $key at rho=$rho")
        level in (0, 1) || error("unexpected rho_level=$level for $key at rho=$rho")
        point = Point(xi, T, rho, mu, residual, level, String(_field(row, :sampling_role, "")))
        push!(get!(groups, key, Point[]), point)
    end
    Set(keys(groups)) == EXPECTED_POINTS || error("source points do not match fixed deep-oracle scope")
    return (; manifest, groups)
end

function _curve(points::Vector{Point})
    points = _deduplicate_points(points)
    length(points) >= MIN_SAMPLES || return nothing
    rho = Float64[point.rho for point in points]
    mu = Float64[point.mu for point in points]
    all(isfinite, rho) && all(isfinite, mu) || return nothing
    return (rho=rho, mu=mu, points=points)
end

function _layer_curves(points::Vector{Point})
    coarse_points = [point for point in points if point.rho_level == 0]
    fine_points = _deduplicate_points(points)
    return (; coarse=_curve(coarse_points), fine=_curve(fine_points))
end

function _evaluate(curve, tol::Float64)
    curve === nothing && return (
        status=:invalid, reason="solver_or_curve_failure", has_s_shape=false,
        sres=MODELS.SShapeResult(), maxwell=MODELS.MaxwellResult(),
        mu_transition=nothing, rho_hadron=nothing, rho_quark=nothing,
        area_residual=Inf, iterations=0,
    )
    sres = MODELS.detect_s_shape(curve.mu, curve.rho)
    if !sres.has_s_shape
        return (
            status=:invalid, reason="no_s_shape", has_s_shape=false, sres=sres,
            maxwell=MODELS.MaxwellResult(), mu_transition=nothing,
            rho_hadron=nothing, rho_quark=nothing, area_residual=Inf, iterations=0,
        )
    end
    mres = MODELS.maxwell_construction(curve.mu, curve.rho;
        spinodal_hint=sres, tol_area=tol, max_iter=MAX_ITER)
    area = mres.area_residual === nothing ? Inf : Float64(mres.area_residual)
    converged = mres.converged && isfinite(area) && area <= tol
    reason = if converged
        "ok"
    elseif !mres.converged
        String(get(mres.details, :reason, "maxwell_failed"))
    else
        "area_residual_above_solver_tol"
    end
    return (
        status=converged ? :valid : :invalid,
        reason=reason,
        has_s_shape=true,
        sres=sres,
        maxwell=mres,
        mu_transition=mres.mu_transition,
        rho_hadron=mres.rho_hadron,
        rho_quark=mres.rho_quark,
        area_residual=area,
        iterations=mres.iterations,
    )
end

function _topology(result)
    result.has_s_shape ? (result.status == :valid ? :first_order_candidate : :s_shape_unresolved) : :no_s_shape
end

function _geometry(coarse, fine)
    MODELS._compare_phase_geometry(coarse, fine, MODELS.PhaseGeometryTolerances(
        position_MeV=POSITION_TOL, density=DENSITY_TOL, maxwell_area=OUTER_AREA_GATE))
end

function _candidate_stable(coarse, fine)
    _topology(coarse) == _topology(fine) &&
        coarse.sres.derivative_sign_changes == fine.sres.derivative_sign_changes &&
        ((coarse.status == :valid && fine.status == :valid) ||
         (coarse.status == :invalid && fine.status == :invalid && !coarse.has_s_shape))
end

function _candidate_stable_across(left, right)
    _topology(left) == _topology(right) &&
        left.sres.derivative_sign_changes == right.sres.derivative_sign_changes &&
        ((left.status == :valid && right.status == :valid) ||
         (left.status == :invalid && right.status == :invalid && !left.has_s_shape))
end

function _point_tolerance_rows(key, curves)
    xi, T = key
    rows = NamedTuple[]
    previous = nothing
    for tol in AREA_TOLS
        coarse = _evaluate(curves.coarse, tol)
        fine = _evaluate(curves.fine, tol)
        geometry = _geometry(coarse, fine)
        pair_stable = _candidate_stable(coarse, fine)
        across_stable = previous === nothing ? true : _candidate_stable_across(previous.fine, fine)
        iteration_ok = coarse.iterations < MAX_ITER && fine.iterations < MAX_ITER
        strict_residual_ok = coarse.area_residual <= STRICT_SOLVER_TOL && fine.area_residual <= STRICT_SOLVER_TOL
        push!(rows, (
            xi=xi, T_MeV=T, solver_tol=tol,
            coarse_status=String(coarse.status), fine_status=String(fine.status),
            coarse_reason=coarse.reason, fine_reason=fine.reason,
            coarse_topology=String(_topology(coarse)), fine_topology=String(_topology(fine)),
            coarse_has_s_shape=coarse.has_s_shape, fine_has_s_shape=fine.has_s_shape,
            coarse_derivative_sign_changes=coarse.sres.derivative_sign_changes,
            fine_derivative_sign_changes=fine.sres.derivative_sign_changes,
            coarse_mu_transition=coarse.mu_transition === nothing ? missing : coarse.mu_transition,
            fine_mu_transition=fine.mu_transition === nothing ? missing : fine.mu_transition,
            coarse_rho_hadron=coarse.rho_hadron === nothing ? missing : coarse.rho_hadron,
            fine_rho_hadron=fine.rho_hadron === nothing ? missing : fine.rho_hadron,
            coarse_rho_quark=coarse.rho_quark === nothing ? missing : coarse.rho_quark,
            fine_rho_quark=fine.rho_quark === nothing ? missing : fine.rho_quark,
            coarse_area_residual=coarse.area_residual,
            fine_area_residual=fine.area_residual,
            coarse_iterations=coarse.iterations, fine_iterations=fine.iterations,
            geometry_converged=geometry.converged,
            position_error_MeV=geometry.position_MeV,
            density_error=geometry.density,
            maxwell_area_gate=geometry.maxwell_area,
            geometry_reason=geometry.reason,
            candidate_stable=pair_stable,
            candidate_stable_vs_previous=across_stable,
            iteration_budget_ok=iteration_ok,
            strict_residual_ok=strict_residual_ok,
            strict_gate=(pair_stable && across_stable && geometry.converged && iteration_ok && strict_residual_ok),
        ))
        previous = (; coarse, fine)
    end
    return rows
end

function _summarize_point(rows)
    strict = only(filter(row -> row.solver_tol == STRICT_SOLVER_TOL, rows))
    fine_tol = only(filter(row -> row.solver_tol == 1e-5, rows))
    return (
        xi=strict.xi, T_MeV=strict.T_MeV,
        strict_solver_tol=STRICT_SOLVER_TOL,
        strict_gate=strict.strict_gate,
        topology_stable=strict.candidate_stable && strict.candidate_stable_vs_previous,
        geometry_gate=strict.geometry_converged,
        area_residual_gate=strict.strict_residual_ok,
        iteration_budget_gate=strict.iteration_budget_ok,
        fine_tol_pair_stable=fine_tol.candidate_stable_vs_previous,
        strict_coarse_area_residual=strict.coarse_area_residual,
        strict_fine_area_residual=strict.fine_area_residual,
        strict_maxwell_area_gate=strict.maxwell_area_gate,
        strict_position_error_MeV=strict.position_error_MeV,
        strict_density_error=strict.density_error,
        strict_geometry_reason=strict.geometry_reason,
        verdict=strict.strict_gate ? "pass" : "inconclusive",
    )
end

function _write_csv(path::String, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _write_docs(output_dir::String, verdict::String, input_dir::String, summary_rows, frontier_rows)
    write(joinpath(output_dir, "README.md"), """# PNJL Maxwell tolerance-contract feasibility v1

verdict: `$(verdict)`。本目录使用 deep-oracle run `$(SOURCE_RUN_ID)` 的五个固定低温
anchor 做 solver-free post-processing replay；没有调用 equilibrium/Newton solver，也没有修改
production、reference 或历史 evidence。

扫描的内部 Maxwell bisection tolerance 为 `$(join(string.(AREA_TOLS), ", "))`。外层
coarse/fine geometry gate 固定为 position `$(POSITION_TOL) MeV`、density `$(DENSITY_TOL)`、
area `$(OUTER_AREA_GATE)`，因此本审计明确区分“求根停止条件”和“跨 rho 层收敛证书”。
严格候选要求 `$(STRICT_SOLVER_TOL)` 下两层 residual 均不超过该值、topology 不切换、
geometry gate 通过且二分未耗尽 `$(MAX_ITER)` 次迭代；`1e-5→5e-6` 只作为跨容差稳定性检查。

本 feasibility 不是 production gate，也不等于放宽 Maxwell area 容差。若 verdict 不是
`feasible_candidate`，不得创建 production tolerance-contract PR；应保留逐点曲线和失败原因
供物理/数值审核。
""")
    write(joinpath(output_dir, "AUDIT.md"), """# Audit boundary

输入文件的 calculation SHA、workflow head SHA、source manifest hash 和本地 producer head
写入 `manifest.json`。曲线只按已有 `(xi,T,rho)` 记录重建；没有插值、补点或 equilibrium
重算。`PhaseCore.maxwell_construction` 使用每个 sweep 的显式 `tol_area`，而
`PhaseGridConvergence._compare_phase_geometry` 使用独立的 outer area gate，并记录两者的
实际值。因此“area residual”与“geometry convergence”均可追溯，且不会因读取旧 slice
summary 而把其 gate 状态冒充为新证书。
""")
    _write_csv(joinpath(output_dir, "claim_ledger.csv"), [
        (claim_id="solver_free_replay", claim="No equilibrium/Newton solver was called", status="pass", boundary="postprocess-only"),
        (claim_id="strict_tolerance_contract", claim="Both rho layers satisfy the strict solver tolerance and outer geometry gate", status=verdict == "feasible_candidate" ? "pass" : "inconclusive", boundary="five deep-oracle anchors only"),
        (claim_id="cross_tolerance_stability", claim="The 1e-5 to 5e-6 transition does not switch the Maxwell candidate", status=all(row -> row.fine_tol_pair_stable, summary_rows) ? "pass" : "inconclusive", boundary="candidate identity, not physics promotion"),
        (claim_id="production_promotion", claim="Tolerance contract is ready for production integration", status="not_claimed", boundary="requires production CI and Actions revalidation"),
    ])
    _write_csv(joinpath(output_dir, "curve_index.csv"), [
        (source="deep_oracle_aggregate", source_run_id=SOURCE_RUN_ID,
         path="curve_points.csv", sha256=_sha256_file(joinpath(input_dir, "curve_points.csv")),
         raw_curve_copy_in_repository=false),
    ])
    open(joinpath(output_dir, "plot_manifest.json"), "w") do io
        JSON3.pretty(io, Dict(
            "schema_version" => "pnjl_maxwell_tolerance_contract_feasibility_v1",
            "figures" => Any[],
            "source_curve_artifact" => "external_actions_or_local_artifact",
        ))
        write(io, '\n')
    end
end

function _write_manifest(output_dir::String, input_dir::String, data, verdict::String, summary_rows, frontier_rows)
    files = Dict{String, String}()
    for (root, _, names) in walkdir(output_dir)
        for name in names
            name == "manifest.json" && continue
            path = joinpath(root, name)
            files[replace(relpath(path, output_dir), '\\' => '/')] = _sha256_file(path)
        end
    end
    manifest = Dict(
        "schema_version" => "pnjl_maxwell_tolerance_contract_feasibility_v1",
        "verdict" => verdict,
        "source_run_id" => SOURCE_RUN_ID,
        "source_calculation_sha" => SOURCE_CALCULATION_SHA,
        "source_workflow_head_sha" => SOURCE_WORKFLOW_HEAD_SHA,
        "source_manifest_sha256" => _sha256_file(joinpath(input_dir, "manifest.json")),
        "producer_head_sha" => try readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`) catch; "unknown" end,
        "solver_called" => false,
        "area_tolerances_swept" => collect(AREA_TOLS),
        "strict_solver_tolerance" => STRICT_SOLVER_TOL,
        "outer_geometry_gate" => Dict("position_MeV" => POSITION_TOL,
            "density" => DENSITY_TOL, "maxwell_area" => OUTER_AREA_GATE),
        "maxwell_max_iter" => MAX_ITER,
        "point_count" => length(summary_rows),
        "summary" => summary_rows,
        "frontier" => frontier_rows,
        "files" => files,
        "promotion_boundary" => "diagnostic_only_until_production_contract_and_actions_revalidation",
    )
    open(joinpath(output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, manifest)
        write(io, '\n')
    end
    return manifest
end

function run(input_dir::String, output_dir::String)
    data = _load_input(input_dir)
    mkpath(output_dir)
    all_tolerance_rows = NamedTuple[]
    summary_rows = NamedTuple[]
    frontier_rows = NamedTuple[]
    for key in sort!(collect(keys(data.groups)); by=x -> (x[1], x[2]))
        curves = _layer_curves(data.groups[key])
        rows = _point_tolerance_rows(key, curves)
        append!(all_tolerance_rows, rows)
        push!(summary_rows, _summarize_point(rows))
        for row in rows
            push!(frontier_rows, row)
        end
    end
    strict_pass = all(row -> row.strict_gate, summary_rows)
    stable_pass = all(row -> row.fine_tol_pair_stable, summary_rows)
    verdict = strict_pass && stable_pass ? "feasible_candidate" : "maxwell_tolerance_inconclusive"
    _write_csv(joinpath(output_dir, "tolerance_frontier.csv"), frontier_rows)
    _write_csv(joinpath(output_dir, "point_results.csv"), summary_rows)
    _write_docs(output_dir, verdict, input_dir, summary_rows, frontier_rows)
    manifest = _write_manifest(output_dir, input_dir, data, verdict, summary_rows, frontier_rows)
    return (; verdict, manifest, summary=summary_rows, frontier=frontier_rows)
end

function _parse_args(args)
    input_dir = nothing
    output_dir = joinpath(PROJECT_ROOT, "docs", "analysis", "pnjl_maxwell_tolerance_contract_feasibility_v1")
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
    options === nothing && (println("Usage: julia --project=. scripts/analysis/pnjl_maxwell_tolerance_feasibility.jl --input-dir=PATH [--output-dir=PATH]"); return 0)
    result = run(options.input_dir, options.output_dir)
    println(JSON3.write(Dict("verdict" => result.verdict,
        "point_count" => length(result.summary), "solver_called" => false)))
    result.verdict == "feasible_candidate" ? 0 : 2
end

end # module MaxwellToleranceFeasibility

if abspath(PROGRAM_FILE) == @__FILE__
    exit(MaxwellToleranceFeasibility.main())
end
