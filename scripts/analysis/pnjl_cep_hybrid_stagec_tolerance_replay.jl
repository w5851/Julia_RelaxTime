#!/usr/bin/env julia

"""Solver-free Stage-C feasibility replay with Maxwell-contract evidence.

The base input is the final 24-anchor hybrid-shadow aggregate.  The second
input is the immutable five-anchor Maxwell-tolerance revalidation aggregate.
For those five anchors, the revalidation rows replace overlapping oracle
points and add the local 0.003125/0.0015625 verification pool.  The base
Stage-A/B and dense artifacts remain the cost reference; no equilibrium
solver is called and no raw curve copy is committed.
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using CSV
using JSON3
using SHA
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const V2_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "analysis",
    "pnjl_cep_hybrid_stagec_certificate_feasibility_v2.jl")
if !isdefined(Main, :StageCCertificateFeasibilityV2)
    include(V2_SCRIPT)
end
const V2 = Main.StageCCertificateFeasibilityV2

const BASE_RUN_ID = "30601857594"
const BASE_CALCULATION_SHA = "fde2b929b60575f1daacb84a1b9b8ff6e3b0a0cc"
const TOLERANCE_RUN_ID = "30730990835"
const TOLERANCE_CALCULATION_SHA = "467be1fce847a9c991ec362c3335be07fccbe604"
const TOLERANCE_COARSE_STEP = 0.003125
const TOLERANCE_FINE_STEP = 0.0015625
const REVALIDATION_ANCHORS = Set([
    (-0.5, 20.0),
    (0.0, 5.0),
    (0.0, 20.0),
    (0.5, 5.0),
    (0.5, 20.0),
])

@inline _finite(value) = value isa Real && isfinite(Float64(value))
_f(value, default=NaN) = V2._float(value, default)
_s(row, name::Symbol, default=nothing) = V2._field(row, name, default)
_key_float(value) = try parse(Float64, string(value)) catch; NaN end

function _sha256_file(path::String)
    bytes2hex(open(sha256, path))
end

function _load_bundle(input_dir::String)
    required = ("manifest.json", "curve_points.csv", "slice_metrics.csv", "method_costs.csv")
    missing = [name for name in required if !isfile(joinpath(input_dir, name))]
    isempty(missing) || error("missing input files in $(input_dir): $(join(missing, ", "))")
    (; manifest=JSON3.read(read(joinpath(input_dir, "manifest.json"), String)),
        curves=collect(CSV.File(joinpath(input_dir, "curve_points.csv"))),
        slices=collect(CSV.File(joinpath(input_dir, "slice_metrics.csv"))),
        costs=collect(CSV.File(joinpath(input_dir, "method_costs.csv"))))
end

function _manifest_value(manifest, name::Symbol, default=nothing)
    V2._manifest_value(manifest, name, default)
end

function _verify_manifest_files(input_dir::String, manifest)
    files = _manifest_value(manifest, :file_sha256, nothing)
    files === nothing && (files = _manifest_value(manifest, :files, nothing))
    files === nothing && error("manifest has no file hash map: $(input_dir)")
    for (name, expected) in pairs(files)
        path = joinpath(input_dir, string(name))
        isfile(path) || error("manifest file is missing: $(path)")
        actual = _sha256_file(path)
        actual == string(expected) || error("manifest hash mismatch: $(path)")
    end
    true
end

function _anchor(row)
    (_f(_s(row, :xi)), _f(_s(row, :T_MeV)))
end

function _is_anchor(row, anchor)
    x, T = _anchor(row)
    isapprox(x, anchor[1]; atol=1e-8, rtol=0.0) &&
        isapprox(T, anchor[2]; atol=1e-8, rtol=0.0)
end

function _validate_curve_rows(curves, label::String)
    seen = Set{Tuple{String, String, String, String}}()
    for row in curves
        key = (
            string(_s(row, :xi, "")),
            string(_s(row, :method, "")),
            string(_s(row, :T_MeV, "")),
            string(_s(row, :rho, "")),
        )
        key in seen && error("duplicate $(label) curve key: $(key)")
        push!(seen, key)
        V2._bool(_s(row, :converged, false)) || error("non-converged $(label) curve: $(key)")
        V2._bool(_s(row, :finite, false)) || error("non-finite $(label) curve: $(key)")
        _finite(_f(_s(row, :rho))) || error("invalid $(label) rho: $(key)")
        _finite(_f(_s(row, :muq_MeV))) || error("invalid $(label) mu: $(key)")
    end
    return length(seen)
end

function _validate_base(base)
    _verify_manifest_files(base.dir, base.manifest)
    string(_manifest_value(base.manifest, :expected_calculation_sha, "")) == BASE_CALCULATION_SHA ||
        error("base calculation SHA mismatch")
    string(_manifest_value(base.manifest, :evidence_state, "")) == "final" ||
        error("base aggregate is not final")
    actions = _manifest_value(base.manifest, :actions, nothing)
    actions === nothing && error("base aggregate has no Actions provenance")
    Bool(_manifest_value(actions, :source_run_completed_success, false)) ||
        error("base source run is not successful")
    string(_manifest_value(actions, :run_id, "")) == BASE_RUN_ID ||
        error("base source run mismatch")
    _validate_curve_rows(base.curves, "base")
end

function _validate_tolerance(tolerance)
    _verify_manifest_files(tolerance.dir, tolerance.manifest)
    string(_manifest_value(tolerance.manifest, :calculation_sha, "")) == TOLERANCE_CALCULATION_SHA ||
        error("tolerance calculation SHA mismatch")
    string(_manifest_value(tolerance.manifest, :workflow_head_sha, "")) == TOLERANCE_CALCULATION_SHA ||
        error("tolerance workflow SHA mismatch")
    string(_manifest_value(tolerance.manifest, :status, "")) == "complete" ||
        error("tolerance aggregate is not complete")
    string(_manifest_value(tolerance.manifest, :physical_verdict, "")) == "author_review_required" ||
        error("unexpected tolerance physical verdict")
    summary_path = joinpath(tolerance.dir, "deep_oracle_summary.json")
    summary = JSON3.read(read(summary_path, String))
    Int(_manifest_value(summary, :job_count, 0)) == length(REVALIDATION_ANCHORS) ||
        error("tolerance job count mismatch")
    _manifest_value(summary, :errors, Any[]) == Any[] || error("tolerance aggregate has errors")
    _validate_curve_rows(tolerance.curves, "tolerance")
    anchors = Set(_anchor(row) for row in tolerance.slices if
        string(_s(row, :method, "")) == "independent_oracle")
    anchors == REVALIDATION_ANCHORS || error("tolerance anchor set mismatch: $(anchors)")
    for row in tolerance.slices
        _is_anchor(row, _anchor(row)) || error("invalid tolerance slice anchor")
    end
    true
end

function _curve_key(row)
    (string(_s(row, :method, "")), string(_s(row, :xi, "")),
        string(_s(row, :T_MeV, "")), string(_s(row, :rho, "")))
end

function _merge_curves(base_curves, tolerance_curves)
    merged = Dict{Tuple{String, String, String, String}, Any}()
    for row in base_curves
        merged[_curve_key(row)] = row
    end
    for row in tolerance_curves
        _s(row, :method, "") == "independent_oracle" || continue
        _is_anchor(row, _anchor(row)) || continue
        merged[_curve_key(row)] = row
    end
    collect(values(merged))
end

function _merge_slices(base_slices, tolerance_slices)
    rows = Any[row for row in base_slices if !(
        string(_s(row, :method, "")) == "independent_oracle" &&
        _anchor(row) in REVALIDATION_ANCHORS)]
    append!(rows, Any[row for row in tolerance_slices if
        string(_s(row, :method, "")) == "independent_oracle"])
    rows
end

function _revalidation_certificate(tolerance_curves, anchor)
    on_rounded_grid = (rho, step) -> isapprox(rho, step * round(rho / step); atol=1.1e-6, rtol=0.0)
    coarse = V2._evaluate(V2._curve(V2._point_rows(tolerance_curves, "independent_oracle",
        anchor[1], anchor[2], rho -> on_rounded_grid(rho, TOLERANCE_COARSE_STEP))))
    fine = V2._evaluate(V2._curve(V2._point_rows(tolerance_curves, "independent_oracle",
        anchor[1], anchor[2], rho -> on_rounded_grid(rho, TOLERANCE_FINE_STEP))))
    geometry = V2._geometry(coarse, fine)
    semantic = V2._semantic(coarse, fine, geometry)
    (
        xi=anchor[1], T_MeV=anchor[2],
        coarse_status=String(coarse.status), fine_status=String(fine.status),
        result_status=String(semantic),
        coarse_area_residual=coarse.area_residual,
        fine_area_residual=fine.area_residual,
        position_error_MeV=geometry.position_MeV,
        density_error=geometry.density,
        maxwell_area_gate=geometry.maxwell_area,
        geometry_converged=geometry.converged,
        coarse_reason=coarse.reason, fine_reason=fine.reason,
    )
end

function _numeric_identity(base_curves, tolerance_curves)
    base = Dict{Tuple{String, String, String, String}, Any}()
    for row in base_curves
        _s(row, :method, "") == "independent_oracle" || continue
        _is_anchor(row, _anchor(row)) || continue
        base[_curve_key(row)] = row
    end
    tolerance = Dict{Tuple{String, String, String, String}, Any}()
    for row in tolerance_curves
        _s(row, :method, "") == "independent_oracle" || continue
        tolerance[_curve_key(row)] = row
    end
    common = intersect(keys(base), keys(tolerance))
    positive_common = [key for key in common if _key_float(key[4]) > 0.0]
    zero_common = [key for key in common if _key_float(key[4]) == 0.0]
    max_mu = 0.0
    max_zero_mu = 0.0
    max_pressure = 0.0
    max_residual = 0.0
    iteration_differences = 0
    for key in common
        a, b = base[key], tolerance[key]
        delta_mu = abs(_f(_s(a, :muq_MeV)) - _f(_s(b, :muq_MeV)))
        if _key_float(key[4]) == 0.0
            max_zero_mu = max(max_zero_mu, delta_mu)
        else
            max_mu = max(max_mu, delta_mu)
        end
        max_pressure = max(max_pressure, abs(_f(_s(a, :pressure_fm4)) - _f(_s(b, :pressure_fm4))))
        max_residual = max(max_residual, abs(_f(_s(a, :residual_norm)) - _f(_s(b, :residual_norm))))
        string(_s(a, :iterations, "")) == string(_s(b, :iterations, "")) || (iteration_differences += 1)
    end
    (
        common_keys=length(common),
        positive_rho_common_keys=length(positive_common),
        rho_zero_common_keys=length(zero_common),
        max_abs_muq_MeV_positive_rho=max_mu,
        max_abs_muq_MeV_rho_zero=max_zero_mu,
        max_abs_pressure_fm4=max_pressure,
        max_abs_residual_norm=max_residual,
        iteration_differences=iteration_differences,
        numeric_identity=(max_mu <= 1e-6 && max_pressure <= 1e-12 && max_residual <= 1e-12),
    )
end

function _write_csv(path::String, rows)
    mkpath(dirname(path))
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _frontier_rows(frontier)
    frontier.frontier
end

function _replay_rows(frontier)
    rows = NamedTuple[]
    for cap in V2.CAPS
        for result in frontier.replays[cap]
            push!(rows, (
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
        end
    end
    rows
end

function _candidate_rows(frontier)
    rows = NamedTuple[]
    for cap in V2.CAPS
        for result in frontier.replays[cap]
            for candidate in result.candidates
                push!(rows, (
                    xi=result.xi, T_MeV=result.T_MeV, cap=cap, level=String(candidate.level),
                    rho_low=candidate.rho_low, rho_high=candidate.rho_high,
                    width=candidate.width, drop_mu=candidate.drop_mu,
                    negative_secants=candidate.negative_secants,
                ))
            end
        end
    end
    rows
end

function _write_text(path::String, text::String)
    mkpath(dirname(path))
    write(path, text)
end

function _write_evidence(output_dir::String, frontier, certificates, identity, base, tolerance,
        base_dir::String, tolerance_dir::String, merged_curve_count::Int)
    mkpath(output_dir)
    selected = frontier.selected
    ambiguous_confirmation = any(result -> result.oracle_status == "ambiguous_near_critical" &&
        result.result_status != :ambiguous_near_critical, frontier.replays[first(V2.CAPS)])
    certificate_pass = all(row -> row.result_status == "confirmed_first_order" &&
        row.geometry_converged && row.maxwell_area_gate <= V2.AREA_TOL &&
        row.position_error_MeV <= V2.POSITION_TOL && row.density_error <= V2.DENSITY_TOL, certificates)
    verdict = if selected === nothing
        ambiguous_confirmation ? "oracle_inconclusive" :
        any(row -> !row.state_gate, frontier.frontier) ? "integration_failed" :
        any(row -> !row.candidate_gate || row.multiple_candidate_anchors > 0, frontier.frontier) ? "maxwell_candidate_inconclusive" :
        any(row -> !row.cost_gate, frontier.frontier) ? "performance_inconclusive" : "integration_failed"
    else
        "feasible_candidate"
    end
    certificate_pass || (verdict = "oracle_inconclusive")
    policy = (
        schema_version="cep_hybrid_stagec_tolerance_replay_v1",
        verdict=verdict,
        selected_cap=selected === nothing ? nothing : frontier.frontier[selected].cap,
        caps=V2.CAPS,
        dense_unique_solves=frontier.dense,
        stage_b_grid_step=V2.STAGE_B_FINE,
        stage_c_local_step=V2.STAGE_C_FINE,
        maxwell_solver_tol=V2.MAXWELL_SOLVER_TOL,
        outer_maxwell_area_tol=V2.AREA_TOL,
        solver_called=false,
        reason=verdict == "feasible_candidate" ? "minimum feasible cap" : "revalidation gate not complete",
    )

    _write_csv(joinpath(output_dir, "anchor_replay.csv"), _replay_rows(frontier))
    _write_csv(joinpath(output_dir, "candidate_runs.csv"), _candidate_rows(frontier))
    _write_csv(joinpath(output_dir, "cost_frontier.csv"), _frontier_rows(frontier))
    _write_csv(joinpath(output_dir, "tolerance_anchor_certificates.csv"), certificates)
    _write_csv(joinpath(output_dir, "curve_identity_audit.csv"), [identity])
    _write_csv(joinpath(output_dir, "selected_points.csv"), NamedTuple[])
    _write_csv(joinpath(output_dir, "claim_ledger.csv"), [
        (claim_id="tolerance_certificate", claim="five revalidated anchors have first-order Maxwell/geometry certificates", status=certificate_pass ? "pass" : "inconclusive", boundary="fixed anchors only"),
        (claim_id="classification", claim="Stage-C replay satisfies all 24-anchor classification gates", status=verdict == "feasible_candidate" ? "pass" : "inconclusive", boundary="solver-free replay"),
        (claim_id="cost", claim="simulated Stage-A/B/C union is below memoized dense", status=all(row -> row.cost_gate, frontier.frontier) ? "pass" : "inconclusive", boundary="runner telemetry requires Actions"),
        (claim_id="promotion", claim="replay authorizes production/reference promotion", status="not_claimed", boundary="targeted/full shadow and author review required"),
    ])
    _write_text(joinpath(output_dir, "README.md"), """# PNJL Hybrid Stage-C tolerance-contract replay v1

verdict: `$(verdict)`。本目录将已完成的 24-anchor hybrid shadow final aggregate 与
Maxwell tolerance-contract 五点 revalidation 做 solver-free 合并重放；没有调用
equilibrium solver，不修改 production、reference 或历史 v1/v2 evidence。

- base source run: `$(BASE_RUN_ID)` / `$(BASE_CALCULATION_SHA)`
- tolerance revalidation run: `$(TOLERANCE_RUN_ID)` / `$(TOLERANCE_CALCULATION_SHA)`
- internal Maxwell solver tolerance: `$(V2.MAXWELL_SOLVER_TOL)`
- outer Maxwell geometry gate: `$(V2.AREA_TOL)`
- tested caps: `$(join(string.(V2.CAPS), ", "))`
- simulated dense reference: `$(frontier.dense)` unique solves
- solver called: `false`

五个 revalidated anchor 的局部 certificate 与 base Stage-A/B 证据分开记录；新深层
`0.0015625` 曲线只作为独立收敛证据，不被免费计入 hybrid cost。该 replay 仍是
candidate/diagnostic evidence，不能直接晋升 production/reference。
""")
    _write_text(joinpath(output_dir, "AUDIT.md"), """# Audit boundary

输入包括 base final aggregate 与五点 tolerance-contract aggregate。两套 manifest、文件
哈希、曲线键、finite/converged 状态和双 SHA provenance 在执行时验证；共同曲线数值字段
另有 `curve_identity_audit.csv` 记录。rho=0 的化学势端点不是唯一可识别的
thermodynamic root，重放保留该差异为 provenance warning，并以 tolerance artifact 覆盖
重叠端点；正密度曲线才用于数值 identity gate。所有 Maxwell/geometry 证书由当前 Julia PhaseCore
与 outer geometry comparison 重算，内部 solver tol 与外层 gate 分离。

本次是 solver-free replay：不估算 residual/Jacobian/runner-minute 成本，不复制 raw
`curve_points.csv` 到仓库。若 verdict 为 `feasible_candidate`，仍必须由 targeted/full
Actions shadow 验证，再由作者审核物理曲线和成本。
""")
    open(joinpath(output_dir, "selected_policy.json"), "w") do io
        JSON3.pretty(io, policy)
        write(io, '\n')
    end
    open(joinpath(output_dir, "provenance.json"), "w") do io
        JSON3.pretty(io, Dict(
            "schema_version" => "cep_hybrid_stagec_tolerance_replay_v1",
            "base_run_id" => BASE_RUN_ID,
            "base_calculation_sha" => BASE_CALCULATION_SHA,
            "base_manifest_sha256" => _sha256_file(joinpath(base_dir, "manifest.json")),
            "tolerance_run_id" => TOLERANCE_RUN_ID,
            "tolerance_calculation_sha" => TOLERANCE_CALCULATION_SHA,
            "tolerance_manifest_sha256" => _sha256_file(joinpath(tolerance_dir, "manifest.json")),
            "producer_head_sha" => try readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`) catch; "unknown" end,
            "merged_curve_count" => merged_curve_count,
            "solver_called" => false,
            "reference_write" => false,
        ))
        write(io, '\n')
    end
    _write_csv(joinpath(output_dir, "curve_index.csv"), [
        (source="base_final_aggregate", path="external:base/curve_points.csv", sha256=_sha256_file(joinpath(base_dir, "curve_points.csv")), raw_curve_copy_in_repository=false),
        (source="tolerance_revalidation_aggregate", path="external:tolerance/curve_points.csv", sha256=_sha256_file(joinpath(tolerance_dir, "curve_points.csv")), raw_curve_copy_in_repository=false),
    ])
    open(joinpath(output_dir, "plot_manifest.json"), "w") do io
        JSON3.pretty(io, Dict("schema_version" => "cep_hybrid_stagec_tolerance_replay_v1", "figures" => Any[], "raw_curves_external" => true))
        write(io, '\n')
    end
    files = Dict{String, String}()
    for (root, _, names) in walkdir(output_dir)
        for name in names
            name == "manifest.json" && continue
            full = joinpath(root, name)
            files[replace(relpath(full, output_dir), '\\' => '/')] = _sha256_file(full)
        end
    end
    manifest = Dict(
        "schema_version" => "cep_hybrid_stagec_tolerance_replay_v1",
        "verdict" => verdict,
        "selected_policy" => policy,
        "base_run_id" => BASE_RUN_ID,
        "base_calculation_sha" => BASE_CALCULATION_SHA,
        "base_manifest_sha256" => _sha256_file(joinpath(base_dir, "manifest.json")),
        "tolerance_run_id" => TOLERANCE_RUN_ID,
        "tolerance_calculation_sha" => TOLERANCE_CALCULATION_SHA,
        "tolerance_manifest_sha256" => _sha256_file(joinpath(tolerance_dir, "manifest.json")),
        "producer_head_sha" => try readchomp(`git -C $PROJECT_ROOT rev-parse HEAD`) catch; "unknown" end,
        "solver_called" => false,
        "reference_write" => false,
        "merged_curve_count" => merged_curve_count,
        "files" => files,
    )
    open(joinpath(output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, manifest)
        write(io, '\n')
    end
    return policy
end

function _parse_tolerance_replay_args(args)
    base_dir = nothing
    tolerance_dir = nothing
    output_dir = joinpath(PROJECT_ROOT, "docs", "analysis", "pnjl",
        "cep_hybrid_stagec_tolerance_replay_v1")
    for arg in args
        startswith(arg, "--base-input-dir=") && (base_dir = abspath(split(arg, "="; limit=2)[2]))
        startswith(arg, "--tolerance-input-dir=") && (tolerance_dir = abspath(split(arg, "="; limit=2)[2]))
        startswith(arg, "--output-dir=") && (output_dir = abspath(split(arg, "="; limit=2)[2]))
        arg in ("-h", "--help") && return nothing
    end
    base_dir === nothing && throw(ArgumentError("--base-input-dir=PATH is required"))
    tolerance_dir === nothing && throw(ArgumentError("--tolerance-input-dir=PATH is required"))
    return (; base_dir, tolerance_dir, output_dir)
end

function main(args=ARGS)
    options = _parse_tolerance_replay_args(args)
    options === nothing && (println("Usage: julia --project=. scripts/analysis/pnjl_cep_hybrid_stagec_tolerance_replay.jl --base-input-dir=PATH --tolerance-input-dir=PATH [--output-dir=PATH]"); return 0)
    base = merge(_load_bundle(options.base_dir), (dir=options.base_dir,))
    tolerance = merge(_load_bundle(options.tolerance_dir), (dir=options.tolerance_dir,))
    _validate_base(base)
    _validate_tolerance(tolerance)
    merged_curves = _merge_curves(base.curves, tolerance.curves)
    merged_slices = _merge_slices(base.slices, tolerance.slices)
    anchors = V2._anchor_set(merged_slices)
    length(anchors) == 24 || error("expected 24 merged anchors, got $(length(anchors))")
    certificates = [_revalidation_certificate(tolerance.curves, anchor) for anchor in sort(collect(REVALIDATION_ANCHORS))]
    identity = _numeric_identity(base.curves, tolerance.curves)
    identity.numeric_identity || error("positive-rho base/tolerance curve values differ beyond replay tolerance")
    frontier = V2._cost_frontier(merged_curves, merged_slices, base.costs, anchors)
    policy = _write_evidence(options.output_dir, frontier, certificates, identity, base, tolerance,
        options.base_dir, options.tolerance_dir, length(merged_curves))
    println(JSON3.write(policy))
    policy.verdict == "feasible_candidate" ? 0 : 2
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
