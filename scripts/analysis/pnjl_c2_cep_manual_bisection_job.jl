#!/usr/bin/env julia

"""Run fixed low/mid/high temperature slices for manual CEP bisection review.

The calculation path is reused from the immutable three-midpoint job. This
wrapper only changes the temperature plan; it never uses oracle labels to
select or omit a slice. Every output remains diagnostic-only until author
review updates the bracket trace.
"""

using CSV
using JSON3
using SHA

# Reuse the existing production-parity slice implementation. The included
# script does not enter its own main block because PROGRAM_FILE is this file.
include(joinpath(@__DIR__, "pnjl_c2_cep_limited_feasibility_job.jl"))

const MANUAL_SCHEMA = "pnjl_c2_cep_manual_bisection_job_v1"
const PLAN_SCHEMA = "pnjl_c2_cep_manual_bisection_plan_v1"
const EXPECTED_XI = (0.125, 0.39375, 0.5)
const EXPECTED_WIDTH = 0.125

function _read_plan(path, xi)
    isfile(path) || error("missing manual bisection plan: $(path)")
    rows = collect(CSV.File(path))
    selected = NamedTuple[]
    for row in rows
        row_xi = _safe_float(_field(row, :xi))
        row_xi == xi || continue
        role = String(_field(row, :role, ""))
        sequence = Int(round(_safe_float(_field(row, :sequence, NaN))))
        T = _safe_float(_field(row, :T_MeV))
        low = _safe_float(_field(row, :bracket_low_MeV))
        high = _safe_float(_field(row, :bracket_high_MeV))
        all(isfinite, (row_xi, T, low, high)) || error("non-finite plan row: $(row)")
        high > low || error("invalid plan bracket: $(row)")
        push!(selected, (xi=xi, sequence=sequence, role=role, T_MeV=T,
            bracket_low_MeV=low, bracket_high_MeV=high))
    end
    sort!(selected; by=row -> row.sequence)
    length(selected) == 3 || error("manual plan must contain three rows for xi=$(xi)")
    [row.role for row in selected] == ["low", "midpoint", "high"] ||
        error("manual plan roles are not low/midpoint/high for xi=$(xi)")
    selected[1].bracket_low_MeV == selected[2].bracket_low_MeV == selected[3].bracket_low_MeV ||
        error("manual plan low bracket mismatch for xi=$(xi)")
    selected[1].bracket_high_MeV == selected[2].bracket_high_MeV == selected[3].bracket_high_MeV ||
        error("manual plan high bracket mismatch for xi=$(xi)")
    isapprox(selected[2].T_MeV, 0.5 * (selected[1].T_MeV + selected[3].T_MeV);
        atol=1e-10, rtol=0.0) || error("manual midpoint mismatch for xi=$(xi)")
    isapprox(selected[3].T_MeV - selected[1].T_MeV, EXPECTED_WIDTH;
        atol=1e-10, rtol=0.0) || error("unexpected initial bracket width for xi=$(xi)")
    selected
end

function _write_manual_csv(path, rows)
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _flatten_cost(cfg, plan_row, costs)
    output = NamedTuple[]
    for item in costs
        snapshot = item.telemetry
        push!(output, (
            xi=cfg.xi, T_MeV=plan_row.T_MeV, role=plan_row.role,
            method=item.method, calculation_sha=cfg.calculation_sha,
            unique_solves=item.cache.unique_solves,
            point_requests=item.cache.point_requests,
            requested_point_calls=item.cache.point_requests,
            cache_hits=item.cache.cache_hits,
            failed_points=item.cache.failed_points,
            equilibrium_requests=snapshot.equilibrium_requests,
            fixedrho_requests=snapshot.fixedrho_requests,
            residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
            jacobian_calls=snapshot.nlsolve_g_calls,
            newton_iterations=snapshot.newton_iterations,
            runner_seconds=item.runner_seconds,
            fallback_count=snapshot.root_fallbacks,
            retry_count=snapshot.scan_retries,
        ))
    end
    output
end

function _run_manual(cfg, plan_path)
    plan = _read_plan(plan_path, cfg.xi)
    mkpath(cfg.output_dir)
    started = time_ns()
    curve_rows = NamedTuple[]
    slice_rows = NamedTuple[]
    cost_rows = NamedTuple[]
    materialization_rows = Any[]
    for plan_row in plan
        slice_dir = joinpath(cfg.output_dir, "slices",
            "T_$(replace(string(plan_row.T_MeV), "." => "p"))")
        rows, hybrid, costs, materialization = _run_slice(cfg, plan_row.T_MeV, slice_dir)
        append!(curve_rows, rows)
        append!(cost_rows, _flatten_cost(cfg, plan_row, costs))
        push!(materialization_rows, Dict(
            "xi" => cfg.xi, "T_MeV" => plan_row.T_MeV, "role" => plan_row.role,
            "path" => replace(relpath(materialization.path, cfg.output_dir), '\\' => '/'),
            "provenance_path" => replace(relpath(materialization.provenance_path,
                cfg.output_dir), '\\' => '/'),
            "sha256" => materialization.sha256, "rows" => materialization.rows,
            "recovered_rows" => materialization.recovered_rows,
            "aggregate_rows" => materialization.aggregate_rows,
            "aggregate_parse_error" => materialization.aggregate_parse_error,
            "source_files" => materialization.source_files,
        ))
        push!(slice_rows, (
            xi=cfg.xi, sequence=plan_row.sequence, role=plan_row.role,
            T_MeV=plan_row.T_MeV, bracket_low_MeV=plan_row.bracket_low_MeV,
            bracket_high_MeV=plan_row.bracket_high_MeV,
            bracket_width_MeV=plan_row.bracket_high_MeV - plan_row.bracket_low_MeV,
            hybrid_status=hybrid.status, hybrid_reason=hybrid.reason,
            manual_decision="manual_pending", oracle_labels_used_for_routing=false,
            oracle_materialized_rows=materialization.rows,
            oracle_recovered_rows=materialization.recovered_rows,
            oracle_materialized_sha256=materialization.sha256,
            oracle_materialization_provenance=replace(
                relpath(materialization.provenance_path, cfg.output_dir), '\\' => '/'),
            rho_rows=length(rows), finite_and_converged=true,
        ))
    end
    sort!(curve_rows; by=row -> (row.T_MeV, row.rho))
    _validate_manual_cost_rows(cost_rows)
    _write_manual_csv(joinpath(cfg.output_dir, "fine_pool.csv"), curve_rows)
    _write_manual_csv(joinpath(cfg.output_dir, "manual_bisection_trace.csv"), slice_rows)
    _write_manual_csv(joinpath(cfg.output_dir, "method_costs.csv"), cost_rows)
    elapsed = (time_ns() - started) / 1e9
    files = Dict(name => _sha(joinpath(cfg.output_dir, name)) for name in
        ("fine_pool.csv", "manual_bisection_trace.csv", "method_costs.csv"))
    plan_sha = bytes2hex(sha256(read(plan_path)))
    summary = Dict(
        "schema_version" => MANUAL_SCHEMA, "plan_schema_version" => PLAN_SCHEMA,
        "scope" => cfg.scope, "xi" => cfg.xi, "method" => "hybrid_and_oracle_fine_pool",
        "tag" => cfg.tag, "source_run_id" => cfg.source_run_id,
        "calculation_sha" => cfg.calculation_sha, "postprocess_sha" => cfg.postprocess_sha,
        "workflow_head_sha" => cfg.postprocess_sha, "rho_fine_step" => RHO_FINE_STEP,
        "rho_max" => RHO_MAX, "manual_decision_required" => true,
        "oracle_labels_used_for_routing" => false, "solver_called" => true,
        "finite_and_converged_final" => true, "initial_bracket_width_MeV" => EXPECTED_WIDTH,
        "plan_sha256" => plan_sha, "temperatures" => [row.T_MeV for row in plan],
        "roles" => [row.role for row in plan], "files" => files,
        "materialization" => materialization_rows, "runner_seconds" => elapsed,
        "cache" => Dict("failed_points" => sum(row.failed_points for row in cost_rows),
            "unique_solves" => sum(row.unique_solves for row in cost_rows),
            "point_requests" => sum(row.point_requests for row in cost_rows),
            "cache_hits" => sum(row.cache_hits for row in cost_rows)),
    )
    open(joinpath(cfg.output_dir, "job_summary.json"), "w") do io
        JSON3.write(io, summary); write(io, '\n')
    end
    open(joinpath(cfg.output_dir, "manifest.json"), "w") do io
        JSON3.pretty(io, merge(summary, Dict("provenance" => Dict(
            "calculation_sha" => cfg.calculation_sha, "postprocess_sha" => cfg.postprocess_sha,
            "reference_write" => false, "manual_route" => true,
            "oracle_labels_used_for_routing" => false,
            "initial_bracket_contract" => "C2 frozen width 0.125 MeV"))))
        write(io, '\n')
    end
    println(JSON3.write(summary))
end

function _validate_manual_cost_rows(cost_rows)
    all(item -> item.failed_points == 0, cost_rows) ||
        error("manual CEP shard has failed solver points")
    all(item -> item.point_requests == item.unique_solves + item.cache_hits, cost_rows) ||
        error("manual CEP shard cache cost does not reconcile")
    true
end

if abspath(PROGRAM_FILE) == @__FILE__
    cfg = _config(ARGS)
    plan_path = abspath(String(_arg(ARGS, "--plan-file", "")))
    _run_manual(cfg, plan_path)
end
