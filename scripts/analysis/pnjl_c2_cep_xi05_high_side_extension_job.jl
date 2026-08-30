#!/usr/bin/env julia

"""Run the fixed high-temperature extension for the xi=0.5 CEP audit.

The underlying slice implementation is included from the frozen limited-
feasibility job. This wrapper changes only the versioned temperature plan and
records the equal-step extension provenance. It never uses oracle labels or
author decisions to select a temperature.
"""

using CSV
using JSON3
using SHA

# The calculation worktree receives this wrapper and the shared slice job from
# the workflow head. Including the shared job keeps the numerical path exactly
# production-parity while this file owns only the extension contract.
include(joinpath(@__DIR__, "pnjl_c2_cep_limited_feasibility_job.jl"))

const HIGH_SIDE_SCHEMA = "pnjl_c2_cep_xi05_high_side_extension_job_v1"
const HIGH_SIDE_PLAN_SCHEMA = "pnjl_c2_cep_xi05_high_side_extension_plan_v1"
const HIGH_SIDE_CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
const HIGH_SIDE_XI = 0.5
const HIGH_SIDE_ANCHOR_T = 107.0625
const HIGH_SIDE_STEP = 0.0625
const HIGH_SIDE_TEMPERATURES = (107.125, 107.1875, 107.25)

@inline function _high_field(row, name::Symbol, default=nothing)
    row === nothing && return default
    try
        if row isa AbstractDict
            haskey(row, name) && return row[name]
            haskey(row, String(name)) && return row[String(name)]
            return default
        end
        hasproperty(row, name) ? getproperty(row, name) : default
    catch
        default
    end
end

function _read_high_side_plan(path)
    isfile(path) || error("missing high-side extension plan: $(path)")
    rows = collect(CSV.File(path))
    length(rows) == length(HIGH_SIDE_TEMPERATURES) ||
        error("high-side extension plan must contain exactly three rows")
    selected = NamedTuple[]
    for row in rows
        xi = _safe_float(_high_field(row, :xi))
        sequence = Int(round(_safe_float(_high_field(row, :sequence, NaN))))
        temperature = _safe_float(_high_field(row, :T_MeV))
        anchor = _safe_float(_high_field(row, :anchor_T_MeV))
        step = _safe_float(_high_field(row, :step_MeV))
        role = String(_high_field(row, :role, ""))
        direction = String(_high_field(row, :direction, ""))
        all(isfinite, (xi, temperature, anchor, step)) ||
            error("non-finite high-side extension plan row: $(row)")
        isapprox(xi, HIGH_SIDE_XI; atol=1e-12, rtol=0.0) ||
            error("high-side plan xi must be $(HIGH_SIDE_XI): $(row)")
        anchor == HIGH_SIDE_ANCHOR_T || error("unexpected high-side anchor: $(row)")
        step == HIGH_SIDE_STEP || error("unexpected high-side step: $(row)")
        direction == "high" || error("high-side plan direction must be high: $(row)")
        push!(selected, (xi=xi, sequence=sequence, role=role,
            T_MeV=temperature, anchor_T_MeV=anchor, step_MeV=step,
            direction=direction))
    end
    sort!(selected; by=row -> row.sequence)
    [row.sequence for row in selected] == [1, 2, 3] ||
        error("high-side plan sequence must be 1,2,3")
    [row.T_MeV for row in selected] == collect(HIGH_SIDE_TEMPERATURES) ||
        error("high-side plan temperatures are not the frozen equal-step list")
    all(row -> startswith(row.role, "high_extension_"), selected) ||
        error("high-side plan roles are not versioned extension roles")
    all(row -> row.T_MeV > row.anchor_T_MeV, selected) ||
        error("high-side extension temperatures must be above the anchor")
    all(isapprox(selected[index].T_MeV - selected[index - 1].T_MeV, HIGH_SIDE_STEP;
        atol=1e-12, rtol=0.0) for index in 2:length(selected)) ||
        error("high-side extension temperatures are not equally spaced")
    selected
end

function _high_write_csv(path, rows)
    isempty(rows) ? write(path, "\n") : CSV.write(path, rows)
end

function _high_cost_rows(cfg, plan_row, costs)
    rows = NamedTuple[]
    for item in costs
        telemetry = item.telemetry
        push!(rows, (
            xi=cfg.xi, T_MeV=plan_row.T_MeV, role=plan_row.role,
            direction=plan_row.direction, method=item.method,
            calculation_sha=cfg.calculation_sha,
            unique_solves=item.cache.unique_solves,
            point_requests=item.cache.point_requests,
            requested_point_calls=item.cache.point_requests,
            cache_hits=item.cache.cache_hits, failed_points=item.cache.failed_points,
            equilibrium_requests=telemetry.equilibrium_requests,
            fixedrho_requests=telemetry.fixedrho_requests,
            residual_calls=telemetry.nlsolve_f_calls + telemetry.postprocess_residual_calls,
            jacobian_calls=telemetry.nlsolve_g_calls,
            newton_iterations=telemetry.newton_iterations,
            runner_seconds=item.runner_seconds,
            fallback_count=telemetry.root_fallbacks,
            retry_count=telemetry.scan_retries,
        ))
    end
    rows
end

function _validate_high_cost_rows(rows)
    all(item -> item.failed_points == 0, rows) ||
        error("high-side extension has failed solver points")
    all(item -> item.point_requests == item.unique_solves + item.cache_hits, rows) ||
        error("high-side extension cache cost does not reconcile")
    true
end

function _run_high_side(cfg, plan_path)
    cfg.calculation_sha == HIGH_SIDE_CALCULATION_SHA ||
        error("high-side extension requires frozen calculation SHA $(HIGH_SIDE_CALCULATION_SHA)")
    isapprox(cfg.xi, HIGH_SIDE_XI; atol=1e-12, rtol=0.0) ||
        error("high-side extension requires xi=$(HIGH_SIDE_XI)")
    plan = _read_high_side_plan(plan_path)
    mkpath(cfg.output_dir)
    started = time_ns()
    curve_rows = NamedTuple[]
    trace_rows = NamedTuple[]
    cost_rows = NamedTuple[]
    materialization_rows = Any[]
    for (index, plan_row) in enumerate(plan)
        slice_dir = joinpath(cfg.output_dir, "slices",
            "T_$(replace(string(plan_row.T_MeV), "." => "p"))")
        rows, hybrid, costs, materialization = _run_slice(cfg, plan_row.T_MeV, slice_dir)
        append!(curve_rows, rows)
        append!(cost_rows, _high_cost_rows(cfg, plan_row, costs))
        push!(materialization_rows, Dict(
            "xi" => cfg.xi, "T_MeV" => plan_row.T_MeV,
            "role" => plan_row.role, "direction" => plan_row.direction,
            "path" => replace(relpath(materialization.path, cfg.output_dir), '\\' => '/'),
            "provenance_path" => replace(relpath(materialization.provenance_path,
                cfg.output_dir), '\\' => '/'),
            "sha256" => materialization.sha256, "rows" => materialization.rows,
            "recovered_rows" => materialization.recovered_rows,
            "aggregate_rows" => materialization.aggregate_rows,
            "aggregate_parse_error" => materialization.aggregate_parse_error,
            "source_files" => materialization.source_files,
        ))
        push!(trace_rows, (
            xi=cfg.xi, sequence=plan_row.sequence, role=plan_row.role,
            direction=plan_row.direction, T_MeV=plan_row.T_MeV,
            anchor_T_MeV=plan_row.anchor_T_MeV, step_MeV=plan_row.step_MeV,
            hybrid_status=hybrid.status, hybrid_reason=hybrid.reason,
            oracle_labels_used_for_routing=false,
            oracle_materialized_rows=materialization.rows,
            oracle_recovered_rows=materialization.recovered_rows,
            oracle_materialized_sha256=materialization.sha256,
            oracle_materialization_provenance=replace(
                relpath(materialization.provenance_path, cfg.output_dir), '\\' => '/'),
            rho_rows=length(rows), finite_and_converged=true, slice_index=index,
        ))
    end
    sort!(curve_rows; by=row -> (row.T_MeV, row.rho))
    _validate_high_cost_rows(cost_rows)
    _high_write_csv(joinpath(cfg.output_dir, "fine_pool.csv"), curve_rows)
    _high_write_csv(joinpath(cfg.output_dir, "high_side_extension_trace.csv"), trace_rows)
    _high_write_csv(joinpath(cfg.output_dir, "method_costs.csv"), cost_rows)
    elapsed = (time_ns() - started) / 1e9
    files = Dict(name => _sha(joinpath(cfg.output_dir, name)) for name in
        ("fine_pool.csv", "high_side_extension_trace.csv", "method_costs.csv"))
    plan_sha = bytes2hex(sha256(read(plan_path)))
    summary = Dict(
        "schema_version" => HIGH_SIDE_SCHEMA,
        "plan_schema_version" => HIGH_SIDE_PLAN_SCHEMA,
        "scope" => cfg.scope, "xi" => cfg.xi,
        "method" => "hybrid_and_oracle_fine_pool",
        "tag" => cfg.tag, "source_run_id" => cfg.source_run_id,
        "calculation_sha" => cfg.calculation_sha, "postprocess_sha" => cfg.postprocess_sha,
        "workflow_head_sha" => cfg.postprocess_sha,
        "rho_fine_step" => RHO_FINE_STEP, "rho_max" => RHO_MAX,
        "manual_decision_required" => true,
        "extension_direction" => "high", "extension_anchor_T_MeV" => HIGH_SIDE_ANCHOR_T,
        "extension_step_MeV" => HIGH_SIDE_STEP,
        "temperatures" => [row.T_MeV for row in plan],
        "oracle_labels_used_for_routing" => false, "solver_called" => true,
        "finite_and_converged_final" => true, "plan_sha256" => plan_sha,
        "files" => files, "materialization" => materialization_rows,
        "runner_seconds" => elapsed,
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
            "calculation_sha" => cfg.calculation_sha,
            "postprocess_sha" => cfg.postprocess_sha,
            "reference_write" => false,
            "manual_route" => true,
            "oracle_labels_used_for_routing" => false,
            "extension_contract" => "equal_high_side_step_v1",
            "anchor_is_prior_evidence_only" => true,
        ))))
        write(io, '\n')
    end
    println(JSON3.write(summary))
end

if abspath(PROGRAM_FILE) == @__FILE__
    cfg = _config(ARGS)
    plan_path = abspath(String(_arg(ARGS, "--plan-file", "")))
    _run_high_side(cfg, plan_path)
end
