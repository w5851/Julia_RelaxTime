using Dates
using JSON3

@inline _json_number(x::Real) = isfinite(x) ? Float64(x) : nothing

function _project_root()
    return normpath(joinpath(@__DIR__, "..", "..", ".."))
end

function _default_run_id(profile)
    return "phase_$(String(profile))_$(Dates.format(now(), "yyyymmdd_HHMMSS"))"
end

function resolve_phase_output_target(model_kind::Symbol; profile=:default, run_id::Union{Nothing, String}=nothing,
        policy::Symbol=:processed_first, project_root::String=_project_root())
    model_tag = lowercase(String(model_kind))
    effective_run_id = isnothing(run_id) ? _default_run_id(profile) : run_id

    processed_root = normpath(joinpath(project_root, "data", "processed", model_tag, "phase_diagram"))
    reference_root = normpath(joinpath(project_root, "data", "reference", model_tag, "phase_diagram"))
    outputs_root = normpath(joinpath(project_root, "data", "outputs", "results", model_tag, "phase"))

    run_dir = if policy == :processed_first
        joinpath(processed_root, effective_run_id)
    elseif policy == :outputs_only
        joinpath(outputs_root, effective_run_id)
    else
        joinpath(processed_root, effective_run_id)
    end

    return (
        processed_root=processed_root,
        reference_root=reference_root,
        outputs_root=outputs_root,
        run_dir=run_dir,
        run_id=effective_run_id,
    )
end

function _write_first_order_boundary(path::String, rows::Vector{NamedTuple})
    open(path, "w") do io
        println(io, "T_MeV,mu_transition_MeV,rho_hadron,rho_quark,area_residual,converged,curve_parameter,plot_order_key")
        for row in rows
            println(io,
                "$(row.T_MeV),$(row.mu_transition_MeV),$(row.rho_hadron),$(row.rho_quark),$(row.area_residual),$(row.converged),$(row.T_MeV),$(row.T_MeV)"
            )
        end
    end
end

function _write_spinodal(path::String, rows::Vector{NamedTuple})
    open(path, "w") do io
        println(io, "T_MeV,mu_spinodal_hadron_MeV,mu_spinodal_quark_MeV,rho_spinodal_hadron,rho_spinodal_quark,curve_parameter,plot_order_key")
        for row in rows
            println(io,
                "$(row.T_MeV),$(row.mu_spinodal_hadron_MeV),$(row.mu_spinodal_quark_MeV),$(row.rho_spinodal_hadron),$(row.rho_spinodal_quark),$(row.T_MeV),$(row.T_MeV)"
            )
        end
    end
end

function _write_crossover_line(path::String, rows::Vector{NamedTuple})
    open(path, "w") do io
        println(io, "mu_MeV,T_crossover_MeV,rho,method,converged,derivative,variable,curve_parameter,plot_order_key")
        for row in rows
            println(io,
                "$(row.mu_MeV),$(row.T_crossover_MeV),$(row.rho),$(row.method),$(row.converged),$(row.derivative),$(row.variable),$(row.mu_MeV),$(row.mu_MeV)"
            )
        end
    end
end

@inline function _phase_record_value(record, key::Symbol, default=nothing)
    return hasproperty(record, key) ? getproperty(record, key) : default
end

@inline function _phase_csv_value(value)
    text = value === nothing ? "" : string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text) || occursin('\r', text)
        return string('"', replace(text, '"' => "\"\""), '"')
    end
    return text
end

function _write_grid_convergence(path::String, result::PhasePipelineResult)
    records = get(result.diagnostics, "grid_convergence_records", NamedTuple[])
    open(path, "w") do io
        println(io, "axis,xi,T_MeV,level,left,right,midpoint,position_error_MeV,density_error,maxwell_area,response_rtol,converged,reason")
        for record in records
            values = (
                _phase_record_value(record, :axis, ""),
                _phase_record_value(record, :xi, result.xi),
                _phase_record_value(record, :T_MeV),
                _phase_record_value(record, :level),
                _phase_record_value(record, :left),
                _phase_record_value(record, :right),
                _phase_record_value(record, :midpoint),
                _phase_record_value(record, :position_error_MeV),
                _phase_record_value(record, :density_error),
                _phase_record_value(record, :maxwell_area),
                _phase_record_value(record, :response_rtol),
                _phase_record_value(record, :converged, false),
                _phase_record_value(record, :reason, ""),
            )
            println(io, join(_phase_csv_value.(values), ','))
        end
    end
end

function _build_conclusion(result::PhasePipelineResult)
    boundary_count = length(result.first_order_boundary)
    crossover_count = length(result.crossover_line)
    phase_structure = if boundary_count > 0
        "first_order_detected"
    elseif result.cep.result_status == :ambiguous
        "ambiguous_near_critical"
    elseif crossover_count > 0
        "crossover_only"
    else
        "no_transition_signal"
    end
    cep_result = String(result.cep.result_status)
    return Dict(
        "phase_structure" => phase_structure,
        "cep_result" => cep_result,
        "boundary_count" => boundary_count,
        "crossover_count" => crossover_count,
    )
end

function _write_phase_report(path::String, result::PhasePipelineResult)
    conclusion = _build_conclusion(result)
    open(path, "w") do io
        println(io, "# Phase Pipeline Report")
        println(io)
        println(io, "- model_kind: $(result.model_kind)")
        println(io, "- run_id: $(result.run_id)")
        println(io, "- xi: $(result.xi)")
        println(io, "- generated_at: $(now())")
        println(io)
        println(io, "## CEP")
        println(io, "- result_status: $(result.cep.result_status)")
        println(io, "- found: $(result.cep.found)")
        println(io, "- T_cep_MeV: $(isfinite(result.cep.T_cep_MeV) ? result.cep.T_cep_MeV : "null")")
        println(io, "- muq_cep_MeV: $(isfinite(result.cep.mu_cep_MeV) ? result.cep.mu_cep_MeV : "null")")
        println(io, "- muB_cep_MeV: $(isfinite(result.cep.mu_cep_MeV) ? 3.0 * result.cep.mu_cep_MeV : "null")")
        println(io, "- mu_cep_MeV: $(isfinite(result.cep.mu_cep_MeV) ? result.cep.mu_cep_MeV : "null")  # compatibility alias for muq_cep_MeV")
        println(io, "- uncertainty_T_MeV: $(isfinite(result.cep.uncertainty_T_MeV) ? result.cep.uncertainty_T_MeV : "null")")
        println(io, "- T_bracket_low_MeV: $(isfinite(result.cep.T_bracket_low_MeV) ? result.cep.T_bracket_low_MeV : "null")")
        println(io, "- T_bracket_high_MeV: $(isfinite(result.cep.T_bracket_high_MeV) ? result.cep.T_bracket_high_MeV : "null")")
        println(io, "- bracket_width_T_MeV: $(isfinite(result.cep.bracket_width_T_MeV) ? result.cep.bracket_width_T_MeV : "null")")
        println(io, "- T_last_first_order_MeV: $(isfinite(result.cep.T_last_first_order_MeV) ? result.cep.T_last_first_order_MeV : "null")")
        println(io, "- mu_last_first_order_MeV: $(isfinite(result.cep.mu_last_first_order_MeV) ? result.cep.mu_last_first_order_MeV : "null")")
        println(io, "- T_first_monotone_MeV: $(isfinite(result.cep.T_first_monotone_MeV) ? result.cep.T_first_monotone_MeV : "null")")
        println(io, "- ambiguity_width_T_MeV: $(isfinite(result.cep.ambiguity_width_T_MeV) ? result.cep.ambiguity_width_T_MeV : "null")")
        println(io, "- temperature_resolution_target_MeV: $(isfinite(result.cep.temperature_resolution_target_MeV) ? result.cep.temperature_resolution_target_MeV : "null")")
        println(io, "- eval_count: $(result.cep.eval_count)")
        println(io, "- unknown_count: $(result.cep.unknown_count)")
        println(io, "- reason: $(isnothing(result.cep.reason) ? "null" : result.cep.reason)")
        println(io, "- method: $(result.cep.method)")
        println(io)
        println(io, "## Summary")
        println(io, "- first_order_boundary_count: $(length(result.first_order_boundary))")
        println(io, "- spinodal_count: $(length(result.spinodal))")
        println(io, "- crossover_count: $(length(result.crossover_line))")
        for (k, v) in sort(collect(result.diagnostics); by=first)
            println(io, "- diag.$k: $v")
        end
        println(io)
        println(io, "## Conclusion")
        println(io, "- phase_structure: $(conclusion["phase_structure"])")
        println(io, "- cep_result: $(conclusion["cep_result"])")
    end
end

function _build_summary(result::PhasePipelineResult)
    conclusion = _build_conclusion(result)
    scan_total = get(result.diagnostics, "scan_total", 0)
    scan_success = get(result.diagnostics, "scan_success", 0)
    scan_failure = get(result.diagnostics, "scan_failure", 0)
    scan_success_rate = scan_total > 0 ? Float64(scan_success) / Float64(scan_total) : nothing
    scan_failure_rate = scan_total > 0 ? Float64(scan_failure) / Float64(scan_total) : nothing

    cep_eval_count = result.cep.eval_count
    cep_unknown_count = result.cep.unknown_count
    cep_unknown_rate = cep_eval_count > 0 ? Float64(cep_unknown_count) / Float64(cep_eval_count) : nothing

    stats_compare = Dict(
        "scan" => Dict(
            "total" => scan_total,
            "success" => scan_success,
            "failure" => scan_failure,
            "success_rate" => scan_success_rate,
            "failure_rate" => scan_failure_rate,
        ),
        "phase" => Dict(
            "curve_count" => get(result.diagnostics, "curve_count", 0),
            "boundary_count" => get(result.diagnostics, "boundary_count", 0),
            "spinodal_count" => get(result.diagnostics, "spinodal_count", 0),
            "crossover_count" => get(result.diagnostics, "crossover_count", 0),
        ),
        "cep" => Dict(
            "method" => String(result.cep.method),
            "result_status" => String(result.cep.result_status),
            "found" => result.cep.found,
            "eval_count" => cep_eval_count,
            "unknown_count" => cep_unknown_count,
            "unknown_rate" => cep_unknown_rate,
            "uncertainty_T_MeV" => _json_number(result.cep.uncertainty_T_MeV),
            "T_bracket_low_MeV" => _json_number(result.cep.T_bracket_low_MeV),
            "T_bracket_high_MeV" => _json_number(result.cep.T_bracket_high_MeV),
            "bracket_width_T_MeV" => _json_number(result.cep.bracket_width_T_MeV),
            "T_last_first_order_MeV" => _json_number(result.cep.T_last_first_order_MeV),
            "mu_last_first_order_MeV" => _json_number(result.cep.mu_last_first_order_MeV),
            "T_first_monotone_MeV" => _json_number(result.cep.T_first_monotone_MeV),
            "ambiguity_width_T_MeV" => _json_number(result.cep.ambiguity_width_T_MeV),
            "temperature_resolution_target_MeV" => _json_number(result.cep.temperature_resolution_target_MeV),
            "reason" => result.cep.reason,
        ),
    )

    return Dict(
        "model_kind" => String(result.model_kind),
        "model_variant" => result.model_variant,
        "schema_version" => "phase-v2",
        "config_hash" => get(result.config_snapshot, "config_hash", ""),
        "run_id" => result.run_id,
        "xi" => result.xi,
        "generated_at" => string(now()),
        "cep" => Dict(
            "result_status" => String(result.cep.result_status),
            "found" => result.cep.found,
            "T_cep_MeV" => _json_number(result.cep.T_cep_MeV),
            "muq_cep_MeV" => _json_number(result.cep.mu_cep_MeV),
            "muB_cep_MeV" => _json_number(3.0 * result.cep.mu_cep_MeV),
            "mu_cep_MeV" => _json_number(result.cep.mu_cep_MeV),
            "uncertainty_T_MeV" => _json_number(result.cep.uncertainty_T_MeV),
            "T_bracket_low_MeV" => _json_number(result.cep.T_bracket_low_MeV),
            "T_bracket_high_MeV" => _json_number(result.cep.T_bracket_high_MeV),
            "bracket_width_T_MeV" => _json_number(result.cep.bracket_width_T_MeV),
            "T_last_first_order_MeV" => _json_number(result.cep.T_last_first_order_MeV),
            "mu_last_first_order_MeV" => _json_number(result.cep.mu_last_first_order_MeV),
            "T_first_monotone_MeV" => _json_number(result.cep.T_first_monotone_MeV),
            "ambiguity_width_T_MeV" => _json_number(result.cep.ambiguity_width_T_MeV),
            "temperature_resolution_target_MeV" => _json_number(result.cep.temperature_resolution_target_MeV),
            "eval_count" => result.cep.eval_count,
            "unknown_count" => result.cep.unknown_count,
            "reason" => result.cep.reason,
            "method" => String(result.cep.method),
        ),
        "stats" => result.diagnostics,
        "stats_compare" => stats_compare,
        "conclusion" => conclusion,
    )
end

function build_phase_artifacts(result::PhasePipelineResult; output_dir::String, format_options...)
    mkpath(output_dir)

    boundary_path = joinpath(output_dir, "first_order_boundary.csv")
    spinodal_path = joinpath(output_dir, "spinodal.csv")
    crossover_path = joinpath(output_dir, "crossover_line.csv")
    grid_convergence_path = joinpath(output_dir, "phase_grid_convergence.csv")
    report_path = joinpath(output_dir, "phase_report.md")
    summary_path = joinpath(output_dir, "phase_summary.json")

    _write_first_order_boundary(boundary_path, result.first_order_boundary)
    _write_spinodal(spinodal_path, result.spinodal)
    _write_crossover_line(crossover_path, result.crossover_line)
    _write_grid_convergence(grid_convergence_path, result)
    _write_phase_report(report_path, result)

    summary = _build_summary(result)
    open(summary_path, "w") do io
        write(io, JSON3.write(summary))
    end

    return Dict(
        "first_order_boundary" => boundary_path,
        "spinodal" => spinodal_path,
        "crossover_line" => crossover_path,
        "phase_grid_convergence" => grid_convergence_path,
        "phase_report" => report_path,
        "phase_summary" => summary_path,
    )
end

function _copy_tree(src::String, dst::String)
    mkpath(dst)
    for (root, dirs, files) in walkdir(src)
        rel = relpath(root, src)
        target_root = rel == "." ? dst : joinpath(dst, rel)
        mkpath(target_root)
        for dir in dirs
            mkpath(joinpath(target_root, dir))
        end
        for file in files
            cp(joinpath(root, file), joinpath(target_root, file); force=true)
        end
    end
end

function _load_summary(path::String)
    isfile(path) || return nothing
    try
        return JSON3.read(read(path, String))
    catch
        return nothing
    end
end

function promote_phase_artifacts(processed_run_dir::String; reference_root::Union{Nothing, String}=nothing,
        gate_options::NamedTuple=(;), write_reference::Bool=true)
    failed = String[]
    required = ["phase_summary.json", "first_order_boundary.csv", "spinodal.csv", "crossover_line.csv", "phase_grid_convergence.csv", "phase_report.md"]

    for file in required
        isfile(joinpath(processed_run_dir, file)) || push!(failed, "missing_file:$file")
    end

    summary = _load_summary(joinpath(processed_run_dir, "phase_summary.json"))
    if summary === nothing
        push!(failed, "invalid_summary_json")
    else
        for key in ("model_kind", "schema_version", "config_hash", "run_id")
            haskey(summary, key) || push!(failed, "missing_summary_field:$key")
        end

        expected_schema = get(gate_options, :expected_schema_version, nothing)
        if expected_schema !== nothing
            actual_schema = get(summary, "schema_version", nothing)
            if actual_schema === nothing
                push!(failed, "schema_mismatch:expected=$(expected_schema),actual=null")
            elseif String(actual_schema) != String(expected_schema)
                push!(failed, "schema_mismatch:expected=$(expected_schema),actual=$(actual_schema)")
            end
        end

        expected_model = get(gate_options, :expected_model_kind, nothing)
        if expected_model !== nothing
            actual_model = get(summary, "model_kind", nothing)
            if actual_model === nothing
                push!(failed, "model_kind_mismatch:expected=$(expected_model),actual=null")
            elseif uppercase(String(actual_model)) != uppercase(String(expected_model))
                push!(failed, "model_kind_mismatch:expected=$(expected_model),actual=$(actual_model)")
            end
        end

        expected_hash = get(gate_options, :expected_config_hash, nothing)
        if expected_hash !== nothing
            actual_hash = get(summary, "config_hash", nothing)
            if actual_hash === nothing
                push!(failed, "config_hash_mismatch")
            elseif String(actual_hash) != String(expected_hash)
                push!(failed, "config_hash_mismatch")
            end
        end
    end

    if !isempty(failed)
        return PromotionResult(passed=false, failed_checks=failed)
    end

    schema_version = String(get(summary, "schema_version", "phase-v1"))
    config_hash = String(get(summary, "config_hash", ""))
    run_id = String(get(summary, "run_id", basename(processed_run_dir)))
    profile = get(gate_options, :profile, "default")
    hash_short = isempty(config_hash) ? "nohash" : first(config_hash, min(12, lastindex(config_hash)))
    baseline_id = "baseline_$(profile)_$(schema_version)_$(hash_short)"

    ref_root = isnothing(reference_root) ? dirname(processed_run_dir) : reference_root
    dst = joinpath(ref_root, baseline_id)

    if write_reference
        _copy_tree(processed_run_dir, dst)
    end

    return PromotionResult(
        passed=true,
        baseline_id=baseline_id,
        failed_checks=String[],
        reference_dir=dst,
    )
end
