using JSON3

@inline _pm_csv_value(x) = isnothing(x) ? "" : string(x)

function _pm_branch_scan_fieldnames()
    return (
        :T_MeV,
        :mu_MeV,
        :branch,
        :branch_status,
        :seed_source,
        :continuity_ok,
        :converged,
        :residual_norm,
        :rho_norm,
        :pressure,
        :omega,
        :endpoint_cause,
    )
end

function _pm_branch_scan_record(row)
    return NamedTuple{_pm_branch_scan_fieldnames()}((
        Float64(row.T_MeV),
        Float64(row.mu_MeV),
        Symbol(row.branch),
        Symbol(row.branch_status),
        Symbol(row.seed_source),
        Bool(row.continuity_ok),
        Bool(row.converged),
        Float64(row.residual_norm),
        Float64(row.rho_norm),
        Float64(row.pressure),
        Float64(row.omega),
        isnothing(row.endpoint_cause) ? nothing : Symbol(row.endpoint_cause),
    ))
end

function _pm_phase_summary_record(row)
    return Dict(
        "T_MeV" => Float64(row.T_MeV),
        "mu_transition_pm_MeV" => row.mu_transition_pm_MeV,
        "hadron_endpoint_mu_MeV" => row.hadron_endpoint_mu_MeV,
        "quark_endpoint_mu_MeV" => row.quark_endpoint_mu_MeV,
        "bistable_window_width_MeV" => Float64(row.bistable_window_width_MeV),
        "comparison_status" => String(row.comparison_status),
        "comparison_mu_tol_MeV" => Float64(row.comparison_mu_tol_MeV),
        "residual_accept_tol" => Float64(row.residual_accept_tol),
        "continuity_x_tol" => Float64(row.continuity_x_tol),
        "continuity_rho_tol" => Float64(row.continuity_rho_tol),
    )
end

function _pm_compare_record(row)
    return (
        T_MeV=Float64(row.T_MeV),
        mu_transition_pm_MeV=row.mu_transition_pm_MeV,
        mu_transition_maxwell_MeV=get(row, :mu_transition_maxwell_MeV, nothing),
        delta_mu_pm_minus_maxwell_MeV=get(row, :delta_mu_pm_minus_maxwell_MeV, nothing),
        comparison_status=Symbol(row.comparison_status),
        branch_disappears_first=Symbol(row.branch_disappears_first),
        hadron_endpoint_mu_MeV=row.hadron_endpoint_mu_MeV,
        quark_endpoint_mu_MeV=row.quark_endpoint_mu_MeV,
        bistable_window_width_MeV=Float64(row.bistable_window_width_MeV),
    )
end

function _pm_write_artifacts(output_dir::String, branch_rows, temperature_summaries, comparison_rows)
    mkpath(output_dir)

    branch_path = joinpath(output_dir, "pm_branch_scan.csv")
    branch_fields = _pm_branch_scan_fieldnames()
    open(branch_path, "w") do io
        println(io, join(string.(branch_fields), ","))
        for row in branch_rows
            record = _pm_branch_scan_record(row)
            println(io, join((_pm_csv_value(getproperty(record, field)) for field in branch_fields), ","))
        end
    end

    summary_path = joinpath(output_dir, "pm_phase_summary.json")
    payload = [_pm_phase_summary_record(row) for row in temperature_summaries]
    open(summary_path, "w") do io
        write(io, JSON3.write(payload))
    end

    compare_path = joinpath(output_dir, "pm_vs_maxwell.csv")
    compare_fields = (
        :T_MeV,
        :mu_transition_pm_MeV,
        :mu_transition_maxwell_MeV,
        :delta_mu_pm_minus_maxwell_MeV,
        :comparison_status,
        :branch_disappears_first,
        :hadron_endpoint_mu_MeV,
        :quark_endpoint_mu_MeV,
        :bistable_window_width_MeV,
    )
    open(compare_path, "w") do io
        println(io, join(string.(compare_fields), ","))
        for row in comparison_rows
            record = _pm_compare_record(row)
            println(io, join((_pm_csv_value(getproperty(record, field)) for field in compare_fields), ","))
        end
    end

    return (
        pm_branch_scan=branch_path,
        pm_phase_summary=summary_path,
        pm_vs_maxwell=compare_path,
    )
end
