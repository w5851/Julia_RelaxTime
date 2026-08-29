module PhaseGuidedTransportScanAssets

using JSON3

function _phase_reference_summary_dict(source)
    summary = Main.PhaseReferenceAdapter.source_summary(source)
    return Dict{String,Any}(String(key) => value for (key, value) in pairs(summary))
end

@inline function _csv_quote(text::AbstractString)
    return "\"" * replace(String(text), "\"" => "\"\"") * "\""
end

function write_plan_csv(path::String, plan)
    open(path, "w") do io
        println(io, "T_MeV,muB_MeV,xi,mode,phase_reference_kind,scan_group,group_label,plot_panel,plot_panel_label,plot_series,plot_series_label,T_phase_base_MeV,alpha_T,phase_anchor_method,coexistence_side,coexistence_certified,coexistence_delta_xi,coexistence_T_lower_MeV,coexistence_T_upper_MeV,anchor_p_num,anchor_t_num,anchor_convergence_delta_MeV,anchor_convergence_certified")
        for point in plan.points
            println(io, join([
                string(point.T_MeV),
                string(point.muB_MeV),
                string(point.xi),
                string(point.mode),
                string(point.phase_reference_kind),
                point.scan_group,
                _csv_quote(point.group_label),
                point.plot_panel,
                _csv_quote(point.plot_panel_label),
                point.plot_series,
                _csv_quote(point.plot_series_label),
                string(point.T_phase_base_MeV),
                string(point.alpha_T),
                string(point.phase_anchor_method),
                string(point.coexistence_side),
                string(point.coexistence_certified),
                string(point.coexistence_delta_xi),
                string(point.coexistence_T_lower_MeV),
                string(point.coexistence_T_upper_MeV),
                string(point.anchor_p_num),
                string(point.anchor_t_num),
                string(point.anchor_convergence_delta_MeV),
                string(point.anchor_convergence_certified),
            ], ','))
        end
    end
end

function write_readme(path::String, opts, plan; result_csv_name::String="phase_guided_transport_scan.csv", figure_dir::Union{Nothing,String}=nothing, phase_reference=nothing)
    mode_dir = opts.mode == :mode_a_fixed_muB_phase_scaled ? "mode_a_fixed_muB_phase_scaled" : "mode_b_fixed_T_sparse_muB"
    mode_summary = opts.mode == :mode_a_fixed_muB_phase_scaled ?
        "固定 muB，沿相变参考温度做 T/T_phase 倍率带，并对每个倍率带连续扫描 xi。" :
        "固定温度、离散 muB、连续 xi 的相图邻域稀疏扫描。"
    open(path, "w") do io
        println(io, "# $(basename(dirname(path)))")
        println(io)
        println(io, "phase-guided transport canonical case for `$(mode_dir)`.")
        println(io)
        println(io, "## Scope")
        println(io, "- mode: `$(opts.mode)`")
        println(io, "- case name: `$(opts.case_name)`")
        println(io, "- summary: $(mode_summary)")
        println(io, "- xi list: `$(join(opts.xi_values, ", "))`")
        println(io, "- muB list (MeV): `$(join(opts.muB_values, ", "))`")
        if opts.mode == :mode_a_fixed_muB_phase_scaled
            println(io, "- alpha_T list: `$(join(opts.alpha_T_values, ", "))`")
        else
            println(io, "- fixed T list (MeV): `$(join(opts.T_values, ", "))`")
        end
        println(io, "- compute bulk viscosity (`zeta`): `$(opts.compute_bulk)`")
        println(io, "- propagator xi policy: `$(opts.propagator_xi_policy)`")
        println(io, "- sigma cache policy: `$(opts.sigma_cache_policy)`")
        println(io, "- thermodynamic nodes: `p_num=$(opts.p_num), t_num=$(opts.t_num)`")
        println(io, "- phase anchor policy: `$(opts.phase_anchor_policy)`")
        phase_root = phase_reference === nothing ? opts.phase_reference_root :
            get(_phase_reference_summary_dict(phase_reference), "source_root", opts.phase_reference_root)
        println(io, "- phase reference root: `$(phase_root)`")
        println(io, "- phase reference layer: `$(opts.phase_reference_layer)`")
        println(io, "- phase reference mode: `$(opts.phase_reference_mode)`")
        if phase_reference !== nothing
            phase_summary = _phase_reference_summary_dict(phase_reference)
            println(io, "- phase reference source: `$(Main.PhaseReferenceAdapter.source_kind(phase_reference))`")
            println(io, "- phase reference primary layer: `$(get(phase_summary, "primary_layer", String(Main.PhaseReferenceAdapter.source_layer(phase_reference))))`")
            println(io, "- phase reference runtime view: `$(get(phase_summary, "runtime_view", ""))`")
        end
        if opts.tau_p_nodes !== nothing || opts.tau_angle_nodes !== nothing || opts.tau_phi_nodes !== nothing ||
           opts.tau_n_sigma_points !== nothing || opts.sigma_grid_n !== nothing
            println(io, "- tau/sigma overrides:")
            println(io, "  - tau_p_nodes: `$(opts.tau_p_nodes)`")
            println(io, "  - tau_angle_nodes: `$(opts.tau_angle_nodes)`")
            println(io, "  - tau_phi_nodes: `$(opts.tau_phi_nodes)`")
            println(io, "  - tau_n_sigma_points: `$(opts.tau_n_sigma_points)`")
            println(io, "  - sigma_grid_n: `$(opts.sigma_grid_n)`")
        end
        println(io, "- channel diagnostics: `$(opts.channel_diagnostics)`")
        println(io, "- total planned points: `$(plan.total)`")
        println(io)
        println(io, "## Key Files")
        println(io, "- `sampling_plan.csv`: phase-guided sampling plan")
        println(io, "- `$(result_csv_name)`: transport scan output CSV")
        println(io, "- `effective_config.json`: effective run config snapshot")
        println(io, "- `run_manifest.json`: provenance metadata")
        if isnothing(figure_dir)
            println(io, "- canonical figures: generated separately under `data/outputs/figures/relaxtime/transport/phase_guided/<mode>/<case_name>/`")
        else
            println(io, "- canonical figures: `$(replace(String(figure_dir), '\\' => '/'))`")
        end
        println(io)
        println(io, "## Interpretation Boundary")
        println(io, "- This directory is a user-facing result asset, not an external validation truth set.")
        println(io, "- Numerical drift should be governed separately by regression coverage.")
        println(io, "- For a directly anchored first-order alpha_T=1 slice, xi=0 is intentionally absent; certified negative/positive near-zero points represent the quark/hadron side limits.")
    end
end

function build_effective_config(opts, result_csv::String, plan_csv::String; figure_dir::Union{Nothing,String}=nothing, phase_reference=nothing)
    resolved_root = phase_reference === nothing ? opts.phase_reference_root :
        get(Main.PhaseReferenceAdapter.source_summary(phase_reference), :source_root, opts.phase_reference_root)
    cfg = Dict(
        "mode" => String(opts.mode),
        "outdir" => opts.outdir,
        "case_name" => opts.case_name,
        "result_csv" => result_csv,
        "plan_csv" => plan_csv,
        "figure_dir" => figure_dir,
        "xi_values" => opts.xi_values,
        "muB_values" => opts.muB_values,
        "alpha_T_values" => opts.alpha_T_values,
        "T_values" => opts.T_values,
        "propagator_xi_policy" => String(opts.propagator_xi_policy),
        "sigma_cache_policy" => String(opts.sigma_cache_policy),
        "tau_p_nodes" => opts.tau_p_nodes,
        "tau_angle_nodes" => opts.tau_angle_nodes,
        "tau_phi_nodes" => opts.tau_phi_nodes,
        "tau_n_sigma_points" => opts.tau_n_sigma_points,
        "sigma_grid_n" => opts.sigma_grid_n,
        "p_num" => opts.p_num,
        "t_num" => opts.t_num,
        "phase_anchor_policy" => String(opts.phase_anchor_policy),
        "phase_reference_root" => resolved_root,
        "phase_reference_layer" => String(opts.phase_reference_layer),
        "phase_reference_mode" => String(opts.phase_reference_mode),
        "phase_reference_source" => phase_reference === nothing ? "none" : String(Main.PhaseReferenceAdapter.source_kind(phase_reference)),
        "phase_reference_primary_layer" => phase_reference === nothing ? "none" : get(_phase_reference_summary_dict(phase_reference), "primary_layer", String(Main.PhaseReferenceAdapter.source_layer(phase_reference))),
        "channel_diagnostics" => opts.channel_diagnostics,
        "compute_bulk" => opts.compute_bulk,
        "dry_run" => opts.dry_run,
        "overwrite" => opts.overwrite,
        "resume" => opts.resume,
    )
    phase_reference !== nothing && (cfg["phase_reference_summary"] = _phase_reference_summary_dict(phase_reference))
    return cfg
end

function write_effective_config(path::String, cfg)
    open(path, "w") do io
        JSON3.pretty(io, cfg)
    end
end

export build_effective_config
export write_effective_config
export write_plan_csv
export write_readme

end # module PhaseGuidedTransportScanAssets
