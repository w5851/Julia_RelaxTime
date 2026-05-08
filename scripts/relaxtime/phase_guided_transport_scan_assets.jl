module PhaseGuidedTransportScanAssets

using JSON3

@inline function _csv_quote(text::AbstractString)
    return "\"" * replace(String(text), "\"" => "\"\"") * "\""
end

function write_plan_csv(path::String, plan)
    open(path, "w") do io
        println(io, "T_MeV,muB_MeV,xi,mode,phase_reference_kind,scan_group,group_label,T_phase_base_MeV,alpha_T")
        for point in plan.points
            println(io, join([
                string(point.T_MeV),
                string(point.muB_MeV),
                string(point.xi),
                string(point.mode),
                string(point.phase_reference_kind),
                point.scan_group,
                _csv_quote(point.group_label),
                string(point.T_phase_base_MeV),
                string(point.alpha_T),
            ], ','))
        end
    end
end

function write_readme(path::String, opts, plan; result_csv_name::String="phase_guided_transport_scan.csv")
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
        println(io, "- total planned points: `$(plan.total)`")
        println(io)
        println(io, "## Key Files")
        println(io, "- `sampling_plan.csv`: phase-guided sampling plan")
        println(io, "- `$(result_csv_name)`: transport scan output CSV")
        println(io, "- `effective_config.json`: effective run config snapshot")
        println(io, "- `run_manifest.json`: provenance metadata")
        println(io, "- `figures/`: reserved for plot-review / canonical figures")
        println(io)
        println(io, "## Interpretation Boundary")
        println(io, "- This directory is a user-facing result asset, not an external validation truth set.")
        println(io, "- Numerical drift should be governed separately by regression coverage.")
    end
end

function build_effective_config(opts, result_csv::String, plan_csv::String)
    return Dict(
        "mode" => String(opts.mode),
        "outdir" => opts.outdir,
        "case_name" => opts.case_name,
        "result_csv" => result_csv,
        "plan_csv" => plan_csv,
        "xi_values" => opts.xi_values,
        "muB_values" => opts.muB_values,
        "alpha_T_values" => opts.alpha_T_values,
        "T_values" => opts.T_values,
        "compute_bulk" => opts.compute_bulk,
        "dry_run" => opts.dry_run,
        "overwrite" => opts.overwrite,
        "resume" => opts.resume,
    )
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
