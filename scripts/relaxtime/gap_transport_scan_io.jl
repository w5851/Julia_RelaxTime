module GapTransportScanIO

using Dates

export write_channel_diagnostics_header_if_needed
export write_failed_points_header_if_needed
export write_failed_point_row!
export ensure_parent_dir
export current_git_commit
export ensure_output_header_compatible
export write_header_if_needed

function write_channel_diagnostics_header_if_needed(io)
    header = join([
        "T_MeV", "muq_MeV", "muB_MeV", "xi",
        "species", "channel", "density_key", "multiplicity",
        "density", "rate", "contribution", "total", "tau_inv_species",
        "equilibrium_backend", "phase_curr", "phase_structure",
    ], ',')
    println(io, header)
end

@inline function _csv_quote(text::AbstractString)
    return "\"" * replace(text, "\"" => "\"\"") * "\""
end

function write_failed_points_header_if_needed(io)
    header = join([
        "T_MeV", "muB_MeV", "xi",
        "seed_source", "phase_prev", "phase_curr_hint",
        "error_type", "error_message", "timestamp",
    ], ',')
    println(io, header)
end

function write_failed_point_row!(io, T_mev::Float64, muB_mev::Float64, xi::Float64, diag, err)
    seed_source = hasproperty(diag, :seed_source) ? string(getproperty(diag, :seed_source)) : "unknown"
    phase_prev = hasproperty(diag, :phase_prev) ? string(getproperty(diag, :phase_prev)) : "unknown"
    phase_curr_hint = hasproperty(diag, :phase_curr) ? string(getproperty(diag, :phase_curr)) : "unknown"
    error_type = string(typeof(err))
    error_message = sprint(showerror, err)
    timestamp = Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS")
    row = join([
        string(T_mev),
        string(muB_mev),
        string(xi),
        _csv_quote(seed_source),
        _csv_quote(phase_prev),
        _csv_quote(phase_curr_hint),
        _csv_quote(error_type),
        _csv_quote(error_message),
        _csv_quote(timestamp),
    ], ',')
    println(io, row)
    flush(io)
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function current_git_commit()
    try
        return readchomp(`git -C $(Main.PROJECT_ROOT) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function ensure_output_header_compatible(path::AbstractString)
    isfile(path) || return
    header = nothing
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, "#")) && continue
            header = s
            break
        end
    end
    header === nothing && return
    required = (
        "mode",
        "phase_reference_kind",
        "scan_group",
        "group_label",
        "plot_panel",
        "plot_panel_label",
        "plot_series",
        "plot_series_label",
        "T_phase_base_MeV",
        "alpha_T",
        "omega_fm4inv",
        "P_fm4inv",
        "epsilon_fm4inv",
        "s_fm3inv",
        "eps_minus_3P_over_T4",
        "eta_over_s",
        "zeta_over_s",
        "sigma_over_T",
        "sigma_over_T_over_eta_over_s",
        "zeta_over_s_over_eta_over_s",
        "quality_flag",
        "quality_reason",
        "quality_metric",
        "run_id",
        "equilibrium_backend",
        "seed_source",
        "phase_prev",
        "phase_curr",
        "phase_structure",
        "phase_boundary_xi_used",
    )
    for c in required
        occursin(c, header) || error("existing output CSV header is incompatible with current script (missing column: $c). Please rerun with --overwrite or choose a new --output path.")
    end
end

function write_header_if_needed(io)
    header = join([
        "T_MeV", "muq_MeV", "muB_MeV", "xi",
        "mode", "phase_reference_kind", "scan_group", "group_label", "plot_panel", "plot_panel_label", "plot_series", "plot_series_label", "T_phase_base_MeV", "alpha_T",
        "T_fm", "muq_fm",
        "converged", "iterations", "residual_norm", "equilibrium_backend", "seed_source", "phase_prev", "phase_curr", "phase_structure", "phase_boundary_xi_used",
        "Phi", "Phibar",
        "m_u", "m_d", "m_s",
        "rho_baryon", "rho_norm",
        "omega_fm4inv", "P_fm4inv", "epsilon_fm4inv", "s_fm3inv",
        "omega_MeV_fm3", "P_MeV_fm3", "epsilon_MeV_fm3",
        "eps_minus_3P_over_T4",
        "n_u", "n_d", "n_s", "n_ubar", "n_dbar", "n_sbar",
        "tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
        "tauinv_u", "tauinv_d", "tauinv_s", "tauinv_ubar", "tauinv_dbar", "tauinv_sbar",
        "eta", "sigma", "zeta", "eta_over_s", "zeta_over_s",
        "sigma_over_T", "sigma_over_T_over_eta_over_s", "zeta_over_s_over_eta_over_s",
        "quality_flag", "quality_reason", "quality_metric", "run_id",
    ], ',')
    println(io, header)
end

end # module GapTransportScanIO
