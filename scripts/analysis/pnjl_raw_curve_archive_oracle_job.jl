#!/usr/bin/env julia

"""Produce one xi shard of full raw PNJL rho-mu curves.

This is intentionally a raw-only producer.  It calls the request-scoped
fixed-rho equilibrium session and materializes the same 28 columns emitted by
TrhoScan.  It does not call the phase pipeline, Maxwell, CEP, boundary,
reference, or transport code.
"""

using Pkg
using CSV
using JSON3
using SHA
using Printf

function _arg(args, name, default=nothing)
    index = findfirst(==(name), args)
    index === nothing && return default
    index == length(args) && throw(ArgumentError("missing value for $(name)"))
    return args[index + 1]
end

function _required_arg(args, name)
    value = _arg(args, name, nothing)
    value === nothing && throw(ArgumentError("missing required argument $(name)"))
    return String(value)
end

const PROJECT_ROOT = abspath(_arg(ARGS, "--project-root", pwd()))
Pkg.activate(PROJECT_ROOT)
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const MODELS = Main.Models

function _parse_sha(value, name)
    text = String(value)
    occursin(r"^[0-9a-fA-F]{40}$", text) ||
        throw(ArgumentError("$(name) must be a 40-character git SHA"))
    return lowercase(text)
end

function _parse_temperatures(args)
    path = _arg(args, "--temperatures-file", nothing)
    raw = if path === nothing
        _required_arg(args, "--temperatures-json")
    else
        read(String(path), String)
    end
    values = JSON3.read(raw)
    temperatures = sort!(unique!(Float64.(collect(values))))
    isempty(temperatures) && throw(ArgumentError("temperature list is empty"))
    all(isfinite, temperatures) || throw(ArgumentError("temperature list contains a non-finite value"))
    return temperatures
end

function _token(value::Real)
    text = replace(@sprintf("%.8f", Float64(value)), "." => "p", "-" => "m")
    text = replace(text, r"0+$" => "")
    endswith(text, "p") && (text *= "0")
    return text
end

function _finite_or_nan(value)
    numeric = Float64(value)
    return isfinite(numeric) ? numeric : NaN
end

const RAW_COLUMNS = (
    :T_MeV,
    :rho,
    :xi,
    :mu_u_MeV,
    :mu_d_MeV,
    :mu_s_MeV,
    :mu_avg_MeV,
    :mu_B_MeV,
    :mu_Q_MeV,
    :mu_S_MeV,
    :pressure_fm4,
    :entropy_fm3,
    :energy_fm4,
    :rho_u_fm3,
    :rho_d_fm3,
    :rho_s_fm3,
    :phi_u,
    :phi_d,
    :phi_s,
    :Phi1,
    :Phi2,
    :M_u_MeV,
    :M_d_MeV,
    :M_s_MeV,
    :iterations,
    :residual_norm,
    :converged,
    :message,
)
const RHO_STEP = 0.003125
const RHO_COUNT = 1281

function _failed_row(T, xi, rho, message)
    return (
        T_MeV=Float64(T), rho=Float64(rho), xi=Float64(xi),
        mu_u_MeV=NaN, mu_d_MeV=NaN, mu_s_MeV=NaN, mu_avg_MeV=NaN,
        mu_B_MeV=NaN, mu_Q_MeV=NaN, mu_S_MeV=NaN,
        pressure_fm4=NaN, entropy_fm3=NaN, energy_fm4=NaN,
        rho_u_fm3=NaN, rho_d_fm3=NaN, rho_s_fm3=NaN,
        phi_u=NaN, phi_d=NaN, phi_s=NaN, Phi1=NaN, Phi2=NaN,
        M_u_MeV=NaN, M_d_MeV=NaN, M_s_MeV=NaN,
        iterations=-1, residual_norm=NaN, converged=false,
        message=String(message),
    )
end

function _materialize_row(model, session_row, T::Float64, xi::Float64, hbar)
    if !session_row.converged || !session_row.finite || session_row.result === nothing
        return _failed_row(T, xi, session_row.rho, session_row.message)
    end
    result = session_row.result
    mu_mev = Float64.(result.mu_vec .* hbar)
    rho_vec = MODELS.model_rho(
        model,
        result.x_state,
        result.mu_vec,
        T / hbar;
        p_num=24,
        t_num=8,
        xi=xi,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8,
        thermo_quadrature_atol=1e-10,
        thermo_quadrature_maxevals=10^7,
    )
    phi = result.x_state[1:3]
    masses_mev = result.masses .* hbar
    return (
        T_MeV=Float64(T), rho=Float64(session_row.rho), xi=Float64(xi),
        mu_u_MeV=mu_mev[1], mu_d_MeV=mu_mev[2], mu_s_MeV=mu_mev[3],
        mu_avg_MeV=sum(mu_mev) / 3,
        mu_B_MeV=mu_mev[1] + 2 * mu_mev[2],
        mu_Q_MeV=mu_mev[1] - mu_mev[2],
        mu_S_MeV=mu_mev[2] - mu_mev[3],
        pressure_fm4=Float64(result.pressure),
        entropy_fm3=Float64(result.entropy),
        energy_fm4=Float64(result.energy),
        rho_u_fm3=Float64(rho_vec[1]),
        rho_d_fm3=Float64(rho_vec[2]),
        rho_s_fm3=Float64(rho_vec[3]),
        phi_u=Float64(phi[1]), phi_d=Float64(phi[2]), phi_s=Float64(phi[3]),
        Phi1=Float64(result.x_state[4]), Phi2=Float64(result.x_state[5]),
        M_u_MeV=Float64(masses_mev[1]), M_d_MeV=Float64(masses_mev[2]),
        M_s_MeV=Float64(masses_mev[3]),
        iterations=Int(result.iterations),
        residual_norm=Float64(result.residual_norm),
        converged=Bool(session_row.converged),
        message=String(session_row.message),
    )
end

function _write_curve(path, rows)
    mkpath(dirname(path))
    CSV.write(path, rows; writeheader=true, header=collect(RAW_COLUMNS))
    return path
end

function _run(args)
    xi = parse(Float64, _required_arg(args, "--xi"))
    isfinite(xi) || throw(ArgumentError("xi must be finite"))
    temperatures = _parse_temperatures(args)
    output_dir = abspath(_required_arg(args, "--output-dir"))
    calculation_sha = _parse_sha(_required_arg(args, "--calculation-sha"), "calculation_sha")
    source_postprocess_sha = _parse_sha(_required_arg(args, "--source-postprocess-sha"), "source_postprocess_sha")
    source_workflow_sha = _parse_sha(_required_arg(args, "--source-workflow-sha"), "source_workflow_sha")
    source_run_id = String(_required_arg(args, "--source-run-id"))
    source_grid_run_id = String(_arg(args, "--source-grid-run-id", "31862752226"))
    source_grid_artifact_name = String(_arg(args, "--source-grid-artifact-name", ""))
    source_grid_manifest_path = String(_arg(args, "--source-grid-manifest-path", ""))
    source_grid_manifest_sha256 = String(_arg(args, "--source-grid-manifest-sha256", ""))

    hbar = Main.Constants_PNJL.ħc_MeV_fm
    rho_grid = Float64.(collect(range(0.0, 4.0; length=RHO_COUNT)))
    isapprox(first(diff(rho_grid)), RHO_STEP; atol=1e-12, rtol=0.0) ||
        error("raw rho grid step is not the frozen 0.003125 contract")
    telemetry = MODELS.SolverWorkTelemetry()
    session = MODELS.TrhoScan.new_rho_point_session(
        model_kind=:PNJL,
        reverse_rho=true,
        seed_policy=:candidates,
        solver_backend=:models,
        p_num=24,
        t_num=8,
        iterations=80,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-8,
        thermo_quadrature_atol=1e-10,
        thermo_quadrature_maxevals=10^7,
        telemetry=telemetry,
    )
    model = MODELS.create_model(:PNJL)
    curve_files = String[]
    failed_points = 0
    for T in temperatures
        for rho in reverse(rho_grid)
            MODELS.TrhoScan.rho_point!(session, T, xi, rho)
        end
        cached_rhos = MODELS.TrhoScan.rho_session_slice_rhos(session, T, xi)
        rows = NamedTuple[]
        for rho in reverse(sort!(copy(cached_rhos)))
            row = session.cache[(T, xi, rho)]
            materialized = _materialize_row(model, row, T, xi, hbar)
            materialized.converged || (failed_points += 1)
            push!(rows, materialized)
        end
        length(rows) == RHO_COUNT || error("raw curve point count mismatch at T=$(T)")
        path = joinpath(output_dir, "production_eval", "prod_eval_xi_$(_token(xi))_T_$(_token(T)).csv")
        _write_curve(path, rows)
        push!(curve_files, relpath(path, output_dir))
    end
    failed_points == 0 || error("raw-only Oracle production had $(failed_points) failed points")

    snapshot = MODELS.TrhoScan.rho_session_snapshot(session)
    summary = Dict(
        "schema_version" => "pnjl_c2_raw_only_oracle_job_v1",
        "scope" => "c2_full_grid_raw_rho_mu",
        "method" => "independent_oracle",
        "xi" => xi,
        "temperatures_MeV" => temperatures,
        "rho_count_per_curve" => length(rho_grid),
        "curve_count" => length(curve_files),
        "source_run_id" => source_run_id,
        "source_grid_run_id" => source_grid_run_id,
        "source_grid_artifact_name" => source_grid_artifact_name,
        "source_grid_manifest_path" => source_grid_manifest_path,
        "source_grid_manifest_sha256" => source_grid_manifest_sha256,
        "calculation_sha" => calculation_sha,
        "source_postprocess_sha" => source_postprocess_sha,
        "source_workflow_sha" => source_workflow_sha,
        "solver_called" => true,
        "phase_pipeline_called" => false,
        "higher_order_postprocess_called" => false,
        "curve_files" => curve_files,
        "session" => Dict(String(field) => getproperty(snapshot, field) for field in propertynames(snapshot)),
    )
    mkpath(output_dir)
    open(joinpath(output_dir, "raw_only_job_summary.json"), "w") do io
        JSON3.pretty(io, summary)
        write(io, '\n')
    end
    println(JSON3.write(summary))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    _run(ARGS)
end
