#!/usr/bin/env julia

"""
Small cross-solver diagnostic for the local pnjl_mag equilibrium replay.

This script intentionally does not create a validation target. It evaluates the
external five-field states through the current Julia MFIR route and runs a
small multi-seed Julia candidate search at the same representative points.
"""

using Printf
using StaticArrays
using LinearAlgebra: norm

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const INPUT_CSV = joinpath(
    PROJECT_ROOT,
    "docs",
    "analysis",
    "historical",
    "legacy",
    "legacy_extraction_v1",
    "pnjl_mag_equilibrium_replay_v1",
    "pnjl_mag_equilibrium_replay_v1.csv",
)
const OUTPUT_DIR = joinpath(
    PROJECT_ROOT,
    "docs",
    "analysis",
    "historical",
    "legacy",
    "legacy_extraction_v1",
    "pnjl_mag_cross_solver_replay_v1",
)

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()

const PNJL = Models.pnjl_module()
const EXTERNAL_HC_MEV_FM = 197.33

function _external_parity_params()
    Λ_inv_fm = 602.3 / EXTERNAL_HC_MEV_FM
    return PNJL.PNJLCore.PNJLParams(
        profile="pnjl_mag_external_replay",
        physics_profile="pnjl_mag_external_replay",
        hbarc_MeV_fm=EXTERNAL_HC_MEV_FM,
        alpha_em=1 / 137.035999084,
        N_color=3,
        N_flavor=3,
        rho0_fm3=0.16,
        Λ_inv_fm=Λ_inv_fm,
        m_ud0_inv_fm=5.5 / EXTERNAL_HC_MEV_FM,
        m_s0_inv_fm=140.7 / EXTERNAL_HC_MEV_FM,
        G_fm2=1.835 / Λ_inv_fm^2,
        K_fm5=12.36 / Λ_inv_fm^5,
        thermal_p_max_inv_fm=40.0,
        T0_inv_fm=210.0 / EXTERNAL_HC_MEV_FM,
        polyakov_scheme=:log,
        a0=3.51,
        a1=-2.47,
        a2=15.2,
        a3=0.0,
        b3=-1.75,
        b4=7.555,
    )
end

function _external_parity_imc()
    return PNJL.MagneticThermodynamics.MagneticIMCParams(
        a=0.0108805,
        b=-1.0133e-4,
        c=0.02228,
        d=1.84558e-4,
        Λ_QCD_MeV=300.0 * Constants_PNJL.ħc_MeV_fm / EXTERNAL_HC_MEV_FM,
    )
end

function _parse_float(value::AbstractString)
    return parse(Float64, strip(value))
end

function _read_external_rows(path::String)
    isfile(path) || error("external replay CSV not found: $path")
    lines = readlines(path)
    length(lines) >= 2 || error("external replay CSV has no data rows: $path")
    header = split(strip(lines[1]), ',')
    expected = [
        "T_MeV", "muB_MeV", "eB_GeV2", "eB_fm_minus2", "phi_u", "phi_d",
        "phi_s", "Phi1", "Phi2", "omega_fm4", "gap_residual_max",
        "committed_state_max_abs_delta", "committed_gap_residual_max",
    ]
    header == expected || error("unexpected external replay header")
    rows = NamedTuple[]
    for line in lines[2:end]
        cols = split(strip(line), ',')
        length(cols) == length(expected) || error("invalid external replay row: $line")
        push!(rows, (
            T_MeV=_parse_float(cols[1]),
            muB_MeV=_parse_float(cols[2]),
            eB_GeV2=_parse_float(cols[3]),
            eB_fm_minus2=_parse_float(cols[4]),
            phi_u=_parse_float(cols[5]),
            phi_d=_parse_float(cols[6]),
            phi_s=_parse_float(cols[7]),
            Phi1=_parse_float(cols[8]),
            Phi2=_parse_float(cols[9]),
            omega_fm4=_parse_float(cols[10]),
            gap_residual_max=_parse_float(cols[11]),
            committed_state_max_abs_delta=_parse_float(cols[12]),
            committed_gap_residual_max=_parse_float(cols[13]),
        ))
    end
    return rows
end

function _parse_args(args::Vector{String})
    scope = :full
    nodes = :screening
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--scope"
            i == length(args) && error("missing value for --scope")
            i += 1
            scope = Symbol(lowercase(args[i]))
            scope in (:smoke, :high, :full) || error("scope must be smoke, high or full")
        elseif arg == "--nodes"
            i == length(args) && error("missing value for --nodes")
            i += 1
            nodes = Symbol(lowercase(args[i]))
            nodes in (:screening, :matched) || error("nodes must be screening or matched")
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/analysis/pnjl/compare_pnjl_mag_julia_equilibrium.jl [--scope smoke|full] [--nodes screening|matched]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return scope, nodes
end

function _controls(nodes::Symbol)
    if nodes === :matched
        return (p_num=128, zeta_num=256, pz_max=40.0, n_max=79, label="matched_external_discretization")
    end
    return (p_num=64, zeta_num=64, pz_max=25.0, n_max=79, label="screening_discretization")
end

function _external_parity_model(eB_fm2::Float64, controls)
    params = _external_parity_params()
    base = PNJL.PNJLModel(params)
    conf = PNJL.MagneticThermodynamics.MagneticConfig(
        eB_fm2=eB_fm2,
        n_max=controls.n_max,
        p_num=controls.p_num,
        pz_max=controls.pz_max,
        imc=_external_parity_imc(),
        route=:mfir,
        zeta_num=controls.zeta_num,
        params=params,
    )
    return PNJL.PNJLMagneticModel(base, conf)
end

function _write_point_header(io)
    println(io, "T_MeV,muB_MeV,eB_GeV2,eB_fm_minus2,external_gap_residual_max,julia_external_residual_norm,julia_external_omega_fm4,external_omega_fm4,omega_delta_julia_minus_external,attempt_count,failed_attempts,candidate_count,status,error")
end

function _write_candidate_header(io)
    println(io, "T_MeV,muB_MeV,eB_GeV2,eB_fm_minus2,candidate_index,seed_index,seed_label,converged,physical,method,iterations,n_max,residual_norm,omega_fm4,state_l2_delta_external,omega_delta_from_julia_external,phi_u,phi_d,phi_s,Phi1,Phi2")
end

function _candidate_seeds(external_state)
    return [
        (label="external_committed", state=external_state),
        (label="external_chiral_shift", state=external_state + SVector{5, Float64}(0.20, 0.20, 0.20, 0.0, 0.0)),
        (label="project_thermal", state=SVector{5, Float64}(-0.03, -0.03, -0.04, 0.50, 0.50)),
    ]
end

function _csv_value(value)
    value isa AbstractString && return replace(value, ',' => ';')
    return string(value)
end

function main(args::Vector{String})
    scope, nodes = _parse_args(args)
    controls = _controls(nodes)
    rows = _read_external_rows(INPUT_CSV)
    if scope === :smoke
        rows = [row for row in rows if row.T_MeV in (50.0, 150.0, 240.0) && row.eB_GeV2 == 0.2]
    elseif scope === :high
        rows = [row for row in rows if row.T_MeV == 240.0]
    end
    isempty(rows) && error("no rows selected for scope=$scope")
    mkpath(OUTPUT_DIR)
    output_points = joinpath(OUTPUT_DIR, "julia_point_summary_$(scope)_$(nodes)_v1.csv")
    output_candidates = joinpath(OUTPUT_DIR, "julia_candidates_$(scope)_$(nodes)_v1.csv")

    open(output_points, "w") do point_io
        open(output_candidates, "w") do candidate_io
            _write_point_header(point_io)
            _write_candidate_header(candidate_io)
            for row in rows
                T_fm = row.T_MeV / EXTERNAL_HC_MEV_FM
                μ_fm = row.muB_MeV / (3.0 * EXTERNAL_HC_MEV_FM)
                mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)
                external_state = SVector{5, Float64}(
                    row.phi_u, row.phi_d, row.phi_s, row.Phi1, row.Phi2,
                )
                model = _external_parity_model(row.eB_fm_minus2, controls)

                external_residual = PNJL.magnetic_gap_residual(
                    model,
                    external_state,
                    T_fm,
                    mu_vec;
                    p_num=controls.p_num,
                    pz_max=controls.pz_max,
                    n_max=controls.n_max,
                )
                external_residual_norm = norm(external_residual)
                julia_external_omega = PNJL.omega(
                    model,
                    external_state,
                    T_fm,
                    mu_vec;
                    p_num=controls.p_num,
                    pz_max=controls.pz_max,
                    n_max=controls.n_max,
                )

                result = nothing
                err_text = ""
                try
                    seeds = _candidate_seeds(external_state)
                    result = PNJL.solve_magnetic_gap(
                        model,
                        T_fm,
                        mu_vec;
                        p_num=controls.p_num,
                        pz_max=controls.pz_max,
                        n_max=controls.n_max,
                        initial_guess=seeds[1].state,
                        seed_candidates=[seed.state for seed in seeds[2:end]],
                        include_default_seeds=false,
                        method=:trust_region,
                        fallback_method=:newton,
                        iterations=120,
                        residual_norm_max=1e-6,
                        root_merge_tol=1e-5,
                    )
                catch err
                    err_text = sprint(showerror, err)
                end

                if result === nothing
                    println(point_io, join(_csv_value.(Any[
                        row.T_MeV, row.muB_MeV, row.eB_GeV2, row.eB_fm_minus2,
                        row.gap_residual_max, external_residual_norm, julia_external_omega,
                        row.omega_fm4, julia_external_omega - row.omega_fm4,
                        0, 0, 0, "solver_error", err_text,
                    ]), ','))
                    continue
                end

                println(point_io, join(_csv_value.(Any[
                    row.T_MeV, row.muB_MeV, row.eB_GeV2, row.eB_fm_minus2,
                    row.gap_residual_max, external_residual_norm, julia_external_omega,
                    row.omega_fm4, julia_external_omega - row.omega_fm4,
                    result.attempt_count, result.failed_attempts, length(result.candidates),
                    "ok", "",
                ]), ','))

                for (candidate_index, candidate) in enumerate(result.candidates)
                    seed_label = candidate.seed_index <= length(_candidate_seeds(external_state)) ?
                        _candidate_seeds(external_state)[candidate.seed_index].label : "unknown"
                    state_delta = norm(candidate.x_state - external_state)
                    println(candidate_io, join(_csv_value.(Any[
                        row.T_MeV, row.muB_MeV, row.eB_GeV2, row.eB_fm_minus2,
                        candidate_index, candidate.seed_index, seed_label,
                        candidate.converged, candidate.physical, candidate.method,
                        candidate.iterations, candidate.n_max, candidate.residual_norm,
                        candidate.omega, state_delta, candidate.omega - julia_external_omega,
                        candidate.x_state[1], candidate.x_state[2], candidate.x_state[3],
                        candidate.x_state[4], candidate.x_state[5],
                    ]), ','))
                end
            end
        end
    end
    println("cross-solver diagnostic exported: $OUTPUT_DIR")
    println("point_summary=$output_points")
    println("candidates=$output_candidates")
    println("scope=$scope nodes=$(controls.label) points=$(length(rows))")
end

main(ARGS)
