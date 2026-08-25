#!/usr/bin/env julia

"""
Fixed-state convergence diagnostic for the pnjl_mag/Julia magnetic comparison.

The source-parity and production-parity outputs are deliberately written to
separate CSV files. This script does not create a validation target.
"""

using Printf
using StaticArrays
using LinearAlgebra: norm

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const INPUT_CSV = joinpath(
    PROJECT_ROOT, "docs", "analysis", "historical", "legacy",
    "legacy_extraction_v1", "pnjl_mag_equilibrium_replay_v1",
    "pnjl_mag_equilibrium_replay_v1.csv",
)
const OUTPUT_DIR = joinpath(
    PROJECT_ROOT, "docs", "analysis", "historical", "legacy",
    "legacy_extraction_v1", "pnjl_mag_convergence_v1",
)
const EXTERNAL_HC_MEV_FM = 197.33

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()
const PNJL = Models.pnjl_module()

function _parse_float(value::AbstractString)
    return parse(Float64, strip(value))
end

function _read_rows(path::String)
    lines = readlines(path)
    header = split(strip(first(lines)), ',')
    expected = [
        "T_MeV", "muB_MeV", "eB_GeV2", "eB_fm_minus2", "phi_u", "phi_d",
        "phi_s", "Phi1", "Phi2", "omega_fm4", "gap_residual_max",
        "committed_state_max_abs_delta", "committed_gap_residual_max",
    ]
    header == expected || error("unexpected replay artifact header")
    rows = NamedTuple[]
    for line in lines[2:end]
        cols = split(strip(line), ',')
        length(cols) == length(expected) || error("invalid replay row")
        push!(rows, (
            T_MeV=_parse_float(cols[1]),
            muB_MeV=_parse_float(cols[2]),
            eB_GeV2=_parse_float(cols[3]),
            eB_fm_minus2=_parse_float(cols[4]),
            state=SVector{5, Float64}(
                _parse_float(cols[5]), _parse_float(cols[6]), _parse_float(cols[7]),
                _parse_float(cols[8]), _parse_float(cols[9]),
            ),
            external_omega=_parse_float(cols[10]),
        ))
    end
    selected = [
        row for row in rows if
        (row.T_MeV == 50.0 && row.eB_GeV2 == 0.2) ||
        (row.T_MeV == 240.0 && row.eB_GeV2 in (0.2, 0.8))
    ]
    length(selected) == 3 || error("expected three convergence anchors, got $(length(selected))")
    sort!(selected; by=row -> (row.T_MeV, row.eB_GeV2))
    return selected
end

function _controls()
    return [
        (name="screening", p_num=64, zeta_num=64, pz_max=25.0, n_max=79),
        (name="mid", p_num=96, zeta_num=128, pz_max=32.0, n_max=79),
        (name="matched", p_num=128, zeta_num=256, pz_max=40.0, n_max=79),
        (name="landau_probe", p_num=128, zeta_num=256, pz_max=40.0, n_max=95),
    ]
end

function _cutoff_controls()
    return [
        (name="nmax31", p_num=128, zeta_num=256, pz_max=40.0, n_max=31),
        (name="nmax47", p_num=128, zeta_num=256, pz_max=40.0, n_max=47),
        (name="nmax63", p_num=128, zeta_num=256, pz_max=40.0, n_max=63),
        (name="nmax79", p_num=128, zeta_num=256, pz_max=40.0, n_max=79),
        (name="nmax95", p_num=128, zeta_num=256, pz_max=40.0, n_max=95),
        (name="nmax127", p_num=128, zeta_num=256, pz_max=40.0, n_max=127),
        (name="nmax159", p_num=128, zeta_num=256, pz_max=40.0, n_max=159),
        (name="nmax255", p_num=128, zeta_num=256, pz_max=40.0, n_max=255),
        (name="nmax383", p_num=128, zeta_num=256, pz_max=40.0, n_max=383),
        (name="nmax511", p_num=128, zeta_num=256, pz_max=40.0, n_max=511),
    ]
end

function _quadrature_controls()
    return [
        (name="p64", p_num=64, zeta_num=256, pz_max=40.0, n_max=383),
        (name="p96", p_num=96, zeta_num=256, pz_max=40.0, n_max=383),
        (name="reference", p_num=128, zeta_num=256, pz_max=40.0, n_max=383),
        (name="p160", p_num=160, zeta_num=256, pz_max=40.0, n_max=383),
        (name="pz25", p_num=128, zeta_num=256, pz_max=25.0, n_max=383),
        (name="pz32", p_num=128, zeta_num=256, pz_max=32.0, n_max=383),
        (name="pz48", p_num=128, zeta_num=256, pz_max=48.0, n_max=383),
        (name="zeta64", p_num=128, zeta_num=64, pz_max=40.0, n_max=383),
        (name="zeta128", p_num=128, zeta_num=128, pz_max=40.0, n_max=383),
        (name="zeta384", p_num=128, zeta_num=384, pz_max=40.0, n_max=383),
        (name="zeta512", p_num=128, zeta_num=512, pz_max=40.0, n_max=383),
    ]
end

function _solver_cutoff_controls()
    return [
        (name="nmax79", p_num=128, zeta_num=384, pz_max=40.0, n_max=79),
        (name="nmax383", p_num=128, zeta_num=384, pz_max=40.0, n_max=383),
    ]
end

function _parse_mode(args::Vector{String})
    mode = :fixed
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--mode"
            i == length(args) && error("missing value for --mode")
            i += 1
            mode = Symbol(lowercase(args[i]))
            mode in (:fixed, :solver, :branch, :cutoff, :quadrature, :solver_cutoff, :default, :both) || error("mode must be fixed, solver, branch, cutoff, quadrature, solver_cutoff, default or both")
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/analysis/pnjl/compare_pnjl_mag_convergence.jl [--mode fixed|solver|branch|cutoff|quadrature|solver_cutoff|default|both]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return mode
end

function _seed_permutations(seeds)
    return [
        (seeds[i], seeds[j], seeds[k])
        for i in eachindex(seeds), j in eachindex(seeds), k in eachindex(seeds)
        if i != j && i != k && j != k
    ]
end

function _run_branch_profile(row, control, profile::Symbol, output_path::String)
    seeds = _solver_seeds(row.state)
    permutations = _seed_permutations(seeds)
    open(output_path, "w") do io
        println(io, "profile,T_MeV,muB_MeV,eB_GeV2,eB_fm_minus2,node,p_num,zeta_num,pz_max_fm_inv,n_max,permutation_index,seed_order,candidate_index,origin_seed_label,method,iterations,residual_norm,omega_fm4,phi_u,phi_d,phi_s,Phi,PhiBar,attempt_count,failed_attempts,candidate_count,status,error")
        for (permutation_index, permutation) in enumerate(permutations)
            model, T_fm, μ_fm = profile === :source ?
                _source_model(row, control) : _production_model(row, control)
            mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)
            seed_order = join((seed.label for seed in permutation), '>')
            result = nothing
            error_text = ""
            try
                result = PNJL.solve_magnetic_gap(
                    model,
                    T_fm,
                    mu_vec;
                    p_num=control.p_num,
                    pz_max=control.pz_max,
                    n_max=control.n_max,
                    initial_guess=permutation[1].state,
                    seed_candidates=[permutation[2].state, permutation[3].state],
                    include_default_seeds=false,
                    method=:trust_region,
                    fallback_method=:newton,
                    iterations=120,
                    residual_norm_max=1e-6,
                    root_merge_tol=1e-5,
                )
            catch err
                error_text = sprint(showerror, err)
            end
            eB_fm2 = model.magnetic.eB_fm2
            if result === nothing
                println(io, join(_csv.(Any[
                    String(profile), row.T_MeV, row.muB_MeV, row.eB_GeV2, eB_fm2,
                    control.name, control.p_num, control.zeta_num, control.pz_max,
                    control.n_max, permutation_index, seed_order, "", "", "", "",
                    "", "", "", "", "", "", "", 0, 0, 0, "solver_error",
                    error_text,
                ]), ','))
                continue
            end
            for (candidate_index, candidate) in enumerate(result.candidates)
                origin_seed_label = permutation[candidate.seed_index].label
                println(io, join(_csv.(Any[
                    String(profile), row.T_MeV, row.muB_MeV, row.eB_GeV2, eB_fm2,
                    control.name, control.p_num, control.zeta_num, control.pz_max,
                    control.n_max, permutation_index, seed_order, candidate_index,
                    origin_seed_label, candidate.method, candidate.iterations,
                    candidate.residual_norm, candidate.omega, candidate.x_state[1],
                    candidate.x_state[2], candidate.x_state[3], candidate.x_state[4],
                    candidate.x_state[5], result.attempt_count, result.failed_attempts,
                    length(result.candidates), "ok", "",
                ]), ','))
            end
        end
    end
end

function _external_params()
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

function _external_imc()
    return PNJL.MagneticThermodynamics.MagneticIMCParams(
        a=0.0108805,
        b=-1.0133e-4,
        c=0.02228,
        d=1.84558e-4,
        Λ_QCD_MeV=300.0 * Constants_PNJL.ħc_MeV_fm / EXTERNAL_HC_MEV_FM,
    )
end

function _source_model(row, control)
    params = _external_params()
    base = PNJL.PNJLModel(params)
    conf = PNJL.MagneticThermodynamics.MagneticConfig(
        eB_fm2=row.eB_fm_minus2,
        n_max=control.n_max,
        p_num=control.p_num,
        pz_max=control.pz_max,
        imc=_external_imc(),
        route=:mfir,
        zeta_num=control.zeta_num,
        params=params,
    )
    return PNJL.PNJLMagneticModel(base, conf), row.T_MeV / EXTERNAL_HC_MEV_FM,
        row.muB_MeV / (3.0 * EXTERNAL_HC_MEV_FM)
end

function _production_model(row, control)
    eB_fm2 = row.eB_GeV2 * (1000.0 / Constants_PNJL.ħc_MeV_fm)^2
    model = PNJL.PNJLMagneticModel(
        eB_fm2=eB_fm2,
        p_num=control.p_num,
        pz_max=control.pz_max,
        n_max=control.n_max,
        zeta_num=control.zeta_num,
    )
    return model, row.T_MeV / Constants_PNJL.ħc_MeV_fm,
        row.muB_MeV / (3.0 * Constants_PNJL.ħc_MeV_fm)
end

function _write_header(io)
    println(io, "profile,T_MeV,muB_MeV,eB_GeV2,eB_fm_minus2,node,p_num,zeta_num,pz_max_fm_inv,n_max,fixed_residual_norm,fixed_omega_fm4,external_omega_fm4,omega_delta_vs_external,omega_comparable,status")
end

function _csv(value)
    value === nothing && return ""
    return value isa AbstractString ? replace(value, ',' => ';') : string(value)
end

function _run_profile(rows, controls, profile::Symbol, output_path::String)
    open(output_path, "w") do io
        _write_header(io)
        for row in rows
            for control in controls
                model, T_fm, μ_fm = profile === :source ?
                    _source_model(row, control) : _production_model(row, control)
                eB_fm2 = model.magnetic.eB_fm2
                mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)
                residual = PNJL.magnetic_gap_residual(
                    model, row.state, T_fm, mu_vec;
                    p_num=control.p_num, pz_max=control.pz_max, n_max=control.n_max,
                )
                omega = PNJL.omega(
                    model, row.state, T_fm, mu_vec;
                    p_num=control.p_num, pz_max=control.pz_max, n_max=control.n_max,
                )
                comparable = profile === :source
                ext_omega = comparable ? row.external_omega : ""
                delta = comparable ? omega - row.external_omega : ""
                println(io, join(_csv.(Any[
                    String(profile), row.T_MeV, row.muB_MeV, row.eB_GeV2, eB_fm2,
                    control.name, control.p_num, control.zeta_num, control.pz_max,
                    control.n_max, norm(residual), omega, ext_omega, delta,
                    comparable, "ok",
                ]), ','))
            end
        end
    end
end

function _solver_seeds(state)
    return [
        (label="external_committed", state=state),
        (label="external_chiral_shift", state=state + SVector{5, Float64}(0.20, 0.20, 0.20, 0.0, 0.0)),
        (label="project_thermal", state=SVector{5, Float64}(-0.03, -0.03, -0.04, 0.50, 0.50)),
    ]
end

function _run_solver_profile(rows, controls, profile::Symbol, output_path::String)
    open(output_path, "w") do io
        println(io, "profile,T_MeV,muB_MeV,eB_GeV2,eB_fm_minus2,node,p_num,zeta_num,pz_max_fm_inv,configured_n_max,resolved_n_max,candidate_index,seed_index,seed_label,converged,physical,method,iterations,residual_norm,omega_fm4,phi_u,phi_d,phi_s,Phi,PhiBar,state_l2_delta_input,status,error,attempt_count,failed_attempts,candidate_count")
        for row in rows
            for control in controls
                model, T_fm, μ_fm = profile === :source ?
                    _source_model(row, control) : _production_model(row, control)
                mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)
                seeds = _solver_seeds(row.state)
                result = nothing
                error_text = ""
                try
                    result = PNJL.solve_magnetic_gap(
                        model,
                        T_fm,
                        mu_vec;
                        p_num=control.p_num,
                        pz_max=control.pz_max,
                        n_max=control.n_max,
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
                    error_text = sprint(showerror, err)
                end
                eB_fm2 = model.magnetic.eB_fm2
                if result === nothing
                    println(io, join(_csv.(Any[
                        String(profile), row.T_MeV, row.muB_MeV, row.eB_GeV2, eB_fm2,
                        control.name, control.p_num, control.zeta_num, control.pz_max,
                        control.n_max, "", "", "", "", "", "", "", "", "", "", "",
                        "", "", "", "", "", "solver_error", error_text, 0, 0, 0,
                    ]), ','))
                    continue
                end
                for (candidate_index, candidate) in enumerate(result.candidates)
                    seed_label = candidate.seed_index <= length(seeds) ? seeds[candidate.seed_index].label : "unknown"
                    println(io, join(_csv.(Any[
                        String(profile), row.T_MeV, row.muB_MeV, row.eB_GeV2, eB_fm2,
                        control.name, control.p_num, control.zeta_num, control.pz_max,
                        control.n_max, candidate.n_max, candidate_index, candidate.seed_index, seed_label,
                        candidate.converged, candidate.physical, candidate.method,
                        candidate.iterations, candidate.residual_norm, candidate.omega,
                        candidate.x_state[1], candidate.x_state[2], candidate.x_state[3],
                        candidate.x_state[4], candidate.x_state[5],
                        norm(candidate.x_state - row.state), "ok", "",
                        result.attempt_count, result.failed_attempts, length(result.candidates),
                    ]), ','))
                end
            end
        end
    end
end

function _run_production_default_probe(rows, output_path::String)
    control = (name="production_default", p_num=96, zeta_num=64, pz_max=25.0, n_max=nothing)
    open(output_path, "w") do io
        println(io, "profile,T_MeV,muB_MeV,eB_GeV2,eB_fm_minus2,p_num,zeta_num,pz_max_fm_inv,configured_n_max,resolved_n_max,fixed_residual_norm,fixed_omega_fm4,status")
        for row in rows
            model, T_fm, μ_fm = _production_model(row, control)
            mu_vec = SVector{3, Float64}(μ_fm, μ_fm, μ_fm)
            comp = PNJL.MagneticThermodynamics.calculate_magnetic_omega_components(
                row.state, mu_vec, T_fm, model.magnetic,
            )
            residual = PNJL.magnetic_gap_residual(
                model, row.state, T_fm, mu_vec;
                p_num=control.p_num, pz_max=control.pz_max, n_max=comp.n_max,
            )
            println(io, join(_csv.(Any[
                "production", row.T_MeV, row.muB_MeV, row.eB_GeV2,
                model.magnetic.eB_fm2, control.p_num, control.zeta_num,
                control.pz_max, control.n_max, comp.n_max, norm(residual),
                comp.omega, "ok",
            ]), ','))
        end
    end
end

function main(args=ARGS)
    mode = _parse_mode(args)
    rows = _read_rows(INPUT_CSV)
    controls = _controls()
    mkpath(OUTPUT_DIR)
    if mode in (:fixed, :both)
        source_path = joinpath(OUTPUT_DIR, "source_parity_fixed_state_v1.csv")
        production_path = joinpath(OUTPUT_DIR, "production_parity_fixed_state_v1.csv")
        _run_profile(rows, controls, :source, source_path)
        _run_profile(rows, controls, :production, production_path)
        println("source_parity_fixed=$source_path")
        println("production_parity_fixed=$production_path")
        println("fixed_rows=$(length(rows) * length(controls)) per profile")
    end
    if mode in (:solver, :both)
        solver_controls = controls[3:4]
        source_path = joinpath(OUTPUT_DIR, "source_parity_solver_v1.csv")
        production_path = joinpath(OUTPUT_DIR, "production_parity_solver_v1.csv")
        _run_solver_profile(rows, solver_controls, :source, source_path)
        _run_solver_profile(rows, solver_controls, :production, production_path)
        println("source_parity_solver=$source_path")
        println("production_parity_solver=$production_path")
        println("solver_rows=$(length(rows) * length(solver_controls)) per profile before candidate expansion")
    end
    if mode in (:branch, :both)
        row = only(row for row in rows if row.T_MeV == 240.0 && row.eB_GeV2 == 0.8)
        control = controls[3]
        source_path = joinpath(OUTPUT_DIR, "source_parity_branch_repeatability_v1.csv")
        production_path = joinpath(OUTPUT_DIR, "production_parity_branch_repeatability_v1.csv")
        _run_branch_profile(row, control, :source, source_path)
        _run_branch_profile(row, control, :production, production_path)
        println("source_parity_branch=$source_path")
        println("production_parity_branch=$production_path")
        println("branch_seed_permutations=6 per profile")
    end
    if mode == :cutoff
        cutoff_controls = _cutoff_controls()
        source_path = joinpath(OUTPUT_DIR, "source_parity_nmax_fixed_state_v1.csv")
        production_path = joinpath(OUTPUT_DIR, "production_parity_nmax_fixed_state_v1.csv")
        _run_profile(rows, cutoff_controls, :source, source_path)
        _run_profile(rows, cutoff_controls, :production, production_path)
        println("source_parity_nmax=$source_path")
        println("production_parity_nmax=$production_path")
        println("nmax_fixed_rows=$(length(rows) * length(cutoff_controls)) per profile")
    end
    if mode == :quadrature
        quadrature_controls = _quadrature_controls()
        source_path = joinpath(OUTPUT_DIR, "source_parity_quadrature_fixed_state_v1.csv")
        production_path = joinpath(OUTPUT_DIR, "production_parity_quadrature_fixed_state_v1.csv")
        _run_profile(rows, quadrature_controls, :source, source_path)
        _run_profile(rows, quadrature_controls, :production, production_path)
        println("source_parity_quadrature=$source_path")
        println("production_parity_quadrature=$production_path")
        println("quadrature_fixed_rows=$(length(rows) * length(quadrature_controls)) per profile")
    end
    if mode == :solver_cutoff
        solver_cutoff_controls = _solver_cutoff_controls()
        source_path = joinpath(OUTPUT_DIR, "source_parity_solver_cutoff_v1.csv")
        production_path = joinpath(OUTPUT_DIR, "production_parity_solver_cutoff_v1.csv")
        _run_solver_profile(rows, solver_cutoff_controls, :source, source_path)
        _run_solver_profile(rows, solver_cutoff_controls, :production, production_path)
        println("source_parity_solver_cutoff=$source_path")
        println("production_parity_solver_cutoff=$production_path")
        println("solver_cutoff_rows=$(length(rows) * length(solver_cutoff_controls)) per profile before candidate expansion")
    end
    if mode == :default
        fixed_path = joinpath(OUTPUT_DIR, "production_default_nmax_fixed_state_v1.csv")
        solver_path = joinpath(OUTPUT_DIR, "production_default_nmax_solver_v1.csv")
        control = (name="production_default", p_num=96, zeta_num=64, pz_max=25.0, n_max=nothing)
        _run_production_default_probe(rows, fixed_path)
        _run_solver_profile(rows, [control], :production, solver_path)
        println("production_default_nmax_fixed=$fixed_path")
        println("production_default_nmax_solver=$solver_path")
        println("production_default_rows=$(length(rows))")
    end
end

main()
