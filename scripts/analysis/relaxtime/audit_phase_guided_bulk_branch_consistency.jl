#!/usr/bin/env julia

using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const DEFAULT_CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1"
const ANALYSIS_TABLE_DIR = joinpath(
    PROJECT_ROOT,
    "docs",
    "analysis",
    "relaxtime",
    "phase_guided_transport",
    "phase_guided_transport_v2_pole_sensitive_rendering",
    "tables",
)
const DEFAULT_OUT = joinpath(ANALYSIS_TABLE_DIR, "bulk_derivative_branch_audit.csv")
const DEFAULT_PHASE_ANCHOR_OUT = joinpath(ANALYSIS_TABLE_DIR, "phase_anchor_coexistence_audit.csv")
const TARGET_XI = (-0.02, -0.01, 0.0)
const BRANCH_MASS_REL_TOL = 0.10
const PRODUCTION_P_NUM = 12
const PRODUCTION_T_NUM = 6
const COEXISTENCE_XI = 0.0
const NEAR_COEXISTENCE_XI = 0.003
const COEXISTENCE_SCAN_OFFSETS_MEV = 0.4:0.1:1.2
const COEXISTENCE_BISECTION_STEPS = 16

function parse_args(args=ARGS)
    case_name = DEFAULT_CASE
    out_path = DEFAULT_OUT
    phase_anchor_out_path = DEFAULT_PHASE_ANCHOR_OUT
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--case-name"
            i += 1
            i <= length(args) || error("missing value for --case-name")
            case_name = args[i]
        elseif arg == "--out"
            i += 1
            i <= length(args) || error("missing value for --out")
            out_path = abspath(args[i])
        elseif arg == "--phase-anchor-out"
            i += 1
            i <= length(args) || error("missing value for --phase-anchor-out")
            phase_anchor_out_path = abspath(args[i])
        elseif arg in ("-h", "--help")
            println("usage: julia --project=. scripts/analysis/relaxtime/audit_phase_guided_bulk_branch_consistency.jl [--case-name NAME] [--out PATH] [--phase-anchor-out PATH]")
            exit(0)
        else
            error("unknown argument: $arg")
        end
        i += 1
    end
    return (; case_name, out_path, phase_anchor_out_path)
end

@inline relative_delta(a::Real, b::Real) =
    abs(Float64(a) - Float64(b)) /
    max(abs(Float64(a)), abs(Float64(b)), eps(Float64))

function production_scan_path(case_name::AbstractString)
    return joinpath(
        PROJECT_ROOT,
        "data",
        "outputs",
        "results",
        "relaxtime",
        "transport",
        "phase_guided",
        "mode_a_fixed_muB_phase_scaled",
        case_name,
        "phase_guided_transport_scan.csv",
    )
end

function select_target(rows, xi::Real)
    matches = [
        row for row in rows if
        row.plot_panel == "muB900.0" &&
        row.plot_series == "alpha1.0" &&
        isapprox(Float64(row.xi), Float64(xi); atol=1.0e-12, rtol=0.0)
    ]
    length(matches) == 1 || error("expected one production row at xi=$xi, got $(length(matches))")
    return only(matches)
end

function solve_seed_candidate(model, T_fm::Real, muq_fm::Real, xi::Real, seed)
    return Main.Models.solve(
        model,
        Main.Models.FixedMu(),
        Float64(T_fm),
        Float64(muq_fm);
        seed_guess=seed,
        xi=Float64(xi),
        p_num=PRODUCTION_P_NUM,
        t_num=PRODUCTION_T_NUM,
        residual_norm_max=1.0e-4,
    )
end

function branch_candidates(model, T_fm::Real, muq_fm::Real, xi::Real)
    candidate_a = solve_seed_candidate(
        model,
        T_fm,
        muq_fm,
        xi,
        Main.Models.HADRON_SEED_5,
    )
    candidate_b = solve_seed_candidate(
        model,
        T_fm,
        muq_fm,
        xi,
        Main.Models.QUARK_SEED_5,
    )
    mass_a = Float64(candidate_a.masses[1])
    mass_b = Float64(candidate_b.masses[1])
    distinct = relative_delta(mass_a, mass_b) > BRANCH_MASS_REL_TOL
    high, low = mass_a >= mass_b ? (candidate_a, candidate_b) : (candidate_b, candidate_a)
    delta_omega = Float64(high.omega) - Float64(low.omega)
    stable_role = if !distinct
        "single_converged_root"
    elseif delta_omega < 0.0
        "high_mass_hadron_candidate"
    elseif delta_omega > 0.0
        "low_mass_quark_candidate"
    else
        "coexisting_candidates"
    end
    return (
        distinct=distinct,
        high=high,
        low=low,
        delta_omega_high_minus_low=delta_omega,
        stable_role=stable_role,
    )
end

function nearest_candidate_role(value::Real, candidates)
    candidates.distinct || return "single_converged_root"
    high_delta = relative_delta(value, candidates.high.masses[1])
    low_delta = relative_delta(value, candidates.low.masses[1])
    return high_delta <= low_delta ?
        "high_mass_hadron_candidate" :
        "low_mass_quark_candidate"
end

function find_coexistence_bracket(model, reference_row)
    reference_T_MeV = Float64(reference_row.T_MeV)
    muq_fm = Float64(reference_row.muq_fm)
    previous = nothing
    for offset_MeV in COEXISTENCE_SCAN_OFFSETS_MEV
        T_MeV = reference_T_MeV + offset_MeV
        candidates = branch_candidates(
            model,
            T_MeV / Main.Models.Constants_PNJL.ħc_MeV_fm,
            muq_fm,
            COEXISTENCE_XI,
        )
        candidates.distinct || continue
        current = (
            T_MeV=T_MeV,
            delta_omega=candidates.delta_omega_high_minus_low,
        )
        if previous !== nothing && previous.delta_omega * current.delta_omega <= 0.0
            return previous, current
        end
        previous = current
    end
    error("failed to bracket the xi=0 coexistence temperature near the production phase anchor")
end

function refine_coexistence_bracket(model, reference_row, left, right)
    muq_fm = Float64(reference_row.muq_fm)
    lower = left
    upper = right
    for _ in 1:COEXISTENCE_BISECTION_STEPS
        midpoint_T_MeV = 0.5 * (lower.T_MeV + upper.T_MeV)
        candidates = branch_candidates(
            model,
            midpoint_T_MeV / Main.Models.Constants_PNJL.ħc_MeV_fm,
            muq_fm,
            COEXISTENCE_XI,
        )
        candidates.distinct || error("coexistence bisection lost the two-root branch bracket")
        midpoint = (
            T_MeV=midpoint_T_MeV,
            delta_omega=candidates.delta_omega_high_minus_low,
        )
        if lower.delta_omega * midpoint.delta_omega <= 0.0
            upper = midpoint
        else
            lower = midpoint
        end
    end
    return lower, upper
end

function near_coexistence_evidence(model, reference_row, lower, upper, xi::Real)
    muq_fm = Float64(reference_row.muq_fm)
    lower_candidates = branch_candidates(
        model,
        lower.T_MeV / Main.Models.Constants_PNJL.ħc_MeV_fm,
        muq_fm,
        xi,
    )
    upper_candidates = branch_candidates(
        model,
        upper.T_MeV / Main.Models.Constants_PNJL.ħc_MeV_fm,
        muq_fm,
        xi,
    )
    return (
        lower_delta_omega=lower_candidates.delta_omega_high_minus_low,
        upper_delta_omega=upper_candidates.delta_omega_high_minus_low,
        lower_distinct=lower_candidates.distinct,
        upper_distinct=upper_candidates.distinct,
    )
end

function main()
    opts = parse_args()
    scan_path = production_scan_path(opts.case_name)
    isfile(scan_path) || error("missing production scan: $scan_path")
    production_rows = collect(CSV.File(scan_path; comment="#"))
    audit_rows = NamedTuple[]
    model = Main.Models.create_model(:PNJL)

    for xi in TARGET_XI
        row = select_target(production_rows, xi)
        bulk = Main.Models.bulk_viscosity_coefficients(
            Float64(row.T_fm),
            Float64(row.muq_fm);
            xi=xi,
            p_num=PRODUCTION_P_NUM,
            t_num=PRODUCTION_T_NUM,
        )
        candidates = branch_candidates(model, Float64(row.T_fm), Float64(row.muq_fm), xi)
        main_m_u = Float64(row.m_u)
        bulk_m_u = Float64(bulk.masses[1])
        mass_rel_delta = relative_delta(main_m_u, bulk_m_u)
        branch_aligned = mass_rel_delta <= BRANCH_MASS_REL_TOL
        main_candidate_role = nearest_candidate_role(main_m_u, candidates)
        bulk_candidate_role = nearest_candidate_role(bulk_m_u, candidates)
        main_is_stable = main_candidate_role == candidates.stable_role || !candidates.distinct
        bulk_is_stable = bulk_candidate_role == candidates.stable_role || !candidates.distinct
        verdict = if branch_aligned
            "main_and_bulk_branch_aligned"
        elseif !main_is_stable && bulk_is_stable
            "main_continuation_metastable_bulk_stable_branch_mismatch"
        else
            "main_and_bulk_unresolved_branch_mismatch"
        end
        note = if branch_aligned
            "bulk derivative mass is compatible with the production equilibrium branch"
        elseif !main_is_stable && bulk_is_stable
            "the production continuation retained the metastable low-mass quark candidate while bulk selected the thermodynamically stable high-mass hadron candidate"
        else
            "main and bulk use different candidate branches; stability attribution requires additional evidence"
        end

        push!(audit_rows, (
            case_name=opts.case_name,
            mode_key="mode_a",
            plot_panel="muB900.0",
            plot_series="alpha1.0",
            xi=xi,
            T_MeV=Float64(row.T_MeV),
            muB_MeV=Float64(row.muB_MeV),
            thermo_p_num=PRODUCTION_P_NUM,
            thermo_t_num=PRODUCTION_T_NUM,
            main_m_u=main_m_u,
            bulk_m_u=bulk_m_u,
            main_m_s=Float64(row.m_s),
            bulk_m_s=Float64(bulk.masses[3]),
            light_mass_rel_delta=mass_rel_delta,
            branch_mass_rel_tol=BRANCH_MASS_REL_TOL,
            branch_aligned=branch_aligned,
            candidate_roots_distinct=candidates.distinct,
            high_mass_candidate_m_u=Float64(candidates.high.masses[1]),
            high_mass_candidate_Phi=Float64(candidates.high.x_state[4]),
            high_mass_candidate_omega=Float64(candidates.high.omega),
            low_mass_candidate_m_u=Float64(candidates.low.masses[1]),
            low_mass_candidate_Phi=Float64(candidates.low.x_state[4]),
            low_mass_candidate_omega=Float64(candidates.low.omega),
            delta_omega_high_minus_low=candidates.delta_omega_high_minus_low,
            stable_candidate_role=candidates.stable_role,
            main_candidate_role=main_candidate_role,
            bulk_candidate_role=bulk_candidate_role,
            main_is_stable=main_is_stable,
            bulk_is_stable=bulk_is_stable,
            production_tau_u=Float64(row.tau_u),
            production_entropy=Float64(row.s_fm3inv),
            production_zeta=Float64(row.zeta),
            production_zeta_over_s=Float64(row.zeta_over_s),
            bulk_v_n_sq=Float64(bulk.v_n_sq),
            bulk_dmuB_dT_sigma=Float64(bulk.dμB_dT_sigma),
            bulk_dM_u_dT=Float64(bulk.dM_dT[1]),
            bulk_dM_u_dmuB=Float64(bulk.dM_dμB[1]),
            verdict=verdict,
            mechanism_note=note,
        ))
    end

    mkpath(dirname(opts.out_path))
    CSV.write(opts.out_path, audit_rows)
    println("wrote $(relpath(opts.out_path, PROJECT_ROOT)) with $(length(audit_rows)) rows")

    reference_row = select_target(production_rows, COEXISTENCE_XI)
    coarse_left, coarse_right = find_coexistence_bracket(model, reference_row)
    lower, upper = refine_coexistence_bracket(model, reference_row, coarse_left, coarse_right)
    minus_evidence = near_coexistence_evidence(
        model,
        reference_row,
        lower,
        upper,
        -NEAR_COEXISTENCE_XI,
    )
    plus_evidence = near_coexistence_evidence(
        model,
        reference_row,
        lower,
        upper,
        NEAR_COEXISTENCE_XI,
    )
    minus_certified = minus_evidence.lower_distinct && minus_evidence.upper_distinct &&
        minus_evidence.lower_delta_omega > 0.0 && minus_evidence.upper_delta_omega > 0.0
    plus_certified = plus_evidence.lower_distinct && plus_evidence.upper_distinct &&
        plus_evidence.lower_delta_omega < 0.0 && plus_evidence.upper_delta_omega < 0.0
    phase_anchor_rows = [(
        case_name=opts.case_name,
        mode_key="mode_a",
        plot_panel="muB900.0",
        plot_series="alpha1.0",
        muB_MeV=Float64(reference_row.muB_MeV),
        reference_xi=COEXISTENCE_XI,
        thermo_p_num=PRODUCTION_P_NUM,
        thermo_t_num=PRODUCTION_T_NUM,
        interpolated_reference_T_MeV=Float64(reference_row.T_MeV),
        coexistence_T_lower_MeV=lower.T_MeV,
        coexistence_T_upper_MeV=upper.T_MeV,
        coexistence_T_mid_MeV=0.5 * (lower.T_MeV + upper.T_MeV),
        coexistence_T_bracket_width_MeV=upper.T_MeV - lower.T_MeV,
        coexistence_minus_interpolated_T_MeV=0.5 * (lower.T_MeV + upper.T_MeV) - Float64(reference_row.T_MeV),
        delta_omega_at_T_lower=lower.delta_omega,
        delta_omega_at_T_upper=upper.delta_omega,
        near_minus_xi=-NEAR_COEXISTENCE_XI,
        near_minus_delta_omega_at_T_lower=minus_evidence.lower_delta_omega,
        near_minus_delta_omega_at_T_upper=minus_evidence.upper_delta_omega,
        near_minus_stable_phase="quark",
        near_minus_certified=minus_certified,
        near_plus_xi=NEAR_COEXISTENCE_XI,
        near_plus_delta_omega_at_T_lower=plus_evidence.lower_delta_omega,
        near_plus_delta_omega_at_T_upper=plus_evidence.upper_delta_omega,
        near_plus_stable_phase="hadron",
        near_plus_certified=plus_certified,
        two_sided_bracket_certified=minus_certified && plus_certified,
        method="direct_two_branch_equal_omega_bisection",
        verdict="interpolated_boundary_anchor_is_not_the_direct_coexistence_temperature",
    )]
    mkpath(dirname(opts.phase_anchor_out_path))
    CSV.write(opts.phase_anchor_out_path, phase_anchor_rows)
    println("wrote $(relpath(opts.phase_anchor_out_path, PROJECT_ROOT)) with $(length(phase_anchor_rows)) row")
end

main()
