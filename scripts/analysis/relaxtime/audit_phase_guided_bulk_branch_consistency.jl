#!/usr/bin/env julia

using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const DEFAULT_CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1"
const DEFAULT_OUT = joinpath(
    PROJECT_ROOT,
    "docs",
    "analysis",
    "relaxtime",
    "phase_guided_transport_v2_pole_sensitive_rendering",
    "tables",
    "bulk_derivative_branch_audit.csv",
)
const TARGET_XI = (-0.02, -0.01, 0.0)
const BRANCH_MASS_REL_TOL = 0.10

function parse_args(args=ARGS)
    case_name = DEFAULT_CASE
    out_path = DEFAULT_OUT
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
        elseif arg in ("-h", "--help")
            println("usage: julia --project=. scripts/analysis/relaxtime/audit_phase_guided_bulk_branch_consistency.jl [--case-name NAME] [--out PATH]")
            exit(0)
        else
            error("unknown argument: $arg")
        end
        i += 1
    end
    return (; case_name, out_path)
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

function main()
    opts = parse_args()
    scan_path = production_scan_path(opts.case_name)
    isfile(scan_path) || error("missing production scan: $scan_path")
    production_rows = collect(CSV.File(scan_path; comment="#"))
    audit_rows = NamedTuple[]

    for xi in TARGET_XI
        row = select_target(production_rows, xi)
        bulk = Main.Models.bulk_viscosity_coefficients(
            Float64(row.T_fm),
            Float64(row.muq_fm);
            xi=xi,
            p_num=Main.Models.default_momentum_count(),
            t_num=Main.Models.default_theta_count(),
        )
        main_m_u = Float64(row.m_u)
        bulk_m_u = Float64(bulk.masses[1])
        mass_rel_delta = relative_delta(main_m_u, bulk_m_u)
        branch_aligned = mass_rel_delta <= BRANCH_MASS_REL_TOL
        verdict = branch_aligned ?
            "main_and_bulk_branch_aligned" :
            "bulk_derivative_branch_mismatch"
        note = branch_aligned ?
            "bulk derivative mass is compatible with the production equilibrium branch" :
            "bulk derivative path has already switched to the high-mass branch while the production equilibrium remains on the low-mass branch"

        push!(audit_rows, (
            case_name=opts.case_name,
            mode_key="mode_a",
            plot_panel="muB900.0",
            plot_series="alpha1.0",
            xi=xi,
            T_MeV=Float64(row.T_MeV),
            muB_MeV=Float64(row.muB_MeV),
            main_m_u=main_m_u,
            bulk_m_u=bulk_m_u,
            main_m_s=Float64(row.m_s),
            bulk_m_s=Float64(bulk.masses[3]),
            light_mass_rel_delta=mass_rel_delta,
            branch_mass_rel_tol=BRANCH_MASS_REL_TOL,
            branch_aligned=branch_aligned,
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
end

main()
