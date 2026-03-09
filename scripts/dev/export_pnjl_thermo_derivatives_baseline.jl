#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_thermo_derivatives_fixedpoints_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

const THERMO_POINTS = [
    (label="mu0_cold", T_fm=0.50, mu_fm=0.00),
    (label="dense_cold", T_fm=0.50, mu_fm=1.50),
    (label="intermediate_lowmu", T_fm=0.66, mu_fm=0.50),
    (label="intermediate_midmu", T_fm=0.75, mu_fm=0.75),
    (label="hot_lowmu", T_fm=0.90, mu_fm=0.15),
]

function parse_args(args::Vector{String})
    output = DEFAULT_OUTPUT
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--output"
            i == length(args) && error("missing value for --output")
            i += 1
            output = args[i]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/export_pnjl_thermo_derivatives_baseline.jl [--output <path>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return output
end

function main(args::Vector{String})
    output = parse_args(args)
    mkpath(dirname(output))

    PNJL = Models.pnjl_module()
    p_num = 24
    t_num = 6
    xi = 0.0

    open(output, "w") do io
        println(io, "label,T_fm,mu_fm,pressure,energy,rho,entropy,dP_dT,dP_dmu,dEpsilon_dT,dEpsilon_dmu,dn_dT,dn_dmu,dP_depsilon_n,dP_dn_epsilon,dM_dT_u,dM_dT_d,dM_dT_s,dM_dmu_u,dM_dmu_d,dM_dmu_s,v_n_sq,dmuB_dT_sigma,n_B")

        for pt in THERMO_POINTS
            td = PNJL.thermo_derivatives(pt.T_fm, pt.mu_fm; xi=xi, p_num=p_num, t_num=t_num)
            bv = PNJL.bulk_viscosity_coefficients(pt.T_fm, pt.mu_fm; xi=xi, p_num=p_num, t_num=t_num)

            @printf(io,
                "%s,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                pt.label,
                pt.T_fm,
                pt.mu_fm,
                td.pressure,
                td.energy,
                td.rho,
                td.entropy,
                td.dP_dT,
                td.dP_dmu,
                td.dEpsilon_dT,
                td.dEpsilon_dmu,
                td.dn_dT,
                td.dn_dmu,
                td.dP_depsilon_n,
                td.dP_dn_epsilon,
                td.dM_dT[1], td.dM_dT[2], td.dM_dT[3],
                td.dM_dmu[1], td.dM_dmu[2], td.dM_dmu[3],
                bv.v_n_sq,
                bv.dμB_dT_sigma,
                bv.n_B,
            )
        end
    end

    println("PNJL thermo derivatives baseline written to: " * output)
    println(@sprintf("points: %d", length(THERMO_POINTS)))
end

main(ARGS)