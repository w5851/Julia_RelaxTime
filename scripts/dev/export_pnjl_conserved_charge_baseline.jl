#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_conserved_charge_fixedpoints_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

const CONSERVED_CHARGE_POINTS = [
    (label="bqs_origin", T_fm=0.50, muB_fm=0.00, muQ_fm=0.00, muS_fm=0.00, V=64.0),
    (label="bqs_asymmetric", T_fm=0.57, muB_fm=0.18, muQ_fm=0.05, muS_fm=0.02, V=64.0),
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
            println("Usage: julia --project=. scripts/dev/export_pnjl_conserved_charge_baseline.jl [--output <path>]")
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

    p_num = 16
    t_num = 6
    xi = 0.0
    PNJL = Models.pnjl_module()

    open(output, "w") do io
        println(io, "label,T_fm,muB_fm,muQ_fm,muS_fm,V,chi1_B,chi2_B,chi3_B,chi4_B,chi1_Q,chi2_Q,chi3_Q,chi4_Q,chi1_S,chi2_S,chi3_S,chi4_S,chi11_BQ,chi11_BS,chi11_QS,C200,C020,C002,C110,C101,C011,Ssigma,kappa_sigma2")

        for pt in CONSERVED_CHARGE_POINTS
            kwargs = (; xi=xi, p_num=p_num, t_num=t_num)
            chi1_B = PNJL.chi1_B(pt.T_fm, pt.muB_fm; kwargs...)
            chi2_B = PNJL.chi2_B(pt.T_fm, pt.muB_fm; kwargs...)
            chi3_B = PNJL.chi3_B(pt.T_fm, pt.muB_fm; kwargs...)
            chi4_B = PNJL.chi4_B(pt.T_fm, pt.muB_fm; kwargs...)
            chi1_Q = PNJL.chi1_Q(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi2_Q = PNJL.chi2_Q(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi3_Q = PNJL.chi3_Q(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi4_Q = PNJL.chi4_Q(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi1_S = PNJL.chi1_S(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi2_S = PNJL.chi2_S(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi3_S = PNJL.chi3_S(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi4_S = PNJL.chi4_S(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi11_BQ = PNJL.chi11_BQ(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi11_BS = PNJL.chi11_BS(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            chi11_QS = PNJL.chi11_QS(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm; kwargs...)
            c200 = PNJL.cumulant_BQS(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm, pt.V; orders=(2, 0, 0), kwargs...)
            c020 = PNJL.cumulant_BQS(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm, pt.V; orders=(0, 2, 0), kwargs...)
            c002 = PNJL.cumulant_BQS(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm, pt.V; orders=(0, 0, 2), kwargs...)
            c110 = PNJL.cumulant_BQS(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm, pt.V; orders=(1, 1, 0), kwargs...)
            c101 = PNJL.cumulant_BQS(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm, pt.V; orders=(1, 0, 1), kwargs...)
            c011 = PNJL.cumulant_BQS(pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm, pt.V; orders=(0, 1, 1), kwargs...)
            ssigma = PNJL.baryon_Ssigma(pt.T_fm, pt.muB_fm; kwargs...)
            kappa_sigma2 = PNJL.baryon_kappa_sigma2(pt.T_fm, pt.muB_fm; kwargs...)

            @printf(io,
                "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                pt.label,
                pt.T_fm, pt.muB_fm, pt.muQ_fm, pt.muS_fm, pt.V,
                chi1_B, chi2_B, chi3_B, chi4_B,
                chi1_Q, chi2_Q, chi3_Q, chi4_Q,
                chi1_S, chi2_S, chi3_S, chi4_S,
                chi11_BQ, chi11_BS, chi11_QS,
                c200, c020, c002, c110, c101, c011,
                ssigma, kappa_sigma2,
            )
        end
    end

    println("PNJL conserved-charge baseline written to: " * output)
    println(@sprintf("points: %d", length(CONSERVED_CHARGE_POINTS)))
end

main(ARGS)