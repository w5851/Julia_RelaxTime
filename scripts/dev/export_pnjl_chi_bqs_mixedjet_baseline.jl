#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_chi_bqs_mixedjet_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

const MIXEDJET_POINTS = [
    (label="bqs_asymmetric_p4_t2", T_fm=0.57, muB_fm=0.18, muQ_fm=0.05, muS_fm=0.02, xi=0.0, p_num=4, t_num=2, orders=(2, 1, 0)),
    (label="bqs_asymmetric_p4_t2", T_fm=0.57, muB_fm=0.18, muQ_fm=0.05, muS_fm=0.02, xi=0.0, p_num=4, t_num=2, orders=(1, 1, 1)),
    (label="bqs_asymmetric_p4_t2", T_fm=0.57, muB_fm=0.18, muQ_fm=0.05, muS_fm=0.02, xi=0.0, p_num=4, t_num=2, orders=(2, 1, 1)),
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
            println("Usage: julia --project=. scripts/dev/export_pnjl_chi_bqs_mixedjet_baseline.jl [--output <path>]")
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

    open(output, "w") do io
        println(io, "label,T_fm,muB_fm,muQ_fm,muS_fm,xi,p_num,t_num,nB,nQ,nS,chi_BQS_mixedjet")
        for pt in MIXEDJET_POINTS
            chi = PNJL.chi_BQS(
                pt.T_fm,
                pt.muB_fm,
                pt.muQ_fm,
                pt.muS_fm;
                orders=pt.orders,
                xi=pt.xi,
                p_num=pt.p_num,
                t_num=pt.t_num,
                derivative_backend=:mixedjet,
            )
            @printf(
                io,
                "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%.16e\n",
                pt.label,
                pt.T_fm,
                pt.muB_fm,
                pt.muQ_fm,
                pt.muS_fm,
                pt.xi,
                pt.p_num,
                pt.t_num,
                pt.orders[1],
                pt.orders[2],
                pt.orders[3],
                chi,
            )
        end
    end

    println("PNJL chi_BQS mixedjet baseline written to: " * output)
    println(@sprintf("rows: %d", length(MIXEDJET_POINTS)))
end

main(ARGS)
