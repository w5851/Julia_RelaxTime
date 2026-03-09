#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_implicit_diff_fixedpoints_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

const IMPLICIT_POINTS = [
    (label="implicit_cold_midmu", T_fm=0.50, mu_fm=1.00),
    (label="implicit_dense_cold", T_fm=0.50, mu_fm=1.50),
    (label="implicit_transition", T_fm=0.75, mu_fm=0.75),
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
            println("Usage: julia --project=. scripts/dev/export_pnjl_implicit_diff_baseline.jl [--output <path>]")
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
        println(io, "label,T_fm,mu_fm,x1,x2,x3,x4,x5,dx_dT_1,dx_dT_2,dx_dT_3,dx_dT_4,dx_dT_5,dx_dmu_1,dx_dmu_2,dx_dmu_3,dx_dmu_4,dx_dmu_5")

        for pt in IMPLICIT_POINTS
            result = PNJL.solve_with_derivatives(pt.T_fm, pt.mu_fm; order=1, xi=xi, p_num=p_num, t_num=t_num)
            x = result.x
            dx_dT = result.dx_dT
            dx_dμ = result.dx_dμ

            @printf(io,
                "%s,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                pt.label,
                pt.T_fm,
                pt.mu_fm,
                x[1], x[2], x[3], x[4], x[5],
                dx_dT[1], dx_dT[2], dx_dT[3], dx_dT[4], dx_dT[5],
                dx_dμ[1], dx_dμ[2], dx_dμ[3], dx_dμ[4], dx_dμ[5],
            )
        end
    end

    println("PNJL implicit diff baseline written to: " * output)
    println(@sprintf("points: %d", length(IMPLICIT_POINTS)))
end

main(ARGS)