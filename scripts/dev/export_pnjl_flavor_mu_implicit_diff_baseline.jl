#!/usr/bin/env julia

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_flavor_mu_implicit_diff_fixedpoints_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

const FLAVOR_IMPLICIT_POINTS = [
    (label="flavor_mu_symmetric", T_fm=0.50, mu_vec=SVector(0.20, 0.20, 0.20)),
    (label="flavor_mu_asymmetric", T_fm=0.52, mu_vec=SVector(0.18, 0.12, 0.06)),
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
            println("Usage: julia --project=. scripts/dev/export_pnjl_flavor_mu_implicit_diff_baseline.jl [--output <path>]")
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

    p_num = 24
    t_num = 6
    xi = 0.0

    open(output, "w") do io
        println(io, "label,T_fm,mu_u_fm,mu_d_fm,mu_s_fm,x1,x2,x3,x4,x5,dx_dT_1,dx_dT_2,dx_dT_3,dx_dT_4,dx_dT_5,dx_dmu_u_1,dx_dmu_u_2,dx_dmu_u_3,dx_dmu_u_4,dx_dmu_u_5,dx_dmu_d_1,dx_dmu_d_2,dx_dmu_d_3,dx_dmu_d_4,dx_dmu_d_5,dx_dmu_s_1,dx_dmu_s_2,dx_dmu_s_3,dx_dmu_s_4,dx_dmu_s_5")

        for pt in FLAVOR_IMPLICIT_POINTS
            result = Models.solve_pnjl_with_flavor_mu_derivatives(pt.T_fm, pt.mu_vec; order=1, xi=xi, p_num=p_num, t_num=t_num)
            x = result.x
            dx_dT = result.dx_dT
            dx_dmu = result.dx_dmu_vec

            @printf(io,
                "%s,%.6f,%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                pt.label,
                pt.T_fm,
                pt.mu_vec[1], pt.mu_vec[2], pt.mu_vec[3],
                x[1], x[2], x[3], x[4], x[5],
                dx_dT[1], dx_dT[2], dx_dT[3], dx_dT[4], dx_dT[5],
                dx_dmu[1, 1], dx_dmu[2, 1], dx_dmu[3, 1], dx_dmu[4, 1], dx_dmu[5, 1],
                dx_dmu[1, 2], dx_dmu[2, 2], dx_dmu[3, 2], dx_dmu[4, 2], dx_dmu[5, 2],
                dx_dmu[1, 3], dx_dmu[2, 3], dx_dmu[3, 3], dx_dmu[4, 3], dx_dmu[5, 3],
            )
        end
    end

    println("PNJL flavor-mu implicit diff baseline written to: " * output)
    println(@sprintf("points: %d", length(FLAVOR_IMPLICIT_POINTS)))
end

main(ARGS)