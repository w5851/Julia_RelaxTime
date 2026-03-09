#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_gap_fixedpoints_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

const GAP_POINTS = [
    (label="hadron_mu0_lowT", T_MeV=50.0, mu_MeV=0.0),
    (label="hadron_mu0_midT", T_MeV=100.0, mu_MeV=0.0),
    (label="crossover_mu0", T_MeV=150.0, mu_MeV=0.0),
    (label="transition_mu0", T_MeV=180.0, mu_MeV=0.0),
    (label="quark_mu0_highT", T_MeV=220.0, mu_MeV=0.0),
    (label="hadron_mu100_lowT", T_MeV=50.0, mu_MeV=100.0),
    (label="hadron_mu100_midT", T_MeV=100.0, mu_MeV=100.0),
    (label="crossover_mu100", T_MeV=150.0, mu_MeV=100.0),
    (label="transition_mu100", T_MeV=180.0, mu_MeV=100.0),
    (label="dense_mu200_midT", T_MeV=100.0, mu_MeV=200.0),
    (label="dense_mu200_highT", T_MeV=150.0, mu_MeV=200.0),
    (label="dense_mu300_midT", T_MeV=100.0, mu_MeV=300.0),
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
            println("Usage: julia --project=. scripts/dev/export_pnjl_gap_fixedpoint_baseline.jl [--output <path>]")
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

    model = Models.create_model(:PNJL; profile="default", physics_profile="default")
    p_num = 12
    t_num = 4
    xi = 0.0
    residual_norm_max = 1e-5

    open(output, "w") do io
        println(io, "label,T_MeV,mu_MeV,phi_u,phi_d,phi_s,Phi,PhiBar,omega,max_abs_residual")

        for pt in GAP_POINTS
            T_fm = pt.T_MeV / Constants_PNJL.ħc_MeV_fm
            mu_fm = pt.mu_MeV / Constants_PNJL.ħc_MeV_fm
            st = Models.solve_gap(model, T_fm, mu_fm;
                p_num=p_num,
                t_num=t_num,
                xi=xi,
                residual_norm_max=residual_norm_max,
            )
            x = Models.state_vector(st)
            residual = Models.gap_residual(model, st, T_fm, mu_fm; p_num=p_num, t_num=t_num, xi=xi)
            omega = Models.omega(model, st, T_fm, mu_fm; p_num=p_num, t_num=t_num, xi=xi)
            max_abs_residual = maximum(abs.(residual))

            @printf(io,
                "%s,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                pt.label,
                pt.T_MeV,
                pt.mu_MeV,
                x[1],
                x[2],
                x[3],
                x[4],
                x[5],
                omega,
                max_abs_residual,
            )
        end
    end

    println("PNJL gap fixedpoint baseline written to: " * output)
    println(@sprintf("points: %d", length(GAP_POINTS)))
end

main(ARGS)