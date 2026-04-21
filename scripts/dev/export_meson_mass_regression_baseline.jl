#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_meson_mass_fixedpoints_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

const POINTS = [
    (label="meson_mu0_T150", T_fm=0.15, mu_fm=0.0, xi=0.0),
    (label="meson_mu015_T150", T_fm=0.15, mu_fm=0.15, xi=0.0),
    (label="meson_mu0_T180", T_fm=0.18, mu_fm=0.0, xi=0.0),
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
            println("Usage: julia --project=. scripts/dev/export_meson_mass_regression_baseline.jl [--output <path>]")
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

    Models.meson_workflow_module()
    mesons = (:pi, :K)

    open(output, "w") do io
        println(io, "label,T_fm,mu_fm,xi,mu_mass,md_mass,ms_mass,pi_converged,pi_mass,pi_threshold,pi_gap,K_converged,K_mass,K_threshold,K_gap")
        for pt in POINTS
            res = MesonMassWorkflow.solve_gap_and_meson_point(
                pt.T_fm,
                pt.mu_fm;
                xi=pt.xi,
                mesons=mesons,
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=30,),
                mass_kwargs=(iterations=30,),
            )
            rpi = res.meson_results[:pi]
            rk = res.meson_results[:K]
            @printf(io, "%s,%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%d,%.16e,%.16e,%.16e,%d,%.16e,%.16e,%.16e\n",
                pt.label, pt.T_fm, pt.mu_fm, pt.xi,
                res.equilibrium.masses[1], res.equilibrium.masses[2], res.equilibrium.masses[3],
                Int(rpi.converged), rpi.mass, rpi.threshold, rpi.gap,
                Int(rk.converged), rk.mass, rk.threshold, rk.gap)
        end
    end

    println("Meson mass regression baseline written to: " * output)
end

main(ARGS)
