#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
using .TransportWorkflow

const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_transport_fixedpoints_v1.csv")

function parse_args(args::Vector{String})
    output = DEFAULT_OUTPUT
    backend = :legacy
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--output"
            i == length(args) && error("missing value for --output")
            i += 1
            output = args[i]
        elseif arg == "--backend"
            i == length(args) && error("missing value for --backend")
            i += 1
            b = Symbol(lowercase(args[i]))
            b in (:legacy, :models) || error("backend must be legacy or models")
            backend = b
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/export_transport_fixedpoint_baseline.jl [--output <path>] [--backend legacy|models]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return output, backend
end

function main(args::Vector{String})
    output, backend = parse_args(args)
    mkpath(dirname(output))

    points = (
        (T=0.75, mu=0.00, xi=0.0),
        (T=0.90, mu=0.00, xi=0.0),
        (T=1.05, mu=0.00, xi=0.0),
        (T=0.90, mu=0.15, xi=0.0),
        (T=0.90, mu=0.00, xi=0.2),
    )

    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
    models_solver = Main.Models.NLsolveGapSolver(
        method=:trust_region,
        jacobian=:finite,
        xtol=1e-10,
        ftol=1e-10,
    )

    open(output, "w") do io
        println(io, "T,mu,xi,eta,sigma,zeta")
        for pt in points
            kwargs = (
                xi=pt.xi,
                thermo_backend=:models,
                tau=tau,
                compute_tau=false,
                compute_bulk=true,
                p_num=8,
                t_num=4,
                transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            )
            result = if backend == :legacy
                TransportWorkflow.solve_gap_and_transport(
                    pt.T,
                    pt.mu;
                    kwargs...,
                    solver_backend=:legacy,
                    solver_kwargs=(iterations=30,),
                )
            else
                TransportWorkflow.solve_gap_and_transport(
                    pt.T,
                    pt.mu;
                    kwargs...,
                    solver_backend=:models,
                    models_solver=models_solver,
                    models_residual_norm_max=1e-4,
                    seed_state=TransportWorkflow.HADRON_SEED_5,
                )
            end

            @printf(io, "%.6f,%.6f,%.6f,%.16e,%.16e,%.16e\n",
                pt.T, pt.mu, pt.xi,
                result.transport.eta,
                result.transport.sigma,
                result.transport.zeta,
            )
        end
    end

    println("baseline exported to: " * output)
    println("backend = " * String(backend))
end

main(ARGS)
