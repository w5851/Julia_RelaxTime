#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "phase", "baseline_phase_pipeline_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

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
            println("Usage: julia --project=. scripts/dev/export_phase_pipeline_regression_baseline.jl [--output <path>]")
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

    tmp = mktempdir()
    result = Models.run_phase_pipeline(
        :PNJL;
        mode=:research,
        T_grid=[120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0],
        rho_grid=collect(0.1:0.1:3.0),
        xi=0.0,
        output_dir=tmp,
        profile=:regression,
        solver_backend=:legacy,
        reverse_rho=true,
        seed_policy=:hybrid_continuity,
        p_num=8,
        t_num=4,
        iterations=40,
        compute_crossover=true,
        crossover_method=:inflection,
        crossover_variable=:phi_u,
        crossover_n_mu=8,
        cep_strategy=:interpolate,
        promote_reference=false,
    )

    open(output, "w") do io
        println(io, "kind,idx,flag,T_MeV,mu_MeV,rho_a,rho_b,aux")
        @printf(io, "cep,0,%d,%.16e,%.16e,%.16e,%.16e,%.16e\n",
            Int(result.cep.found),
            result.cep.T_cep_MeV,
            result.cep.mu_cep_MeV,
            NaN, NaN, NaN)

        for (idx, row) in enumerate(result.first_order_boundary)
            @printf(io, "boundary,%d,%d,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                idx, Int(row.converged), row.T_MeV, row.mu_transition_MeV, row.rho_hadron, row.rho_quark, row.area_residual)
        end

        for (idx, row) in enumerate(result.crossover_line)
            @printf(io, "crossover,%d,%d,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                idx, Int(row.converged), row.T_crossover_MeV, row.mu_MeV, row.rho, row.derivative, NaN)
        end
    end

    println("Phase pipeline regression baseline written to: " * output)
end

main(ARGS)
