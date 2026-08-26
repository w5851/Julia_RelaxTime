#!/usr/bin/env julia

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_SMOKE_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_magnetic_fixedpoints_smoke_v1.csv")
const DEFAULT_NIGHTLY_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_magnetic_fixedpoints_nightly_v1.csv")

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()

using .Models: calculate_magnetic_number_densities,
               calculate_magnetic_omega_components,
               calculate_magnetic_rho,
               default_magnetic_config,
               magnetic_nmax_convergence_report

const PNJL = Models.pnjl_module()

function parse_args(args::Vector{String})
    scope = :smoke
    output = ""
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--scope"
            i == length(args) && error("missing value for --scope")
            i += 1
            scope = Symbol(lowercase(args[i]))
            scope in (:smoke, :nightly) || error("scope must be smoke or nightly")
        elseif arg == "--output"
            i == length(args) && error("missing value for --output")
            i += 1
            output = args[i]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/export_pnjl_magnetic_fixedpoint_baseline.jl [--scope smoke|nightly] [--output <path>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    if isempty(output)
        output = scope === :smoke ? DEFAULT_SMOKE_OUTPUT : DEFAULT_NIGHTLY_OUTPUT
    end
    return scope, output
end

function _points(scope::Symbol)
    if scope === :smoke
        return [
            (T=140.0, mu=0.0, eB=1.0e4, xi=0.0),
            (T=150.0, mu=60.0, eB=2.0e4, xi=0.0),
            (T=165.0, mu=120.0, eB=3.0e4, xi=0.0),
        ]
    end
    return [
        (T=130.0, mu=0.0, eB=8.0e3, xi=0.0),
        (T=140.0, mu=60.0, eB=1.5e4, xi=0.0),
        (T=150.0, mu=90.0, eB=2.0e4, xi=0.0),
        (T=160.0, mu=120.0, eB=2.5e4, xi=0.0),
        (T=170.0, mu=180.0, eB=3.0e4, xi=0.0),
        (T=180.0, mu=240.0, eB=3.5e4, xi=0.0),
    ]
end

function main(args::Vector{String})
    scope, output = parse_args(args)
    mkpath(dirname(output))

    points = _points(scope)
    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)

    open(output, "w") do io
        println(io, "scope,T_MeV,mu_MeV,eB_MeV2,xi,omega_fm4,pressure_fm4,rho_u,rho_d,rho_s,baryon_rho0,n_max,G_B,nmax_rel_diff,nmax_converged")
        for pt in points
            T_fm = pt.T / Constants_PNJL.ħc_MeV_fm
            mu_fm = pt.mu / Constants_PNJL.ħc_MeV_fm
            eB_fm2 = pt.eB / (Constants_PNJL.ħc_MeV_fm^2)
            mu_vec = SVector{3, Float64}(mu_fm, mu_fm, mu_fm)

            conf = default_magnetic_config(eB_fm2=eB_fm2)
            comp = calculate_magnetic_omega_components(x_state, mu_vec, T_fm, conf)
            rho = calculate_magnetic_rho(x_state, mu_vec, T_fm, conf)
            nd = calculate_magnetic_number_densities(x_state, mu_vec, T_fm, conf)
            conv = magnetic_nmax_convergence_report(x_state, mu_vec, T_fm, conf; delta_n=6, rtol=3e-2)

            @printf(io, "%s,%.6f,%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%d,%.16e,%.16e,%s\n",
                String(scope), pt.T, pt.mu, pt.eB, pt.xi,
                comp.omega, -comp.omega,
                rho[1], rho[2], rho[3], nd.baryon,
                comp.n_max, comp.G_B,
                conv.rel_diff,
                string(conv.converged),
            )
        end
    end

    println("magnetic baseline exported: " * output)
end

main(ARGS)
