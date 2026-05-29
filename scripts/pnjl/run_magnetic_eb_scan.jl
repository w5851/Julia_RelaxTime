#!/usr/bin/env julia

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()

using .Models: calculate_magnetic_number_densities,
               calculate_magnetic_omega_components,
               calculate_magnetic_rho,
               default_magnetic_config,
               magnetic_nmax_convergence_report

const PNJL = Models.pnjl_module()

function main(; T_MeV::Float64=150.0, mu_MeV::Float64=60.0, eB_min_MeV2::Float64=5.0e3, eB_max_MeV2::Float64=4.0e4, points::Int=8, output::String=joinpath(PROJECT_ROOT, "data", "outputs", "results", "pnjl_magnetic", "scan", "pnjl_magnetic_eb_scan.csv"))
    points >= 2 || error("points must be >= 2")
    mkpath(dirname(output))

    T_fm = T_MeV / Constants_PNJL.ħc_MeV_fm
    mu_fm = mu_MeV / Constants_PNJL.ħc_MeV_fm
    mu_vec = SVector{3, Float64}(mu_fm, mu_fm, mu_fm)
    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)

    eB_values = range(eB_min_MeV2, eB_max_MeV2; length=points)

    open(output, "w") do io
        println(io, "T_MeV,mu_MeV,eB_MeV2,omega_fm4,pressure_fm4,rho_u,rho_d,rho_s,baryon_rho0,n_max,G_B,nmax_rel_diff,nmax_converged")
        for eB_MeV2 in eB_values
            eB_fm2 = eB_MeV2 / (Constants_PNJL.ħc_MeV_fm^2)
            conf = default_magnetic_config(eB_fm2=eB_fm2)

            comp = calculate_magnetic_omega_components(x_state, mu_vec, T_fm, conf)
            rho = calculate_magnetic_rho(x_state, mu_vec, T_fm, conf)
            nd = calculate_magnetic_number_densities(x_state, mu_vec, T_fm, conf)
            conv = magnetic_nmax_convergence_report(x_state, mu_vec, T_fm, conf; delta_n=6, rtol=5e-3)

            @printf(io, "%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%d,%.16e,%.16e,%s\n",
                T_MeV, mu_MeV, eB_MeV2,
                comp.omega, -comp.omega,
                rho[1], rho[2], rho[3], nd.baryon,
                comp.n_max, comp.G_B,
                conv.rel_diff,
                string(conv.converged),
            )
        end
    end

    println("magnetic eB scan saved to: $output")
end

function _option_value(args::Vector{String}, i::Int)
    arg = args[i]
    if occursin("=", arg)
        return split(arg, "="; limit=2)[2], i + 1
    end
    i < length(args) || error("missing value for $arg")
    return args[i + 1], i + 2
end

function _parse_cli_args(args::Vector{String})
    kwargs = Dict{Symbol, Any}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--output",) || startswith(arg, "--output=")
            value, i = _option_value(args, i)
            kwargs[:output] = value
        elseif arg in ("--T", "--T-MeV", "--T_MeV") || startswith(arg, "--T=") || startswith(arg, "--T-MeV=") || startswith(arg, "--T_MeV=")
            value, i = _option_value(args, i)
            kwargs[:T_MeV] = parse(Float64, value)
        elseif arg in ("--mu", "--mu-MeV", "--mu_MeV") || startswith(arg, "--mu=") || startswith(arg, "--mu-MeV=") || startswith(arg, "--mu_MeV=")
            value, i = _option_value(args, i)
            kwargs[:mu_MeV] = parse(Float64, value)
        elseif arg in ("--eB-min", "--eB_min", "--eB_min_MeV2") || startswith(arg, "--eB-min=") || startswith(arg, "--eB_min=") || startswith(arg, "--eB_min_MeV2=")
            value, i = _option_value(args, i)
            kwargs[:eB_min_MeV2] = parse(Float64, value)
        elseif arg in ("--eB-max", "--eB_max", "--eB_max_MeV2") || startswith(arg, "--eB-max=") || startswith(arg, "--eB_max=") || startswith(arg, "--eB_max_MeV2=")
            value, i = _option_value(args, i)
            kwargs[:eB_max_MeV2] = parse(Float64, value)
        elseif arg in ("--points",) || startswith(arg, "--points=")
            value, i = _option_value(args, i)
            kwargs[:points] = parse(Int, value)
        else
            error("unknown argument: $arg")
        end
    end
    return kwargs
end

main(; _parse_cli_args(ARGS)...)
