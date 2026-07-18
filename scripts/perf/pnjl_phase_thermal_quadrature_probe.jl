#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using BenchmarkTools
using Printf
using StaticArrays

include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
using .Models

function _parse_values(args, prefix, default)
    match = findfirst(arg -> startswith(arg, prefix), args)
    match === nothing && return default
    raw = split(args[match], "="; limit=2)[2]
    return parse.(Float64, filter(!isempty, strip.(split(raw, ','))))
end

function _parse_int(args, prefix, default)
    match = findfirst(arg -> startswith(arg, prefix), args)
    match === nothing && return default
    return parse(Int, split(args[match], "="; limit=2)[2])
end

function main(args=ARGS)
    temperatures_MeV = _parse_values(args, "--temperatures=", [1.0, 5.0, 20.0, 60.0, 150.0, 240.0])
    xi_values = _parse_values(args, "--xi=", [-0.5, 0.0, 0.5])
    samples = _parse_int(args, "--samples=", 10)

    model = Models.PNJLModel()
    state = Models.MeanFieldState(SVector(-1.8, -1.8, -2.2); Phi=0.1, PhiBar=0.1)
    x_state = Models.state_vector(state)
    mu_vec = SVector(1.8, 1.8, 1.8)

    println("T_MeV,xi,policy,p_num,t_num,value,min_time_ms,alloc_bytes")
    for T_MeV in temperatures_MeV, xi in xi_values
        T_fm = T_MeV / 197.327
        for (policy, p_num, t_num) in (
            (:tensor_gauss, 24, 8),
            (:tensor_gauss, 32, 10),
            (:rs_reduced_adaptive, 0, 0),
        )
            kwargs = if policy === :tensor_gauss
                (; p_num=p_num, t_num=t_num, xi=xi, thermo_quadrature_policy=policy)
            else
                (;
                    xi=xi,
                    thermo_quadrature_policy=policy,
                    thermo_quadrature_rtol=1e-8,
                    thermo_quadrature_atol=1e-10,
                )
            end
            value = Models.omega(model, x_state, T_fm, mu_vec; kwargs...)
            evaluate = () -> Models.omega(model, x_state, T_fm, mu_vec; kwargs...)
            trial = @benchmark $(evaluate)() samples=samples evals=1
            estimate = minimum(trial)
            @printf(
                "%.8g,%.8g,%s,%d,%d,%.16g,%.8g,%d\n",
                T_MeV,
                xi,
                String(policy),
                p_num,
                t_num,
                value,
                estimate.time / 1e6,
                estimate.memory,
            )
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
