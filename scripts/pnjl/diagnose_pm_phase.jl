#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
using .Models

Base.@kwdef mutable struct PMDiagnosticCliConfig
    T_values::Vector{Float64} = [130.9]
    mu_start::Float64 = 289.5
    mu_stop::Float64 = 292.5
    mu_step::Float64 = 0.1
    xi::Float64 = 0.0
    solver_backend::Symbol = :models
    p_num::Int = 24
    t_num::Int = 8
    output_dir::String = mktempdir()
end

function _usage()
    println("用法: julia scripts/pnjl/diagnose_pm_phase.jl [options]")
    println("选项:")
    println("  --T_values=130.9,131.0      温度列表 (MeV)")
    println("  --mu_start=289.5            mu 网格起点 (MeV)")
    println("  --mu_stop=292.5             mu 网格终点 (MeV)")
    println("  --mu_step=0.1               mu 网格步长 (MeV)")
    println("  --xi=0.0                    各向异性参数")
    println("  --solver_backend=models     models|auto")
    println("  --p_num=24                  动量积分点数")
    println("  --t_num=8                   角度积分点数")
    println("  --output_dir=...            输出目录")
    println("  -h, --help                  显示帮助")
end

function _parse_float_list(raw::AbstractString)
    isempty(strip(String(raw))) && return Float64[]
    return [parse(Float64, strip(String(part))) for part in split(String(raw), ',') if !isempty(strip(String(part)))]
end

function parse_args(args)
    cfg = PMDiagnosticCliConfig()
    for arg in args
        if arg in ("-h", "--help")
            _usage()
            exit(0)
        elseif startswith(arg, "--T_values=")
            cfg.T_values = _parse_float_list(split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--mu_start=")
            cfg.mu_start = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--mu_stop=")
            cfg.mu_stop = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--mu_step=")
            cfg.mu_step = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--xi=")
            cfg.xi = parse(Float64, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--solver_backend=")
            backend = Symbol(lowercase(split(arg, "="; limit=2)[2]))
            backend == :models || throw(ArgumentError("--solver_backend only supports models in current PM diagnostic path"))
            cfg.solver_backend = backend
        elseif startswith(arg, "--p_num=")
            cfg.p_num = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--t_num=")
            cfg.t_num = parse(Int, split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--output_dir=")
            cfg.output_dir = split(arg, "="; limit=2)[2]
        else
            error("Unknown option: $arg")
        end
    end
    return cfg
end

function main(args=ARGS)
    cfg = parse_args(args)
    mu_grid = collect(cfg.mu_start:cfg.mu_step:cfg.mu_stop)

    result = Models.analyze_pm_branch_competition(
        T_values=cfg.T_values,
        mu_grid=mu_grid,
        xi=cfg.xi,
        solver_backend=cfg.solver_backend,
        p_num=cfg.p_num,
        t_num=cfg.t_num,
        output_dir=cfg.output_dir,
    )

    println("PM diagnostic complete")
    println("  output_dir = $(cfg.output_dir)")
    println("  branch_rows = $(length(result.branch_rows))")
    println("  summaries = $(length(result.temperature_summaries))")
    return result
end

main(ARGS)
