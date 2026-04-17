#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflows", "MesonMassWorkflow.jl"))

using Main.Constants_PNJL: ħc_MeV_fm
using Main.MesonMassWorkflow: solve_gap_and_meson_point, _mixed_branch_scores

struct Opts
    out_csv::String
    T_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    meson::Symbol
    m_min::Float64
    m_max::Float64
    m_step::Float64
    g_min::Float64
    g_max::Float64
    g_step::Float64
    p_num::Int
    t_num::Int
    max_iter::Int
end

function _usage()
    println("Usage: julia --project=. scripts/analysis/relaxtime/export_direct_eval_scatter.jl [options]")
    println("Options: --out --T --muB --xi --meson --m-min --m-max --m-step --g-min --g-max --g-step --p-num --t-num --max-iter")
end

function _parse(args::Vector{String})
    out_csv = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "mott_phase", "residual_surface", "edge255_xi-0p3_etaplus", "direct_eval_scatter.csv")
    T_MeV = 255.0
    muB_MeV = 0.0
    xi = -0.3
    meson = :eta_prime
    m_min = 2.2
    m_max = 3.3
    m_step = 0.005
    g_min = 0.0
    g_max = 3.0
    g_step = 0.005
    p_num = 8
    t_num = 4
    max_iter = 25

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $arg"))
            i += 1
            return args[i]
        end
        if arg == "--out"
            out_csv = require_value()
        elseif arg == "--T"
            T_MeV = parse(Float64, require_value())
        elseif arg == "--muB"
            muB_MeV = parse(Float64, require_value())
        elseif arg == "--xi"
            xi = parse(Float64, require_value())
        elseif arg == "--meson"
            meson = Symbol(require_value())
        elseif arg == "--m-min"
            m_min = parse(Float64, require_value())
        elseif arg == "--m-max"
            m_max = parse(Float64, require_value())
        elseif arg == "--m-step"
            m_step = parse(Float64, require_value())
        elseif arg == "--g-min"
            g_min = parse(Float64, require_value())
        elseif arg == "--g-max"
            g_max = parse(Float64, require_value())
        elseif arg == "--g-step"
            g_step = parse(Float64, require_value())
        elseif arg == "--p-num"
            p_num = parse(Int, require_value())
        elseif arg == "--t-num"
            t_num = parse(Int, require_value())
        elseif arg == "--max-iter"
            max_iter = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            _usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end

    return Opts(out_csv, T_MeV, muB_MeV, xi, meson, m_min, m_max, m_step, g_min, g_max, g_step, p_num, t_num, max_iter)
end

function main()
    opts = _parse(ARGS)
    T_fm = opts.T_MeV / ħc_MeV_fm
    mu_fm = (opts.muB_MeV / ħc_MeV_fm) / 3.0

    base = solve_gap_and_meson_point(
        T_fm,
        mu_fm;
        xi=opts.xi,
        mesons=(opts.meson,),
        p_num=opts.p_num,
        t_num=opts.t_num,
        solver_kwargs=(iterations=opts.max_iter,),
        mass_kwargs=(iterations=opts.max_iter,),
    )

    m_vals = collect(opts.m_min:opts.m_step:opts.m_max)
    g_vals = collect(opts.g_min:opts.g_step:opts.g_max)

    mkpath(dirname(opts.out_csv))
    open(opts.out_csv, "w") do io
        println(io, "m_fm,g_fm,m_MeV,g_MeV,residual_norm,log10_residual_norm")
        for g in g_vals
            for m in m_vals
                splus, sminus = _mixed_branch_scores(opts.meson, m, g, base.quark_params, base.thermo_params, 0.0)
                n = if opts.meson == :eta || opts.meson == :sigma
                    sminus
                else
                    splus
                end
                l = n > 0 ? log10(n) : -Inf
                println(io, string(m, ",", g, ",", m * ħc_MeV_fm, ",", g * ħc_MeV_fm, ",", n, ",", l))
            end
        end
    end

    println("Wrote direct-eval scatter CSV: ", opts.out_csv)
end

main()
