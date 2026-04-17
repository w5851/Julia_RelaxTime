#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflows", "MesonMassWorkflow.jl"))

using Main.Constants_PNJL: ħc_MeV_fm
using Main.MesonMassWorkflow: solve_gap_and_meson_point
using Main.MesonMass: solve_meson_mass

struct Opts
    out_csv::String
    T_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    m0::Float64
    g0::Float64
    method::Symbol
    p_num::Int
    t_num::Int
    max_iter::Int
    nlsolve_iter::Int
end

function _usage()
    println("Usage: julia --project=. scripts/analysis/relaxtime/export_eta_prime_solver_trace.jl [options]")
    println("Options:")
    println("  --out <path>         Output trace CSV")
    println("  --T <MeV>            Temperature in MeV (default 255)")
    println("  --muB <MeV>          Baryon chemical potential in MeV (default 0)")
    println("  --xi <value>         Anisotropy xi (default -0.3)")
    println("  --m0 <fm^-1>         Initial mass seed")
    println("  --g0 <fm^-1>         Initial gamma seed")
    println("  --method <symbol>    newton|trust_region (default trust_region)")
    println("  --p-num <int>        Momentum nodes (default 8)")
    println("  --t-num <int>        Angle nodes (default 4)")
    println("  --max-iter <int>     Workflow iteration cap (default 25)")
    println("  --nlsolve-iter <int> NLsolve iterations (default 200)")
end

function _parse(args::Vector{String})
    out_csv = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "mott_phase", "residual_surface", "edge255_xi-0p3_etaplus", "solver_trace.csv")
    T_MeV = 255.0
    muB_MeV = 0.0
    xi = -0.3
    m0 = 2.7099381339372197
    g0 = 1.1231259359255894
    method = :trust_region
    p_num = 8
    t_num = 4
    max_iter = 25
    nlsolve_iter = 200

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
        elseif arg == "--m0"
            m0 = parse(Float64, require_value())
        elseif arg == "--g0"
            g0 = parse(Float64, require_value())
        elseif arg == "--method"
            method = Symbol(require_value())
        elseif arg == "--p-num"
            p_num = parse(Int, require_value())
        elseif arg == "--t-num"
            t_num = parse(Int, require_value())
        elseif arg == "--max-iter"
            max_iter = parse(Int, require_value())
        elseif arg == "--nlsolve-iter"
            nlsolve_iter = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            _usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end

    method in (:newton, :trust_region) || throw(ArgumentError("method must be newton|trust_region"))
    return Opts(out_csv, T_MeV, muB_MeV, xi, m0, g0, method, p_num, t_num, max_iter, nlsolve_iter)
end

function main()
    opts = _parse(ARGS)

    T_fm = opts.T_MeV / ħc_MeV_fm
    mu_fm = (opts.muB_MeV / ħc_MeV_fm) / 3.0

    base = solve_gap_and_meson_point(
        T_fm,
        mu_fm;
        xi=opts.xi,
        mesons=(:eta_prime,),
        p_num=opts.p_num,
        t_num=opts.t_num,
        solver_kwargs=(iterations=opts.max_iter,),
        mass_kwargs=(iterations=opts.max_iter,),
    )

    qp = (
        m=(u=Float64(base.quark_params.m.u), d=Float64(base.quark_params.m.d), s=Float64(base.quark_params.m.s)),
        μ=(u=Float64(base.quark_params.μ.u), d=Float64(base.quark_params.μ.d), s=Float64(base.quark_params.μ.s)),
    )
    tp = (
        T=Float64(base.thermo_params.T),
        Φ=Float64(base.thermo_params.Φ),
        Φbar=Float64(base.thermo_params.Φbar),
        ξ=Float64(base.thermo_params.ξ),
    )

    result = solve_meson_mass(
        :eta_prime,
        qp,
        tp;
        k_norm=0.0,
        initial_mass=opts.m0,
        initial_gamma=opts.g0,
        method=opts.method,
        iterations=opts.nlsolve_iter,
        store_trace=true,
        extended_trace=true,
    )

    states = result.solution.trace.states

    mkpath(dirname(opts.out_csv))
    open(opts.out_csv, "w") do io
        println(io, "iter,m_fm,g_fm,residual_norm")
        for st in states
            m = st.metadata["x"][1]
            g = st.metadata["x"][2]
            rn = st.fnorm
            println(io, string(st.iteration, ",", m, ",", g, ",", rn))
        end
    end

    println("Wrote solver trace CSV: ", opts.out_csv)
end

main()
