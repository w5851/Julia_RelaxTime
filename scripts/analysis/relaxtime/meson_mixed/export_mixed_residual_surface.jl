#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "workflow_apps", "MesonMassWorkflow.jl"))

using Main.Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using Main.MesonMassWorkflow: solve_gap_and_meson_point
using Main.MesonMass: ensure_quark_params_has_A
using Main.PolarizationAniso: polarization_with_width
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings, mixing_matrix_elements

struct Opts
    outdir::String
    T_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    meson::Symbol
    branch::Symbol
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
    println("Usage: julia --project=. scripts/analysis/relaxtime/export_mixed_residual_surface.jl [options]")
    println("Options:")
    println("  --outdir <path>      Output directory")
    println("  --T <MeV>            Temperature in MeV (default 255)")
    println("  --muB <MeV>          Baryon chemical potential in MeV (default 0)")
    println("  --xi <value>         Anisotropy parameter xi (default -0.3)")
    println("  --meson <symbol>     eta|eta_prime|sigma|sigma_prime (default eta_prime)")
    println("  --branch <symbol>    plus|minus (default inferred from meson)")
    println("  --m-min <fm^-1>      Mass grid lower bound (default 2.2)")
    println("  --m-max <fm^-1>      Mass grid upper bound (default 3.3)")
    println("  --m-step <fm^-1>     Mass grid step (default 0.01)")
    println("  --g-min <fm^-1>      Gamma grid lower bound (default 0.0)")
    println("  --g-max <fm^-1>      Gamma grid upper bound (default 3.0)")
    println("  --g-step <fm^-1>     Gamma grid step (default 0.01)")
    println("  --p-num <int>        Momentum nodes (default 8)")
    println("  --t-num <int>        Angle nodes (default 4)")
    println("  --max-iter <int>     Iteration cap (default 25)")
    println("  -h, --help           Show help")
end

function _parse(args::Vector{String})
    outdir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "mott_phase", "residual_surface", "default")
    T_MeV = 255.0
    muB_MeV = 0.0
    xi = -0.3
    meson = :eta_prime
    branch = nothing
    m_min = 2.2
    m_max = 3.3
    m_step = 0.01
    g_min = 0.0
    g_max = 3.0
    g_step = 0.01
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
        if arg == "--outdir"
            outdir = require_value()
        elseif arg == "--T"
            T_MeV = parse(Float64, require_value())
        elseif arg == "--muB"
            muB_MeV = parse(Float64, require_value())
        elseif arg == "--xi"
            xi = parse(Float64, require_value())
        elseif arg == "--meson"
            meson = Symbol(require_value())
        elseif arg == "--branch"
            branch = Symbol(require_value())
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

    if branch === nothing
        branch = (meson == :eta_prime || meson == :sigma_prime) ? :plus : :minus
    end
    branch in (:plus, :minus) || throw(ArgumentError("branch must be plus or minus"))
    meson in (:eta, :eta_prime, :sigma, :sigma_prime) || throw(ArgumentError("meson must be eta|eta_prime|sigma|sigma_prime"))

    return Opts(outdir, T_MeV, muB_MeV, xi, meson, branch, m_min, m_max, m_step, g_min, g_max, g_step, p_num, t_num, max_iter)
end

@inline function _channel_symbol(meson::Symbol)
    return (meson == :eta || meson == :eta_prime) ? :P : :S
end

function _residual_complex(meson::Symbol, branch::Symbol, mass::Float64, gamma::Float64, qp, tp)
    chan = _channel_symbol(meson)
    qpA = ensure_quark_params_has_A(qp, tp)

    Gu = calculate_G_from_A(qpA.A.u, qpA.m.u)
    Gs = calculate_G_from_A(qpA.A.s, qpA.m.s)
    Kc = calculate_effective_couplings(G_fm2, K_fm5, Gu, Gs)

    μq = Float64(qp.μ.u)
    mu = Float64(qp.m.u)
    ms = Float64(qp.m.s)

    Puu_re, Puu_im = polarization_with_width(
        chan,
        mass,
        gamma,
        0.0,
        mu,
        mu,
        μq,
        μq,
        Float64(tp.T),
        Float64(tp.Φ),
        Float64(tp.Φbar),
        Float64(tp.ξ),
        qpA.A.u,
        qpA.A.u,
        0,
    )
    Pss_re, Pss_im = polarization_with_width(
        chan,
        mass,
        gamma,
        0.0,
        ms,
        ms,
        μq,
        μq,
        Float64(tp.T),
        Float64(tp.Φ),
        Float64(tp.Φbar),
        Float64(tp.ξ),
        qpA.A.s,
        qpA.A.s,
        2,
    )

    elems = mixing_matrix_elements(ComplexF64(Puu_re, Puu_im), ComplexF64(Pss_re, Pss_im), Kc, chan)
    root = sqrt((elems.M00 - elems.M88)^2 + 4.0 * elems.M08^2)
    return branch == :plus ? (elems.M00 + elems.M88) + root : (elems.M00 + elems.M88) - root
end

function _write_surface(path::String, m_vals::Vector{Float64}, g_vals::Vector{Float64}, qp, tp, meson::Symbol, branch::Symbol)
    open(path, "w") do io
        println(io, "m_fm,g_fm,m_MeV,g_MeV,residual_re,residual_im,residual_norm,log10_residual_norm")
        for g in g_vals
            for m in m_vals
                z = _residual_complex(meson, branch, m, g, qp, tp)
                n = hypot(real(z), imag(z))
                ln = n > 0 ? log10(n) : -Inf
                println(io, string(m, ",", g, ",", m * ħc_MeV_fm, ",", g * ħc_MeV_fm, ",", real(z), ",", imag(z), ",", n, ",", ln))
            end
        end
    end
end

function _write_slice(path::String, vary::Symbol, fixed::Float64, vals::Vector{Float64}, qp, tp, meson::Symbol, branch::Symbol)
    open(path, "w") do io
        if vary == :m
            println(io, "m_fm,g_fixed_fm,m_MeV,residual_re,residual_im,residual_norm,log10_residual_norm")
            for m in vals
                z = _residual_complex(meson, branch, m, fixed, qp, tp)
                n = hypot(real(z), imag(z))
                ln = n > 0 ? log10(n) : -Inf
                println(io, string(m, ",", fixed, ",", m * ħc_MeV_fm, ",", real(z), ",", imag(z), ",", n, ",", ln))
            end
        else
            println(io, "g_fm,m_fixed_fm,g_MeV,residual_re,residual_im,residual_norm,log10_residual_norm")
            for g in vals
                z = _residual_complex(meson, branch, fixed, g, qp, tp)
                n = hypot(real(z), imag(z))
                ln = n > 0 ? log10(n) : -Inf
                println(io, string(g, ",", fixed, ",", g * ħc_MeV_fm, ",", real(z), ",", imag(z), ",", n, ",", ln))
            end
        end
    end
end

function main()
    opts = _parse(ARGS)
    mkpath(opts.outdir)

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

    center_m = Float64(base.meson_results[opts.meson].mass)
    center_g = Float64(base.meson_results[opts.meson].gamma)

    m_vals = collect(opts.m_min:opts.m_step:opts.m_max)
    g_vals = collect(opts.g_min:opts.g_step:opts.g_max)

    surface_csv = joinpath(opts.outdir, "residual_surface.csv")
    slice_m_csv = joinpath(opts.outdir, "residual_slice_mass.csv")
    slice_g_csv = joinpath(opts.outdir, "residual_slice_gamma.csv")
    meta_txt = joinpath(opts.outdir, "residual_surface_meta.txt")

    _write_surface(surface_csv, m_vals, g_vals, qp, tp, opts.meson, opts.branch)
    _write_slice(slice_m_csv, :m, center_g, m_vals, qp, tp, opts.meson, opts.branch)
    _write_slice(slice_g_csv, :g, center_m, g_vals, qp, tp, opts.meson, opts.branch)

    open(meta_txt, "w") do io
        println(io, "T_MeV=", opts.T_MeV)
        println(io, "muB_MeV=", opts.muB_MeV)
        println(io, "xi=", opts.xi)
        println(io, "meson=", opts.meson)
        println(io, "branch=", opts.branch)
        println(io, "center_mass_fm=", center_m)
        println(io, "center_gamma_fm=", center_g)
        println(io, "center_mass_MeV=", center_m * ħc_MeV_fm)
        println(io, "center_gamma_MeV=", center_g * ħc_MeV_fm)
        println(io, "grid_mass=", opts.m_min, ":", opts.m_step, ":", opts.m_max)
        println(io, "grid_gamma=", opts.g_min, ":", opts.g_step, ":", opts.g_max)
    end

    println("Wrote residual surface CSV: ", surface_csv)
    println("Wrote residual mass slice CSV: ", slice_m_csv)
    println("Wrote residual gamma slice CSV: ", slice_g_csv)
    println("Wrote residual metadata: ", meta_txt)
end

main()
