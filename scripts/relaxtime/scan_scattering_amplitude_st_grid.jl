#!/usr/bin/env julia

raw"""
按散射过程在物理允许的 (s,t) 区域内扫描散射矩阵元 |M|^2。

核心要求：
1) 先确定每个过程的阈值 s_th 以及 s 扫描范围；
2) 对每个 s 计算该过程的物理 t 边界 (t_min, t_max)，并只在 [t_min, t_max] 内采样 t。

本脚本提供一种常用的“阈值附近加密”的 s 网格：
  δ = √s - √s_th  (单位 MeV)
  δ 用 logspace 采样，覆盖 [delta_sqrt_mev_min, delta_sqrt_mev_max]
  √s = √s_th + δ

对每个 s 的 t 网格：
  在线性区间 [t_min, t_max] 上等距采样 n_t 个点。

输出 CSV（每行一个 (s,t) 点）：
  process,T_MeV,muB_MeV,xi,s_th,sqrt_s_th_MeV,delta_sqrt_MeV,s,sqrt_s_MeV,t_min,t_max,t,M2

示例：
  julia --threads 8 --project=. scripts/relaxtime/scan_scattering_amplitude_st_grid.jl \
    --T-MeV 150 --muB-MeV 800 --xi 0 \
    --processes ssbar_to_uubar,uubar_to_ssbar,ud_to_ud,udbar_to_udbar \
    --delta-sqrt-mev-min 1e-3 --delta-sqrt-mev-max 60 --n-s 120 --n-t 60 \
    --out data/processed/results/relaxtime/M2_st_grid_T150_muB800.csv --overwrite
"""

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "ScatteringAmplitude.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .PNJL: solve, FixedMu, calculate_mass_vec, HADRON_SEED_5
using .PNJL: DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT
using .PNJL.Integrals: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .ScatteringAmplitude: scattering_amplitude_squared
using .TotalCrossSection: calculate_t_bounds
const ASR = AverageScatteringRate

struct Options
    T_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    max_iter::Int
    processes::Vector{Symbol}
    delta_sqrt_mev_min::Float64
    delta_sqrt_mev_max::Float64
    n_s::Int
    n_t::Int
    out::String
    overwrite::Bool
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/scan_scattering_amplitude_st_grid.jl [options]\n")
    println("Options:")
    println("  --T-MeV <float>              Temperature (default 150)")
    println("  --muB-MeV <float>            Baryon chemical potential μ_B (default 800)")
    println("  --xi <float>                 Momentum anisotropy ξ (default 0)")
    println("  --p-num <int>                Gap-solver momentum nodes (default 12)")
    println("  --t-num <int>                Gap-solver angle nodes (default 6)")
    println("  --max-iter <int>             Gap-solver max iterations (default 40)")
    println("  --processes p1,p2,...         Process list")
    println("  --delta-sqrt-mev-min <float>  Threshold window min δ=√s-√s_th in MeV (default 1e-3)")
    println("  --delta-sqrt-mev-max <float>  Threshold window max δ in MeV (default 60)")
    println("  --n-s <int>                   #s points in threshold window (default 240)")
    println("  --n-t <int>                   #t points per s (default 60)")
    println("  --out <path>                  Output CSV path")
    println("  --overwrite                   Overwrite output")
    println("  -h, --help                    Show this help")
end

function default_processes()
    return Symbol[:ssbar_to_uubar, :uubar_to_ssbar, :ud_to_ud, :us_to_us, :udbar_to_udbar, :uu_to_uu, :uubar_to_uubar]
end

function parse_args(args::Vector{String})::Options
    opts = Dict{Symbol,Any}(
        :T_MeV => 150.0,
        :muB_MeV => 800.0,
        :xi => 0.0,
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :processes => default_processes(),
        :delta_sqrt_mev_min => 1e-3,
        :delta_sqrt_mev_max => 60.0,
        :n_s => 240,
        :n_t => 60,
        :out => joinpath("data", "processed", "results", "relaxtime", "M2_st_grid.csv"),
        :overwrite => false,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            val = args[i + 1]
            i += 1
            return val
        end

        if arg == "--T-MeV"
            opts[:T_MeV] = parse(Float64, require_value())
        elseif arg == "--muB-MeV"
            opts[:muB_MeV] = parse(Float64, require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, require_value())
        elseif arg == "--processes"
            raw = split(require_value(), ',')
            vals = Symbol[Symbol(strip(x)) for x in raw if !isempty(strip(x))]
            isempty(vals) && error("empty processes")
            opts[:processes] = vals
        elseif arg == "--delta-sqrt-mev-min"
            opts[:delta_sqrt_mev_min] = parse(Float64, require_value())
        elseif arg == "--delta-sqrt-mev-max"
            opts[:delta_sqrt_mev_max] = parse(Float64, require_value())
        elseif arg == "--n-s"
            opts[:n_s] = parse(Int, require_value())
        elseif arg == "--n-t"
            opts[:n_t] = parse(Int, require_value())
        elseif arg == "--out"
            opts[:out] = require_value()
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg in ("-h", "--help")
            print_usage(); exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return Options(
        Float64(opts[:T_MeV]),
        Float64(opts[:muB_MeV]),
        Float64(opts[:xi]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Vector{Symbol}(opts[:processes]),
        Float64(opts[:delta_sqrt_mev_min]),
        Float64(opts[:delta_sqrt_mev_max]),
        Int(opts[:n_s]),
        Int(opts[:n_t]),
        String(opts[:out]),
        Bool(opts[:overwrite]),
    )
end

@inline function logspace(a::Float64, b::Float64, n::Int)
    n <= 1 && return Float64[b]
    a > 0 || error("logspace requires a>0")
    b > a || error("logspace requires b>a")
    la = log10(a)
    lb = log10(b)
    return [10.0^(la + (lb - la) * (i - 1) / (n - 1)) for i in 1:n]
end

@inline function linspace(a::Float64, b::Float64, n::Int)
    n <= 1 && return Float64[0.5 * (a + b)]
    return [a + (b - a) * (i - 1) / (n - 1) for i in 1:n]
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function build_equilibrium_params(opts::Options)
    T_fm = opts.T_MeV / ħc_MeV_fm
    muq_mev = opts.muB_MeV / 3.0
    muq_fm = muq_mev / ħc_MeV_fm

    base = solve(FixedMu(), T_fm, muq_fm;
        xi=opts.xi,
        p_num=opts.p_num,
        t_num=opts.t_num,
    )
    Bool(base.converged) || error("equilibrium not converged: iters=$(base.iterations) residual=$(base.residual_norm)")

    Φ = Float64(base.x_state[4])
    Φbar = Float64(base.x_state[5])
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

    A_u = A(masses.u, muq_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(masses.s, muq_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    quark_params = (m=masses, μ=(u=Float64(muq_fm), d=Float64(muq_fm), s=Float64(muq_fm)), A=(u=A_u, d=A_u, s=A_s))
    thermo_params = (T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(opts.xi))

    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

function process_threshold(process::Symbol, quark_params)
    i, j, c, d = ASR.parse_particles_from_process(process)
    mi = ASR.get_mass(i, quark_params)
    mj = ASR.get_mass(j, quark_params)
    mc = ASR.get_mass(c, quark_params)
    md = ASR.get_mass(d, quark_params)
    s_in = (mi + mj)^2
    s_out = (mc + md)^2
    s_th = max(s_in, s_out)
    return (s_th=s_th, mi=mi, mj=mj, mc=mc, md=md)
end

function safe_M2(process::Symbol, s::Float64, t::Float64, params)
    try
        v = scattering_amplitude_squared(process, s, t, params.quark_params, params.thermo_params, params.K_coeffs)
        return isfinite(v) ? v : NaN
    catch
        return NaN
    end
end

function run_scan(opts::Options)
    ensure_parent_dir(opts.out)
    if isfile(opts.out) && opts.overwrite
        rm(opts.out)
    elseif isfile(opts.out) && !opts.overwrite
        error("output exists; use --overwrite: $(opts.out)")
    end

    params = build_equilibrium_params(opts)

    io = open(opts.out, "w")
    try
        println(io, "# schema: M2_st_grid_v1")
        println(io, "# T_MeV: $(opts.T_MeV)")
        println(io, "# muB_MeV: $(opts.muB_MeV)")
        println(io, "# xi: $(opts.xi)")
        println(io, "# delta_sqrt_mev_min: $(opts.delta_sqrt_mev_min)")
        println(io, "# delta_sqrt_mev_max: $(opts.delta_sqrt_mev_max)")
        println(io, "# n_s: $(opts.n_s)")
        println(io, "# n_t: $(opts.n_t)")
        println(io, "process,T_MeV,muB_MeV,xi,s_th,sqrt_s_th_MeV,delta_sqrt_MeV,s,sqrt_s_MeV,t_min,t_max,t,M2")

        δs = logspace(max(opts.delta_sqrt_mev_min, 1e-12), max(opts.delta_sqrt_mev_max, opts.delta_sqrt_mev_min * 10), opts.n_s)

        for process in opts.processes
            th = process_threshold(process, params.quark_params)
            sqrt_s_th = sqrt(th.s_th)
            sqrt_s_th_mev = sqrt_s_th * ħc_MeV_fm

            @printf("process=%s  sqrt_s_th=%.3f MeV  s_points=%d  t_points=%d\n", string(process), sqrt_s_th_mev, length(δs), opts.n_t)

            for δ_mev in δs
                sqrt_s_mev = sqrt_s_th_mev + δ_mev
                s = (sqrt_s_mev / ħc_MeV_fm)^2

                tb = try
                    calculate_t_bounds(Float64(s), Float64(th.mi), Float64(th.mj), Float64(th.mc), Float64(th.md))
                catch
                    (t_min=NaN, t_max=NaN)
                end

                if !(isfinite(tb.t_min) && isfinite(tb.t_max))
                    row = (string(process), opts.T_MeV, opts.muB_MeV, opts.xi, th.s_th, sqrt_s_th_mev, δ_mev, s, sqrt_s_mev, tb.t_min, tb.t_max, NaN, NaN)
                    println(io, join(row, ','))
                    continue
                end

                tmin = Float64(tb.t_min)
                tmax = Float64(tb.t_max)
                if tmin > tmax
                    tmin, tmax = tmax, tmin
                end

                ts = linspace(tmin, tmax, opts.n_t)
                for t in ts
                    M2 = safe_M2(process, s, t, params)
                    row = (
                        string(process),
                        opts.T_MeV, opts.muB_MeV, opts.xi,
                        th.s_th,
                        sqrt_s_th_mev,
                        δ_mev,
                        s,
                        sqrt_s_mev,
                        tmin,
                        tmax,
                        t,
                        M2,
                    )
                    println(io, join(row, ','))
                end
            end
            flush(io)
        end
    finally
        close(io)
    end

    @printf("Done. Wrote %s\n", opts.out)
end

if abspath(PROGRAM_FILE) == @__FILE__
    opts = parse_args(ARGS)
    run_scan(opts)
end
