#!/usr/bin/env julia

raw"""
按散射过程扫描总截面 σ(s)（固定 T, μ_B, ξ），并自动计算每个过程的物理阈值。

关注点：特别输出阈值附近（s→s_th^+）的 σ(s) 行为，用于验证你在文档中分析的
\sigma(s) ~ A / sqrt(s - s0) 等可积奇异/光滑起步特征。

特点：
- 先解 PNJL 平衡（得到 m_u/m_d/m_s 与 Φ, Φbar），保证与弛豫时间计算一致
- 对每个 process 自动计算 s_th = max((m_i+m_j)^2, (m_c+m_d)^2)
- 在 √s 变量下做“阈值附近对数加密”的采样：δ = √s-√s_th 用 logspace
- 可选再做一段线性 √s 扫描（更远离阈值）

输出：CSV（每行一个 s 点），包含：process, s_th, √s_th, √s, s, σ, 1/σ^2 等。

示例：
  julia --project=. scripts/relaxtime/scan_cross_section_vs_s_by_process.jl \
    --T-MeV 150 --muB-MeV 800 --xi 0 \
    --processes ssbar_to_uubar,uubar_to_ssbar \
    --delta-sqrt-mev-min 1e-3 --delta-sqrt-mev-max 60 --n-threshold 240 \
    --n-points 64 \
    --out data/processed/results/relaxtime/xs_vs_s_threshold_T150_muB800.csv
"""

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .PNJL: solve, FixedMu, calculate_mass_vec, HADRON_SEED_5
using .PNJL: DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT
using .PNJL.Integrals: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: total_cross_section
const ASR = AverageScatteringRate

struct Options
    T_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    max_iter::Int
    n_points::Int
    processes::Vector{Symbol}
    # threshold sampling in δ=√s-√s_th [MeV]
    delta_sqrt_mev_min::Float64
    delta_sqrt_mev_max::Float64
    n_threshold::Int
    # optional linear tail in √s [MeV]
    sqrt_tail_max_mev::Float64
    n_tail::Int
    out::String
    overwrite::Bool
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/scan_cross_section_vs_s_by_process.jl [options]\n")
    println("Options:")
    println("  --T-MeV <float>              Temperature (default 150)")
    println("  --muB-MeV <float>            Baryon chemical potential μ_B (default 800)")
    println("  --xi <float>                 Momentum anisotropy ξ (default 0)")
    println("  --p-num <int>                Gap-solver momentum nodes (default 12)")
    println("  --t-num <int>                Gap-solver angle nodes (default 6)")
    println("  --max-iter <int>             Gap-solver max iterations (default 40)")
    println("  --n-points <int>             t-integration Gauss points for σ(s) (default 64)")
    println("  --processes p1,p2,...         Process list (default: common 11 processes)")
    println("  --delta-sqrt-mev-min <float>  Threshold window min δ=√s-√s_th in MeV (default 1e-3)")
    println("  --delta-sqrt-mev-max <float>  Threshold window max δ in MeV (default 60)")
    println("  --n-threshold <int>           #points in threshold window (default 240)")
    println("  --sqrt-tail-max-mev <float>   Optional extra scan to √s_th+tail (MeV, default 0=off)")
    println("  --n-tail <int>                #points in tail scan (default 0)")
    println("  --out <path>                  Output CSV path")
    println("  --overwrite                   Overwrite output")
    println("  -h, --help                    Show this help")
end

function default_processes()
    return Symbol[:uu_to_uu, :ud_to_ud, :us_to_us, :usbar_to_usbar, :uubar_to_uubar, :uubar_to_ddbar, :uubar_to_ssbar,
                  :udbar_to_udbar, :ss_to_ss, :ssbar_to_ssbar, :ssbar_to_uubar]
end

function parse_args(args::Vector{String})::Options
    opts = Dict{Symbol,Any}(
        :T_MeV => 150.0,
        :muB_MeV => 800.0,
        :xi => 0.0,
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :n_points => 64,
        :processes => default_processes(),
        :delta_sqrt_mev_min => 1e-3,
        :delta_sqrt_mev_max => 60.0,
        :n_threshold => 240,
        :sqrt_tail_max_mev => 0.0,
        :n_tail => 0,
        :out => joinpath("data", "processed", "results", "relaxtime", "xs_vs_s_by_process.csv"),
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
        elseif arg == "--n-points"
            opts[:n_points] = parse(Int, require_value())
        elseif arg == "--processes"
            raw = split(require_value(), ',')
            vals = Symbol[Symbol(strip(x)) for x in raw if !isempty(strip(x))]
            isempty(vals) && error("empty processes")
            opts[:processes] = vals
        elseif arg == "--delta-sqrt-mev-min"
            opts[:delta_sqrt_mev_min] = parse(Float64, require_value())
        elseif arg == "--delta-sqrt-mev-max"
            opts[:delta_sqrt_mev_max] = parse(Float64, require_value())
        elseif arg == "--n-threshold"
            opts[:n_threshold] = parse(Int, require_value())
        elseif arg == "--sqrt-tail-max-mev"
            opts[:sqrt_tail_max_mev] = parse(Float64, require_value())
        elseif arg == "--n-tail"
            opts[:n_tail] = parse(Int, require_value())
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
        Int(opts[:n_points]),
        Vector{Symbol}(opts[:processes]),
        Float64(opts[:delta_sqrt_mev_min]),
        Float64(opts[:delta_sqrt_mev_max]),
        Int(opts[:n_threshold]),
        Float64(opts[:sqrt_tail_max_mev]),
        Int(opts[:n_tail]),
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
    return (s_th=s_th, mi=mi, mj=mj, mc=mc, md=md, i=i, j=j, c=c, d=d)
end

function safe_sigma(process::Symbol, s::Float64, params; n_points::Int)
    try
        σ = total_cross_section(process, s, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_points)
        return isfinite(σ) ? σ : NaN
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
        # minimal metadata (comment header)
        println(io, "# schema: xs_vs_s_by_process_v1")
        println(io, "# T_MeV: $(opts.T_MeV)")
        println(io, "# muB_MeV: $(opts.muB_MeV)")
        println(io, "# xi: $(opts.xi)")
        println(io, "# n_points: $(opts.n_points)")
        println(io, "# delta_sqrt_mev_min: $(opts.delta_sqrt_mev_min)")
        println(io, "# delta_sqrt_mev_max: $(opts.delta_sqrt_mev_max)")
        println(io, "# n_threshold: $(opts.n_threshold)")
        println(io, "# sqrt_tail_max_mev: $(opts.sqrt_tail_max_mev)")
        println(io, "# n_tail: $(opts.n_tail)")
        println(io, "process,T_MeV,muB_MeV,xi,i,j,c,d,mi_fm,mj_fm,mc_fm,md_fm,s_th,sqrt_s_th_MeV,delta_sqrt_MeV,s,sqrt_s_MeV,sigma,inv_sigma2")

        for process in opts.processes
            th = process_threshold(process, params.quark_params)
            sqrt_s_th = sqrt(th.s_th)
            sqrt_s_th_mev = sqrt_s_th * ħc_MeV_fm

            # threshold window: δ in MeV (log spaced)
            δ_min = max(opts.delta_sqrt_mev_min, 1e-12)
            δ_max = max(opts.delta_sqrt_mev_max, δ_min * 10.0)
            δs = logspace(δ_min, δ_max, opts.n_threshold)

            # optional tail: linear in √s up to √s_th + sqrt_tail_max_mev
            tail = Float64[]
            if opts.n_tail > 0 && opts.sqrt_tail_max_mev > 0
                tail = collect(range(δ_max, opts.sqrt_tail_max_mev; length=opts.n_tail))
            end

            all_δ = vcat(δs, tail)

            @printf("process=%s  sqrt_s_th=%.3f MeV  points=%d\n", string(process), sqrt_s_th_mev, length(all_δ))

            for δ_mev in all_δ
                sqrt_s_mev = sqrt_s_th_mev + δ_mev
                s = (sqrt_s_mev / ħc_MeV_fm)^2
                σ = safe_sigma(process, s, params; n_points=opts.n_points)
                invσ2 = (isfinite(σ) && σ > 0) ? 1.0 / (σ * σ) : NaN

                row = (
                    string(process),
                    opts.T_MeV, opts.muB_MeV, opts.xi,
                    string(th.i), string(th.j), string(th.c), string(th.d),
                    th.mi, th.mj, th.mc, th.md,
                    th.s_th,
                    sqrt_s_th_mev,
                    δ_mev,
                    s,
                    sqrt_s_mev,
                    σ,
                    invσ2,
                )
                println(io, join(row, ','))
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
