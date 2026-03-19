#!/usr/bin/env julia

"""
T200 usbar_to_usbar 单通道缓存网格对照实验：
- 基线：adaptive=false + 原始 w0cdf 缓存
- 参考：adaptive=true
- 对照：adaptive=false + 在 s_peak 附近手工加密缓存点

输出 CSV 用于检验是否可连续收敛到 adaptive=true 结果。
"""

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))

using .Constants_PNJL: G_fm2, K_fm5, Λ_inv_fm
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS

const RT_ASR = Main.AverageScatteringRate
const TotalCS = Main.TotalCrossSection

Base.@kwdef struct Options
    scan_csv::String = joinpath("data", "outputs", "results", "relaxtime", "plan_b", "transport_vs_xi_T200_muB0.csv")
    xi_list::Vector{Float64} = Float64[-0.22, -0.20, -0.18]
    process::Symbol = :usbar_to_usbar
    sigma_grid_n::Int = 64
    tau_p_nodes::Int = 28
    tau_angle_nodes::Int = 8
    tau_phi_nodes::Int = 8
    tau_n_sigma_points::Int = 6
    local_halfwidth::Float64 = 0.02
    local_points_list::Vector{Int} = Int[0, 8, 16, 32, 64, 128, 256]
    out_csv::String = joinpath("data", "outputs", "results", "relaxtime", "scan", "_t200_usbar_cache_grid_study.csv")
end

function usage()
    println("Usage: julia --project=. scripts/analysis/relaxtime_t200_usbar_cache_grid_study.jl [options]\n")
    println("Options:")
    println("  --scan-csv <path>               输入 plan_b CSV")
    println("  --xi-list <csv-float>           xi 列表，默认 -0.22,-0.2,-0.18")
    println("  --process <symbol>              默认 usbar_to_usbar")
    println("  --sigma-grid-n <int>            缓存网格基准点数")
    println("  --tau-p-nodes <int>             τ 动量积分节点")
    println("  --tau-angle-nodes <int>         τ 角度积分节点")
    println("  --tau-phi-nodes <int>           τ 方位角积分节点")
    println("  --tau-n-sigma <int>             σ(s) t 积分点")
    println("  --local-halfwidth <float>       s_peak 局部加密半宽")
    println("  --local-points-list <csv-int>   局部加密点数列表（含 0 表示不加密）")
    println("  --out-csv <path>                输出 CSV")
    println("  -h, --help                      帮助")
end

function parse_float_list(raw::AbstractString)
    vals = Float64[]
    for p in split(raw, ',')
        s = strip(p)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && error("empty float list")
    return vals
end

function parse_int_list(raw::AbstractString)
    vals = Int[]
    for p in split(raw, ',')
        s = strip(p)
        isempty(s) && continue
        push!(vals, parse(Int, s))
    end
    isempty(vals) && error("empty int list")
    return vals
end

function parse_args(args::Vector{String})
    opts = Options()
    scan_csv = opts.scan_csv
    xi_list = opts.xi_list
    process = opts.process
    sigma_grid_n = opts.sigma_grid_n
    tau_p_nodes = opts.tau_p_nodes
    tau_angle_nodes = opts.tau_angle_nodes
    tau_phi_nodes = opts.tau_phi_nodes
    tau_n_sigma_points = opts.tau_n_sigma_points
    local_halfwidth = opts.local_halfwidth
    local_points_list = opts.local_points_list
    out_csv = opts.out_csv

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            value = args[i + 1]
            i += 1
            return value
        end

        if arg == "--scan-csv"
            scan_csv = require_value()
        elseif arg == "--xi-list"
            xi_list = parse_float_list(require_value())
        elseif arg == "--process"
            process = Symbol(require_value())
        elseif arg == "--sigma-grid-n"
            sigma_grid_n = parse(Int, require_value())
        elseif arg == "--tau-p-nodes"
            tau_p_nodes = parse(Int, require_value())
        elseif arg == "--tau-angle-nodes"
            tau_angle_nodes = parse(Int, require_value())
        elseif arg == "--tau-phi-nodes"
            tau_phi_nodes = parse(Int, require_value())
        elseif arg == "--tau-n-sigma"
            tau_n_sigma_points = parse(Int, require_value())
        elseif arg == "--local-halfwidth"
            local_halfwidth = parse(Float64, require_value())
        elseif arg == "--local-points-list"
            local_points_list = parse_int_list(require_value())
        elseif arg == "--out-csv"
            out_csv = require_value()
        elseif arg in ("-h", "--help")
            usage(); exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return Options(
        scan_csv=scan_csv,
        xi_list=sort(unique(xi_list)),
        process=process,
        sigma_grid_n=sigma_grid_n,
        tau_p_nodes=tau_p_nodes,
        tau_angle_nodes=tau_angle_nodes,
        tau_phi_nodes=tau_phi_nodes,
        tau_n_sigma_points=tau_n_sigma_points,
        local_halfwidth=local_halfwidth,
        local_points_list=sort(unique(local_points_list)),
        out_csv=out_csv,
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function parse_csv_line(line::String)
    return split(chomp(line), ',')
end

function read_comment_csv(path::AbstractString)
    rows = Vector{Dict{String, String}}()
    open(path, "r") do io
        header = String[]
        for raw in eachline(io)
            line = strip(raw)
            isempty(line) && continue
            startswith(line, "#") && continue
            if isempty(header)
                header = parse_csv_line(raw)
                continue
            end
            values = parse_csv_line(raw)
            push!(rows, Dict(header[i] => values[i] for i in eachindex(header)))
        end
    end
    return rows
end

@inline to_f(row::Dict{String, String}, key::String) = parse(Float64, row[key])

function build_K_data(T_fm::Float64, muq_fm::Float64, masses::NamedTuple, Phi::Float64, Phibar::Float64)
    A_u = A(masses.u, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(masses.s, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return (
        K_coeffs=calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s),
        A_vals=(u=A_u, d=A_u, s=A_s),
    )
end

function row_to_state(row::Dict{String, String})
    T_fm = to_f(row, "T_fm")
    muq_fm = to_f(row, "muq_fm")
    xi = to_f(row, "xi")
    Phi = to_f(row, "Phi")
    Phibar = to_f(row, "Phibar")
    masses = (u=to_f(row, "m_u"), d=to_f(row, "m_d"), s=to_f(row, "m_s"))
    ktmp = build_K_data(T_fm, muq_fm, masses, Phi, Phibar)
    quark_params = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm), A=ktmp.A_vals)
    thermo_params = (T=T_fm, Φ=Phi, Φbar=Phibar, ξ=xi)
    return (
        xi=xi,
        quark_params=quark_params,
        thermo_params=thermo_params,
        K_coeffs=ktmp.K_coeffs,
    )
end

function process_window(process::Symbol, state, sigma_cutoff::Float64)
    pi_sym, pj_sym, pc_sym, pd_sym = RT_ASR.parse_particles_from_process(process)
    mi = RT_ASR.get_mass(pi_sym, state.quark_params)
    mj = RT_ASR.get_mass(pj_sym, state.quark_params)
    mc = RT_ASR.get_mass(pc_sym, state.quark_params)
    md = RT_ASR.get_mass(pd_sym, state.quark_params)
    s_bo = max((mi + mj)^2, (mc + md)^2)
    Λ = sigma_cutoff
    s_up = min(
        (sqrt(mi^2 + Λ^2) + sqrt(mj^2 + Λ^2))^2,
        (sqrt(mc^2 + Λ^2) + sqrt(md^2 + Λ^2))^2,
    )
    return (s_bo=s_bo, s_up=s_up)
end

function eval_rate(process::Symbol, state, cache;
    tau_p_nodes::Int,
    tau_angle_nodes::Int,
    tau_phi_nodes::Int,
    tau_n_sigma_points::Int,
    sigma_grid_n::Int,
    sigma_cutoff::Float64,
)
    return RT_ASR.average_scattering_rate(
        process,
        state.quark_params,
        state.thermo_params,
        state.K_coeffs;
        p_nodes=tau_p_nodes,
        angle_nodes=tau_angle_nodes,
        phi_nodes=tau_phi_nodes,
        cs_cache=cache,
        n_sigma_points=tau_n_sigma_points,
        sigma_grid_n=sigma_grid_n,
        sigma_cutoff=sigma_cutoff,
        apply_s_domain_cut=true,
    )
end

function refine_cache_local(process::Symbol, state, base_cache;
    n_local::Int,
    halfwidth::Float64,
    n_sigma_points::Int,
)
    n_local <= 0 && return base_cache
    length(base_cache.s_vals) == 0 && return base_cache

    peak_idx = argmax(base_cache.sigma_vals)
    s_peak = base_cache.s_vals[peak_idx]
    s_min = base_cache.s_vals[1]
    s_max = base_cache.s_vals[end]
    lo = max(s_min, s_peak - halfwidth)
    hi = min(s_max, s_peak + halfwidth)
    hi <= lo && return base_cache

    local_grid = collect(range(lo, hi, length=n_local))
    merged = Dict{Float64, Float64}()
    for (s, sigma) in zip(base_cache.s_vals, base_cache.sigma_vals)
        merged[s] = sigma
    end
    for s in local_grid
        sigma = TotalCS.total_cross_section(process, s, state.quark_params, state.thermo_params, state.K_coeffs; n_points=n_sigma_points)
        merged[s] = sigma
    end

    s_all = sort(collect(keys(merged)))
    sigma_all = [merged[s] for s in s_all]
    return RT_ASR.CrossSectionCache(process, s_all, sigma_all)
end

function main(args=ARGS)
    opts = parse_args(args)
    ensure_parent_dir(opts.out_csv)

    rows = read_comment_csv(opts.scan_csv)
    row_by_xi = Dict{Float64, Dict{String, String}}(to_f(r, "xi") => r for r in rows)

    io = open(opts.out_csv, "w")
    try
        println(io, join([
            "xi", "process", "config", "adaptive_refinement", "local_points", "local_halfwidth",
            "cache_n_points", "cache_sigma_max", "cache_sigma_mean",
            "s_bo", "s_up", "s_peak", "s_peak_minus_bo", "rate",
        ], ','))

        for xi in opts.xi_list
            haskey(row_by_xi, xi) || error("xi not found in scan csv: $xi")
            state = row_to_state(row_by_xi[xi])

            base_cache = RT_ASR.build_w0cdf_pchip_cache(
                opts.process,
                state.quark_params,
                state.thermo_params,
                state.K_coeffs;
                N=opts.sigma_grid_n,
                p_cutoff=Λ_inv_fm,
                n_sigma_points=opts.tau_n_sigma_points,
                adaptive_refinement=false,
            )
            base_rate = eval_rate(opts.process, state, base_cache;
                tau_p_nodes=opts.tau_p_nodes,
                tau_angle_nodes=opts.tau_angle_nodes,
                tau_phi_nodes=opts.tau_phi_nodes,
                tau_n_sigma_points=opts.tau_n_sigma_points,
                sigma_grid_n=opts.sigma_grid_n,
                sigma_cutoff=Λ_inv_fm,
            )

            adapt_cache = RT_ASR.build_w0cdf_pchip_cache(
                opts.process,
                state.quark_params,
                state.thermo_params,
                state.K_coeffs;
                N=opts.sigma_grid_n,
                p_cutoff=Λ_inv_fm,
                n_sigma_points=opts.tau_n_sigma_points,
                adaptive_refinement=true,
            )
            adapt_rate = eval_rate(opts.process, state, adapt_cache;
                tau_p_nodes=opts.tau_p_nodes,
                tau_angle_nodes=opts.tau_angle_nodes,
                tau_phi_nodes=opts.tau_phi_nodes,
                tau_n_sigma_points=opts.tau_n_sigma_points,
                sigma_grid_n=opts.sigma_grid_n,
                sigma_cutoff=Λ_inv_fm,
            )

            for (cfg_name, adaptive, local_points, cache, rate) in (
                ("base", false, 0, base_cache, base_rate),
                ("adaptive_true", true, 0, adapt_cache, adapt_rate),
            )
                win = process_window(opts.process, state, Λ_inv_fm)
                peak_idx = argmax(cache.sigma_vals)
                s_peak = cache.s_vals[peak_idx]
                println(io, join([
                    string(xi), string(opts.process), cfg_name, string(adaptive), string(local_points), string(opts.local_halfwidth),
                    string(length(cache.s_vals)), string(maximum(cache.sigma_vals)), string(sum(cache.sigma_vals) / length(cache.sigma_vals)),
                    string(win.s_bo), string(win.s_up), string(s_peak), string(s_peak - win.s_bo), string(rate),
                ], ','))
            end

            for local_points in opts.local_points_list
                local_points <= 0 && continue
                cache_refined = refine_cache_local(opts.process, state, base_cache;
                    n_local=local_points,
                    halfwidth=opts.local_halfwidth,
                    n_sigma_points=opts.tau_n_sigma_points,
                )
                rate_refined = eval_rate(opts.process, state, cache_refined;
                    tau_p_nodes=opts.tau_p_nodes,
                    tau_angle_nodes=opts.tau_angle_nodes,
                    tau_phi_nodes=opts.tau_phi_nodes,
                    tau_n_sigma_points=opts.tau_n_sigma_points,
                    sigma_grid_n=opts.sigma_grid_n,
                    sigma_cutoff=Λ_inv_fm,
                )
                win = process_window(opts.process, state, Λ_inv_fm)
                peak_idx = argmax(cache_refined.sigma_vals)
                s_peak = cache_refined.s_vals[peak_idx]
                println(io, join([
                    string(xi), string(opts.process), "manual_local_refine", "false", string(local_points), string(opts.local_halfwidth),
                    string(length(cache_refined.s_vals)), string(maximum(cache_refined.sigma_vals)), string(sum(cache_refined.sigma_vals) / length(cache_refined.sigma_vals)),
                    string(win.s_bo), string(win.s_up), string(s_peak), string(s_peak - win.s_bo), string(rate_refined),
                ], ','))
            end
        end

        @printf("saved: %s\n", opts.out_csv)
    finally
        close(io)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
