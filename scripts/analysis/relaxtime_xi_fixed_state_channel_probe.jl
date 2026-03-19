#!/usr/bin/env julia

"""
固定态（冻结平衡态）下的 ξ 邻域 τ^-1 通道诊断：
1) 扫描 sigma 网格/自适应参数，检查 τ^-1 是否收敛；
2) 对可疑通道临时关闭 s 硬切（apply_s_domain_cut=false）做对照，比较尖点是否收缩。

输入依赖：run_gap_transport_scan.jl 生成的 plan_b CSV（包含 T_fm, muq_fm, xi, Phi, Phibar, m_*, n_*）。
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

const RT = Main.RelaxationTime
const RT_ASR = Main.AverageScatteringRate
const REQUIRED_PROCESSES = RT.REQUIRED_PROCESSES

function required_processes_list()
    out = Symbol[]
    for item in REQUIRED_PROCESSES
        if item isa Symbol
            push!(out, item)
        elseif item isa Tuple
            for sub in item
                sub isa Symbol && push!(out, sub)
            end
        else
            error("unsupported REQUIRED_PROCESSES entry type: $(typeof(item))")
        end
    end
    return unique(out)
end

Base.@kwdef struct Options
    scan_csv::String = joinpath("data", "outputs", "results", "relaxtime", "plan_b", "transport_vs_xi_T190_muB0.csv")
    centers::Vector{Float64} = Float64[-0.44]
    species::Symbol = :u
    sigma_grid_list::Vector{Int} = Int[64, 128, 256, 512]
    adaptive_modes::Vector{Bool} = Bool[true, false]
    disable_s_cut_channels::Vector{Symbol} = Symbol[]
    processes::Vector{Symbol} = Symbol[]
    tau_p_nodes::Int = 28
    tau_angle_nodes::Int = 8
    tau_phi_nodes::Int = 8
    tau_n_sigma_points::Int = 6
    fine_xi_min::Union{Nothing, Float64} = nothing
    fine_xi_max::Union{Nothing, Float64} = nothing
    fine_xi_step::Union{Nothing, Float64} = nothing
    anchor_xi::Union{Nothing, Float64} = nothing
    summary_out::String = joinpath("data", "outputs", "results", "relaxtime", "scan", "_xi_fixed_state_probe_summary.csv")
    channels_out::String = joinpath("data", "outputs", "results", "relaxtime", "scan", "_xi_fixed_state_probe_channels.csv")
    switch_log_out::Union{Nothing, String} = nothing
end

function usage()
    println("Usage: julia --project=. scripts/analysis/relaxtime_xi_fixed_state_channel_probe.jl [options]\n")
    println("Options:")
    println("  --scan-csv <path>               输入 plan_b CSV")
    println("  --centers <csv>                 中心 xi 列表，如 -0.44,-0.2")
    println("  --species <u|s|ubar|sbar>       关注物种（default u）")
    println("  --sigma-grid-list <csv-int>     sigma 网格列表（default 64,128,256,512）")
    println("  --adaptive-modes <csv-bool>     自适应模式列表（default true,false）")
    println("  --disable-s-cut-channels <csv>  对照场景里关闭硬切的 channel 列表")
    println("  --processes <csv>               仅计算给定 process 列表；默认全部 REQUIRED_PROCESSES")
    println("  --tau-p-nodes <int>             τ 动量积分节点")
    println("  --tau-angle-nodes <int>         τ 角度积分节点")
    println("  --tau-phi-nodes <int>           τ 方位角积分节点")
    println("  --tau-n-sigma <int>             σ(s) t 积分点")
    println("  --fine-xi-min <float>           fixed-state 细扫 xi 最小值")
    println("  --fine-xi-max <float>           fixed-state 细扫 xi 最大值")
    println("  --fine-xi-step <float>          fixed-state 细扫 xi 步长")
    println("  --anchor-xi <float>             fixed-state 锚点 xi（默认用 centers[1]）")
    println("  --summary-out <path>            摘要输出")
    println("  --channels-out <path>           通道明细输出")
    println("  --switch-log-out <path>         通道开关命中日志输出")
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

function parse_bool_list(raw::AbstractString)
    vals = Bool[]
    for p in split(raw, ',')
        s = lowercase(strip(p))
        isempty(s) && continue
        if s in ("1", "true", "yes")
            push!(vals, true)
        elseif s in ("0", "false", "no")
            push!(vals, false)
        else
            error("invalid bool token: $p")
        end
    end
    isempty(vals) && error("empty bool list")
    return vals
end

function parse_symbol_list(raw::AbstractString)
    syms = Symbol[]
    for p in split(raw, ',')
        s = strip(p)
        isempty(s) && continue
        push!(syms, Symbol(s))
    end
    return syms
end

function parse_args(args::Vector{String})
    scan_csv = joinpath("data", "outputs", "results", "relaxtime", "plan_b", "transport_vs_xi_T190_muB0.csv")
    centers = Float64[-0.44]
    species = :u
    sigma_grid_list = Int[64, 128, 256, 512]
    adaptive_modes = Bool[true, false]
    disable_s_cut_channels = Symbol[]
    processes = Symbol[]
    tau_p_nodes = 28
    tau_angle_nodes = 8
    tau_phi_nodes = 8
    tau_n_sigma_points = 6
    fine_xi_min = nothing
    fine_xi_max = nothing
    fine_xi_step = nothing
    anchor_xi = nothing
    summary_out = joinpath("data", "outputs", "results", "relaxtime", "scan", "_xi_fixed_state_probe_summary.csv")
    channels_out = joinpath("data", "outputs", "results", "relaxtime", "scan", "_xi_fixed_state_probe_channels.csv")
    switch_log_out = nothing

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
        elseif arg == "--centers"
            centers = parse_float_list(require_value())
        elseif arg == "--species"
            species = Symbol(require_value())
        elseif arg == "--sigma-grid-list"
            sigma_grid_list = parse_int_list(require_value())
        elseif arg == "--adaptive-modes"
            adaptive_modes = parse_bool_list(require_value())
        elseif arg == "--disable-s-cut-channels"
            disable_s_cut_channels = parse_symbol_list(require_value())
        elseif arg == "--processes"
            processes = parse_symbol_list(require_value())
        elseif arg == "--tau-p-nodes"
            tau_p_nodes = parse(Int, require_value())
        elseif arg == "--tau-angle-nodes"
            tau_angle_nodes = parse(Int, require_value())
        elseif arg == "--tau-phi-nodes"
            tau_phi_nodes = parse(Int, require_value())
        elseif arg == "--tau-n-sigma"
            tau_n_sigma_points = parse(Int, require_value())
        elseif arg == "--fine-xi-min"
            fine_xi_min = parse(Float64, require_value())
        elseif arg == "--fine-xi-max"
            fine_xi_max = parse(Float64, require_value())
        elseif arg == "--fine-xi-step"
            fine_xi_step = parse(Float64, require_value())
        elseif arg == "--anchor-xi"
            anchor_xi = parse(Float64, require_value())
        elseif arg == "--summary-out"
            summary_out = require_value()
        elseif arg == "--channels-out"
            channels_out = require_value()
        elseif arg == "--switch-log-out"
            switch_log_out = require_value()
        elseif arg in ("-h", "--help")
            usage(); exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return Options(
        scan_csv=scan_csv,
        centers=centers,
        species=species,
        sigma_grid_list=sigma_grid_list,
        adaptive_modes=adaptive_modes,
        disable_s_cut_channels=disable_s_cut_channels,
        processes=processes,
        tau_p_nodes=tau_p_nodes,
        tau_angle_nodes=tau_angle_nodes,
        tau_phi_nodes=tau_phi_nodes,
        tau_n_sigma_points=tau_n_sigma_points,
        fine_xi_min=fine_xi_min,
        fine_xi_max=fine_xi_max,
        fine_xi_step=fine_xi_step,
        anchor_xi=anchor_xi,
        summary_out=summary_out,
        channels_out=channels_out,
        switch_log_out=switch_log_out,
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
    densities = (
        u=to_f(row, "n_u"), d=to_f(row, "n_d"), s=to_f(row, "n_s"),
        ubar=to_f(row, "n_ubar"), dbar=to_f(row, "n_dbar"), sbar=to_f(row, "n_sbar"),
    )
    ktmp = build_K_data(T_fm, muq_fm, masses, Phi, Phibar)
    quark_params = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm), A=ktmp.A_vals)
    thermo_params = (T=T_fm, Φ=Phi, Φbar=Phibar, ξ=xi)
    return (
        T_fm=T_fm,
        muq_fm=muq_fm,
        xi=xi,
        quark_params=quark_params,
        thermo_params=thermo_params,
        densities=densities,
        K_coeffs=ktmp.K_coeffs,
    )
end

function state_with_xi(anchor_state, xi::Float64)
    T_fm = anchor_state.T_fm
    muq_fm = anchor_state.muq_fm
    masses = anchor_state.quark_params.m
    Phi = anchor_state.thermo_params.Φ
    Phibar = anchor_state.thermo_params.Φbar
    ktmp = build_K_data(T_fm, muq_fm, masses, Phi, Phibar)
    quark_params = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm), A=ktmp.A_vals)
    thermo_params = (T=T_fm, Φ=Phi, Φbar=Phibar, ξ=xi)
    return (
        T_fm=T_fm,
        muq_fm=muq_fm,
        xi=xi,
        quark_params=quark_params,
        thermo_params=thermo_params,
        densities=anchor_state.densities,
        K_coeffs=ktmp.K_coeffs,
    )
end

function rates_dict_to_namedtuple(d::Dict{Symbol, Float64})
    keys_sorted = sort(collect(keys(d)); by=String)
    return (; (k => d[k] for k in keys_sorted)...)
end

function complete_rates_namedtuple(d::Dict{Symbol, Float64})
    req = required_processes_list()
    all_keys = sort(collect(Set(vcat(req, collect(keys(d))))); by=String)
    return (; (k => get(d, k, 0.0) for k in all_keys)...)
end

function contribution_rows_from_rates(densities, rates_nt)
    return NamedTuple[]
end

function process_window_diagnostics(process::Symbol, state, sigma_cutoff::Float64, apply_requested::Bool)
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
    apply_effective = apply_requested && !isnan(Λ)
    return (s_bo=s_bo, s_up=s_up, apply_effective=apply_effective)
end

function compute_state_rates(state;
    sigma_grid_n::Int,
    adaptive_refinement::Bool,
    tau_p_nodes::Int,
    tau_angle_nodes::Int,
    tau_phi_nodes::Int,
    tau_n_sigma_points::Int,
    processes::Vector{Symbol}=required_processes_list(),
    sigma_cutoff::Float64=Λ_inv_fm,
    disable_s_cut_channels::Set{Symbol}=Set{Symbol}(),
)
    rates = Dict{Symbol, Float64}()
    switch_hits = NamedTuple[]

    for process in processes
        cache = RT_ASR.build_w0cdf_pchip_cache(
            process,
            state.quark_params,
            state.thermo_params,
            state.K_coeffs;
            N=sigma_grid_n,
            design_p_nodes=tau_p_nodes,
            design_angle_nodes=tau_angle_nodes,
            design_phi_nodes=tau_phi_nodes,
            p_cutoff=sigma_cutoff,
            n_sigma_points=tau_n_sigma_points,
        )

        apply_s_cut = !(process in disable_s_cut_channels)
        win_diag = process_window_diagnostics(process, state, sigma_cutoff, apply_s_cut)
        rates[process] = RT_ASR.average_scattering_rate(
            process,
            state.quark_params,
            state.thermo_params,
            state.K_coeffs;
            p_nodes=tau_p_nodes,
            angle_nodes=tau_angle_nodes,
            phi_nodes=tau_phi_nodes,
            cs_cache=cache,
            n_sigma_points=tau_n_sigma_points,
            sigma_cutoff=sigma_cutoff,
            apply_s_domain_cut=apply_s_cut,
        )

        cache_n = length(cache.s_vals)
        cache_s_min = cache_n > 0 ? cache.s_vals[1] : NaN
        cache_s_max = cache_n > 0 ? cache.s_vals[end] : NaN
        cache_sigma_min = cache_n > 0 ? minimum(cache.sigma_vals) : NaN
        cache_sigma_max = cache_n > 0 ? maximum(cache.sigma_vals) : NaN
        cache_sigma_mean = cache_n > 0 ? (sum(cache.sigma_vals) / cache_n) : NaN
        peak_idx = cache_n > 0 ? argmax(cache.sigma_vals) : 0
        s_peak = peak_idx > 0 ? cache.s_vals[peak_idx] : NaN
        s_peak_minus_bo = peak_idx > 0 ? (s_peak - win_diag.s_bo) : NaN
        push!(switch_hits, (
            process=process,
            xi=state.xi,
            apply_s_cut_requested=apply_s_cut,
            apply_s_cut_effective=win_diag.apply_effective,
            disabled_by_scenario=(process in disable_s_cut_channels),
            sigma_cutoff=sigma_cutoff,
            cache_n_points=cache_n,
            cache_s_min=cache_s_min,
            cache_s_max=cache_s_max,
            cache_sigma_min=cache_sigma_min,
            cache_sigma_max=cache_sigma_max,
            cache_sigma_mean=cache_sigma_mean,
            s_peak=s_peak,
            s_peak_minus_bo=s_peak_minus_bo,
            cache_asym_enabled=cache.asym_enabled,
            cache_asym_s0=cache.asym_s0,
            cache_asym_A=cache.asym_A,
            s_bo=win_diag.s_bo,
            s_up=win_diag.s_up,
            rate=rates[process],
        ))
    end

    rates_nt = complete_rates_namedtuple(rates)
    tau_inv = RT.relaxation_rates(state.densities, rates_nt)
    tau = (
        u = 1.0 / tau_inv.u,
        d = 1.0 / tau_inv.d,
        s = 1.0 / tau_inv.s,
        ubar = 1.0 / tau_inv.ubar,
        dbar = 1.0 / tau_inv.dbar,
        sbar = 1.0 / tau_inv.sbar,
    )
    contrib = if isdefined(RT, :relaxation_rate_contribution_rows)
        RT.relaxation_rate_contribution_rows(state.densities, rates_nt)
    else
        contribution_rows_from_rates(state.densities, rates_nt)
    end
    return (tau=tau, tau_inv=tau_inv, rates=rates_nt, contributions=contrib, switch_hits=switch_hits)
end

@inline function kink_metric(prev::Float64, curr::Float64, next::Float64)
    return abs(curr - 0.5 * (prev + next)) / max(abs(curr), 1e-12)
end

function main(args=ARGS)
    opts = parse_args(args)
    ensure_parent_dir(opts.summary_out)
    ensure_parent_dir(opts.channels_out)
    if opts.switch_log_out !== nothing
        ensure_parent_dir(opts.switch_log_out)
    end

    rows = read_comment_csv(opts.scan_csv)
    sort!(rows, by=row -> to_f(row, "xi"))
    xis = [to_f(r, "xi") for r in rows]

    row_by_xi = Dict{Float64, Dict{String, String}}(to_f(r, "xi") => r for r in rows)
    process_list = isempty(opts.processes) ? required_processes_list() : opts.processes
    process_set = Set{Symbol}(process_list)

    scenarios = [
        (name="baseline", disable=Set{Symbol}()),
        (name="softcut_off", disable=Set{Symbol}(opts.disable_s_cut_channels)),
    ]
    if isempty(opts.disable_s_cut_channels)
        scenarios = [(name="baseline", disable=Set{Symbol}())]
    end

    summary_io = open(opts.summary_out, "w")
    channels_io = open(opts.channels_out, "w")
    switch_io = opts.switch_log_out === nothing ? nothing : open(opts.switch_log_out, "w")
    try
        println(summary_io, join([
            "center_xi", "scenario", "sigma_grid_n", "adaptive_refinement", "species",
            "xi_prev", "xi_curr", "xi_next",
            "tauinv_prev", "tauinv_curr", "tauinv_next", "kink_metric_tauinv",
            "tau_prev", "tau_curr", "tau_next", "kink_metric_tau",
        ], ','))

        println(channels_io, join([
            "center_xi", "scenario", "sigma_grid_n", "adaptive_refinement", "species",
            "xi", "channel", "density_key", "multiplicity", "density", "rate", "contribution", "total",
        ], ','))
        if switch_io !== nothing
            println(switch_io, join([
                "center_xi", "scenario", "sigma_grid_n", "adaptive_refinement", "species",
                "xi", "process", "apply_s_cut_requested", "apply_s_cut_effective", "disabled_by_scenario",
                "sigma_cutoff", "cache_n_points", "cache_s_min", "cache_s_max",
                "cache_sigma_min", "cache_sigma_max", "cache_sigma_mean",
                "s_peak", "s_peak_minus_bo",
                "cache_asym_enabled", "cache_asym_s0", "cache_asym_A", "s_bo", "s_up", "rate",
            ], ','))
        end

        species = opts.species

        state_map = Dict{Float64, Any}()
        eval_centers = Float64[]
        if opts.fine_xi_min !== nothing || opts.fine_xi_max !== nothing || opts.fine_xi_step !== nothing
            opts.fine_xi_min === nothing && error("fine scan requires --fine-xi-min")
            opts.fine_xi_max === nothing && error("fine scan requires --fine-xi-max")
            opts.fine_xi_step === nothing && error("fine scan requires --fine-xi-step")
            opts.fine_xi_step <= 0 && error("--fine-xi-step must be positive")

            anchor_xi = opts.anchor_xi === nothing ? opts.centers[1] : opts.anchor_xi
            haskey(row_by_xi, anchor_xi) || error("anchor xi not found in scan csv: $anchor_xi")
            anchor_state = row_to_state(row_by_xi[anchor_xi])
            xi_grid = collect(opts.fine_xi_min:opts.fine_xi_step:opts.fine_xi_max)
            length(xi_grid) >= 3 || error("fine xi grid must have at least 3 points")
            for x in xi_grid
                state_map[x] = state_with_xi(anchor_state, x)
            end
            eval_centers = xi_grid[2:(end - 1)]
        else
            for x in xis
                state_map[x] = row_to_state(row_by_xi[x])
            end
            eval_centers = copy(opts.centers)
        end

        for center in eval_centers
            idx = findfirst(x -> abs(x - center) <= 1e-9, sort(collect(keys(state_map))))
            idx === nothing && error("center xi not found in state grid: $center")
            grid_sorted = sort(collect(keys(state_map)))
            (idx == 1 || idx == length(grid_sorted)) && continue

            xi_prev = grid_sorted[idx - 1]
            xi_curr = grid_sorted[idx]
            xi_next = grid_sorted[idx + 1]

            state_prev = state_map[xi_prev]
            state_curr = state_map[xi_curr]
            state_next = state_map[xi_next]

            for sigma_n in opts.sigma_grid_list
                for adaptive in opts.adaptive_modes
                    for sc in scenarios
                        rp = compute_state_rates(state_prev;
                            sigma_grid_n=sigma_n,
                            adaptive_refinement=adaptive,
                            tau_p_nodes=opts.tau_p_nodes,
                            tau_angle_nodes=opts.tau_angle_nodes,
                            tau_phi_nodes=opts.tau_phi_nodes,
                            tau_n_sigma_points=opts.tau_n_sigma_points,
                            processes=process_list,
                            sigma_cutoff=Λ_inv_fm,
                            disable_s_cut_channels=sc.disable,
                        )
                        rc = compute_state_rates(state_curr;
                            sigma_grid_n=sigma_n,
                            adaptive_refinement=adaptive,
                            tau_p_nodes=opts.tau_p_nodes,
                            tau_angle_nodes=opts.tau_angle_nodes,
                            tau_phi_nodes=opts.tau_phi_nodes,
                            tau_n_sigma_points=opts.tau_n_sigma_points,
                            processes=process_list,
                            sigma_cutoff=Λ_inv_fm,
                            disable_s_cut_channels=sc.disable,
                        )
                        rn = compute_state_rates(state_next;
                            sigma_grid_n=sigma_n,
                            adaptive_refinement=adaptive,
                            tau_p_nodes=opts.tau_p_nodes,
                            tau_angle_nodes=opts.tau_angle_nodes,
                            tau_phi_nodes=opts.tau_phi_nodes,
                            tau_n_sigma_points=opts.tau_n_sigma_points,
                            processes=process_list,
                            sigma_cutoff=Λ_inv_fm,
                            disable_s_cut_channels=sc.disable,
                        )

                        tauinv_prev = getproperty(rp.tau_inv, species)
                        tauinv_curr = getproperty(rc.tau_inv, species)
                        tauinv_next = getproperty(rn.tau_inv, species)

                        tau_prev = getproperty(rp.tau, species)
                        tau_curr = getproperty(rc.tau, species)
                        tau_next = getproperty(rn.tau, species)

                        k_tauinv = kink_metric(tauinv_prev, tauinv_curr, tauinv_next)
                        k_tau = kink_metric(tau_prev, tau_curr, tau_next)

                        println(summary_io, join([
                            string(center), sc.name, string(sigma_n), string(adaptive), string(species),
                            string(xi_prev), string(xi_curr), string(xi_next),
                            string(tauinv_prev), string(tauinv_curr), string(tauinv_next), string(k_tauinv),
                            string(tau_prev), string(tau_curr), string(tau_next), string(k_tau),
                        ], ','))

                        for (xi_val, bundle) in ((xi_prev, rp), (xi_curr, rc), (xi_next, rn))
                            for row in bundle.contributions
                                row.species == species || continue
                                (row.channel in process_set) || continue
                                println(channels_io, join([
                                    string(center), sc.name, string(sigma_n), string(adaptive), string(species),
                                    string(xi_val), string(row.channel), string(row.density_key), string(row.multiplicity),
                                    string(row.density), string(row.rate), string(row.contribution), string(row.total),
                                ], ','))
                            end
                            if switch_io !== nothing
                                for sw in bundle.switch_hits
                                    println(switch_io, join([
                                        string(center), sc.name, string(sigma_n), string(adaptive), string(species),
                                        string(sw.xi), string(sw.process), string(sw.apply_s_cut_requested),
                                        string(sw.apply_s_cut_effective), string(sw.disabled_by_scenario),
                                        string(sw.sigma_cutoff), string(sw.cache_n_points), string(sw.cache_s_min),
                                        string(sw.cache_s_max), string(sw.cache_sigma_min), string(sw.cache_sigma_max),
                                        string(sw.cache_sigma_mean), string(sw.s_peak), string(sw.s_peak_minus_bo),
                                        string(sw.cache_asym_enabled), string(sw.cache_asym_s0), string(sw.cache_asym_A),
                                        string(sw.s_bo), string(sw.s_up), string(sw.rate),
                                    ], ','))
                                end
                            end
                        end
                    end
                end
            end
        end

        @printf("saved summary: %s\n", opts.summary_out)
        @printf("saved channels: %s\n", opts.channels_out)
        if opts.switch_log_out !== nothing
            @printf("saved switch log: %s\n", opts.switch_log_out)
        end
    finally
        close(summary_io)
        close(channels_io)
        if switch_io !== nothing
            close(switch_io)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
