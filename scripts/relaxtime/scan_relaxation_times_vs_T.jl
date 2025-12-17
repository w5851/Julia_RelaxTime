"""
扫描弛豫时间 τ(T) 在给定 μ_B 下的温度依赖。

需求场景：输出 u, s, ubar, sbar 的 τ(T)（及 τ^{-1}），并支持 μ_B=0/800 MeV 两条曲线。

实现要点：
- 先解各向同性/各向异性 PNJL 平衡（能隙+Polyakov 环）
- 构造 TotalCrossSection 所需的 quark_params(A 字段) 与有效耦合 K_coeffs
- 用 RelaxationTime.relaxation_times 计算 τ
- 为了性能：对每个 T/μ_B，给每个散射过程预计算一个粗网格 σ(s)，并在 ω 积分中仅做线性插值
  （CrossSectionCache(process; compute_missing=false)）。

单位约定：
- CLI 温度/化学势使用 MeV
- 内部计算使用 fm⁻¹
- 输出 τ 单位为 fm

示例：
  julia --project=. scripts/relaxtime/scan_relaxation_times_vs_T.jl \
    --tmin 80 --tmax 220 --tstep 10 \
    --mub-list 0,800 \
    --tau-phi-nodes 6 \
    --out data/outputs/results/relaxtime/relaxation_times_vs_T.csv
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "solvers", "AnisoGapSolver.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "analysis", "ThermoDerivatives.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxationTime.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using Printf
using StaticArrays

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm
using .GaussLegendre: gauleg, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .ThermoDerivatives: solve_equilibrium_mu
using .AnisoGapSolver: cached_nodes, calculate_mass_vec, calculate_number_densities, DEFAULT_MU_GUESS
using .RelaxationTime: relaxation_times, REQUIRED_PROCESSES
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .ScanCSV: ScanCSV

const RT_ASR = RelaxationTime.AverageScatteringRate
const RT_TCS = RelaxationTime.TotalCrossSection

struct Options
    out::String
    overwrite::Bool
    resume::Bool
    xi::Float64
    tmin_mev::Float64
    tmax_mev::Float64
    tstep_mev::Float64
    mub_list_mev::Vector{Float64}
    p_num::Int
    t_num::Int
    max_iter::Int
    # tau integration nodes
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
    # sigma precompute grid
    sigma_grid_n::Int
    p_cut_factor::Float64
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/scan_relaxation_times_vs_T.jl [options]\n")
    println("Options:")
    println("  --out <path>                 输出 CSV (default data/outputs/results/relaxtime/relaxation_times_vs_T.csv)")
    println("  --overwrite                  覆盖输出文件")
    println("  --no-resume                  禁用跳过逻辑，强制重算")
    println("  --xi <value>                 各向异性参数 ξ (default 0.0)")
    println("  --tmin/--tmax/--tstep <MeV>  温度范围与步长")
    println("  --mub-list v1,v2,...         μ_B 列表 (MeV), default 0,800")
    println("  --p-num <int>                能隙求解/密度的动量节点数 (default 12)")
    println("  --t-num <int>                能隙求解/密度的角度节点数 (default 6)")
    println("  --max-iter <int>             NLsolve iterations 上限 (default 40)")
    println("  --tau-p-nodes <int>          τ 平均散射率动量节点 (default 6)")
    println("  --tau-angle-nodes <int>      τ 平均散射率 cosθ 节点 (default 2)")
    println("  --tau-phi-nodes <int>        τ 平均散射率 φ 节点 (default 6)")
    println("  --tau-n-sigma <int>          σ(s) 的 t 积分点数 (default 6)")
    println("  --sigma-grid-n <int>         每个过程预计算 σ(s) 的网格点数 (default 12)")
    println("  --p-cut-factor <float>       σ(s) 预计算 s_max 使用的 p_cut = min(Λ, factor*T) (default 8.0)")
    println("  -h, --help                   显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :out => joinpath("data", "outputs", "results", "relaxtime", "relaxation_times_vs_T.csv"),
        :overwrite => false,
        :resume => true,
        :xi => 0.0,
        :tmin => 80.0,
        :tmax => 220.0,
        :tstep => 10.0,
        :mub_list => Float64[0.0, 800.0],
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :tau_p_nodes => 6,
        :tau_angle_nodes => 2,
        :tau_phi_nodes => 6,
        :tau_n_sigma => 6,
        :sigma_grid_n => 12,
        :p_cut_factor => 8.0,
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

        if arg == "--out"
            opts[:out] = require_value()
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--no-resume"
            opts[:resume] = false
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--tmin"
            opts[:tmin] = parse(Float64, require_value())
        elseif arg == "--tmax"
            opts[:tmax] = parse(Float64, require_value())
        elseif arg == "--tstep"
            opts[:tstep] = parse(Float64, require_value())
        elseif arg == "--mub-list"
            raw = split(require_value(), ',')
            vals = Float64[parse(Float64, strip(v)) for v in raw if !isempty(strip(v))]
            isempty(vals) && error("empty mub-list")
            opts[:mub_list] = vals
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, require_value())
        elseif arg == "--tau-p-nodes"
            opts[:tau_p_nodes] = parse(Int, require_value())
        elseif arg == "--tau-angle-nodes"
            opts[:tau_angle_nodes] = parse(Int, require_value())
        elseif arg == "--tau-phi-nodes"
            opts[:tau_phi_nodes] = parse(Int, require_value())
        elseif arg == "--tau-n-sigma"
            opts[:tau_n_sigma] = parse(Int, require_value())
        elseif arg == "--sigma-grid-n"
            opts[:sigma_grid_n] = parse(Int, require_value())
        elseif arg == "--p-cut-factor"
            opts[:p_cut_factor] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage(); exit(0)
        else
            error("unknown option: $arg")
        end

        i += 1
    end

    tstep = Float64(opts[:tstep]); tstep > 0 || error("tstep must be positive")
    return Options(
        String(opts[:out]),
        Bool(opts[:overwrite]),
        Bool(opts[:resume]),
        Float64(opts[:xi]),
        Float64(opts[:tmin]),
        Float64(opts[:tmax]),
        Float64(opts[:tstep]),
        Float64.(opts[:mub_list]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Int(opts[:tau_p_nodes]),
        Int(opts[:tau_angle_nodes]),
        Int(opts[:tau_phi_nodes]),
        Int(opts[:tau_n_sigma]),
        Int(opts[:sigma_grid_n]),
        Float64(opts[:p_cut_factor]),
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function write_header(io)
    cols = [
        "T_MeV", "muB_MeV", "muq_MeV", "xi",
        "T_fm", "muq_fm",
        "converged", "iterations", "residual_norm",
        "Phi", "Phibar",
        "m_u", "m_s",
        "n_u", "n_s", "n_ubar", "n_sbar",
        "tau_u", "tau_s", "tau_ubar", "tau_sbar",
        "tauinv_u", "tauinv_s", "tauinv_ubar", "tauinv_sbar",
        "tau_p_nodes", "tau_angle_nodes", "tau_phi_nodes", "tau_n_sigma_points",
        "sigma_grid_n", "p_cut_factor",
    ]
    println(io, join(cols, ','))
end

@inline function csv_bool(x::Bool)
    return x ? "true" : "false"
end

function build_K_coeffs(T_fm::Float64, muq_fm::Float64, masses::NamedTuple, Φ::Float64, Φbar::Float64)
    # Need A values for TotalPropagator/EffectiveCouplings.
    A_u = A(masses.u, muq_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(masses.s, muq_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return (K_coeffs=calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s), A_vals=(u=A_u, d=A_u, s=A_s))
end

function densities_from_equilibrium(x_state, mu_vec, T_fm, thermal_nodes, xi)
    nd = calculate_number_densities(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return (
        u=Float64(nd.quark[1]),
        d=Float64(nd.quark[2]),
        s=Float64(nd.quark[3]),
        ubar=Float64(nd.antiquark[1]),
        dbar=Float64(nd.antiquark[2]),
        sbar=Float64(nd.antiquark[3]),
    )
end

@inline function s_bounds_for_process(process::Symbol, quark_params::NamedTuple, thermo_params::NamedTuple; p_cut::Float64)
    pi_sym, pj_sym, _, _ = RT_ASR.parse_particles_from_process(process)
    mi = RT_ASR.get_mass(pi_sym, quark_params)
    mj = RT_ASR.get_mass(pj_sym, quark_params)

    s_thr = (mi + mj)^2
    # Stay away from exact threshold to avoid integrand singularities (B0/Π) causing NaNs.
    s_min = s_thr + max(1e-6, 1e-4 * s_thr)

    Ei = sqrt(p_cut^2 + mi^2)
    Ej = sqrt(p_cut^2 + mj^2)
    # Max when cosΘ=-1 -> p_dot = Ei*Ej + p_cut^2
    s_max = mi^2 + mj^2 + 2.0 * (Ei * Ej + p_cut^2)
    return (s_min=s_min, s_max=s_max)
end

function safe_total_cross_section(process::Symbol, s::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int, max_tries::Int=4)

    s_try = s
    last_err = nothing
    for _ in 1:max_tries
        try
            σ = RT_TCS.total_cross_section(process, s_try, quark_params, thermo_params, K_coeffs; n_points=n_points)
            if isfinite(σ)
                return σ
            end
        catch err
            last_err = err
        end
        # Nudge upward and retry
        s_try = s_try * (1.0 + 1e-6) + 1e-10
    end

    @warn "failed to compute sigma; returning NaN and continuing" process=process s=s last_error=last_err
    return NaN
end

function build_sigma_caches(processes::Tuple, quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_sigma_points::Int, sigma_grid_n::Int, p_cut_factor::Float64)

    p_cut = min(Λ_inv_fm, p_cut_factor * thermo_params.T)
    cs_caches = Dict{Symbol,RT_ASR.CrossSectionCache}()

    for process in processes
        bounds = s_bounds_for_process(process, quark_params, thermo_params; p_cut=p_cut)
        s_min = bounds.s_min
        s_max = bounds.s_max
        s_min < s_max || error("invalid s bounds for $process")

        # Use uniform grid in sqrt(s) for better resolution near threshold.
        sqrt_s_grid = range(sqrt(s_min), sqrt(s_max); length=sigma_grid_n)
        s_grid = Float64[(x * x) for x in sqrt_s_grid]

        cache = RT_ASR.CrossSectionCache(process; compute_missing=true)
        n_ok = 0
        n_bad = 0
        for s in s_grid
            σ = safe_total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
            if !isfinite(σ)
                n_bad += 1
                continue
            end
            RT_ASR.insert_sigma!(cache, s, σ)
            n_ok += 1
        end
        if n_bad > 0
            @warn "sigma grid had non-finite points (skipped)" process=process n_ok=n_ok n_bad=n_bad
        end
        if n_ok < 2
            @warn "sigma grid is sparse; will rely on on-demand sigma computation" process=process n_ok=n_ok n_bad=n_bad
        end
        cs_caches[process] = cache
    end

    return cs_caches
end

function run_scan(opts::Options)
    ensure_parent_dir(opts.out)
    existing = opts.resume && isfile(opts.out) && !opts.overwrite ? ScanCSV.read_existing_keys(opts.out, ["T_MeV", "muB_MeV", "xi"]) : Set{Tuple{Float64,Float64,Float64}}()
    if opts.overwrite && isfile(opts.out)
        rm(opts.out)
    end

    new_file = !isfile(opts.out)
    io = open(opts.out, "a")
    try
        if new_file
            ScanCSV.write_metadata(io, Dict(
                "schema" => "scan_csv_v1",
                "title" => "relaxation_times_vs_T",
                "script" => "scripts/relaxtime/scan_relaxation_times_vs_T.jl",
                # plot hints (consumed by scripts/plot_scan_csv.py)
                "x_label" => "T",
                "x_unit" => "MeV",
                "y_unit.tau_u" => "fm",
                "y_unit.tau_s" => "fm",
                "y_unit.tau_ubar" => "fm",
                "y_unit.tau_sbar" => "fm",
                "y_unit.tauinv_u" => "1/fm",
                "y_unit.tauinv_s" => "1/fm",
                "y_unit.tauinv_ubar" => "1/fm",
                "y_unit.tauinv_sbar" => "1/fm",
                # relaxation times span many orders; log often reads better
                "y_scale.tau_u" => "log",
                "y_scale.tau_s" => "log",
                "y_scale.tau_ubar" => "log",
                "y_scale.tau_sbar" => "log",
            ))
            write_header(io)
        end

        T_vals = collect(range(opts.tmin_mev, opts.tmax_mev; step=opts.tstep_mev))
        for muB_mev in opts.mub_list_mev
            muq_mev = muB_mev / 3.0
            muq_fm = muq_mev / ħc_MeV_fm

            seed_state = DEFAULT_MU_GUESS
            for T_mev in T_vals
                if opts.resume && !opts.overwrite
                    if (T_mev, muB_mev, opts.xi) in existing
                        continue
                    end
                end
                T_fm = T_mev / ħc_MeV_fm

                base = solve_equilibrium_mu(
                    T_fm,
                    muq_fm;
                    xi=opts.xi,
                    seed_state=seed_state,
                    p_num=opts.p_num,
                    t_num=opts.t_num,
                    iterations=opts.max_iter,
                )
                seed_state = base.x_state

                Φ = Float64(base.x_state[4])
                Φbar = Float64(base.x_state[5])
                masses_vec = calculate_mass_vec(SVector{3}(base.x_state[1], base.x_state[2], base.x_state[3]))
                masses = (u=Float64(masses_vec[1]), d=Float64(masses_vec[2]), s=Float64(masses_vec[3]))

                thermo_params = (T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(opts.xi))
                quark_params_basic = (m=masses, μ=(u=Float64(muq_fm), d=Float64(muq_fm), s=Float64(muq_fm)))

                ktmp = build_K_coeffs(Float64(T_fm), Float64(muq_fm), masses, Φ, Φbar)
                quark_params = (m=quark_params_basic.m, μ=quark_params_basic.μ, A=ktmp.A_vals)
                K_coeffs = ktmp.K_coeffs

                thermal_nodes = cached_nodes(opts.p_num, opts.t_num)
                densities = densities_from_equilibrium(base.x_state, base.mu_vec, T_fm, thermal_nodes, Float64(opts.xi))

                tau = (u=NaN, s=NaN, ubar=NaN, sbar=NaN)
                tauinv = (u=NaN, s=NaN, ubar=NaN, sbar=NaN)
                if Bool(base.converged)
                    cs_caches = build_sigma_caches(
                        REQUIRED_PROCESSES,
                        quark_params,
                        thermo_params,
                        K_coeffs;
                        n_sigma_points=opts.tau_n_sigma_points,
                        sigma_grid_n=opts.sigma_grid_n,
                        p_cut_factor=opts.p_cut_factor,
                    )

                    tau_res = relaxation_times(
                        quark_params,
                        thermo_params,
                        K_coeffs;
                        densities=densities,
                        cs_caches=cs_caches,
                        p_nodes=opts.tau_p_nodes,
                        angle_nodes=opts.tau_angle_nodes,
                        phi_nodes=opts.tau_phi_nodes,
                        n_sigma_points=opts.tau_n_sigma_points,
                    )

                    tau = tau_res.tau
                    tauinv = tau_res.tau_inv
                end

                row = Any[
                    T_mev, muB_mev, muq_mev, opts.xi,
                    T_fm, muq_fm,
                    csv_bool(Bool(base.converged)), Int(base.iterations), Float64(base.residual_norm),
                    Φ, Φbar,
                    masses.u, masses.s,
                    densities.u, densities.s, densities.ubar, densities.sbar,
                    tau.u, tau.s, tau.ubar, tau.sbar,
                    tauinv.u, tauinv.s, tauinv.ubar, tauinv.sbar,
                    opts.tau_p_nodes, opts.tau_angle_nodes, opts.tau_phi_nodes, opts.tau_n_sigma_points,
                    opts.sigma_grid_n, opts.p_cut_factor,
                ]
                println(io, join(row, ','))
                flush(io)

                if Bool(base.converged)
                    @printf("T=%6.1f MeV  muB=%6.1f MeV | tau_u=%.3e tau_s=%.3e tau_ubar=%.3e tau_sbar=%.3e\n",
                        T_mev, muB_mev, tau.u, tau.s, tau.ubar, tau.sbar)
                else
                    @printf("T=%6.1f MeV  muB=%6.1f MeV | equilibrium NOT converged; tau=NaN\n", T_mev, muB_mev)
                end
            end
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
