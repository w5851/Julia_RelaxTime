"""
批量扫描 (T, μ) 网格，串联：各向异性 PNJL 平衡求解 → τ → RTA 输运系数。

输出：CSV（默认 data/outputs/results/relaxtime/gap_transport_scan.csv）

单位约定：
- CLI 的温度/化学势默认使用 MeV（更符合扫描习惯）；脚本内部会换算到 fm⁻¹。
- 输出同时包含 MeV 与 fm⁻¹ 的关键量。

示例：
  julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --tmin 120 --tmax 400 --tstep 10 --mubmin 0 --mubmax 800 --mubstep 800 --xi 0.0 --mode finite_15 --compute-bulk --overwrite

注意：
- compute_bulk 默认关闭（体粘滞需要多次自动微分+求解，扫描会很慢）。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm, ρ0_inv_fm3
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS, gauleg
using StaticArrays
using .ScanCSV: ScanCSV

const RT_ASR = TransportWorkflow.RelaxationTime.AverageScatteringRate
const RT_TCS = TransportWorkflow.RelaxationTime.TotalCrossSection
const REQUIRED_PROCESSES = TransportWorkflow.RelaxationTime.REQUIRED_PROCESSES

const MODULE_DEFAULT_P_NODES = RT_ASR.DEFAULT_P_NODES           # 20
const MODULE_DEFAULT_ANGLE_NODES = RT_ASR.DEFAULT_ANGLE_NODES   # 4
const MODULE_DEFAULT_PHI_NODES = RT_ASR.DEFAULT_PHI_NODES       # 8
const MODULE_DEFAULT_SIGMA_GRID_N = RT_ASR.DEFAULT_SIGMA_GRID_N # 60
struct ScanOptions
    output::String
    xi_values::Vector{Float64}
    tmin_mev::Float64
    tmax_mev::Float64
    tstep_mev::Float64
    mubmin_mev::Float64
    mubmax_mev::Float64
    mubstep_mev::Float64
    overwrite::Bool
    resume::Bool
    compute_bulk::Bool
    p_num::Int
    t_num::Int
    max_iter::Int
    # tau / cross-section settings
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
    sigma_grid_n::Int
    integration_mode::Symbol
    gc_every_n::Int
    # transport settings
    tr_p_nodes::Int
    tr_p_max_fm::Float64
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_gap_transport_scan.jl [options]\n")
    println("Options:")
    println("  --output <path>             输出 CSV (default data/outputs/results/relaxtime/gap_transport_scan.csv)")
    println("  --xi <value>                追加一个 ξ 值（可多次传入）")
    println("  --xi-list v1,v2,...         用逗号分隔的 ξ 列表替换")
    println("  --tmin/--tmax/--tstep <MeV> 温度范围与步长")
    println("  --mubmin/--mubmax/--mubstep <MeV> 重子化学势 μ_B 范围与步长")
    println("  --mumin/--mumax/--mustep <MeV> (兼容旧参数) 夸克化学势 μ_q 范围与步长")
    println("  --overwrite                 覆盖输出文件")
    println("  --no-resume                 禁用跳过逻辑，强制重算")
    println("  --compute-bulk              计算体粘滞 ζ（很慢；默认关闭）")
    println("  --p-num <int>               能隙/密度的动量节点数 (default 12)")
    println("  --t-num <int>               能隙/密度的角度节点数 (default 6)")
    println("  --max-iter <int>            NLsolve iterations 上限 (default 40)")
    println("  --tau-p-nodes <int>         τ 平均散射率动量节点 (default $(MODULE_DEFAULT_P_NODES))")
    println("  --tau-angle-nodes <int>     τ 平均散射率 cosθ 节点 (default $(MODULE_DEFAULT_ANGLE_NODES))")
    println("  --tau-phi-nodes <int>       τ 平均散射率 φ 节点 (default $(MODULE_DEFAULT_PHI_NODES))")
    println("  --tau-n-sigma <int>         σ(s) 的 t 积分点数 (default $(RT_TCS.DEFAULT_T_INTEGRAL_POINTS))")
    println("  --sigma-grid-n <int>        σ(s) 预计算网格点数 (default $(MODULE_DEFAULT_SIGMA_GRID_N))")
    println("  --mode <mode>               τ 积分模式: semi_infinite | finite_15 | finite_lambda (default semi_infinite)")
    println("  --gc-every-n <int>          每 N 个点触发一次 GC (default 5; 0 表示关闭)")
    println("  --tr-p-nodes <int>          输运积分动量节点 (default 24)")
    println("  --tr-p-max <fm^-1>          输运积分 p 上限 (default 8.0)")
    println("  -h, --help                  显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :output => joinpath("data", "outputs", "results", "relaxtime", "gap_transport_scan.csv"),
        :xi_values => Float64[0.0],
        :tmin => 50.0,
        :tmax => 200.0,
        :tstep => 10.0,
        :mubmin => 0.0,
        :mubmax => 1200.0,
        :mubstep => 60.0,
        :overwrite => false,
        :resume => true,
        :compute_bulk => false,
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :tau_p_nodes => MODULE_DEFAULT_P_NODES,
        :tau_angle_nodes => MODULE_DEFAULT_ANGLE_NODES,
        :tau_phi_nodes => MODULE_DEFAULT_PHI_NODES,
        :tau_n_sigma => RT_TCS.DEFAULT_T_INTEGRAL_POINTS,
        :sigma_grid_n => MODULE_DEFAULT_SIGMA_GRID_N,
        :integration_mode => :semi_infinite,
        :gc_every_n => 5,
        :tr_p_nodes => 24,
        :tr_p_max => 8.0,
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
        if arg == "--output"
            opts[:output] = require_value()
        elseif arg == "--xi"
            val = parse(Float64, require_value())
            if opts[:xi_values] == Float64[0.0]
                opts[:xi_values] = Float64[]
            end
            push!(opts[:xi_values], val)
        elseif arg == "--xi-list"
            raw = split(require_value(), ',')
            vals = Float64[parse(Float64, strip(v)) for v in raw if !isempty(strip(v))]
            opts[:xi_values] = vals
        elseif arg == "--tmin"
            opts[:tmin] = parse(Float64, require_value())
        elseif arg == "--tmax"
            opts[:tmax] = parse(Float64, require_value())
        elseif arg == "--tstep"
            opts[:tstep] = parse(Float64, require_value())
        elseif arg == "--mubmin"
            opts[:mubmin] = parse(Float64, require_value())
        elseif arg == "--mubmax"
            opts[:mubmax] = parse(Float64, require_value())
        elseif arg == "--mubstep"
            opts[:mubstep] = parse(Float64, require_value())
        # 兼容旧参数：把 μ_q 转成 μ_B
        elseif arg == "--mumin"
            opts[:mubmin] = 3.0 * parse(Float64, require_value())
        elseif arg == "--mumax"
            opts[:mubmax] = 3.0 * parse(Float64, require_value())
        elseif arg == "--mustep"
            opts[:mubstep] = 3.0 * parse(Float64, require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--no-resume"
            opts[:resume] = false
        elseif arg == "--compute-bulk"
            opts[:compute_bulk] = true
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
        elseif arg == "--mode"
            mode_str = require_value()
            mode_sym = Symbol(mode_str)
            mode_sym in (:semi_infinite, :finite_15, :finite_lambda) || error("unknown mode: $mode_str (expected: semi_infinite, finite_15, finite_lambda)")
            opts[:integration_mode] = mode_sym
        elseif arg == "--gc-every-n"
            opts[:gc_every_n] = parse(Int, require_value())
        elseif arg == "--tr-p-nodes"
            opts[:tr_p_nodes] = parse(Int, require_value())
        elseif arg == "--tr-p-max"
            opts[:tr_p_max] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage(); exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    xi_vals = opts[:xi_values]
    isempty(xi_vals) && (xi_vals = Float64[0.0])
    xi_vals = unique(sort(Float64.(xi_vals)))

    tstep = Float64(opts[:tstep]); tstep > 0 || error("tstep must be positive")
    mubstep = Float64(opts[:mubstep]); mubstep > 0 || error("mubstep must be positive")

    return ScanOptions(
        String(opts[:output]),
        xi_vals,
        Float64(opts[:tmin]),
        Float64(opts[:tmax]),
        tstep,
        Float64(opts[:mubmin]),
        Float64(opts[:mubmax]),
        mubstep,
        Bool(opts[:overwrite]),
        Bool(opts[:resume]),
        Bool(opts[:compute_bulk]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Int(opts[:tau_p_nodes]),
        Int(opts[:tau_angle_nodes]),
        Int(opts[:tau_phi_nodes]),
        Int(opts[:tau_n_sigma]),
        Int(opts[:sigma_grid_n]),
        Symbol(opts[:integration_mode]),
        Int(opts[:gc_every_n]),
        Int(opts[:tr_p_nodes]),
        Float64(opts[:tr_p_max]),
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function read_existing_keys(path::AbstractString)
    return ScanCSV.read_existing_keys(path, ["T_MeV", "muB_MeV", "xi"])
end

function ensure_output_header_compatible(path::AbstractString)
    isfile(path) || return
    header = nothing
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, "#")) && continue
            header = s
            break
        end
    end
    header === nothing && return
    required = (
        "omega_fm4inv",
        "P_fm4inv",
        "epsilon_fm4inv",
        "s_fm3inv",
        "eps_minus_3P_over_T4",
        "eta_over_s",
        "zeta_over_s",
    )
    for c in required
        occursin(c, header) || error("existing output CSV header is incompatible with current script (missing column: $c). Please rerun with --overwrite or choose a new --output path.")
    end
end

function write_header_if_needed(io)
    header = join([
        "T_MeV", "muq_MeV", "muB_MeV", "xi",
        "T_fm", "muq_fm",
        "converged", "iterations", "residual_norm",
        "Phi", "Phibar",
        "m_u", "m_d", "m_s",
        "rho_baryon", "rho_norm",
        "omega_fm4inv", "P_fm4inv", "epsilon_fm4inv", "s_fm3inv",
        "omega_MeV_fm3", "P_MeV_fm3", "epsilon_MeV_fm3",
        "eps_minus_3P_over_T4",
        "n_u", "n_d", "n_s", "n_ubar", "n_dbar", "n_sbar",
        "tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
        "tauinv_u", "tauinv_d", "tauinv_s", "tauinv_ubar", "tauinv_dbar", "tauinv_sbar",
        "eta", "sigma", "zeta", "eta_over_s", "zeta_over_s",
    ], ',')
    println(io, header)
end

function csv_bool(x::Bool)
    return x ? "true" : "false"
end

function build_K_data(T_fm::Float64, mu_fm::Float64, masses::NamedTuple, Φ::Float64, Φbar::Float64)
    nodes = DEFAULT_MOMENTUM_NODES
    weights = DEFAULT_MOMENTUM_WEIGHTS
    A_u = A(masses.u, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    A_s = A(masses.s, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return (K_coeffs=calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s), A_vals=(u=A_u, d=A_u, s=A_s))
end

function integration_grids(opts::ScanOptions)
    if opts.integration_mode == :finite_15
        pg, pw = gauleg(0.0, 15.0, opts.tau_p_nodes)
        return (pg, pw, Λ_inv_fm)
    elseif opts.integration_mode == :finite_lambda
        pg, pw = gauleg(0.0, Λ_inv_fm, opts.tau_p_nodes)
        return (pg, pw, Λ_inv_fm)
    else
        return (nothing, nothing, nothing)
    end
end

function safe_total_cross_section(process::Symbol, s::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int, max_tries::Int=4)
    s_try = s
    last_err = nothing
    for _ in 1:max_tries
        try
            σ = RT_TCS.total_cross_section(process, s_try, quark_params, thermo_params, K_coeffs; n_points=n_points)
            isfinite(σ) && return σ
        catch err
            last_err = err
        end
        s_try = s_try * (1.0 + 1e-6) + 1e-10
    end
    @warn "failed to compute sigma; returning NaN" process=process s=s last_error=last_err
    return NaN
end

function build_sigma_caches(processes::Tuple, quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_sigma_points::Int, sigma_grid_n::Int)
    caches = Dict{Symbol,RT_ASR.CrossSectionCache}()
    for process in processes
        s_grid = RT_ASR.design_w0cdf_s_grid(process, quark_params, thermo_params;
            N=sigma_grid_n, p_cutoff=Λ_inv_fm)

        cache = RT_ASR.CrossSectionCache(process)
        n_ok, n_bad = 0, 0
        for s in s_grid
            σ = safe_total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
            if !isfinite(σ)
                n_bad += 1
                continue
            end
            RT_ASR.insert_sigma!(cache, s, σ)
            n_ok += 1
        end
        n_bad > 0 && @warn "sigma grid had non-finite points" process=process n_ok=n_ok n_bad=n_bad
        n_ok >= 2 || error("sigma cache has too few valid points for $process (n_ok=$n_ok)")
        caches[process] = cache
    end
    return caches
end

function run_scan(opts::ScanOptions)
    ensure_parent_dir(opts.output)

    if opts.resume && isfile(opts.output) && !opts.overwrite
        ensure_output_header_compatible(opts.output)
    end

    existing = opts.resume && isfile(opts.output) && !opts.overwrite ? read_existing_keys(opts.output) : Set{Tuple{Float64,Float64,Float64}}()

    if opts.overwrite && isfile(opts.output)
        rm(opts.output)
    end

    new_file = !isfile(opts.output) || filesize(opts.output) == 0
    p_grid, p_w, sigma_cutoff = integration_grids(opts)
    cos_grid, cos_w = gauleg(-1.0, 1.0, opts.tau_angle_nodes)
    phi_grid, phi_w = gauleg(0.0, 2 * pi, opts.tau_phi_nodes)
    solver_kwargs = (iterations=opts.max_iter,)
    io = open(opts.output, "a")
    try
        if new_file
            ScanCSV.write_metadata(io, Dict(
                "schema" => "scan_csv_v1",
                "title" => "gap_transport_scan",
                "script" => "scripts/relaxtime/run_gap_transport_scan.jl",
                "sigma_grid_n" => string(opts.sigma_grid_n),
                "integration_mode" => string(opts.integration_mode),
                "gc_every_n" => string(opts.gc_every_n),
                "tau_p_nodes" => string(opts.tau_p_nodes),
                "tau_angle_nodes" => string(opts.tau_angle_nodes),
                "tau_phi_nodes" => string(opts.tau_phi_nodes),
                "tau_n_sigma_points" => string(opts.tau_n_sigma_points),
                "tr_p_nodes" => string(opts.tr_p_nodes),
                "tr_p_max_fm" => string(opts.tr_p_max_fm),
            ))
            write_header_if_needed(io)
        end

        T_values = collect(range(opts.tmin_mev; stop=opts.tmax_mev, step=opts.tstep_mev))
        muB_values = collect(range(opts.mubmin_mev; stop=opts.mubmax_mev, step=opts.mubstep_mev))
        muB_values = unique(sort(Float64.(muB_values)))

        total = length(opts.xi_values) * length(T_values) * length(muB_values)
        done = 0

        for xi in opts.xi_values
            for muB_mev in muB_values
                # 固定同一个 μ_B：按温度递增扫描，输出自然按 T 排序
                seed_state = nothing
                # 初值策略：首点用 MultiSeed 选最低 Ω；后续点用相变感知的连续性跟踪，避免跨一阶相变线时跳到错误分支。
                phase_tracker = try
                    TransportWorkflow.PNJL.PhaseAwareContinuitySeed(xi)
                catch
                    TransportWorkflow.PNJL.PhaseAwareContinuitySeed()
                end
                for T_mev in T_values
                    done += 1
                    key = (T_mev, muB_mev, xi)
                    if opts.resume && (key in existing)
                        continue
                    end

                    T_fm = T_mev / ħc_MeV_fm
                    muq_mev = muB_mev / 3.0
                    muq_fm = muq_mev / ħc_MeV_fm

                    # 先求一次平衡解（用于 K_coeffs，并作为后续 workflow 的 equilibrium 复用）
                    # 首点：MultiSeed；后续：相变感知的连续性跟踪（phase_tracker 会用上一次解作为种子，并在跨相变线时切换默认相的种子）。
                    seed_strategy = (phase_tracker.previous_solution === nothing) ? TransportWorkflow.PNJL.MultiSeed() : phase_tracker
                    base = TransportWorkflow.PNJL.solve(TransportWorkflow.PNJL.FixedMu(), T_fm, muq_fm;
                        xi=xi,
                        p_num=opts.p_num,
                        t_num=opts.t_num,
                        seed_strategy=seed_strategy,
                        solver_kwargs...,
                    )

                    if base.converged
                        # PhaseAwareContinuitySeed 的 update! 需要 MeV 单位（与其内部 boundary.csv 数据一致）
                        TransportWorkflow.PNJL.update!(phase_tracker, base.solution, T_mev, muq_mev)
                    end

                    seed_state = Vector(base.solution)

                    Φ = Float64(base.x_state[4])
                    Φbar = Float64(base.x_state[5])
                    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

                    ktmp = build_K_data(Float64(T_fm), Float64(muq_fm), masses, Φ, Φbar)
                    thermo_params = (T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(xi))
                    quark_params = (m=masses, μ=(u=Float64(muq_fm), d=Float64(muq_fm), s=Float64(muq_fm)), A=ktmp.A_vals)

                    cs_caches = nothing
                    if Bool(base.converged)
                        cs_caches = build_sigma_caches(REQUIRED_PROCESSES, quark_params, thermo_params, ktmp.K_coeffs;
                            n_sigma_points=opts.tau_n_sigma_points, sigma_grid_n=opts.sigma_grid_n)
                    end

                    # 正式计算：τ + 输运（η/σ）; ζ 默认关
                    res = solve_gap_and_transport(
                        T_fm,
                        muq_fm;
                        xi=xi,
                        equilibrium=base,
                        compute_tau=true,
                        K_coeffs=ktmp.K_coeffs,
                        tau=nothing,
                        compute_bulk=opts.compute_bulk,
                        p_num=opts.p_num,
                        t_num=opts.t_num,
                        seed_state=seed_state,
                        solver_kwargs=solver_kwargs,
                        tau_kwargs=(
                            p_nodes=opts.tau_p_nodes,
                            angle_nodes=opts.tau_angle_nodes,
                            phi_nodes=opts.tau_phi_nodes,
                            n_sigma_points=opts.tau_n_sigma_points,
                            cs_caches=cs_caches,
                            p_grid=p_grid,
                            p_w=p_w,
                            cos_grid=cos_grid,
                            cos_w=cos_w,
                            phi_grid=phi_grid,
                            phi_w=phi_w,
                            sigma_cutoff=sigma_cutoff,
                        ),
                        transport_kwargs=(p_nodes=opts.tr_p_nodes, p_max=opts.tr_p_max_fm,),
                    )

                    eq = res.equilibrium
                    dens = res.densities
                    tau = res.tau
                    tauinv = res.tau_inv
                    tr = res.transport

                    # 重建重子数密度（旧版 eq.rho/eq.rho_norm 已移除）
                    rho_quark_net = (dens.u - dens.ubar) + (dens.d - dens.dbar) + (dens.s - dens.sbar)
                    rho_baryon = rho_quark_net / 3.0
                    rho_norm = rho_baryon / ρ0_inv_fm3

                    omega_fm4inv = eq.omega
                    P_fm4inv = eq.pressure
                    epsilon_fm4inv = eq.energy
                    s_fm3inv = eq.entropy
                    omega_MeV_fm3 = omega_fm4inv * ħc_MeV_fm
                    P_MeV_fm3 = P_fm4inv * ħc_MeV_fm
                    epsilon_MeV_fm3 = epsilon_fm4inv * ħc_MeV_fm
                    eps_minus_3P_over_T4 = (isfinite(epsilon_fm4inv) && isfinite(P_fm4inv) && isfinite(T_fm) && T_fm != 0.0) ? ((epsilon_fm4inv - 3.0 * P_fm4inv) / T_fm^4) : NaN
                    eta_over_s = (isfinite(tr.eta) && isfinite(s_fm3inv) && s_fm3inv != 0.0) ? (tr.eta / s_fm3inv) : NaN
                    zeta_over_s = (isfinite(tr.zeta) && isfinite(s_fm3inv) && s_fm3inv != 0.0) ? (tr.zeta / s_fm3inv) : NaN

                    row = join([
                        string(T_mev), string(muq_mev), string(muB_mev), string(xi),
                        string(T_fm), string(muq_fm),
                        csv_bool(eq.converged), string(eq.iterations), string(eq.residual_norm),
                        string(Φ), string(Φbar),
                        string(masses.u), string(masses.d), string(masses.s),
                        string(rho_baryon), string(rho_norm),
                        string(omega_fm4inv), string(P_fm4inv), string(epsilon_fm4inv), string(s_fm3inv),
                        string(omega_MeV_fm3), string(P_MeV_fm3), string(epsilon_MeV_fm3),
                        string(eps_minus_3P_over_T4),
                        string(dens.u), string(dens.d), string(dens.s), string(dens.ubar), string(dens.dbar), string(dens.sbar),
                        string(tau.u), string(tau.d), string(tau.s), string(tau.ubar), string(tau.dbar), string(tau.sbar),
                        string(tauinv.u), string(tauinv.d), string(tauinv.s), string(tauinv.ubar), string(tauinv.dbar), string(tauinv.sbar),
                        string(tr.eta), string(tr.sigma), string(tr.zeta), string(eta_over_s), string(zeta_over_s),
                    ], ',')
                    println(io, row)
                    flush(io)

                    if eq.converged
                        seed_state = eq.x_state
                    end

                    if opts.gc_every_n > 0 && done % opts.gc_every_n == 0
                        GC.gc()
                    end

                    if done % 10 == 0
                        println("progress: $(done)/$(total) (last muB=$(muB_mev) MeV, T=$(T_mev) MeV, xi=$(xi))")
                    end
                end
            end
        end
    finally
        close(io)
    end

    println("Scan finished. Output: $(opts.output)")
end

function main()
    opts = parse_args(copy(ARGS))
    run_scan(opts)
end

main()
