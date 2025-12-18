#!/usr/bin/env julia

raw"""
固定 T=150 MeV、μ_B=800 MeV，扫描各向异性参数 ξ 下的 PNJL 平衡解 + τ + RTA 输运系数。

输出：scan_csv_v1（带 # 元数据头）的 CSV。

默认输出列与 scripts/relaxtime/run_gap_transport_scan.jl 保持一致，便于复用绘图/分析脚本。

示例：
  julia --threads 8 --project=. scripts/relaxtime/scan_transport_vs_xi_T150_muB800.jl \
    --xi-min 0 --xi-max 1.0 --xi-step 0.1 --overwrite

  # 指定离散 ξ
  julia --threads 8 --project=. scripts/relaxtime/scan_transport_vs_xi_T150_muB800.jl \
    --xi-list 0,0.2,0.4,0.6,0.8,1.0 --overwrite

建议绘图：
  python scripts/plot_scan_csv.py --mode lines \
    --csv data/processed/results/relaxtime/transport_vs_xi_T150_muB800.csv \
    --x xi --ys eta,sigma,zeta,tau_u,tau_s --out-dir data/processed/figures/relaxtime/transport_vs_xi
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using StaticArrays
using Dates
using .ScanCSV: ScanCSV
using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS

struct Options
    output::String
    T_mev::Float64
    muB_mev::Float64
    xi_values::Vector{Float64}
    overwrite::Bool
    resume::Bool
    ignore_lock::Bool

    # equilibrium solver
    p_num::Int
    t_num::Int
    max_iter::Int

    # tau / cross-section settings
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int

    # transport settings
    tr_p_nodes::Int
    tr_p_max_fm::Float64

    # bulk viscosity (slow)
    compute_bulk::Bool
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/scan_transport_vs_xi_T150_muB800.jl [options]\n")
    println("Options:")
    println("  --output <path>             输出 CSV (default data/processed/results/relaxtime/transport_vs_xi_T150_muB800.csv)")
    println("  --out <path>                同 --output")
    println("  --T-MeV <value>             Temperature in MeV (default 150)")
    println("  --muB-MeV <value>           Baryon chemical potential μ_B in MeV (default 800)")
    println("  --xi <value>                追加一个 ξ 值（可多次传入）")
    println("  --xi-list v1,v2,...         用逗号分隔的 ξ 列表替换")
    println("  --xi-min/--xi-max/--xi-step 扫描 ξ 的区间与步长（默认 0..1 step 0.1）")
    println("  --overwrite                 覆盖输出文件")
    println("  --no-resume                 禁用跳过逻辑，强制重算")
    println("  --ignore-lock               忽略输出锁文件（不建议；用于确认没有重复进程时手动解锁）")
    println("  --p-num <int>               能隙/密度的动量节点数 (default 12)")
    println("  --t-num <int>               能隙/密度的角度节点数 (default 6)")
    println("  --max-iter <int>            NLsolve iterations 上限 (default 40)")
    println("  --tau-p-nodes <int>         τ 平均散射率动量节点 (default 6)")
    println("  --tau-angle-nodes <int>     τ 平均散射率 cosθ 节点 (default 2)")
    println("  --tau-phi-nodes <int>       τ 平均散射率 φ 节点 (default 4)")
    println("  --tau-n-sigma <int>         截面 t 积分点数 (default 8)")
    println("  --tr-p-nodes <int>          输运积分动量节点 (default 24)")
    println("  --tr-p-max <fm^-1>          输运积分 p 上限 (default 8.0)")
    println("  --compute-bulk              计算体粘滞 ζ（很慢；默认关闭）")
    println("  -h, --help                  显示帮助")
end

function parse_args(args::Vector{String})::Options
    opts = Dict{Symbol,Any}(
        :output => joinpath("data", "processed", "results", "relaxtime", "transport_vs_xi_T150_muB800.csv"),
        :T_mev => 150.0,
        :muB_mev => 800.0,
        :xi_values => Float64[],
        :xi_min => 0.0,
        :xi_max => 1.0,
        :xi_step => 0.1,
        :overwrite => false,
        :resume => true,
        :ignore_lock => false,
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :tau_p_nodes => 6,
        :tau_angle_nodes => 2,
        :tau_phi_nodes => 4,
        :tau_n_sigma => 8,
        :tr_p_nodes => 24,
        :tr_p_max => 8.0,
        :compute_bulk => false,
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

        if arg == "--output" || arg == "--out"
            opts[:output] = require_value()
        elseif arg == "--T-MeV"
            opts[:T_mev] = parse(Float64, require_value())
        elseif arg == "--muB-MeV"
            opts[:muB_mev] = parse(Float64, require_value())
        elseif arg == "--xi"
            push!(opts[:xi_values], parse(Float64, require_value()))
        elseif arg == "--xi-list"
            raw = split(require_value(), ',')
            vals = Float64[parse(Float64, strip(v)) for v in raw if !isempty(strip(v))]
            opts[:xi_values] = vals
        elseif arg == "--xi-min"
            opts[:xi_min] = parse(Float64, require_value())
        elseif arg == "--xi-max"
            opts[:xi_max] = parse(Float64, require_value())
        elseif arg == "--xi-step"
            opts[:xi_step] = parse(Float64, require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--no-resume"
            opts[:resume] = false
        elseif arg == "--ignore-lock"
            opts[:ignore_lock] = true
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
        elseif arg == "--tr-p-nodes"
            opts[:tr_p_nodes] = parse(Int, require_value())
        elseif arg == "--tr-p-max"
            opts[:tr_p_max] = parse(Float64, require_value())
        elseif arg == "--compute-bulk"
            opts[:compute_bulk] = true
        elseif arg in ("-h", "--help")
            print_usage(); exit(0)
        else
            error("unknown option: $arg")
        end

        i += 1
    end

    xi_vals = Float64.(opts[:xi_values])
    if isempty(xi_vals)
        xmin = Float64(opts[:xi_min])
        xmax = Float64(opts[:xi_max])
        xstep = Float64(opts[:xi_step])
        xstep > 0 || error("xi-step must be positive")
        xmax >= xmin || error("xi-max must be >= xi-min")
        xi_vals = collect(range(xmin; stop=xmax, step=xstep))
    end
    xi_vals = unique(sort(xi_vals))

    return Options(
        String(opts[:output]),
        Float64(opts[:T_mev]),
        Float64(opts[:muB_mev]),
        xi_vals,
        Bool(opts[:overwrite]),
        Bool(opts[:resume]),
        Bool(opts[:ignore_lock]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Int(opts[:tau_p_nodes]),
        Int(opts[:tau_angle_nodes]),
        Int(opts[:tau_phi_nodes]),
        Int(opts[:tau_n_sigma]),
        Int(opts[:tr_p_nodes]),
        Float64(opts[:tr_p_max]),
        Bool(opts[:compute_bulk]),
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function acquire_output_lock(output_path::AbstractString; ignore_lock::Bool)
    lock_path = output_path * ".lock"
    if isfile(lock_path) && !ignore_lock
        error("Output appears to be locked (another run may be in progress): $(lock_path). " *
              "If you are sure no other process is running, delete the lock file or pass --ignore-lock.")
    end
    open(lock_path, "w") do io
        println(io, "pid=$(getpid())")
        println(io, "time=$(Dates.now())")
        println(io, "output=$(abspath(output_path))")
    end
    return lock_path
end

function write_header(io)
    ScanCSV.write_header(io, [
        "T_MeV", "muq_MeV", "muB_MeV", "xi",
        "T_fm", "muq_fm",
        "converged", "iterations", "residual_norm",
        "Phi", "Phibar",
        "m_u", "m_d", "m_s",
        "rho_baryon", "rho_norm",
        "n_u", "n_d", "n_s", "n_ubar", "n_dbar", "n_sbar",
        "tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
        "tauinv_u", "tauinv_d", "tauinv_s", "tauinv_ubar", "tauinv_dbar", "tauinv_sbar",
        "eta", "sigma", "zeta",
    ])
end

function csv_bool(x::Bool)
    return x ? "true" : "false"
end

function compute_K_coeffs(T_fm::Float64, mu_fm::Float64, masses::NamedTuple, Φ::Float64, Φbar::Float64)
    nodes = DEFAULT_MOMENTUM_NODES
    weights = DEFAULT_MOMENTUM_WEIGHTS
    A_u = A(masses.u, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    A_s = A(masses.s, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
end

function run_scan(opts::Options)
    ensure_parent_dir(opts.output)

    lock_path = acquire_output_lock(opts.output; ignore_lock=opts.ignore_lock)

    if opts.overwrite && isfile(opts.output)
        rm(opts.output)
    end

    existing = (opts.resume && isfile(opts.output) && !opts.overwrite) ? ScanCSV.read_existing_keys(opts.output, ["T_MeV", "muB_MeV", "xi"]) : Set{Tuple{Float64,Float64,Float64}}()

    new_file = !isfile(opts.output)
    io = open(opts.output, "a")
    try
        if new_file
            ScanCSV.write_metadata(io, Dict(
                "schema" => "scan_csv_v1",
                "title" => "transport_vs_xi",
                "script" => "scripts/relaxtime/scan_transport_vs_xi_T150_muB800.jl",
                "x_label" => "ξ",
                "x_label.xi" => "ξ",
                "y_label.eta" => "η",
                "y_label.sigma" => "σ",
                "y_label.zeta" => "ζ",
                "y_label.tau_u" => "τ_u",
                "y_label.tau_s" => "τ_s",
                "y_scale.tau_u" => "log",
                "y_scale.tau_s" => "log",
            ))
            write_header(io)
        end

        T_mev = opts.T_mev
        muB_mev = opts.muB_mev
        muq_mev = muB_mev / 3.0

        T_fm = T_mev / ħc_MeV_fm
        muq_fm = muq_mev / ħc_MeV_fm

        seed_state = nothing

        total = length(opts.xi_values)
        done = 0

        for xi in opts.xi_values
            done += 1
            key = (T_mev, muB_mev, xi)
            if opts.resume && (key in existing)
                continue
            end

            base = TransportWorkflow.ThermoDerivatives.solve_equilibrium_mu(
                T_fm,
                muq_fm;
                xi=xi,
                seed_state=(seed_state === nothing ? TransportWorkflow.AnisoGapSolver.DEFAULT_MU_GUESS : seed_state),
                p_num=opts.p_num,
                t_num=opts.t_num,
                iterations=opts.max_iter,
            )
            seed_state = base.x_state

            Φ = Float64(base.x_state[4])
            Φbar = Float64(base.x_state[5])
            φ = SVector{3}(base.x_state[1], base.x_state[2], base.x_state[3])
            mvec = TransportWorkflow.AnisoGapSolver.calculate_mass_vec(φ)
            masses = (u=Float64(mvec[1]), d=Float64(mvec[2]), s=Float64(mvec[3]))

            K_coeffs = compute_K_coeffs(Float64(T_fm), Float64(muq_fm), masses, Φ, Φbar)

            res = solve_gap_and_transport(
                T_fm,
                muq_fm;
                xi=xi,
                equilibrium=base,
                compute_tau=true,
                K_coeffs=K_coeffs,
                tau=nothing,
                compute_bulk=opts.compute_bulk,
                p_num=opts.p_num,
                t_num=opts.t_num,
                seed_state=seed_state,
                solver_kwargs=(iterations=opts.max_iter,),
                tau_kwargs=(
                    p_nodes=opts.tau_p_nodes,
                    angle_nodes=opts.tau_angle_nodes,
                    phi_nodes=opts.tau_phi_nodes,
                    n_sigma_points=opts.tau_n_sigma_points,
                ),
                transport_kwargs=(p_nodes=opts.tr_p_nodes, p_max=opts.tr_p_max_fm,),
            )

            eq = res.equilibrium
            dens = res.densities
            tau = res.tau
            tauinv = res.tau_inv
            tr = res.transport

            row = join([
                string(T_mev), string(muq_mev), string(muB_mev), string(xi),
                string(T_fm), string(muq_fm),
                csv_bool(eq.converged), string(eq.iterations), string(eq.residual_norm),
                string(Φ), string(Φbar),
                string(masses.u), string(masses.d), string(masses.s),
                string(eq.rho), string(eq.rho_norm),
                string(dens.u), string(dens.d), string(dens.s), string(dens.ubar), string(dens.dbar), string(dens.sbar),
                string(tau.u), string(tau.d), string(tau.s), string(tau.ubar), string(tau.dbar), string(tau.sbar),
                string(tauinv.u), string(tauinv.d), string(tauinv.s), string(tauinv.ubar), string(tauinv.dbar), string(tauinv.sbar),
                string(tr.eta), string(tr.sigma), string(tr.zeta),
            ], ',')
            println(io, row)
            flush(io)

            if eq.converged
                seed_state = eq.x_state
            end

            if done % 5 == 0
                println("progress: $(done)/$(total) (xi=$(xi))")
            end
        end
    finally
        close(io)
        try
            isfile(lock_path) && rm(lock_path; force=true)
        catch err
            @warn "failed to remove lock file" lock_path=lock_path err=err
        end
    end

    println("Scan finished. Output: $(opts.output)")
end

function main()
    opts = parse_args(copy(ARGS))
    run_scan(opts)
end

main()
