"""
扫描弛豫时间 τ(T) 在给定 μ_B 下的温度依赖。

需求场景：输出 u, s, ubar, sbar 的 τ(T)（及 τ^{-1}），并支持 μ_B=0/800 MeV 两条曲线。

实现要点：
- 先解各向同性/各向异性 PNJL 平衡（能隙+Polyakov 环）
- 构造 TotalCrossSection 所需的 quark_params(A 字段) 与有效耦合 K_coeffs
- 用 RelaxationTime.relaxation_times 计算 τ
- 默认使用模块内部的默认参数（p_nodes=20, angle_nodes=4, phi_nodes=8, sigma_grid_n=60）
- 支持三种动量积分模式：
  - semi_infinite: 半无穷积分 [0, ∞)，默认行为
  - finite_15: 有限截断积分 [0, 15] fm⁻¹ + Λ截断的σ(s)范围（与C++一致）
  - finite_lambda: 有限截断积分 [0, Λ] fm⁻¹

单位约定：
- CLI 温度/化学势使用 MeV
- 内部计算使用 fm⁻¹
- 输出 τ 单位为 fm

示例：
  # 半无穷积分（默认）
  julia --project=. scripts/relaxtime/scan_relaxation_times_vs_T.jl --mode semi_infinite --tmin 100 --tmax 400 --tstep 10 --mub-list 0,800

  # 有限截断到15 fm⁻¹（与C++一致）
  julia --project=. scripts/relaxtime/scan_relaxation_times_vs_T.jl --mode finite_15 --tmin 100 --tmax 400 --tstep 10 --mub-list 0,800

  # 有限截断到Λ
  julia --project=. scripts/relaxtime/scan_relaxation_times_vs_T.jl --mode finite_lambda --tmin 100 --tmax 400 --tstep 10 --mub-list 0,800
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxationTime.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))

using Printf
using StaticArrays

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm
using .PNJL: solve, FixedMu, cached_nodes, calculate_number_densities
using .PNJL: HADRON_SEED_5
using .PNJL.Integrals: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .RelaxationTime: relaxation_times, REQUIRED_PROCESSES
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .ScanCSV: ScanCSV
using .GaussLegendre: gauleg

const RT_ASR = RelaxationTime.AverageScatteringRate
const RT_TCS = RelaxationTime.TotalCrossSection

# 从模块获取默认参数
const MODULE_DEFAULT_P_NODES = RT_ASR.DEFAULT_P_NODES           # 20
const MODULE_DEFAULT_ANGLE_NODES = RT_ASR.DEFAULT_ANGLE_NODES   # 4
const MODULE_DEFAULT_PHI_NODES = RT_ASR.DEFAULT_PHI_NODES       # 8
const MODULE_DEFAULT_SIGMA_GRID_N = RT_ASR.DEFAULT_SIGMA_GRID_N # 60

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
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
    sigma_grid_n::Int
    integration_mode::Symbol  # :semi_infinite, :finite_15, :finite_lambda
    gc_every_n::Int
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
    println("  --tau-p-nodes <int>          τ 平均散射率动量节点 (default $(MODULE_DEFAULT_P_NODES))")
    println("  --tau-angle-nodes <int>      τ 平均散射率 cosθ 节点 (default $(MODULE_DEFAULT_ANGLE_NODES))")
    println("  --tau-phi-nodes <int>        τ 平均散射率 φ 节点 (default $(MODULE_DEFAULT_PHI_NODES))")
    println("  --tau-n-sigma <int>          σ(s) 的 t 积分点数 (default $(RT_TCS.DEFAULT_T_INTEGRAL_POINTS))")
    println("  --sigma-grid-n <int>         每个过程预计算 σ(s) 的网格点数 (default $(MODULE_DEFAULT_SIGMA_GRID_N))")
    println("  --mode <mode>                积分模式: semi_infinite | finite_15 | finite_lambda (default semi_infinite)")
    println("  --gc-every-n <int>           每 N 个点触发一次 GC (default 5; 0 表示关闭)")
    println("  -h, --help                   显示帮助")
    println()
    println("积分模式说明:")
    println("  semi_infinite  : 半无穷积分 [0, ∞)，使用变换 p = scale*t/(1-t)")
    println("  finite_15      : 有限截断积分 [0, 15] fm⁻¹，σ(s)范围基于Λ截断（与C++一致）")
    println("  finite_lambda  : 有限截断积分 [0, Λ] fm⁻¹，Λ ≈ 3.05 fm⁻¹")
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
        :tau_p_nodes => MODULE_DEFAULT_P_NODES,
        :tau_angle_nodes => MODULE_DEFAULT_ANGLE_NODES,
        :tau_phi_nodes => MODULE_DEFAULT_PHI_NODES,
        :tau_n_sigma => RT_TCS.DEFAULT_T_INTEGRAL_POINTS,
        :sigma_grid_n => MODULE_DEFAULT_SIGMA_GRID_N,
        :integration_mode => :semi_infinite,
        :gc_every_n => 5,
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
        elseif arg == "--mode"
            mode_str = require_value()
            mode_sym = Symbol(mode_str)
            mode_sym in (:semi_infinite, :finite_15, :finite_lambda) || error("unknown mode: $mode_str (expected: semi_infinite, finite_15, finite_lambda)")
            opts[:integration_mode] = mode_sym
        elseif arg == "--gc-every-n"
            opts[:gc_every_n] = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            print_usage(); exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    tstep = Float64(opts[:tstep]); tstep > 0 || error("tstep must be positive")
    return Options(
        String(opts[:out]), Bool(opts[:overwrite]), Bool(opts[:resume]),
        Float64(opts[:xi]), Float64(opts[:tmin]), Float64(opts[:tmax]), Float64(opts[:tstep]),
        Float64.(opts[:mub_list]), Int(opts[:p_num]), Int(opts[:t_num]), Int(opts[:max_iter]),
        Int(opts[:tau_p_nodes]), Int(opts[:tau_angle_nodes]), Int(opts[:tau_phi_nodes]),
        Int(opts[:tau_n_sigma]), Int(opts[:sigma_grid_n]),
        Symbol(opts[:integration_mode]), Int(opts[:gc_every_n]),
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function write_header(io, opts::Options)
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
        "sigma_grid_n", "integration_mode",
        "gc_every_n",
    ]
    println(io, join(cols, ','))
end

@inline csv_bool(x::Bool) = x ? "true" : "false"

function build_K_coeffs(T_fm::Float64, muq_fm::Float64, masses::NamedTuple, Φ::Float64, Φbar::Float64)
    A_u = A(masses.u, muq_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(masses.s, muq_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return (K_coeffs=calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s), A_vals=(u=A_u, d=A_u, s=A_s))
end

function densities_from_equilibrium(x_state, mu_vec, T_fm, thermal_nodes, xi)
    nd = calculate_number_densities(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return (
        u=Float64(nd.quark[1]), d=Float64(nd.quark[2]), s=Float64(nd.quark[3]),
        ubar=Float64(nd.antiquark[1]), dbar=Float64(nd.antiquark[2]), sbar=Float64(nd.antiquark[3]),
    )
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
    n_sigma_points::Int, sigma_grid_n::Int, integration_mode::Symbol)
    cs_caches = Dict{Symbol,RT_ASR.CrossSectionCache}()
    
    # 所有模式都使用基于 Λ 截断的 σ(s) 范围（确保一致性）
    # 这是物理上正确的行为：PNJL 模型的动量截断 Λ 决定了 σ(s) 的有效范围
    for process in processes
        # 使用 w0cdf 设计，但限制在 Λ 截断范围内
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
        cs_caches[process] = cache
    end
    return cs_caches
end

function run_scan(opts::Options)
    ensure_parent_dir(opts.out)
    existing = opts.resume && isfile(opts.out) && !opts.overwrite ? 
        ScanCSV.read_existing_keys(opts.out, ["T_MeV", "muB_MeV", "xi"]) : Set{Tuple{Float64,Float64,Float64}}()
    opts.overwrite && isfile(opts.out) && rm(opts.out)

    new_file = !isfile(opts.out) || filesize(opts.out) == 0
    io = open(opts.out, "a")
    try
        if new_file
            ScanCSV.write_metadata(io, Dict(
                "schema" => "scan_csv_v1",
                "title" => "relaxation_times_vs_T",
                "script" => "scripts/relaxtime/scan_relaxation_times_vs_T.jl",
                "x_label" => "T", "x_unit" => "MeV",
                "y_unit.tau_u" => "fm", "y_unit.tau_s" => "fm",
                "y_unit.tau_ubar" => "fm", "y_unit.tau_sbar" => "fm",
                "y_scale.tau_u" => "log", "y_scale.tau_s" => "log",
                "y_scale.tau_ubar" => "log", "y_scale.tau_sbar" => "log",
                "sigma_grid_n" => string(opts.sigma_grid_n),
                "integration_mode" => string(opts.integration_mode),
                "gc_every_n" => string(opts.gc_every_n),
            ))
            write_header(io, opts)
        end

        # 根据积分模式准备动量网格
        p_grid, p_w, sigma_cutoff = if opts.integration_mode == :finite_15
            # 有限截断到15 fm⁻¹，σ(s)范围基于Λ
            pg, pw = gauleg(0.0, 15.0, opts.tau_p_nodes)
            (pg, pw, Λ_inv_fm)
        elseif opts.integration_mode == :finite_lambda
            # 有限截断到Λ
            pg, pw = gauleg(0.0, Λ_inv_fm, opts.tau_p_nodes)
            (pg, pw, Λ_inv_fm)
        else  # :semi_infinite
            # 半无穷积分，使用模块默认行为
            (nothing, nothing, nothing)
        end
        
        cos_grid, cos_w = gauleg(-1.0, 1.0, opts.tau_angle_nodes)
        phi_grid, phi_w = gauleg(0.0, 2π, opts.tau_phi_nodes)

        @printf("Integration mode: %s\n", opts.integration_mode)
        if opts.integration_mode == :finite_15
            @printf("  p_grid: [0, 15] fm⁻¹, sigma_cutoff: Λ = %.4f fm⁻¹\n", Λ_inv_fm)
        elseif opts.integration_mode == :finite_lambda
            @printf("  p_grid: [0, Λ] = [0, %.4f] fm⁻¹\n", Λ_inv_fm)
        else
            @printf("  p_grid: [0, ∞) via semi-infinite transform\n")
        end

        T_vals = collect(range(opts.tmin_mev, opts.tmax_mev; step=opts.tstep_mev))
        for muB_mev in opts.mub_list_mev
            muq_mev = muB_mev / 3.0
            muq_fm = muq_mev / ħc_MeV_fm
            local_iter = 0

            for T_mev in T_vals
                opts.resume && !opts.overwrite && (T_mev, muB_mev, opts.xi) in existing && continue
                T_fm = T_mev / ħc_MeV_fm

                base = solve(FixedMu(), T_fm, muq_fm; xi=opts.xi, p_num=opts.p_num, t_num=opts.t_num)
                Φ, Φbar = Float64(base.x_state[4]), Float64(base.x_state[5])
                masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

                thermo_params = (T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(opts.xi))
                ktmp = build_K_coeffs(Float64(T_fm), Float64(muq_fm), masses, Φ, Φbar)
                quark_params = (m=masses, μ=(u=Float64(muq_fm), d=Float64(muq_fm), s=Float64(muq_fm)), A=ktmp.A_vals)
                K_coeffs = ktmp.K_coeffs

                thermal_nodes = cached_nodes(opts.p_num, opts.t_num)
                densities = densities_from_equilibrium(base.x_state, base.mu_vec, T_fm, thermal_nodes, Float64(opts.xi))

                tau = (u=NaN, s=NaN, ubar=NaN, sbar=NaN)
                tauinv = (u=NaN, s=NaN, ubar=NaN, sbar=NaN)
                if Bool(base.converged)
                    cs_caches = build_sigma_caches(REQUIRED_PROCESSES, quark_params, thermo_params, K_coeffs;
                        n_sigma_points=opts.tau_n_sigma_points, sigma_grid_n=opts.sigma_grid_n,
                        integration_mode=opts.integration_mode)

                    tau_res = relaxation_times(quark_params, thermo_params, K_coeffs;
                        densities=densities, cs_caches=cs_caches,
                        p_nodes=opts.tau_p_nodes, angle_nodes=opts.tau_angle_nodes, phi_nodes=opts.tau_phi_nodes,
                        p_grid=p_grid, p_w=p_w,
                        cos_grid=cos_grid, cos_w=cos_w, phi_grid=phi_grid, phi_w=phi_w,
                        n_sigma_points=opts.tau_n_sigma_points,
                        sigma_cutoff=sigma_cutoff,
                    )
                    tau, tauinv = tau_res.tau, tau_res.tau_inv
                end

                row = Any[
                    T_mev, muB_mev, muq_mev, opts.xi, T_fm, muq_fm,
                    csv_bool(Bool(base.converged)), Int(base.iterations), Float64(base.residual_norm),
                    Φ, Φbar, masses.u, masses.s,
                    densities.u, densities.s, densities.ubar, densities.sbar,
                    tau.u, tau.s, tau.ubar, tau.sbar,
                    tauinv.u, tauinv.s, tauinv.ubar, tauinv.sbar,
                    opts.tau_p_nodes, opts.tau_angle_nodes, opts.tau_phi_nodes, opts.tau_n_sigma_points,
                    opts.sigma_grid_n, string(opts.integration_mode),
                    opts.gc_every_n,
                ]
                println(io, join(row, ',')); flush(io)

                if Bool(base.converged)
                    @printf("T=%6.1f MeV  muB=%6.1f MeV | tau_u=%.3e tau_s=%.3e tau_ubar=%.3e tau_sbar=%.3e\n",
                        T_mev, muB_mev, tau.u, tau.s, tau.ubar, tau.sbar)
                else
                    @printf("T=%6.1f MeV  muB=%6.1f MeV | equilibrium NOT converged; tau=NaN\n", T_mev, muB_mev)
                end

                local_iter += 1
                opts.gc_every_n > 0 && local_iter % opts.gc_every_n == 0 && GC.gc()
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
