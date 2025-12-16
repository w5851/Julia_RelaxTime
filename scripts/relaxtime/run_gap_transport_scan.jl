"""
批量扫描 (T, μ) 网格，串联：各向异性 PNJL 平衡求解 → τ → RTA 输运系数。

输出：CSV（默认 data/outputs/results/relaxtime/gap_transport_scan.csv）

单位约定：
- CLI 的温度/化学势默认使用 MeV（更符合扫描习惯）；脚本内部会换算到 fm⁻¹。
- 输出同时包含 MeV 与 fm⁻¹ 的关键量。

示例：
  julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --tmin 50 --tmax 200 --tstep 10 --mumin 0 --mumax 400 --mustep 20 --xi 0.2 --overwrite

注意：
- compute_bulk 默认关闭（体粘滞需要多次自动微分+求解，扫描会很慢）。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .TransportWorkflow: solve_gap_and_transport
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .EffectiveCouplings.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using StaticArrays

struct ScanOptions
    output::String
    xi_values::Vector{Float64}
    tmin_mev::Float64
    tmax_mev::Float64
    tstep_mev::Float64
    mumin_mev::Float64
    mumax_mev::Float64
    mustep_mev::Float64
    overwrite::Bool
    resume::Bool
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
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_gap_transport_scan.jl [options]\n")
    println("Options:")
    println("  --output <path>             输出 CSV (default data/outputs/results/relaxtime/gap_transport_scan.csv)")
    println("  --xi <value>                追加一个 ξ 值（可多次传入）")
    println("  --xi-list v1,v2,...         用逗号分隔的 ξ 列表替换")
    println("  --tmin/--tmax/--tstep <MeV> 温度范围与步长")
    println("  --mumin/--mumax/--mustep <MeV> 化学势范围与步长")
    println("  --overwrite                 覆盖输出文件")
    println("  --no-resume                 禁用跳过逻辑，强制重算")
    println("  --p-num <int>               能隙/密度的动量节点数 (default 12)")
    println("  --t-num <int>               能隙/密度的角度节点数 (default 6)")
    println("  --max-iter <int>            NLsolve iterations 上限 (default 40)")
    println("  --tau-p-nodes <int>         τ 平均散射率动量节点 (default 6)")
    println("  --tau-angle-nodes <int>     τ 平均散射率 cosθ 节点 (default 2)")
    println("  --tau-phi-nodes <int>       τ 平均散射率 φ 节点 (default 4)")
    println("  --tau-n-sigma <int>         截面 t 积分点数 (default 8)")
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
        :mumin => 0.0,
        :mumax => 400.0,
        :mustep => 20.0,
        :overwrite => false,
        :resume => true,
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :tau_p_nodes => 6,
        :tau_angle_nodes => 2,
        :tau_phi_nodes => 4,
        :tau_n_sigma => 8,
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
        elseif arg == "--mumin"
            opts[:mumin] = parse(Float64, require_value())
        elseif arg == "--mumax"
            opts[:mumax] = parse(Float64, require_value())
        elseif arg == "--mustep"
            opts[:mustep] = parse(Float64, require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--no-resume"
            opts[:resume] = false
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
    mustep = Float64(opts[:mustep]); mustep > 0 || error("mustep must be positive")

    return ScanOptions(
        String(opts[:output]),
        xi_vals,
        Float64(opts[:tmin]),
        Float64(opts[:tmax]),
        tstep,
        Float64(opts[:mumin]),
        Float64(opts[:mumax]),
        mustep,
        Bool(opts[:overwrite]),
        Bool(opts[:resume]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Int(opts[:tau_p_nodes]),
        Int(opts[:tau_angle_nodes]),
        Int(opts[:tau_phi_nodes]),
        Int(opts[:tau_n_sigma]),
        Int(opts[:tr_p_nodes]),
        Float64(opts[:tr_p_max]),
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function read_existing_keys(path::AbstractString)
    keys = Set{Tuple{Float64,Float64,Float64}}()
    isfile(path) || return keys
    open(path, "r") do io
        header = readline(io; keep=true)
        isempty(header) && return keys
        # assume columns contain T_MeV, mu_MeV, xi
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(line, ',')
            length(parts) < 3 && continue
            T = tryparse(Float64, parts[1]); μ = tryparse(Float64, parts[2]); xi = tryparse(Float64, parts[3])
            (T === nothing || μ === nothing || xi === nothing) && continue
            push!(keys, (T, μ, xi))
        end
    end
    return keys
end

function write_header_if_needed(io)
    header = join([
        "T_MeV", "mu_MeV", "xi",
        "T_fm", "mu_fm",
        "converged", "iterations", "residual_norm",
        "Phi", "Phibar",
        "m_u", "m_d", "m_s",
        "rho_baryon", "rho_norm",
        "n_u", "n_d", "n_s", "n_ubar", "n_dbar", "n_sbar",
        "tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
        "tauinv_u", "tauinv_d", "tauinv_s", "tauinv_ubar", "tauinv_dbar", "tauinv_sbar",
        "eta", "sigma", "zeta",
    ], ',')
    println(io, header)
end

function csv_bool(x::Bool)
    return x ? "true" : "false"
end

function compute_K_coeffs(T_fm::Float64, mu_fm::Float64, masses::NamedTuple, Φ::Float64, Φbar::Float64)
    nodes = DEFAULT_MOMENTUM_NODES
    weights = DEFAULT_MOMENTUM_WEIGHTS
    A_u = A(masses.u, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    A_s = A(masses.s, mu_fm, T_fm, Φ, Φbar, nodes, weights)
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    return calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
end

function run_scan(opts::ScanOptions)
    ensure_parent_dir(opts.output)

    existing = opts.resume && isfile(opts.output) && !opts.overwrite ? read_existing_keys(opts.output) : Set{Tuple{Float64,Float64,Float64}}()

    if opts.overwrite && isfile(opts.output)
        rm(opts.output)
    end

    new_file = !isfile(opts.output)
    io = open(opts.output, "a")
    try
        if new_file
            write_header_if_needed(io)
        end

        T_values = collect(range(opts.tmin_mev; stop=opts.tmax_mev, step=opts.tstep_mev))
        mu_values = collect(range(opts.mumin_mev; stop=opts.mumax_mev, step=opts.mustep_mev))

        total = length(opts.xi_values) * length(T_values) * length(mu_values)
        done = 0

        for xi in opts.xi_values
            for T_mev in T_values
                seed_state = nothing
                for mu_mev in mu_values
                    done += 1
                    key = (T_mev, mu_mev, xi)
                    if opts.resume && (key in existing)
                        continue
                    end

                    T_fm = T_mev / ħc_MeV_fm
                    mu_fm = mu_mev / ħc_MeV_fm

                    # 先求一次平衡解（用于 K_coeffs，并作为后续 workflow 的 equilibrium 复用）
                    base = TransportWorkflow.ThermoDerivatives.solve_equilibrium_mu(
                        T_fm,
                        mu_fm;
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

                    K_coeffs = compute_K_coeffs(Float64(T_fm), Float64(mu_fm), masses, Φ, Φbar)

                    # 正式计算：τ + 输运（η/σ）; ζ 默认关
                    res = solve_gap_and_transport(
                        T_fm,
                        mu_fm;
                        xi=xi,
                        equilibrium=base,
                        compute_tau=true,
                        K_coeffs=K_coeffs,
                        tau=nothing,
                        compute_bulk=false,
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
                        string(T_mev), string(mu_mev), string(xi),
                        string(T_fm), string(mu_fm),
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

                    if done % 10 == 0
                        println("progress: $(done)/$(total) (last T=$(T_mev) MeV, mu=$(mu_mev) MeV, xi=$(xi))")
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
