"""
批量扫描 Phase-E3 最小相移双积分 meson density workflow，输出 `π/K` 数密度与 `K/π` 比值。

主链约束：
- 统一通过 `Models.solve_gap_and_phase_shift_meson_density_point` 进入；
- 不在脚本层手工拼平衡态 / 介子质量 / 相移双积分流程；
- 当前仅支持 `xi = 0`。

输出：CSV（默认 `data/outputs/results/relaxtime/scan/phase_shift_meson_density_scan.csv`）
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))
include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .ScanCSV: ScanCSV
using .Constants_PNJL: ħc_MeV_fm
using .Models: solve_gap_and_phase_shift_meson_density_point

struct ScanOptions
    output::String
    xi_values::Vector{Float64}
    tmin_mev::Float64
    tmax_mev::Float64
    tstep_mev::Float64
    mub_mev::Float64
    overwrite::Bool
    resume::Bool
    p_num::Int
    t_num::Int
    max_iter::Int
    qmax::Float64
    q_nodes::Int
    omega_min::Float64
    omega_max::Float64
    omega_nodes::Int
    eta::Float64
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_phase_shift_meson_density_scan.jl [options]\n")
    println("Options:")
    println("  --output <path>             输出 CSV (default data/outputs/results/relaxtime/scan/phase_shift_meson_density_scan.csv)")
    println("  --xi <value>                追加一个 ξ 值（当前仅支持 0，可多次传入）")
    println("  --xi-list v1,v2,...         用逗号分隔的 ξ 列表替换（当前仅支持 0）")
    println("  --tmin/--tmax/--tstep <MeV> 温度范围与步长")
    println("  --mub <MeV>                 固定重子化学势 μ_B (default 0)")
    println("  --overwrite                 覆盖输出文件")
    println("  --no-resume                 禁用跳过逻辑，强制重算")
    println("  --p-num <int>               平衡态求解动量节点数 (default 12)")
    println("  --t-num <int>               平衡态求解角度节点数 (default 6)")
    println("  --max-iter <int>            平衡态 / 介子求解迭代上限 (default 40)")
    println("  --qmax <fm^-1>              外层 q 硬截断 (default 12)")
    println("  --q-nodes <int>             外层 q 积分节点数 (default 48)")
    println("  --omega-min <fm^-1>         内层 omega 下限 (default 0.05)")
    println("  --omega-max <fm^-1>         内层 omega 上限 (default 10)")
    println("  --omega-nodes <int>         内层 omega 节点数 (default 48)")
    println("  --eta <float>               相位诊断中的小虚部参数 (default 1e-6)")
    println("  -h, --help                  显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :output => joinpath("data", "outputs", "results", "relaxtime", "scan", "phase_shift_meson_density_scan.csv"),
        :xi_values => Float64[0.0],
        :tmin => 208.0,
        :tmax => 212.0,
        :tstep => 2.0,
        :mub => 0.0,
        :overwrite => false,
        :resume => true,
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :qmax => 12.0,
        :q_nodes => 48,
        :omega_min => 0.05,
        :omega_max => 10.0,
        :omega_nodes => 48,
        :eta => 1e-6,
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
        elseif arg == "--mub"
            opts[:mub] = parse(Float64, require_value())
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
        elseif arg == "--qmax"
            opts[:qmax] = parse(Float64, require_value())
        elseif arg == "--q-nodes"
            opts[:q_nodes] = parse(Int, require_value())
        elseif arg == "--omega-min"
            opts[:omega_min] = parse(Float64, require_value())
        elseif arg == "--omega-max"
            opts[:omega_max] = parse(Float64, require_value())
        elseif arg == "--omega-nodes"
            opts[:omega_nodes] = parse(Int, require_value())
        elseif arg == "--eta"
            opts[:eta] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage(); exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    xi_vals = unique(sort(Float64.(opts[:xi_values])))
    isempty(xi_vals) && (xi_vals = Float64[0.0])
    any(!iszero, xi_vals) && error("current Phase-E3 phase-shift meson-density scan only supports xi = 0")

    tstep = Float64(opts[:tstep]); tstep > 0 || error("tstep must be positive")
    q_nodes = Int(opts[:q_nodes]); q_nodes > 1 || error("q-nodes must be > 1")
    omega_nodes = Int(opts[:omega_nodes]); omega_nodes > 1 || error("omega-nodes must be > 1")
    qmax = Float64(opts[:qmax]); qmax > 0.0 || error("qmax must be positive")
    omega_min = Float64(opts[:omega_min])
    omega_max = Float64(opts[:omega_max]); omega_max > omega_min || error("omega-max must exceed omega-min")

    return ScanOptions(
        String(opts[:output]),
        xi_vals,
        Float64(opts[:tmin]),
        Float64(opts[:tmax]),
        tstep,
        Float64(opts[:mub]),
        Bool(opts[:overwrite]),
        Bool(opts[:resume]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        qmax,
        q_nodes,
        omega_min,
        omega_max,
        omega_nodes,
        Float64(opts[:eta]),
    )
end

@inline function _format(x)
    x isa Bool && return (x ? "true" : "false")
    return string(x)
end

const OUTPUT_COLUMNS = [
    "T_MeV", "muB_MeV", "xi",
    "T_fm", "muB_fm", "mu_fm",
    "Phi", "Phibar",
    "m_u", "m_d", "m_s",
    "m_pi", "m_K",
    "n_pi", "n_K", "kpi_ratio",
    "d_pi", "d_K",
    "qmax", "q_nodes",
    "omega_min", "omega_max", "omega_nodes",
    "eta",
    "pi_q_integral_estimate", "pi_omega_shell_at_qmax",
    "K_q_integral_estimate", "K_omega_shell_at_qmax",
]

function _row_to_values(cols::Vector{String}, row::Dict{String,Any})
    return [_format(get(row, c, "")) for c in cols]
end

function main()
    opts = parse_args(ARGS)
    outpath = opts.output
    mkpath(dirname(outpath))

    key_cols = ["T_MeV", "muB_MeV", "xi"]
    existing = (isfile(outpath) && opts.resume && !opts.overwrite) ? ScanCSV.read_existing_keys(outpath, key_cols) : Set{Tuple{Vararg{Float64}}}()

    if opts.overwrite && isfile(outpath)
        rm(outpath)
        empty!(existing)
    end

    is_new = !isfile(outpath)
    if !is_new
        ScanCSV.assert_required_columns(outpath, OUTPUT_COLUMNS)
    end

    open(outpath, is_new ? "w" : "a") do io
        if is_new
            ScanCSV.write_metadata(io, Dict(
                "format" => "scan_csv_v1",
                "script" => "scripts/relaxtime/run_phase_shift_meson_density_scan.jl",
                "workflow_entry" => "Models.solve_gap_and_phase_shift_meson_density_point",
                "mesons" => "pi,K",
                "continuation" => "MesonMassWorkflow.continuation_state",
                "phase_shift_scheme" => "GL_hard_cutoff_phase_e3",
                "note" => "mu_fm denotes quark chemical potential (muB/3); current helper only supports xi=0"
            ))
            ScanCSV.write_header(io, OUTPUT_COLUMNS)
        end

        muB = opts.mub_mev
        muB_fm = muB / ħc_MeV_fm
        mu_fm = muB_fm / 3.0

        for xi in opts.xi_values
            continuation_state = nothing
            T = opts.tmin_mev
            while T <= opts.tmax_mev + 1e-9
                key = (Float64(T), Float64(muB), Float64(xi))
                if key in existing
                    T += opts.tstep_mev
                    continue
                end

                T_fm = T / ħc_MeV_fm
                res = solve_gap_and_phase_shift_meson_density_point(
                    T_fm,
                    mu_fm;
                    xi=xi,
                    mesons=(:pi, :K),
                    continuation_state=continuation_state,
                    mixed_branch_align=:strict_sign_binding,
                    p_num=opts.p_num,
                    t_num=opts.t_num,
                    solver_kwargs=(; iterations=opts.max_iter),
                    mass_kwargs=(; iterations=opts.max_iter),
                    density_kwargs=(;
                        qmax=opts.qmax,
                        q_nodes=opts.q_nodes,
                        omega_min=opts.omega_min,
                        omega_max=opts.omega_max,
                        omega_nodes=opts.omega_nodes,
                        eta=opts.eta,
                    ),
                )

                qp = res.quark_params
                tp = res.thermo_params
                md = res.phase_shift_meson_density

                row = Dict{String,Any}(
                    "T_MeV" => T,
                    "muB_MeV" => muB,
                    "xi" => xi,
                    "T_fm" => T_fm,
                    "muB_fm" => muB_fm,
                    "mu_fm" => mu_fm,
                    "Phi" => tp.Φ,
                    "Phibar" => tp.Φbar,
                    "m_u" => qp.m.u,
                    "m_d" => qp.m.d,
                    "m_s" => qp.m.s,
                    "m_pi" => md.m_pi,
                    "m_K" => md.m_K,
                    "n_pi" => md.n_pi,
                    "n_K" => md.n_K,
                    "kpi_ratio" => md.kpi_ratio,
                    "d_pi" => md.d_pi,
                    "d_K" => md.d_K,
                    "qmax" => md.qmax,
                    "q_nodes" => md.q_nodes,
                    "omega_min" => md.omega_min,
                    "omega_max" => md.omega_max,
                    "omega_nodes" => md.omega_nodes,
                    "eta" => md.eta,
                    "pi_q_integral_estimate" => md.pi_density.q_integral_estimate,
                    "pi_omega_shell_at_qmax" => md.pi_density.omega_shell_at_qmax,
                    "K_q_integral_estimate" => md.k_density.q_integral_estimate,
                    "K_omega_shell_at_qmax" => md.k_density.omega_shell_at_qmax,
                )

                println(io, join(_row_to_values(OUTPUT_COLUMNS, row), ','))
                flush(io)

                continuation_state = res.continuation_state
                T += opts.tstep_mev
            end
        end
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
