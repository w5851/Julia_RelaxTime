"""
沿 freeze-out baseline 路径批量运行介子数密度 workflow。

主链约束：
- 统一通过 `Models.run_freezeout_meson_density_scan` 进入；
- 不在脚本层手工拼 freeze-out 映射、介子化学势 profile、平衡态与数密度流程。

输出：CSV（默认 `data/outputs/results/relaxtime/meson_density/freezeout/freezeout_meson_density_scan.csv`）
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Models: run_freezeout_meson_density_scan

struct ScanOptions
    output::String
    sqrt_s_values::Vector{Float64}
    xi_values::Vector{Float64}
    freezeout_profile::String
    path_profile::String
    flavor_chemical_profile::String
    meson_chemical_profile::String
    regime::Symbol
    traversal::Symbol
    overwrite::Bool
    resume::Bool
    p_num::Int
    t_num::Int
    max_iter::Int
    stable_q_nodes::Int
    qmax::Float64
    q_nodes::Int
    omega_min::Float64
    omega_max::Float64
    omega_nodes::Int
    eta::Float64
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_freezeout_meson_density_scan.jl [options]\n")
    println("Options:")
    println("  --output <path>                     输出 CSV")
    println("  --sqrts <value>                     追加一个 sqrt(s_NN) 点（GeV，可多次）")
    println("  --sqrts-list v1,v2,...              用逗号分隔的 sqrt(s_NN) 列表替换")
    println("  --xi <value>                        追加一个 ξ 值（可多次）")
    println("  --xi-list v1,v2,...                 用逗号分隔的 ξ 列表替换")
    println("  --freezeout-profile <name>          freeze-out profile (default default)")
    println("  --path-profile <name>               path profile (default baseline_freezeout)")
    println("  --flavor-chemical-profile <name>    flavor chemical profile (default default)")
    println("  --meson-chemical-profile <name>     meson chemical profile (default default)")
    println("  --regime <stable|strict_bw_stage1|strict_bw_stage2|phase_shift_current|phase_shift_gbu>")
    println("  --traversal <as_given|sqrts_ascending|sqrts_descending|muB_descending|muB_ascending>")
    println("  --overwrite                         覆盖输出文件")
    println("  --no-resume                         禁用续算跳过")
    println("  --p-num <int>                       平衡态动量节点数 (default 24)")
    println("  --t-num <int>                       平衡态角度节点数 (default 8)")
    println("  --max-iter <int>                    solver/mass iterations 上限 (default 40)")
    println("  --stable-q-nodes <int>              stable 数密度 q 节点数 (default 256)")
    println("  --qmax <float>                      strict BW / phase-shift q 上限 (default 12)")
    println("  --q-nodes <int>                     strict BW / phase-shift q 节点数 (default 48)")
    println("  --omega-min <float>                 phase-shift omega 下限 (default 0.05)")
    println("  --omega-max <float>                 strict BW / phase-shift omega 上限 (default 10)")
    println("  --omega-nodes <int>                 strict BW / phase-shift omega 节点数 (default 48)")
    println("  --eta <float>                       phase-shift eta (default 1e-6)")
    println("  -h, --help                          显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :output => joinpath("data", "outputs", "results", "relaxtime", "meson_density", "freezeout", "freezeout_meson_density_scan.csv"),
        :sqrts_values => Float64[20.0],
        :xi_values => Float64[0.0],
        :freezeout_profile => "default",
        :path_profile => "baseline_freezeout",
        :flavor_chemical_profile => "default",
        :meson_chemical_profile => "default",
        :regime => :stable,
        :traversal => :sqrts_ascending,
        :overwrite => false,
        :resume => true,
        :p_num => 24,
        :t_num => 8,
        :max_iter => 40,
        :stable_q_nodes => 256,
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
        elseif arg == "--sqrts"
            val = parse(Float64, require_value())
            if opts[:sqrts_values] == Float64[20.0]
                opts[:sqrts_values] = Float64[]
            end
            push!(opts[:sqrts_values], val)
        elseif arg == "--sqrts-list"
            raw = split(require_value(), ',')
            vals = Float64[parse(Float64, strip(v)) for v in raw if !isempty(strip(v))]
            opts[:sqrts_values] = vals
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
        elseif arg == "--freezeout-profile"
            opts[:freezeout_profile] = require_value()
        elseif arg == "--path-profile"
            opts[:path_profile] = require_value()
        elseif arg == "--flavor-chemical-profile"
            opts[:flavor_chemical_profile] = require_value()
        elseif arg == "--meson-chemical-profile"
            opts[:meson_chemical_profile] = require_value()
        elseif arg == "--regime"
            opts[:regime] = Symbol(require_value())
        elseif arg == "--traversal"
            opts[:traversal] = Symbol(require_value())
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
        elseif arg == "--stable-q-nodes"
            opts[:stable_q_nodes] = parse(Int, require_value())
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
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    sqrt_s_vals = unique(Float64.(opts[:sqrts_values]))
    xi_vals = unique(Float64.(opts[:xi_values]))
    isempty(sqrt_s_vals) && error("sqrt(s_NN) grid must not be empty")
    isempty(xi_vals) && error("xi grid must not be empty")

    return ScanOptions(
        String(opts[:output]),
        sqrt_s_vals,
        xi_vals,
        String(opts[:freezeout_profile]),
        String(opts[:path_profile]),
        String(opts[:flavor_chemical_profile]),
        String(opts[:meson_chemical_profile]),
        Symbol(opts[:regime]),
        Symbol(opts[:traversal]),
        Bool(opts[:overwrite]),
        Bool(opts[:resume]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Int(opts[:stable_q_nodes]),
        Float64(opts[:qmax]),
        Int(opts[:q_nodes]),
        Float64(opts[:omega_min]),
        Float64(opts[:omega_max]),
        Int(opts[:omega_nodes]),
        Float64(opts[:eta]),
    )
end

function main()
    opts = parse_args(ARGS)
    mkpath(dirname(opts.output))

    out = run_freezeout_meson_density_scan(
        sqrt_s_NN_values=opts.sqrt_s_values,
        xi_values=opts.xi_values,
        freezeout_profile_name=opts.freezeout_profile,
        path_profile_name=opts.path_profile,
        flavor_chemical_profile_name=opts.flavor_chemical_profile,
        meson_chemical_profile_name=opts.meson_chemical_profile,
        regime=opts.regime,
        output_path=opts.output,
        traversal=opts.traversal,
        overwrite=opts.overwrite,
        resume=opts.resume,
        p_num=opts.p_num,
        t_num=opts.t_num,
        stable_num_q_nodes=opts.stable_q_nodes,
        strict_bw_qmax=opts.qmax,
        strict_bw_q_nodes=opts.q_nodes,
        strict_bw_omega_max=opts.omega_max,
        strict_bw_omega_nodes=opts.omega_nodes,
        phase_shift_qmax=opts.qmax,
        phase_shift_q_nodes=opts.q_nodes,
        phase_shift_omega_min=opts.omega_min,
        phase_shift_omega_max=opts.omega_max,
        phase_shift_omega_nodes=opts.omega_nodes,
        phase_shift_eta=opts.eta,
        solver_kwargs=(; iterations=opts.max_iter),
        mass_kwargs=(; iterations=opts.max_iter),
    )

    println(out)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
