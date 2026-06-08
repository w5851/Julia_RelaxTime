module GapTransportScanCLI

export ScanOptions
export parse_args
export print_usage

struct ScanOptions
    output::String
    channel_diagnostics_output::Union{Nothing, String}
    failed_points_output::Union{Nothing, String}
    provenance_dir::Union{Nothing, String}
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
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
    tau_threshold_subtraction::Bool
    tau_asym_window::Float64
    tau_asym_fit_min_points::Int
    tau_asym_extra_points::Int
    tau_interpolation_mode::Symbol
    propagator_xi_policy::Symbol
    sigma_grid_n::Int
    integration_mode::Symbol
    gc_every_n::Int
    tr_p_nodes::Int
    tr_p_max_fm::Float64
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_gap_transport_scan.jl [options]\n")
    println("Options:")
    println("  --output <path>             输出 CSV (default data/outputs/results/relaxtime/scan/gap_transport_scan.csv)")
    println("  --channel-diagnostics-output <path>  可选：输出每点每通道的 τ^-1 贡献明细 CSV")
    println("  --failed-points-output <path>        可选：输出跳过失败点 sidecar CSV")
    println("  --provenance-dir <path>              可选：写入 run_manifest/effective_config 的目录（默认 output 同目录）")
    println("  --xi <value>                追加一个 ξ 值（可多次传入）")
    println("  --xi-list v1,v2,...         用逗号分隔的 ξ 列表替换")
    println("  --tmin/--tmax/--tstep <MeV> 温度范围与步长")
    println("  --mubmin/--mubmax/--mubstep <MeV> 重子化学势 μ_B 范围与步长")
    println("  --mumin/--mumax/--mustep <MeV> (兼容旧参数) 夸克化学势 μ_q 范围与步长")
    println("  --overwrite                 覆盖输出文件")
    println("  --no-resume                 禁用跳过逻辑，强制重算")
    println("  --compute-bulk              显式开启体粘滞 ζ 计算（默认开启）")
    println("  --no-compute-bulk           显式关闭体粘滞 ζ 计算")
    println("  --p-num <int>               能隙/密度的动量节点数 (default 12)")
    println("  --t-num <int>               能隙/密度的角度节点数 (default 6)")
    println("  --max-iter <int>            NLsolve iterations 上限 (default 40)")
    println("  --tau-p-nodes <int>         τ 平均散射率动量节点 (default $(Main.MODULE_DEFAULT_P_NODES))")
    println("  --tau-angle-nodes <int>     τ 平均散射率 cosθ 节点 (default $(Main.MODULE_DEFAULT_ANGLE_NODES))")
    println("  --tau-phi-nodes <int>       τ 平均散射率 φ 节点 (default $(Main.MODULE_DEFAULT_PHI_NODES))")
    println("  --tau-n-sigma <int>         σ(s) 的 t 积分点数 (default $(Main.RT_TCS.DEFAULT_T_INTEGRAL_POINTS))")
    println("  --tau-threshold-subtraction / --tau-no-threshold-subtraction  是否启用阈值减法缓存 (default true)")
    println("  --tau-asym-window <float>   阈值减法拟合窗口 Δs (default 0.6)")
    println("  --tau-asym-fit-min-points <int> 阈值减法最小拟合点数 (default 8)")
    println("  --tau-asym-extra-points <int> 阈值减法补点数量 (default 10)")
    println("  --tau-interpolation-mode <pchip|linear|direct|hybrid_threshold>  τ 中 σ(s) 取值模式 (default linear)")
    println("  --propagator-xi-policy <match_thermo|isotropic>  σ(s)/propagator 的 ξ 口径；isotropic 表示传播子用 ξ=0，外层分布仍用真实 ξ (default match_thermo)")
    println("  --sigma-grid-n <int>        σ(s) 预计算网格点数 (default $(Main.MODULE_DEFAULT_SIGMA_GRID_N))")
    println("  --mode <mode>               τ 积分模式: semi_infinite | finite_15 | finite_lambda (default semi_infinite)")
    println("  --gc-every-n <int>          每 N 个点触发一次 GC (default 5; 0 表示关闭)")
    println("  --tr-p-nodes <int>          输运积分动量节点 (default 24)")
    println("  --tr-p-max <fm^-1>          输运积分 p 上限 (default 8.0)")
    println("  -h, --help                  显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :output => joinpath("data", "outputs", "results", "relaxtime", "scan", "gap_transport_scan.csv"),
        :channel_diagnostics_output => nothing,
        :failed_points_output => nothing,
        :provenance_dir => nothing,
        :xi_values => Float64[0.0],
        :tmin => 50.0,
        :tmax => 200.0,
        :tstep => 10.0,
        :mubmin => 0.0,
        :mubmax => 1200.0,
        :mubstep => 60.0,
        :overwrite => false,
        :resume => true,
        :compute_bulk => true,
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :tau_p_nodes => Main.MODULE_DEFAULT_P_NODES,
        :tau_angle_nodes => Main.MODULE_DEFAULT_ANGLE_NODES,
        :tau_phi_nodes => Main.MODULE_DEFAULT_PHI_NODES,
        :tau_n_sigma => Main.RT_TCS.DEFAULT_T_INTEGRAL_POINTS,
        :tau_threshold_subtraction => true,
        :tau_asym_window => 0.6,
        :tau_asym_fit_min_points => 8,
        :tau_asym_extra_points => 10,
        :tau_interpolation_mode => :linear,
        :propagator_xi_policy => :match_thermo,
        :sigma_grid_n => Main.MODULE_DEFAULT_SIGMA_GRID_N,
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
        elseif arg == "--channel-diagnostics-output"
            opts[:channel_diagnostics_output] = require_value()
        elseif arg == "--failed-points-output"
            opts[:failed_points_output] = require_value()
        elseif arg == "--provenance-dir"
            opts[:provenance_dir] = require_value()
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
        elseif arg == "--no-compute-bulk"
            opts[:compute_bulk] = false
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
        elseif arg == "--tau-threshold-subtraction"
            opts[:tau_threshold_subtraction] = true
        elseif arg == "--tau-no-threshold-subtraction"
            opts[:tau_threshold_subtraction] = false
        elseif arg == "--tau-asym-window"
            opts[:tau_asym_window] = parse(Float64, require_value())
        elseif arg == "--tau-asym-fit-min-points"
            opts[:tau_asym_fit_min_points] = parse(Int, require_value())
        elseif arg == "--tau-asym-extra-points"
            opts[:tau_asym_extra_points] = parse(Int, require_value())
        elseif arg == "--tau-interpolation-mode"
            mode = Symbol(require_value())
            mode in (:pchip, :linear, :direct, :hybrid_threshold) || error("unknown tau interpolation mode: $(mode) (expected: pchip|linear|direct|hybrid_threshold)")
            opts[:tau_interpolation_mode] = mode
        elseif arg == "--propagator-xi-policy"
            policy = Symbol(require_value())
            policy in (:match_thermo, :isotropic) || error("unknown propagator xi policy: $(policy) (expected: match_thermo|isotropic)")
            opts[:propagator_xi_policy] = policy
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
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    xi_vals = opts[:xi_values]
    isempty(xi_vals) && (xi_vals = Float64[0.0])
    xi_vals = unique(sort(Float64.(xi_vals)))

    tstep = Float64(opts[:tstep])
    tstep > 0 || error("tstep must be positive")
    mubstep = Float64(opts[:mubstep])
    mubstep > 0 || error("mubstep must be positive")

    return ScanOptions(
        String(opts[:output]),
        isnothing(opts[:channel_diagnostics_output]) ? nothing : String(opts[:channel_diagnostics_output]),
        isnothing(opts[:failed_points_output]) ? nothing : String(opts[:failed_points_output]),
        isnothing(opts[:provenance_dir]) ? nothing : String(opts[:provenance_dir]),
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
        Bool(opts[:tau_threshold_subtraction]),
        Float64(opts[:tau_asym_window]),
        Int(opts[:tau_asym_fit_min_points]),
        Int(opts[:tau_asym_extra_points]),
        Symbol(opts[:tau_interpolation_mode]),
        Symbol(opts[:propagator_xi_policy]),
        Int(opts[:sigma_grid_n]),
        Symbol(opts[:integration_mode]),
        Int(opts[:gc_every_n]),
        Int(opts[:tr_p_nodes]),
        Float64(opts[:tr_p_max]),
    )
end

end # module GapTransportScanCLI
