"""
批量扫描 (T, μ) 网格，串联：各向异性 PNJL 平衡求解 → τ → RTA 输运系数。

输出：CSV（默认 data/outputs/results/relaxtime/scan/gap_transport_scan.csv）

单位约定：
- CLI 的温度/化学势默认使用 MeV（更符合扫描习惯）；脚本内部会换算到 fm⁻¹。
- 输出同时包含 MeV 与 fm⁻¹ 的关键量。

示例：
    julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --tmin 120 --tmax 400 --tstep 10 --mubmin 0 --mubmax 800 --mubstep 800 --xi-list -0.6,-0.4,-0.2,0.0,0.2,0.4,0.6 --mode finite_15 --compute-bulk --overwrite --output data/outputs/results/relaxtime/scan/gap_transport_scan_xi-0p6to0p6.csv

注意：
- compute_bulk 默认关闭（体粘滞需要多次自动微分+求解，扫描会很慢）。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "scan_quality.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "provenance_metadata.jl"))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm, ρ0_inv_fm3
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.OneLoopIntegrals: A
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS, gauleg
using StaticArrays
using Dates
using .ScanCSV: ScanCSV
using .ScanQuality: assess_point_quality
using .ProvenanceMetadata: ProvenanceMetadata

const PNJL_MODEL = Models.create_model(:PNJL)
const TransportWorkflow = Models.transport_workflow_module()
const RT_ASR = Main.AverageScatteringRate
const RT_TCS = Main.TotalCrossSection
const REQUIRED_PROCESSES = TransportWorkflow.RelaxationTime.REQUIRED_PROCESSES

function preferred_phase_reference_path(dense_name::String, legacy_name::String)
    dense_path = joinpath(PROJECT_ROOT, "data", "reference", "pnjl", dense_name)
    return isfile(dense_path) ? dense_path : joinpath(PROJECT_ROOT, "data", "reference", "pnjl", legacy_name)
end

const DEFAULT_PHASE_BOUNDARY_PATH = preferred_phase_reference_path("boundary_dense.csv", "boundary.csv")
const DEFAULT_PHASE_CEP_PATH = preferred_phase_reference_path("cep_dense.csv", "cep.csv")
const DEFAULT_PHASE_CROSSOVER_PATH = preferred_phase_reference_path("crossover_dense.csv", "crossover.csv")

const MODULE_DEFAULT_P_NODES = RT_ASR.DEFAULT_P_NODES           # 20
const MODULE_DEFAULT_ANGLE_NODES = RT_ASR.DEFAULT_ANGLE_NODES   # 4
const MODULE_DEFAULT_PHI_NODES = RT_ASR.DEFAULT_PHI_NODES       # 8
const MODULE_DEFAULT_SIGMA_GRID_N = RT_ASR.DEFAULT_SIGMA_GRID_N # 60
const PHASE_BOUNDARY_XI_CACHE = Ref{Union{Nothing, Vector{Float64}}}(nothing)
const PHASE_CROSSOVER_XI_CACHE = Ref{Union{Nothing, Vector{Float64}}}(nothing)
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
    # tau / cross-section settings
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
    tau_threshold_subtraction::Bool
    tau_asym_window::Float64
    tau_asym_fit_min_points::Int
    tau_asym_extra_points::Int
    tau_interpolation_mode::Symbol
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
    println("  --compute-bulk              计算体粘滞 ζ（很慢；默认关闭）")
    println("  --p-num <int>               能隙/密度的动量节点数 (default 12)")
    println("  --t-num <int>               能隙/密度的角度节点数 (default 6)")
    println("  --max-iter <int>            NLsolve iterations 上限 (default 40)")
    println("  --tau-p-nodes <int>         τ 平均散射率动量节点 (default $(MODULE_DEFAULT_P_NODES))")
    println("  --tau-angle-nodes <int>     τ 平均散射率 cosθ 节点 (default $(MODULE_DEFAULT_ANGLE_NODES))")
    println("  --tau-phi-nodes <int>       τ 平均散射率 φ 节点 (default $(MODULE_DEFAULT_PHI_NODES))")
    println("  --tau-n-sigma <int>         σ(s) 的 t 积分点数 (default $(RT_TCS.DEFAULT_T_INTEGRAL_POINTS))")
    println("  --tau-threshold-subtraction / --tau-no-threshold-subtraction  是否启用阈值减法缓存 (default true)")
    println("  --tau-asym-window <float>   阈值减法拟合窗口 Δs (default 0.6)")
    println("  --tau-asym-fit-min-points <int> 阈值减法最小拟合点数 (default 8)")
    println("  --tau-asym-extra-points <int> 阈值减法补点数量 (default 10)")
    println("  --tau-interpolation-mode <pchip|linear|direct|hybrid_threshold>  τ 中 σ(s) 取值模式 (default linear)")
    println("  --sigma-grid-n <int>        σ(s) 预计算网格点数 (default $(MODULE_DEFAULT_SIGMA_GRID_N))")
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
        :compute_bulk => false,
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :tau_p_nodes => MODULE_DEFAULT_P_NODES,
        :tau_angle_nodes => MODULE_DEFAULT_ANGLE_NODES,
        :tau_phi_nodes => MODULE_DEFAULT_PHI_NODES,
        :tau_n_sigma => RT_TCS.DEFAULT_T_INTEGRAL_POINTS,
        :tau_threshold_subtraction => true,
        :tau_asym_window => 0.6,
        :tau_asym_fit_min_points => 8,
        :tau_asym_extra_points => 10,
        :tau_interpolation_mode => :linear,
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
        Int(opts[:sigma_grid_n]),
        Symbol(opts[:integration_mode]),
        Int(opts[:gc_every_n]),
        Int(opts[:tr_p_nodes]),
        Float64(opts[:tr_p_max]),
    )
end

function write_channel_diagnostics_header_if_needed(io)
    header = join([
        "T_MeV", "muq_MeV", "muB_MeV", "xi",
        "species", "channel", "density_key", "multiplicity",
        "density", "rate", "contribution", "total", "tau_inv_species",
        "equilibrium_backend", "phase_curr", "phase_structure",
    ], ',')
    println(io, header)
end

@inline function _csv_quote(text::AbstractString)
    return "\"" * replace(text, "\"" => "\"\"") * "\""
end

function write_failed_points_header_if_needed(io)
    header = join([
        "T_MeV", "muB_MeV", "xi",
        "seed_source", "phase_prev", "phase_curr_hint",
        "error_type", "error_message", "timestamp",
    ], ',')
    println(io, header)
end

function write_failed_point_row!(io, T_mev::Float64, muB_mev::Float64, xi::Float64, diag, err)
    seed_source = hasproperty(diag, :seed_source) ? string(getproperty(diag, :seed_source)) : "unknown"
    phase_prev = hasproperty(diag, :phase_prev) ? string(getproperty(diag, :phase_prev)) : "unknown"
    phase_curr_hint = hasproperty(diag, :phase_curr) ? string(getproperty(diag, :phase_curr)) : "unknown"
    error_type = string(typeof(err))
    error_message = sprint(showerror, err)
    timestamp = Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS")
    row = join([
        string(T_mev),
        string(muB_mev),
        string(xi),
        _csv_quote(seed_source),
        _csv_quote(phase_prev),
        _csv_quote(phase_curr_hint),
        _csv_quote(error_type),
        _csv_quote(error_message),
        _csv_quote(timestamp),
    ], ',')
    println(io, row)
    flush(io)
end

@inline function _rate_with_alias(rates, key::Symbol)
    if hasproperty(rates, key)
        return getproperty(rates, key)
    end
    if key == :dubar_to_dubar && hasproperty(rates, :udbar_to_udbar)
        return getproperty(rates, :udbar_to_udbar)
    elseif key == :subar_to_subar && hasproperty(rates, :usbar_to_usbar)
        return getproperty(rates, :usbar_to_usbar)
    elseif key == :ubardbar_to_ubardbar && hasproperty(rates, :ud_to_ud)
        return getproperty(rates, :ud_to_ud)
    elseif key == :ubarubar_to_ubarubar && hasproperty(rates, :uu_to_uu)
        return getproperty(rates, :uu_to_uu)
    elseif key == :ubarsbar_to_ubarsbar && hasproperty(rates, :us_to_us)
        return getproperty(rates, :us_to_us)
    elseif key == :sbarsbar_to_sbarsbar && hasproperty(rates, :ss_to_ss)
        return getproperty(rates, :ss_to_ss)
    end
    return 0.0
end

function _fallback_relaxation_rate_contribution_rows(dens, rates)
    rows = NamedTuple[]

    function add_row(species::Symbol, channel::Symbol, density_key::Symbol, multiplicity::Float64)
        density = getproperty(dens, density_key)
        rate = Float64(_rate_with_alias(rates, channel))
        contribution = multiplicity * density * rate
        push!(rows, (
            species=species,
            channel=channel,
            density_key=density_key,
            multiplicity=multiplicity,
            density=density,
            rate=rate,
            contribution=contribution,
            total=0.0,
        ))
    end

    # u / d (isospin)
    add_row(:u, :uu_to_uu, :u, 1.0)
    add_row(:u, :ud_to_ud, :u, 1.0)
    add_row(:u, :uubar_to_uubar, :ubar, 1.0)
    add_row(:u, :uubar_to_ddbar, :ubar, 1.0)
    add_row(:u, :uubar_to_ssbar, :ubar, 1.0)
    add_row(:u, :udbar_to_udbar, :ubar, 1.0)
    add_row(:u, :us_to_us, :s, 1.0)
    add_row(:u, :usbar_to_usbar, :sbar, 1.0)

    add_row(:d, :uu_to_uu, :d, 1.0)
    add_row(:d, :ud_to_ud, :d, 1.0)
    add_row(:d, :uubar_to_uubar, :dbar, 1.0)
    add_row(:d, :uubar_to_ddbar, :dbar, 1.0)
    add_row(:d, :uubar_to_ssbar, :dbar, 1.0)
    add_row(:d, :udbar_to_udbar, :dbar, 1.0)
    add_row(:d, :us_to_us, :s, 1.0)
    add_row(:d, :usbar_to_usbar, :sbar, 1.0)

    # s
    add_row(:s, :us_to_us, :u, 2.0)
    add_row(:s, :usbar_to_usbar, :ubar, 2.0)
    add_row(:s, :ss_to_ss, :s, 1.0)
    add_row(:s, :ssbar_to_ssbar, :sbar, 1.0)
    add_row(:s, :ssbar_to_uubar, :sbar, 2.0)

    # anti-u / anti-d (isospin)
    add_row(:ubar, :uubar_to_uubar, :u, 1.0)
    add_row(:ubar, :uubar_to_ddbar, :u, 1.0)
    add_row(:ubar, :uubar_to_ssbar, :u, 1.0)
    add_row(:ubar, :dubar_to_dubar, :u, 1.0)
    add_row(:ubar, :ubardbar_to_ubardbar, :ubar, 1.0)
    add_row(:ubar, :ubarubar_to_ubarubar, :ubar, 1.0)
    add_row(:ubar, :subar_to_subar, :s, 1.0)
    add_row(:ubar, :ubarsbar_to_ubarsbar, :sbar, 1.0)

    add_row(:dbar, :uubar_to_uubar, :d, 1.0)
    add_row(:dbar, :uubar_to_ddbar, :d, 1.0)
    add_row(:dbar, :uubar_to_ssbar, :d, 1.0)
    add_row(:dbar, :dubar_to_dubar, :d, 1.0)
    add_row(:dbar, :ubardbar_to_ubardbar, :dbar, 1.0)
    add_row(:dbar, :ubarubar_to_ubarubar, :dbar, 1.0)
    add_row(:dbar, :subar_to_subar, :s, 1.0)
    add_row(:dbar, :ubarsbar_to_ubarsbar, :sbar, 1.0)

    # anti-s
    add_row(:sbar, :usbar_to_usbar, :u, 2.0)
    add_row(:sbar, :ubarsbar_to_ubarsbar, :ubar, 2.0)
    add_row(:sbar, :sbarsbar_to_sbarsbar, :sbar, 1.0)
    add_row(:sbar, :ssbar_to_ssbar, :s, 1.0)
    add_row(:sbar, :ssbar_to_uubar, :s, 2.0)

    totals = Dict{Symbol, Float64}()
    for row in rows
        totals[row.species] = get(totals, row.species, 0.0) + row.contribution
    end

    out = NamedTuple[]
    for row in rows
        push!(out, (
            species=row.species,
            channel=row.channel,
            density_key=row.density_key,
            multiplicity=row.multiplicity,
            density=row.density,
            rate=row.rate,
            contribution=row.contribution,
            total=get(totals, row.species, 0.0),
        ))
    end
    return out
end

function write_channel_diagnostics_rows!(io, T_mev::Float64, muq_mev::Float64, muB_mev::Float64, xi::Float64,
    dens, rates, tauinv, eq_backend, diag)
    rows = if isdefined(TransportWorkflow.RelaxationTime, :relaxation_rate_contribution_rows)
        TransportWorkflow.RelaxationTime.relaxation_rate_contribution_rows(dens, rates)
    else
        _fallback_relaxation_rate_contribution_rows(dens, rates)
    end
    for row in rows
        species = row.species
        tauinv_species = getproperty(tauinv, species)
        line = join([
            string(T_mev), string(muq_mev), string(muB_mev), string(xi),
            string(species), string(row.channel), string(row.density_key), string(row.multiplicity),
            string(row.density), string(row.rate), string(row.contribution), string(row.total), string(tauinv_species),
            string(eq_backend), string(diag.phase_curr), string(diag.phase_structure),
        ], ',')
        println(io, line)
    end
    flush(io)
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function current_git_commit()
    try
        return readchomp(`git -C $(PROJECT_ROOT) rev-parse HEAD`)
    catch
        return "unknown"
    end
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
        "sigma_over_T",
        "sigma_over_T_over_eta_over_s",
        "zeta_over_s_over_eta_over_s",
        "quality_flag",
        "quality_reason",
        "quality_metric",
        "run_id",
        "equilibrium_backend",
        "seed_source",
        "phase_prev",
        "phase_curr",
        "phase_structure",
        "phase_boundary_xi_used",
    )
    for c in required
        occursin(c, header) || error("existing output CSV header is incompatible with current script (missing column: $c). Please rerun with --overwrite or choose a new --output path.")
    end
end

function write_header_if_needed(io)
    header = join([
        "T_MeV", "muq_MeV", "muB_MeV", "xi",
        "T_fm", "muq_fm",
        "converged", "iterations", "residual_norm", "equilibrium_backend", "seed_source", "phase_prev", "phase_curr", "phase_structure", "phase_boundary_xi_used",
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
        "sigma_over_T", "sigma_over_T_over_eta_over_s", "zeta_over_s_over_eta_over_s",
        "quality_flag", "quality_reason", "quality_metric", "run_id",
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

mutable struct LocalPhaseTracker
    boundary_data::Union{Nothing, NamedTuple{(:T_values, :mu_values, :T_CEP, :mu_CEP, :xi), Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64, Float64}}}
    previous_solution::Union{Nothing, Vector{Float64}}
    previous_phase::Symbol
    hadron_seed::Vector{Float64}
    quark_seed::Vector{Float64}
end

function load_phase_boundary_data(xi::Float64; boundary_path::String=DEFAULT_PHASE_BOUNDARY_PATH, cep_path::String=DEFAULT_PHASE_CEP_PATH)
    T_CEP = NaN
    mu_CEP = NaN
    if isfile(cep_path)
        for line in eachline(cep_path)
            startswith(line, "xi") && continue
            parts = split(line, ',')
            length(parts) >= 3 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            abs(xi_val - xi) > 1e-6 && continue
            T_CEP = tryparse(Float64, parts[2])
            mu_CEP = tryparse(Float64, parts[3])
            break
        end
    end

    T_values = Float64[]
    mu_values = Float64[]
    if isfile(boundary_path)
        for line in eachline(boundary_path)
            startswith(line, "xi") && continue
            parts = split(line, ',')
            length(parts) >= 3 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            abs(xi_val - xi) > 1e-6 && continue
            T_val = tryparse(Float64, parts[2])
            mu_val = tryparse(Float64, parts[3])
            (T_val === nothing || mu_val === nothing) && continue
            push!(T_values, T_val)
            push!(mu_values, mu_val)
        end
    end

    if !isempty(T_values)
        order = sortperm(T_values)
        T_values = T_values[order]
        mu_values = mu_values[order]
    end

    return (T_values=T_values, mu_values=mu_values, T_CEP=T_CEP, mu_CEP=mu_CEP, xi=Float64(xi))
end

function interpolate_boundary_mu_c(data, T_mev::Float64)
    if !isnan(data.T_CEP) && T_mev > data.T_CEP
        return NaN
    end
    isempty(data.T_values) && return NaN
    Ts = data.T_values
    μs = data.mu_values
    if T_mev <= Ts[1]
        return μs[1]
    elseif T_mev >= Ts[end]
        return μs[end]
    end
    for i in 1:(length(Ts)-1)
        if Ts[i] <= T_mev <= Ts[i+1]
            t = (T_mev - Ts[i]) / (Ts[i+1] - Ts[i])
            return μs[i] + t * (μs[i+1] - μs[i])
        end
    end
    return NaN
end

function current_phase_hint(tracker::LocalPhaseTracker, T_mev::Float64, muq_mev::Float64)
    data = tracker.boundary_data
    data === nothing && return :unknown
    if !isnan(data.T_CEP) && T_mev > data.T_CEP
        return :crossover
    end
    μ_c_mev = interpolate_boundary_mu_c(data, T_mev)
    isnan(μ_c_mev) && return :unknown
    return muq_mev < μ_c_mev ? :hadron : :quark
end

@inline function local_is_phase_transition(prev_phase::Symbol, current_phase::Symbol)
    return (prev_phase == :hadron && current_phase == :quark) ||
           (prev_phase == :quark && current_phase == :hadron)
end

function tracker_seed(tracker::LocalPhaseTracker, T_fm::Float64, muq_fm::Float64)
    mode = Models.FixedMu()
    if tracker.previous_solution !== nothing
        if length(tracker.previous_solution) == Models.state_dim(mode)
            return copy(tracker.previous_solution)
        elseif length(tracker.previous_solution) >= 5
            return Models.extend_seed(tracker.previous_solution, mode)
        end
    end
    hint = Models.auto_phase_hint(T_fm, muq_fm)
    base = (hint === :quark) ? tracker.quark_seed : tracker.hadron_seed
    return Models.extend_seed(base, mode)
end

function build_phase_tracker(xi::Float64, previous_solution, previous_phase::Symbol)
    boundary_xi_used = xi
    boundary_data = try
        load_phase_boundary_data(xi; boundary_path=DEFAULT_PHASE_BOUNDARY_PATH, cep_path=DEFAULT_PHASE_CEP_PATH)
    catch
        nothing
    end
    tracker = LocalPhaseTracker(boundary_data, nothing, previous_phase, copy(TransportWorkflow.PNJL.HADRON_SEED_5), copy(TransportWorkflow.PNJL.QUARK_SEED_5))
    if tracker.boundary_data !== nothing && isempty(tracker.boundary_data.T_values)
        nearest_xi = nearest_phase_boundary_xi(xi)
        if nearest_xi !== nothing
            nearest_data = try
                load_phase_boundary_data(nearest_xi; boundary_path=DEFAULT_PHASE_BOUNDARY_PATH, cep_path=DEFAULT_PHASE_CEP_PATH)
            catch
                tracker.boundary_data
            end
            tracker.boundary_data = nearest_data
            boundary_xi_used = nearest_xi
        end
    end
    if previous_solution !== nothing
        tracker.previous_solution = collect(Float64, previous_solution)
        tracker.previous_phase = previous_phase
    end
    return tracker, boundary_xi_used
end

function available_phase_boundary_xis()
    if PHASE_BOUNDARY_XI_CACHE[] !== nothing
        return PHASE_BOUNDARY_XI_CACHE[]
    end
    xis = Float64[]
    if isfile(DEFAULT_PHASE_BOUNDARY_PATH)
        for line in eachline(DEFAULT_PHASE_BOUNDARY_PATH)
            startswith(line, "xi") && continue
            parts = split(line, ',')
            length(parts) >= 1 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            push!(xis, xi_val)
        end
    end
    PHASE_BOUNDARY_XI_CACHE[] = unique(sort(xis))
    return PHASE_BOUNDARY_XI_CACHE[]
end

function nearest_phase_boundary_xi(xi::Float64)
    xis = available_phase_boundary_xis()
    isempty(xis) && return nothing
    distances = abs.(xis .- xi)
    return xis[argmin(distances)]
end

function available_phase_crossover_xis()
    if PHASE_CROSSOVER_XI_CACHE[] !== nothing
        return PHASE_CROSSOVER_XI_CACHE[]
    end
    xis = Float64[]
    if isfile(DEFAULT_PHASE_CROSSOVER_PATH)
        for line in eachline(DEFAULT_PHASE_CROSSOVER_PATH)
            startswith(line, "xi") && continue
            parts = split(line, ',')
            length(parts) >= 1 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            push!(xis, xi_val)
        end
    end
    PHASE_CROSSOVER_XI_CACHE[] = unique(sort(xis))
    return PHASE_CROSSOVER_XI_CACHE[]
end

function nearest_phase_crossover_xi(xi::Float64)
    xis = available_phase_crossover_xis()
    isempty(xis) && return nothing
    distances = abs.(xis .- xi)
    return xis[argmin(distances)]
end

function load_crossover_reference(xi::Float64)
    isfile(DEFAULT_PHASE_CROSSOVER_PATH) || return nothing, nothing

    header = nothing
    rows = NamedTuple{(:mu_MeV, :T_crossover_MeV), Tuple{Float64, Float64}}[]
    xi_used = xi
    exact_match_found = false
    nearest_xi = nearest_phase_crossover_xi(xi)

    open(DEFAULT_PHASE_CROSSOVER_PATH, "r") do io
        for raw_line in eachline(io)
            line = strip(raw_line)
            isempty(line) && continue
            if header === nothing
                header = split(line, ',')
                continue
            end

            parts = split(line, ',')
            length(parts) == length(header) || continue
            idx_xi = findfirst(==("xi"), header)
            idx_mu = findfirst(==("mu_MeV"), header)
            idx_T = findfirst(==("T_crossover_MeV"), header)
            if idx_T === nothing
                idx_T = findfirst(==("T_crossover_chiral_MeV"), header)
            end
            (idx_xi === nothing || idx_mu === nothing || idx_T === nothing) && return nothing, nothing

            row_xi = tryparse(Float64, parts[idx_xi])
            row_mu = tryparse(Float64, parts[idx_mu])
            row_T = tryparse(Float64, parts[idx_T])
            (row_xi === nothing || row_mu === nothing || row_T === nothing) && continue

            if abs(row_xi - xi) <= 1e-6
                exact_match_found = true
                push!(rows, (mu_MeV=Float64(row_mu), T_crossover_MeV=Float64(row_T)))
                xi_used = Float64(row_xi)
            elseif !exact_match_found && nearest_xi !== nothing && abs(row_xi - nearest_xi) <= 1e-6
                push!(rows, (mu_MeV=Float64(row_mu), T_crossover_MeV=Float64(row_T)))
                xi_used = Float64(row_xi)
            end
        end
    end

    filtered = filter(row -> isfinite(row.mu_MeV) && isfinite(row.T_crossover_MeV), rows)
    isempty(filtered) && return nothing, nothing
    sort!(filtered, by=row -> row.mu_MeV)
    return filtered, xi_used
end

function interpolate_crossover_temperature(xi::Float64, muq_mev::Float64)
    data, xi_used = load_crossover_reference(xi)
    data === nothing && return NaN, xi_used
    length(data) == 1 && return data[1].T_crossover_MeV, xi_used

    if muq_mev <= data[1].mu_MeV
        return data[1].T_crossover_MeV, xi_used
    elseif muq_mev >= data[end].mu_MeV
        return data[end].T_crossover_MeV, xi_used
    end

    for i in 1:(length(data) - 1)
        left = data[i]
        right = data[i + 1]
        if left.mu_MeV <= muq_mev <= right.mu_MeV
            weight = (muq_mev - left.mu_MeV) / (right.mu_MeV - left.mu_MeV)
            return left.T_crossover_MeV + weight * (right.T_crossover_MeV - left.T_crossover_MeV), xi_used
        end
    end

    return NaN, xi_used
end

function tracker_phase(tracker::LocalPhaseTracker, T_mev::Float64, muq_mev::Float64, xi::Float64)
    boundary_phase = current_phase_hint(tracker, T_mev, muq_mev)

    if boundary_phase in (:hadron, :quark)
        return boundary_phase
    end

    Tc_mev, _ = interpolate_crossover_temperature(xi, muq_mev)
    if isfinite(Tc_mev)
        if abs(T_mev - Tc_mev) <= 2.0
            return :crossover
        elseif T_mev < Tc_mev
            return :hadron
        else
            return :quark
        end
    end

    if boundary_phase == :crossover
        return :quark
    end

    return boundary_phase
end

function is_phase_transition(prev_phase::Symbol, current_phase::Symbol)
    return local_is_phase_transition(prev_phase, current_phase)
end

function describe_seed_source(tracker, current_phase::Symbol)
    if tracker.previous_solution === nothing
        return "phase_aware_default_$(String(current_phase))"
    elseif is_phase_transition(tracker.previous_phase, current_phase)
        return "phase_aware_phase_switch_$(String(tracker.previous_phase))_to_$(String(current_phase))"
    else
        return "phase_aware_continuity_$(String(current_phase))"
    end
end

function phase_structure(tracker, T_mev::Float64, muq_mev::Float64, xi::Float64)
    data = tracker.boundary_data
    if data !== nothing && !isempty(data.T_values) && !isnan(data.T_CEP) && T_mev <= data.T_CEP
        return :first_order_possible
    end

    Tc_mev, _ = interpolate_crossover_temperature(xi, muq_mev)
    if isfinite(Tc_mev)
        return T_mev <= Tc_mev ? :crossover_possible : :no_transition
    end

    return :unknown
end

@inline function _seed_state_5(seed_state)
    return Float64.(seed_state[1:min(5, length(seed_state))])
end

function _normalize_equilibrium_result(raw; solver_backend::Symbol=:models)
    Bool(raw.converged) || return nothing
    x_state = SVector{5}(Tuple(Float64.(raw.x_state[1:5])))
    mu_vec = SVector{3}(Tuple(Float64.(raw.mu_vec[1:3])))
    masses = SVector{3}(Tuple(Float64.(raw.masses[1:3])))
    (all(isfinite, masses) && all(>(0.0), masses)) || return nothing
    return (
        converged=true,
        x_state=x_state,
        mu_vec=mu_vec,
        masses=masses,
        iterations=Int(raw.iterations),
        residual_norm=Float64(raw.residual_norm),
        solver_backend=solver_backend,
        omega=Float64(raw.omega),
    )
end

function _solve_fixedmu_via_models_solve(T_fm::Float64, muq_fm::Float64, xi::Float64, seed_state, opts::ScanOptions)
    return Models.solve(
        PNJL_MODEL,
        Models.FixedMu(),
        T_fm,
        muq_fm;
        seed_guess=_seed_state_5(seed_state),
        xi=xi,
        p_num=opts.p_num,
        t_num=opts.t_num,
        residual_norm_max=1e-4,
    )
end

function _solve_fixedmu_via_models_constraint(T_fm::Float64, muq_fm::Float64, xi::Float64, seed_state, opts::ScanOptions)
    return Models.solve_constraint(
        PNJL_MODEL,
        Models.FixedMu(),
        T_fm;
        μ_fm=muq_fm,
        seed_guess=_seed_state_5(seed_state),
        xi=xi,
        p_num=opts.p_num,
        t_num=opts.t_num,
        residual_norm_max=1e-4,
        physicality_check=Models.is_physical_solution,
    )
end

function solve_models_equilibrium(T_fm::Float64, muq_fm::Float64, xi::Float64, seed_state, opts::ScanOptions)
    models_err = nothing
    try
        raw = _solve_fixedmu_via_models_solve(T_fm, muq_fm, xi, seed_state, opts)
        eq = _normalize_equilibrium_result(raw; solver_backend=:models)
        eq !== nothing && return eq
    catch err
        models_err = err
    end

    try
        raw = _solve_fixedmu_via_models_constraint(T_fm, muq_fm, xi, seed_state, opts)
        return _normalize_equilibrium_result(raw; solver_backend=:models)
    catch
        if models_err !== nothing
            rethrow(models_err)
        end
        rethrow()
    end
end

# 多初值池治理求解：使用 6 种子池 + governance 选优，不依赖 boundary 数据
# 与旧的 2 种子 omega 比较不同，solve_multi 使用 physicality → residual → pressure 优先级
function solve_with_multiseed_governance(T_fm::Float64, muq_fm::Float64, xi::Float64, opts::ScanOptions)
    result = try
        Models.solve_multi(
            PNJL_MODEL,
            Models.FixedMu(),
            T_fm,
            muq_fm;
            xi=xi,
            p_num=opts.p_num,
            t_num=opts.t_num,
            residual_norm_max=1e-4,
            evaluate_all_attempts=true,
        )
    catch
        nothing
    end
    result === nothing && return nothing
    return _normalize_equilibrium_result(result; solver_backend=:models_multi)
end

# 从解的状态推断物相（强子/夸克），作为 boundary 数据的备选
function classify_phase_from_solution(eq)
    m_u = eq.masses[1]
    Phi = eq.x_state[4]
    # 夸克相判据：质量较轻（手征部分恢复）或有明显的 Polyakov loop（退禁闭）
    if m_u < 0.8 || Phi > 0.1
        return :quark
    else
        return :hadron
    end
end

function solve_equilibrium_with_diagnostics(T_mev::Float64, muB_mev::Float64, xi::Float64, opts::ScanOptions;
    previous_solution=nothing,
    previous_phase::Symbol=:unknown,
)
    T_fm = T_mev / ħc_MeV_fm
    muq_mev = muB_mev / 3.0
    muq_fm = muq_mev / ħc_MeV_fm

    tracker, boundary_xi_used = build_phase_tracker(xi, previous_solution, previous_phase)
    phase_prev = tracker.previous_phase
    phase_curr_hint = tracker_phase(tracker, T_mev, muq_mev, xi)
    structure = phase_structure(tracker, T_mev, muq_mev, xi)
    seed_source = describe_seed_source(tracker, phase_curr_hint)
    seed_state = tracker_seed(tracker, T_fm, muq_fm)

    eq = nothing
    phase_curr = phase_curr_hint
    models_err = nothing

    if tracker.previous_solution !== nothing && is_phase_transition(phase_prev, phase_curr_hint) && phase_curr_hint in (:hadron, :quark)
        eq_multi = solve_with_multiseed_governance(T_fm, muq_fm, xi, opts)
        if eq_multi !== nothing
            eq = eq_multi
            # 物相标签：优先使用 boundary 判定（现已准确），仅在其不可用时回退到解的分类
            phase_curr = if phase_curr_hint in (:hadron, :quark)
                phase_curr_hint
            else
                classify_phase_from_solution(eq)
            end
            seed_source = "phase_aware_multiseed_governance_$(String(phase_curr))"
        end
    end

    if eq === nothing
        try
            eq = solve_models_equilibrium(T_fm, muq_fm, xi, seed_state, opts)
        catch err
            models_err = err
        end
    end

    if eq === nothing
        if models_err !== nothing
            rethrow(models_err)
        end
        throw(ArgumentError("models equilibrium solver returned no valid candidate"))
    end

    next_solution = collect(Float64, eq.x_state)
    next_phase = phase_curr
    return eq, next_solution, next_phase, (
        seed_source=seed_source,
        phase_prev=phase_prev,
        phase_curr=next_phase,
        phase_structure=structure,
        phase_boundary_xi_used=boundary_xi_used,
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

function run_scan(opts::ScanOptions, ctx::ProvenanceMetadata.RunContext)
    ensure_parent_dir(opts.output)
    if opts.channel_diagnostics_output !== nothing
        ensure_parent_dir(opts.channel_diagnostics_output)
    end
    if opts.failed_points_output !== nothing
        ensure_parent_dir(opts.failed_points_output)
    end
    provenance_dir = isnothing(opts.provenance_dir) ? dirname(opts.output) : opts.provenance_dir
    mkpath(provenance_dir)

    if opts.resume && isfile(opts.output) && !opts.overwrite
        ensure_output_header_compatible(opts.output)
    end

    existing = opts.resume && isfile(opts.output) && !opts.overwrite ? read_existing_keys(opts.output) : Set{Tuple{Float64,Float64,Float64}}()

    if opts.overwrite && isfile(opts.output)
        rm(opts.output)
    end
    if opts.channel_diagnostics_output !== nothing && opts.overwrite && isfile(opts.channel_diagnostics_output)
        rm(opts.channel_diagnostics_output)
    end
    if opts.failed_points_output !== nothing && opts.overwrite && isfile(opts.failed_points_output)
        rm(opts.failed_points_output)
    end

    new_file = !isfile(opts.output) || filesize(opts.output) == 0
    p_grid, p_w, sigma_cutoff = integration_grids(opts)
    cos_grid, cos_w = gauleg(-1.0, 1.0, opts.tau_angle_nodes)
    phi_grid, phi_w = gauleg(0.0, 2 * pi, opts.tau_phi_nodes)
    solver_kwargs = (iterations=opts.max_iter,)
    io = open(opts.output, "a")
    channel_io = opts.channel_diagnostics_output === nothing ? nothing : open(opts.channel_diagnostics_output, "a")
    failed_io = opts.failed_points_output === nothing ? nothing : open(opts.failed_points_output, "a")
    stats_success = 0
    stats_error = 0
    stats_skipped = 0
    try
        if new_file
            ScanCSV.write_metadata(io, Dict(
                "schema" => "scan_csv_v1",
                "title" => "gap_transport_scan",
                "script" => "scripts/relaxtime/run_gap_transport_scan.jl",
                "git_commit" => current_git_commit(),
                "provenance.entrypoint" => "workflow",
                "provenance.equilibrium_backend" => "models.solve_constraint(FixedMu) with phase-aware xi/T continuity",
                "provenance.tau_path" => "TransportWorkflow.solve_gap_and_transport",
                "provenance.integration_mode" => string(opts.integration_mode),
                "sigma_grid_n" => string(opts.sigma_grid_n),
                "integration_mode" => string(opts.integration_mode),
                "gc_every_n" => string(opts.gc_every_n),
                "tau_p_nodes" => string(opts.tau_p_nodes),
                "tau_angle_nodes" => string(opts.tau_angle_nodes),
                "tau_phi_nodes" => string(opts.tau_phi_nodes),
                "tau_n_sigma_points" => string(opts.tau_n_sigma_points),
                "tau_threshold_subtraction" => string(opts.tau_threshold_subtraction),
                "tau_asym_window" => string(opts.tau_asym_window),
                "tau_asym_fit_min_points" => string(opts.tau_asym_fit_min_points),
                "tau_asym_extra_points" => string(opts.tau_asym_extra_points),
                "tau_interpolation_mode" => string(opts.tau_interpolation_mode),
                "note.tau_threshold_hint" => "for near-threshold sharp channels, linear+threshold_subtraction often more robust than pchip",
                "tr_p_nodes" => string(opts.tr_p_nodes),
                "tr_p_max_fm" => string(opts.tr_p_max_fm),

                # labels for plotting convenience
                "y_label.sigma_over_T" => "σ/T",
                "y_label.sigma_over_T_over_eta_over_s" => "(σ/T)/(η/s)",
                "y_label.zeta_over_s_over_eta_over_s" => "(ζ/s)/(η/s)",
            ))
            write_header_if_needed(io)
        end
        if channel_io !== nothing
            channel_new_file = !isfile(opts.channel_diagnostics_output) || filesize(opts.channel_diagnostics_output) == 0
            if channel_new_file
                ScanCSV.write_metadata(channel_io, Dict(
                    "schema" => "scan_csv_v1",
                    "title" => "gap_transport_scan_channel_diagnostics",
                    "script" => "scripts/relaxtime/run_gap_transport_scan.jl",
                    "git_commit" => current_git_commit(),
                    "source_csv" => opts.output,
                ))
                write_channel_diagnostics_header_if_needed(channel_io)
            end
        end
        if failed_io !== nothing
            failed_new_file = !isfile(opts.failed_points_output) || filesize(opts.failed_points_output) == 0
            if failed_new_file
                ScanCSV.write_metadata(failed_io, Dict(
                    "schema" => "scan_csv_v1",
                    "title" => "gap_transport_scan_failed_points",
                    "script" => "scripts/relaxtime/run_gap_transport_scan.jl",
                    "git_commit" => current_git_commit(),
                    "source_csv" => opts.output,
                ))
                write_failed_points_header_if_needed(failed_io)
            end
        end

        T_values = collect(range(opts.tmin_mev; stop=opts.tmax_mev, step=opts.tstep_mev))
        muB_values = collect(range(opts.mubmin_mev; stop=opts.mubmax_mev, step=opts.mubstep_mev))
        muB_values = unique(sort(Float64.(muB_values)))
        xi_continuity_mode = length(T_values) == 1 && length(muB_values) == 1 && length(opts.xi_values) > 1

        total = length(opts.xi_values) * length(T_values) * length(muB_values)
        done = 0

        if xi_continuity_mode
            T_mev = T_values[1]
            muB_mev = muB_values[1]
            previous_solution = nothing
            previous_phase = :unknown
            for xi in opts.xi_values
                done += 1
                key = (T_mev, muB_mev, xi)
                if opts.resume && (key in existing)
                    stats_skipped += 1
                    continue
                end

                try
                    eq, previous_solution, previous_phase, diag = solve_equilibrium_with_diagnostics(
                        T_mev,
                        muB_mev,
                        xi,
                        opts;
                        previous_solution=previous_solution,
                        previous_phase=previous_phase,
                    )

                    T_fm = T_mev / ħc_MeV_fm
                    muq_mev = muB_mev / 3.0
                    muq_fm = muq_mev / ħc_MeV_fm
                    seed_state = Vector(eq.x_state)

                    Φ = Float64(eq.x_state[4])
                    Φbar = Float64(eq.x_state[5])
                    masses = (u=Float64(eq.masses[1]), d=Float64(eq.masses[2]), s=Float64(eq.masses[3]))

                    ktmp = build_K_data(Float64(T_fm), Float64(muq_fm), masses, Φ, Φbar)

                    _compute_bulk_this_point = opts.compute_bulk
                    res = try
                        TransportWorkflow.solve_gap_and_transport(
                            T_fm,
                            muq_fm;
                            xi=xi,
                            equilibrium=eq,
                            compute_tau=true,
                            K_coeffs=ktmp.K_coeffs,
                            tau=nothing,
                            compute_bulk=_compute_bulk_this_point,
                            p_num=opts.p_num,
                            t_num=opts.t_num,
                            seed_state=seed_state,
                            solver_kwargs=solver_kwargs,
                            tau_kwargs=(
                                p_nodes=opts.tau_p_nodes,
                                angle_nodes=opts.tau_angle_nodes,
                                phi_nodes=opts.tau_phi_nodes,
                                n_sigma_points=opts.tau_n_sigma_points,
                                p_grid=p_grid,
                                p_w=p_w,
                                cos_grid=cos_grid,
                                cos_w=cos_w,
                                phi_grid=phi_grid,
                                phi_w=phi_w,
                                sigma_cutoff=sigma_cutoff,
                                threshold_subtraction=opts.tau_threshold_subtraction,
                                asym_window=opts.tau_asym_window,
                                asym_fit_min_points=opts.tau_asym_fit_min_points,
                                asym_extra_points=opts.tau_asym_extra_points,
                                interpolation_mode=opts.tau_interpolation_mode,
                            ),
                            transport_config=TransportWorkflow.TransportIntegrationConfig(p_nodes=opts.tr_p_nodes, p_max=opts.tr_p_max_fm),
                        )
                    catch bulk_err
                        if _compute_bulk_this_point
                            @warn "transport with bulk failed, retrying without bulk" T_mev=T_mev muB_mev=muB_mev xi=xi err=bulk_err
                            _compute_bulk_this_point = false
                            TransportWorkflow.solve_gap_and_transport(
                                T_fm,
                                muq_fm;
                                xi=xi,
                                equilibrium=eq,
                                compute_tau=true,
                                K_coeffs=ktmp.K_coeffs,
                                tau=nothing,
                                compute_bulk=false,
                                p_num=opts.p_num,
                                t_num=opts.t_num,
                                seed_state=seed_state,
                                solver_kwargs=solver_kwargs,
                                tau_kwargs=(
                                    p_nodes=opts.tau_p_nodes,
                                    angle_nodes=opts.tau_angle_nodes,
                                    phi_nodes=opts.tau_phi_nodes,
                                    n_sigma_points=opts.tau_n_sigma_points,
                                    p_grid=p_grid,
                                    p_w=p_w,
                                    cos_grid=cos_grid,
                                    cos_w=cos_w,
                                    phi_grid=phi_grid,
                                    phi_w=phi_w,
                                    sigma_cutoff=sigma_cutoff,
                                    threshold_subtraction=opts.tau_threshold_subtraction,
                                    asym_window=opts.tau_asym_window,
                                    asym_fit_min_points=opts.tau_asym_fit_min_points,
                                    asym_extra_points=opts.tau_asym_extra_points,
                                    interpolation_mode=opts.tau_interpolation_mode,
                                ),
                                transport_config=TransportWorkflow.TransportIntegrationConfig(p_nodes=opts.tr_p_nodes, p_max=opts.tr_p_max_fm),
                            )
                        else
                            rethrow(bulk_err)
                        end
                    end

                    dens = res.densities
                    tau = res.tau
                    tauinv = res.tau_inv
                    tr = res.transport
                    rates = res.rates

                    P_fm4inv, _, s_fm3inv, epsilon_fm4inv = Models.model_thermo(
                        PNJL_MODEL,
                        eq.x_state,
                        eq.mu_vec,
                        T_fm,
                        p_num=opts.p_num,
                        t_num=opts.t_num,
                        xi=xi,
                    )

                    rho_quark_net = (dens.u - dens.ubar) + (dens.d - dens.dbar) + (dens.s - dens.sbar)
                    rho_baryon = rho_quark_net / 3.0
                    rho_norm = rho_baryon / ρ0_inv_fm3

                    omega_fm4inv = Models.omega(
                        PNJL_MODEL,
                        eq.x_state,
                        T_fm,
                        eq.mu_vec;
                        p_num=opts.p_num,
                        t_num=opts.t_num,
                        xi=xi,
                    )
                    omega_MeV_fm3 = omega_fm4inv * ħc_MeV_fm
                    P_MeV_fm3 = P_fm4inv * ħc_MeV_fm
                    epsilon_MeV_fm3 = epsilon_fm4inv * ħc_MeV_fm
                    eps_minus_3P_over_T4 = (isfinite(epsilon_fm4inv) && isfinite(P_fm4inv) && isfinite(T_fm) && T_fm != 0.0) ? ((epsilon_fm4inv - 3.0 * P_fm4inv) / T_fm^4) : NaN
                    eta_over_s = (isfinite(tr.eta) && isfinite(s_fm3inv) && s_fm3inv != 0.0) ? (tr.eta / s_fm3inv) : NaN
                    zeta_over_s = (isfinite(tr.zeta) && isfinite(s_fm3inv) && s_fm3inv != 0.0) ? (tr.zeta / s_fm3inv) : NaN
                    sigma_over_T = (isfinite(tr.sigma) && isfinite(T_fm) && T_fm != 0.0) ? (tr.sigma / T_fm) : NaN
                    sigma_over_T_over_eta_over_s = (isfinite(sigma_over_T) && isfinite(eta_over_s) && eta_over_s != 0.0) ? (sigma_over_T / eta_over_s) : NaN
                    zeta_over_s_over_eta_over_s = (isfinite(zeta_over_s) && isfinite(eta_over_s) && eta_over_s != 0.0) ? (zeta_over_s / eta_over_s) : NaN
                    quality_flag, quality_reason, quality_metric = assess_point_quality(tau, eta_over_s, sigma_over_T)

                    row = join([
                        string(T_mev), string(muq_mev), string(muB_mev), string(xi),
                        string(T_fm), string(muq_fm),
                        csv_bool(eq.converged), string(eq.iterations), string(eq.residual_norm), string(eq.solver_backend), diag.seed_source, string(diag.phase_prev), string(diag.phase_curr), string(diag.phase_structure), string(diag.phase_boundary_xi_used),
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
                        string(sigma_over_T), string(sigma_over_T_over_eta_over_s), string(zeta_over_s_over_eta_over_s),
                        csv_bool(quality_flag), quality_reason, string(quality_metric), string(ctx.run_id),
                    ], ',')
                    println(io, row)
                    flush(io)

                    if channel_io !== nothing && rates !== nothing
                        write_channel_diagnostics_rows!(channel_io, T_mev, muq_mev, muB_mev, xi,
                            dens, rates, tauinv, eq.solver_backend, diag)
                    end
                    stats_success += 1
                catch point_err
                    @warn "SCAN POINT FAILED — skipped" T_mev=T_mev muB_mev=muB_mev xi=xi err=point_err
                    stats_error += 1
                    if failed_io !== nothing
                        diag_or_hint = (seed_source="unknown", phase_prev=previous_phase, phase_curr=:unknown)
                        write_failed_point_row!(failed_io, T_mev, muB_mev, xi, diag_or_hint, point_err)
                    end
                end

                if opts.gc_every_n > 0 && done % opts.gc_every_n == 0
                    GC.gc()
                end

                if done % 10 == 0
                    println("progress: $(done)/$(total) (last muB=$(muB_mev) MeV, T=$(T_mev) MeV, xi=$(xi))")
                end
            end
        else
            for xi in opts.xi_values
                for muB_mev in muB_values
                    previous_solution = nothing
                    previous_phase = :unknown
                    for T_mev in T_values
                    done += 1
                    key = (T_mev, muB_mev, xi)
                    if opts.resume && (key in existing)
                        stats_skipped += 1
                        continue
                    end

                    try  # 单点容错：失败不中断后续扫描

                    eq, previous_solution, previous_phase, diag = solve_equilibrium_with_diagnostics(
                        T_mev,
                        muB_mev,
                        xi,
                        opts;
                        previous_solution=previous_solution,
                        previous_phase=previous_phase,
                    )

                    T_fm = T_mev / ħc_MeV_fm
                    muq_mev = muB_mev / 3.0
                    muq_fm = muq_mev / ħc_MeV_fm

                    seed_state = Vector(eq.x_state)

                    Φ = Float64(eq.x_state[4])
                    Φbar = Float64(eq.x_state[5])
                    masses = (u=Float64(eq.masses[1]), d=Float64(eq.masses[2]), s=Float64(eq.masses[3]))

                    ktmp = build_K_data(Float64(T_fm), Float64(muq_fm), masses, Φ, Φbar)
                    thermo_params = (T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(xi))
                    quark_params = (m=masses, μ=(u=Float64(muq_fm), d=Float64(muq_fm), s=Float64(muq_fm)), A=ktmp.A_vals)

                    cs_caches = nothing
                    """
                    if Bool(base.converged)
                        cs_caches = build_sigma_caches(REQUIRED_PROCESSES, quark_params, thermo_params, ktmp.K_coeffs;
                            n_sigma_points=opts.tau_n_sigma_points, sigma_grid_n=opts.sigma_grid_n)
                    end
                    """
                    # 正式计算：τ + 输运（η/σ）; ζ 默认关
                    # 容错策略：先尝试含 bulk，若 bulk 计算抛错则回退到只算 η/σ
                    _compute_bulk_this_point = opts.compute_bulk
                    res = try
                        TransportWorkflow.solve_gap_and_transport(
                            T_fm,
                            muq_fm;
                            xi=xi,
                            equilibrium=eq,
                            compute_tau=true,
                            K_coeffs=ktmp.K_coeffs,
                            tau=nothing,
                            compute_bulk=_compute_bulk_this_point,
                            p_num=opts.p_num,
                            t_num=opts.t_num,
                            seed_state=seed_state,
                            solver_kwargs=solver_kwargs,
                            tau_kwargs=(
                                p_nodes=opts.tau_p_nodes,
                                angle_nodes=opts.tau_angle_nodes,
                                phi_nodes=opts.tau_phi_nodes,
                                n_sigma_points=opts.tau_n_sigma_points,
                                #cs_caches=cs_caches,
                                p_grid=p_grid,
                                p_w=p_w,
                                cos_grid=cos_grid,
                                cos_w=cos_w,
                                phi_grid=phi_grid,
                                phi_w=phi_w,
                                sigma_cutoff=sigma_cutoff,
                                threshold_subtraction=opts.tau_threshold_subtraction,
                                asym_window=opts.tau_asym_window,
                                asym_fit_min_points=opts.tau_asym_fit_min_points,
                                asym_extra_points=opts.tau_asym_extra_points,
                                interpolation_mode=opts.tau_interpolation_mode,
                            ),
                            transport_config=TransportWorkflow.TransportIntegrationConfig(p_nodes=opts.tr_p_nodes, p_max=opts.tr_p_max_fm),
                        )
                    catch bulk_err
                        if _compute_bulk_this_point
                            @warn "transport with bulk failed, retrying without bulk" T_mev=T_mev muB_mev=muB_mev xi=xi err=bulk_err
                            _compute_bulk_this_point = false
                            TransportWorkflow.solve_gap_and_transport(
                                T_fm,
                                muq_fm;
                                xi=xi,
                                equilibrium=eq,
                                compute_tau=true,
                                K_coeffs=ktmp.K_coeffs,
                                tau=nothing,
                                compute_bulk=false,
                                p_num=opts.p_num,
                                t_num=opts.t_num,
                                seed_state=seed_state,
                                solver_kwargs=solver_kwargs,
                                tau_kwargs=(
                                    p_nodes=opts.tau_p_nodes,
                                    angle_nodes=opts.tau_angle_nodes,
                                    phi_nodes=opts.tau_phi_nodes,
                                    n_sigma_points=opts.tau_n_sigma_points,
                                    p_grid=p_grid,
                                    p_w=p_w,
                                    cos_grid=cos_grid,
                                    cos_w=cos_w,
                                    phi_grid=phi_grid,
                                    phi_w=phi_w,
                                    sigma_cutoff=sigma_cutoff,
                                    threshold_subtraction=opts.tau_threshold_subtraction,
                                    asym_window=opts.tau_asym_window,
                                    asym_fit_min_points=opts.tau_asym_fit_min_points,
                                    asym_extra_points=opts.tau_asym_extra_points,
                                    interpolation_mode=opts.tau_interpolation_mode,
                                ),
                                transport_config=TransportWorkflow.TransportIntegrationConfig(p_nodes=opts.tr_p_nodes, p_max=opts.tr_p_max_fm),
                            )
                        else
                            rethrow(bulk_err)
                        end
                    end

                    eq = res.equilibrium
                    dens = res.densities
                    tau = res.tau
                    tauinv = res.tau_inv
                    tr = res.transport
                    rates = res.rates

                    P_fm4inv, _, s_fm3inv, epsilon_fm4inv = Models.model_thermo(
                        PNJL_MODEL,
                        eq.x_state,
                        eq.mu_vec,
                        T_fm,
                        p_num=opts.p_num,
                        t_num=opts.t_num,
                        xi=xi,
                    )

                    # 重建重子数密度（旧版 eq.rho/eq.rho_norm 已移除）
                    rho_quark_net = (dens.u - dens.ubar) + (dens.d - dens.dbar) + (dens.s - dens.sbar)
                    rho_baryon = rho_quark_net / 3.0
                    rho_norm = rho_baryon / ρ0_inv_fm3

                    omega_fm4inv = Models.omega(
                        PNJL_MODEL,
                        eq.x_state,
                        T_fm,
                        eq.mu_vec;
                        p_num=opts.p_num,
                        t_num=opts.t_num,
                        xi=xi,
                    )
                    omega_MeV_fm3 = omega_fm4inv * ħc_MeV_fm
                    P_MeV_fm3 = P_fm4inv * ħc_MeV_fm
                    epsilon_MeV_fm3 = epsilon_fm4inv * ħc_MeV_fm
                    eps_minus_3P_over_T4 = (isfinite(epsilon_fm4inv) && isfinite(P_fm4inv) && isfinite(T_fm) && T_fm != 0.0) ? ((epsilon_fm4inv - 3.0 * P_fm4inv) / T_fm^4) : NaN
                    eta_over_s = (isfinite(tr.eta) && isfinite(s_fm3inv) && s_fm3inv != 0.0) ? (tr.eta / s_fm3inv) : NaN
                    zeta_over_s = (isfinite(tr.zeta) && isfinite(s_fm3inv) && s_fm3inv != 0.0) ? (tr.zeta / s_fm3inv) : NaN

                    sigma_over_T = (isfinite(tr.sigma) && isfinite(T_fm) && T_fm != 0.0) ? (tr.sigma / T_fm) : NaN
                    sigma_over_T_over_eta_over_s = (isfinite(sigma_over_T) && isfinite(eta_over_s) && eta_over_s != 0.0) ? (sigma_over_T / eta_over_s) : NaN
                    zeta_over_s_over_eta_over_s = (isfinite(zeta_over_s) && isfinite(eta_over_s) && eta_over_s != 0.0) ? (zeta_over_s / eta_over_s) : NaN
                    quality_flag, quality_reason, quality_metric = assess_point_quality(tau, eta_over_s, sigma_over_T)

                    row = join([
                        string(T_mev), string(muq_mev), string(muB_mev), string(xi),
                        string(T_fm), string(muq_fm),
                        csv_bool(eq.converged), string(eq.iterations), string(eq.residual_norm), string(eq.solver_backend), diag.seed_source, string(diag.phase_prev), string(diag.phase_curr), string(diag.phase_structure), string(diag.phase_boundary_xi_used),
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
                        string(sigma_over_T), string(sigma_over_T_over_eta_over_s), string(zeta_over_s_over_eta_over_s),
                        csv_bool(quality_flag), quality_reason, string(quality_metric), string(ctx.run_id),
                    ], ',')
                    println(io, row)
                    flush(io)

                    if channel_io !== nothing && rates !== nothing
                        write_channel_diagnostics_rows!(channel_io, T_mev, muq_mev, muB_mev, xi,
                            dens, rates, tauinv, eq.solver_backend, diag)
                    end
                    stats_success += 1

                    catch point_err
                        @warn "SCAN POINT FAILED — skipped" T_mev=T_mev muB_mev=muB_mev xi=xi err=point_err
                        if failed_io !== nothing
                            diag_or_hint = (seed_source="unknown", phase_prev=previous_phase, phase_curr=:unknown)
                            write_failed_point_row!(failed_io, T_mev, muB_mev, xi, diag_or_hint, point_err)
                        end
                        stats_error += 1
                    end  # try 单点容错

                    if opts.gc_every_n > 0 && done % opts.gc_every_n == 0
                        GC.gc()
                    end

                    if done % 10 == 0
                        println("progress: $(done)/$(total) (last muB=$(muB_mev) MeV, T=$(T_mev) MeV, xi=$(xi))")
                    end
                end
            end
        end
        end
    finally
        close(io)
        if channel_io !== nothing
            close(channel_io)
        end
        if failed_io !== nothing
            close(failed_io)
        end
    end

    effective_config = Dict{String,Any}(
        "output" => opts.output,
        "channel_diagnostics_output" => opts.channel_diagnostics_output,
        "failed_points_output" => opts.failed_points_output,
        "xi_values" => opts.xi_values,
        "tmin_mev" => opts.tmin_mev,
        "tmax_mev" => opts.tmax_mev,
        "tstep_mev" => opts.tstep_mev,
        "mubmin_mev" => opts.mubmin_mev,
        "mubmax_mev" => opts.mubmax_mev,
        "mubstep_mev" => opts.mubstep_mev,
        "overwrite" => opts.overwrite,
        "resume" => opts.resume,
        "compute_bulk" => opts.compute_bulk,
        "p_num" => opts.p_num,
        "t_num" => opts.t_num,
        "max_iter" => opts.max_iter,
        "tau_p_nodes" => opts.tau_p_nodes,
        "tau_angle_nodes" => opts.tau_angle_nodes,
        "tau_phi_nodes" => opts.tau_phi_nodes,
        "tau_n_sigma_points" => opts.tau_n_sigma_points,
        "tau_threshold_subtraction" => opts.tau_threshold_subtraction,
        "tau_asym_window" => opts.tau_asym_window,
        "tau_asym_fit_min_points" => opts.tau_asym_fit_min_points,
        "tau_asym_extra_points" => opts.tau_asym_extra_points,
        "tau_interpolation_mode" => String(opts.tau_interpolation_mode),
        "sigma_grid_n" => opts.sigma_grid_n,
        "integration_mode" => String(opts.integration_mode),
        "gc_every_n" => opts.gc_every_n,
        "tr_p_nodes" => opts.tr_p_nodes,
        "tr_p_max_fm" => opts.tr_p_max_fm,
    )

    summary = Dict{String,Any}(
        "points_total" => stats_success + stats_error,
        "success_count" => stats_success,
        "error_count" => stats_error,
        "skipped_count" => stats_skipped,
    )

    artifacts = String[opts.output]
    if opts.channel_diagnostics_output !== nothing
        push!(artifacts, opts.channel_diagnostics_output)
    end
    if opts.failed_points_output !== nothing
        push!(artifacts, opts.failed_points_output)
    end

    ProvenanceMetadata.write_run_sidecars(
        provenance_dir;
        ctx=ctx,
        effective_config=effective_config,
        artifacts=artifacts,
        summary=summary,
    )

    println("Scan finished. Output: $(opts.output)")
end

function main()
    opts = parse_args(copy(ARGS))
    ctx = ProvenanceMetadata.new_run_context("scripts/relaxtime/run_gap_transport_scan.jl", copy(ARGS))
    run_scan(opts, ctx)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
