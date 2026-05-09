"""
固定 `mu_B = 0` 的 canonical meson thermo 温度扫描，输出 phase-shift EOS 结果与最小结果资产。

主链约束：
- 统一通过 `Models.solve_gap_and_phase_shift_meson_thermo_point` 进入；
- CSV 平铺字段统一复用 `Models.build_meson_thermo_contract_row`；
- 默认 canonical case 收敛到 `pi/sigma_pi` + `phase_shift_current`；
- 扫描层只负责编排、落盘与 provenance sidecar，不在脚本层手工拼总热力学派生量。

结果资产（默认 outdir）：
- `scan.csv`
- `README.md`
- `effective_config.json`
- `run_manifest.json`
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "provenance_metadata.jl"))
include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .ScanCSV: ScanCSV
using .ProvenanceMetadata: ProvenanceMetadata
using .Constants_PNJL: ħc_MeV_fm
using .Models: build_meson_thermo_contract_row, solve_gap_and_phase_shift_meson_thermo_point

Base.@kwdef struct ScanOptions
    outdir::String=joinpath(
        PROJECT_ROOT,
        "data", "outputs", "results", "relaxtime", "meson_thermo",
        "canonical_muB0_phase_shift_current_pi_sigma_pi",
    )
    scheme::Symbol=:current
    primary_channel::Symbol=:pi
    secondary_channel::Symbol=:sigma_pi
    xi::Float64=0.0
    tmin_mev::Float64=140.0
    tmax_mev::Float64=260.0
    tstep_mev::Float64=5.0
    overwrite::Bool=false
    resume::Bool=true
    p_num::Int=12
    t_num::Int=6
    max_iter::Int=40
    qmax::Float64=12.0
    q_nodes::Int=48
    omega_min::Float64=0.05
    omega_max::Float64=10.0
    omega_nodes::Int=48
    eta::Float64=1e-6
    ld_cutoff::Union{Nothing,Float64}=nothing
    ld_cutoff_mode::Symbol=:match_qmax
    ld_threshold_mode::Symbol=:omega_lt_q
    allow_legacy_fd_fallback::Bool=false
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl [options]\n")
    println("Options:")
    println("  --outdir <path>               输出目录（默认 canonical_muB0_phase_shift_current_pi_sigma_pi）")
    println("  --scheme <current|gbu>        相移口径 (default current)")
    println("  --mesons <a,b>                双通道介子集合 (default pi,sigma_pi)")
    println("  --xi <value>                  各向异性参数 ξ (default 0)")
    println("  --tmin/--tmax/--tstep <MeV>   温度范围与步长")
    println("  --overwrite                   覆盖现有 case 输出")
    println("  --no-resume                   禁用按 key 跳过逻辑")
    println("  --p-num <int>                 平衡态求解动量节点数 (default 12)")
    println("  --t-num <int>                 平衡态求解角度节点数 (default 6)")
    println("  --max-iter <int>              平衡态/介子求解迭代上限 (default 40)")
    println("  --qmax <fm^-1>                外层 q 硬截断 (default 12)")
    println("  --q-nodes <int>               外层 q 积分节点数 (default 48)")
    println("  --omega-min <fm^-1>           内层 omega 下限 (default 0.05)")
    println("  --omega-max <fm^-1>           内层 omega 上限 (default 10)")
    println("  --omega-nodes <int>           内层 omega 节点数 (default 48)")
    println("  --eta <float>                 相移诊断小虚部参数 (default 1e-6)")
    println("  --ld-cutoff <fm^-1>           显式 LD cutoff；未给出时沿 workflow 默认治理")
    println("  --ld-cutoff-mode <symbol>     LD cutoff 治理口径 (default match_qmax)")
    println("  --ld-threshold-mode <symbol>  QP/LD 分区口径 (default omega_lt_q)")
    println("  --allow-legacy-fd-fallback    允许 AD 失败时回落 workflow 私有差分")
    println("  -h, --help                    显示帮助")
end

function _parse_channels(raw::AbstractString)
    parts = [Symbol(strip(seg)) for seg in split(raw, ',') if !isempty(strip(seg))]
    length(parts) == 2 || error("--mesons must contain exactly two channels, got: $(raw)")
    return Tuple(parts)
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}()

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            val = args[i + 1]
            i += 1
            return val
        end

        if arg == "--outdir"
            opts[:outdir] = require_value()
        elseif arg == "--scheme"
            raw = lowercase(strip(require_value()))
            opts[:scheme] = raw == "gbu" ? :gbu_reference : :current
        elseif arg == "--mesons"
            opts[:channels] = _parse_channels(require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--tmin"
            opts[:tmin_mev] = parse(Float64, require_value())
        elseif arg == "--tmax"
            opts[:tmax_mev] = parse(Float64, require_value())
        elseif arg == "--tstep"
            opts[:tstep_mev] = parse(Float64, require_value())
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
        elseif arg == "--ld-cutoff"
            opts[:ld_cutoff] = parse(Float64, require_value())
        elseif arg == "--ld-cutoff-mode"
            opts[:ld_cutoff_mode] = Symbol(strip(require_value()))
        elseif arg == "--ld-threshold-mode"
            opts[:ld_threshold_mode] = Symbol(strip(require_value()))
        elseif arg == "--allow-legacy-fd-fallback"
            opts[:allow_legacy_fd_fallback] = true
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    base = ScanOptions()
    channels = get(opts, :channels, (base.primary_channel, base.secondary_channel))
    resolved = ScanOptions(
        outdir=String(get(opts, :outdir, base.outdir)),
        scheme=Symbol(get(opts, :scheme, base.scheme)),
        primary_channel=channels[1],
        secondary_channel=channels[2],
        xi=Float64(get(opts, :xi, base.xi)),
        tmin_mev=Float64(get(opts, :tmin_mev, base.tmin_mev)),
        tmax_mev=Float64(get(opts, :tmax_mev, base.tmax_mev)),
        tstep_mev=Float64(get(opts, :tstep_mev, base.tstep_mev)),
        overwrite=Bool(get(opts, :overwrite, base.overwrite)),
        resume=Bool(get(opts, :resume, base.resume)),
        p_num=Int(get(opts, :p_num, base.p_num)),
        t_num=Int(get(opts, :t_num, base.t_num)),
        max_iter=Int(get(opts, :max_iter, base.max_iter)),
        qmax=Float64(get(opts, :qmax, base.qmax)),
        q_nodes=Int(get(opts, :q_nodes, base.q_nodes)),
        omega_min=Float64(get(opts, :omega_min, base.omega_min)),
        omega_max=Float64(get(opts, :omega_max, base.omega_max)),
        omega_nodes=Int(get(opts, :omega_nodes, base.omega_nodes)),
        eta=Float64(get(opts, :eta, base.eta)),
        ld_cutoff=get(opts, :ld_cutoff, base.ld_cutoff),
        ld_cutoff_mode=Symbol(get(opts, :ld_cutoff_mode, base.ld_cutoff_mode)),
        ld_threshold_mode=Symbol(get(opts, :ld_threshold_mode, base.ld_threshold_mode)),
        allow_legacy_fd_fallback=Bool(get(opts, :allow_legacy_fd_fallback, base.allow_legacy_fd_fallback)),
    )

    resolved.tstep_mev > 0.0 || error("tstep must be positive")
    resolved.tmax_mev >= resolved.tmin_mev || error("tmax must be >= tmin")
    resolved.qmax > 0.0 || error("qmax must be positive")
    resolved.q_nodes > 1 || error("q-nodes must be > 1")
    resolved.omega_nodes > 1 || error("omega-nodes must be > 1")
    resolved.omega_max > resolved.omega_min || error("omega-max must exceed omega-min")
    return resolved
end

@inline function _format(x)
    x isa Bool && return x ? "true" : "false"
    x isa Symbol && return String(x)
    return string(x)
end

function _row_to_values(cols::Vector{String}, row::Dict{String,Any})
    return [_format(get(row, c, "")) for c in cols]
end

function _output_columns()
    return [
        "T_MeV", "muB_MeV", "xi",
        "T_fm", "muB_fm", "mu_fm",
        "workflow", "channel_set", "primary_channel", "secondary_channel",
        "Phi", "Phibar",
        "m_u", "m_d", "m_s",
        "m_primary", "m_secondary",
        "P_meson", "P_meson_qp", "P_meson_ld",
        "P_quark_meanfield", "P_total",
        "entropy", "epsilon", "trace_anomaly",
        "P_meson_over_P_total", "delta_P_vs_no_meson",
        "P_pi_qp", "P_pi_ld", "P_K_qp", "P_K_ld",
        "P_primary", "P_secondary", "P_primary_qp", "P_primary_ld", "P_secondary_qp", "P_secondary_ld",
        "equilibrium_converged", "phase_structure",
        "phase_shift_variant", "thermo_derivation_mode",
        "ld_cutoff", "ld_cutoff_mode", "ld_threshold_mode",
        "scheme",
        "qmax", "q_nodes", "omega_min", "omega_max", "omega_nodes", "eta",
    ]
end

function _effective_config(opts::ScanOptions, out_csv::String)
    return Dict{String,Any}(
        "schema_version" => "v1",
        "profile_name" => "canonical_muB0_phase_shift_meson_thermo",
        "workflow_entry" => "Models.solve_gap_and_phase_shift_meson_thermo_point",
        "output_csv" => replace(relpath(out_csv, PROJECT_ROOT), '\\' => '/'),
        "muB_MeV" => 0.0,
        "xi" => opts.xi,
        "scheme" => String(opts.scheme),
        "channels" => [String(opts.primary_channel), String(opts.secondary_channel)],
        "temperature_scan" => Dict{String,Any}(
            "T_min_MeV" => opts.tmin_mev,
            "T_max_MeV" => opts.tmax_mev,
            "T_step_MeV" => opts.tstep_mev,
        ),
        "solver" => Dict{String,Any}(
            "p_num" => opts.p_num,
            "t_num" => opts.t_num,
            "max_iter" => opts.max_iter,
            "allow_legacy_fd_fallback" => opts.allow_legacy_fd_fallback,
        ),
        "phase_shift" => Dict{String,Any}(
            "qmax" => opts.qmax,
            "q_nodes" => opts.q_nodes,
            "omega_min" => opts.omega_min,
            "omega_max" => opts.omega_max,
            "omega_nodes" => opts.omega_nodes,
            "eta" => opts.eta,
            "ld_cutoff" => opts.ld_cutoff,
            "ld_cutoff_mode" => String(opts.ld_cutoff_mode),
            "ld_threshold_mode" => String(opts.ld_threshold_mode),
        ),
    )
end

function _write_readme(path::String, opts::ScanOptions, out_csv::String, stats)
    open(path, "w") do io
        println(io, "# canonical mu_B = 0 phase-shift meson thermo case")
        println(io)
        println(io, "- workflow_entry: `Models.solve_gap_and_phase_shift_meson_thermo_point`")
        println(io, "- scheme: `$(opts.scheme)`")
        println(io, "- channel_set: `$(opts.primary_channel),$(opts.secondary_channel)`")
        println(io, "- mu_B: `0 MeV`")
        println(io, "- xi: `$(opts.xi)`")
        println(io, "- temperature_window: `$(opts.tmin_mev)` to `$(opts.tmax_mev)` MeV, step `$(opts.tstep_mev)` MeV")
        println(io, "- solver_grid: `p_num=$(opts.p_num)`, `t_num=$(opts.t_num)`, `max_iter=$(opts.max_iter)`")
        println(io, "- phase_shift_grid: `qmax=$(opts.qmax) fm^-1`, `q_nodes=$(opts.q_nodes)`, `omega=[$(opts.omega_min), $(opts.omega_max)] fm^-1`, `omega_nodes=$(opts.omega_nodes)`, `eta=$(opts.eta)`")
        println(io, "- AD_fallback_allowed: `$(opts.allow_legacy_fd_fallback)`")
        println(io, "- points_total: `$(stats.points_total)`")
        println(io, "- success_count: `$(stats.success_count)`")
        println(io, "- error_count: `$(stats.error_count)`")
        println(io)
        println(io, "## Artifacts")
        println(io)
        println(io, "- `scan.csv`: canonical temperature scan output")
        println(io, "- `effective_config.json`: effective run config snapshot")
        println(io, "- `run_manifest.json`: provenance metadata")
        println(io)
        println(io, "## Notes")
        println(io)
        println(io, "- 当前 `primary/secondary` 字段是解释双通道结果的首选字段。")
        println(io, "- 当第二通道为 `sigma_pi` 时，兼容列 `P_K_*` 仅表示“第二通道 legacy 槽位”，不应按 kaon 物理解读。")
        println(io, "- 本 case 当前只收口 canonical scan/README/provenance；图资产未在本脚本内生成。")
        println(io, "- 主结果文件：`$(replace(relpath(out_csv, dirname(path)), '\\' => '/'))`")
    end
end

function main()
    opts = parse_args(ARGS)
    outdir = opts.outdir
    out_csv = joinpath(outdir, "scan.csv")
    readme_path = joinpath(outdir, "README.md")

    mkpath(outdir)

    if opts.overwrite
        for path in (out_csv, readme_path, joinpath(outdir, "effective_config.json"), joinpath(outdir, "run_manifest.json"))
            isfile(path) && rm(path)
        end
    end

    output_columns = _output_columns()
    key_cols = ["T_MeV", "muB_MeV", "xi"]
    existing = (isfile(out_csv) && opts.resume && !opts.overwrite) ? ScanCSV.read_existing_keys(out_csv, key_cols) : Set{Tuple{Vararg{Float64}}}()

    if isfile(out_csv)
        ScanCSV.assert_required_columns(out_csv, output_columns)
    end

    is_new = !isfile(out_csv)
    open(out_csv, is_new ? "w" : "a") do io
        if is_new
            ScanCSV.write_metadata(io, Dict(
                "format" => "scan_csv_v1",
                "script" => "scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl",
                "workflow_entry" => "Models.solve_gap_and_phase_shift_meson_thermo_point",
                "canonical_case" => "muB=0",
                "phase_shift_scheme" => String(opts.scheme),
                "channel_set" => string(opts.primary_channel, ",", opts.secondary_channel),
                "thermo_derivation" => "Omega_total -> Models.model_thermo -> ForwardDiff",
                "note" => "Prefer primary/secondary columns when second channel is sigma_pi; legacy P_K_* columns map to the secondary slot.",
            ))
            ScanCSV.write_header(io, output_columns)
        end

        muB = 0.0
        muB_fm = 0.0
        mu_fm = 0.0
        continuation_state = nothing
        T = opts.tmin_mev

        while T <= opts.tmax_mev + 1e-9
            key = (Float64(T), muB, Float64(opts.xi))
            if key in existing
                T += opts.tstep_mev
                continue
            end

            T_fm = T / ħc_MeV_fm
            point = solve_gap_and_phase_shift_meson_thermo_point(
                T_fm,
                mu_fm;
                xi=opts.xi,
                mesons=(opts.primary_channel, opts.secondary_channel),
                continuation_state=continuation_state,
                mixed_branch_align=:strict_sign_binding,
                p_num=opts.p_num,
                t_num=opts.t_num,
                solver_kwargs=(; iterations=opts.max_iter),
                mass_kwargs=(; iterations=opts.max_iter),
                thermo_kwargs=(;
                    pi_channel=opts.primary_channel,
                    k_channel=opts.secondary_channel,
                    scheme=opts.scheme,
                    qmax=opts.qmax,
                    q_nodes=opts.q_nodes,
                    omega_min=opts.omega_min,
                    omega_max=opts.omega_max,
                    omega_nodes=opts.omega_nodes,
                    eta=opts.eta,
                    ld_cutoff=opts.ld_cutoff,
                    ld_cutoff_mode=opts.ld_cutoff_mode,
                    ld_threshold_mode=opts.ld_threshold_mode,
                    p_num=opts.p_num,
                    t_num=opts.t_num,
                ),
                allow_legacy_fd_fallback=opts.allow_legacy_fd_fallback,
            )

            row = build_meson_thermo_contract_row(point)
            qp = point.quark_params
            tp = point.thermo_params
            thermo = point.phase_shift_meson_thermo
            output_row = Dict{String,Any}(
                "T_MeV" => row.T_MeV,
                "muB_MeV" => row.muB_MeV,
                "xi" => opts.xi,
                "T_fm" => T_fm,
                "muB_fm" => muB_fm,
                "mu_fm" => mu_fm,
                "workflow" => row.workflow,
                "channel_set" => replace(row.channel_set, "," => "|"),
                "primary_channel" => row.primary_channel,
                "secondary_channel" => row.secondary_channel,
                "Phi" => tp.Φ,
                "Phibar" => tp.Φbar,
                "m_u" => qp.m.u,
                "m_d" => qp.m.d,
                "m_s" => qp.m.s,
                "m_primary" => hasproperty(thermo, :m_pi) ? Float64(thermo.m_pi) : NaN,
                "m_secondary" => hasproperty(thermo, :m_K) ? Float64(thermo.m_K) : NaN,
                "P_meson" => row.P_meson,
                "P_meson_qp" => row.P_meson_qp,
                "P_meson_ld" => row.P_meson_ld,
                "P_quark_meanfield" => row.P_quark_meanfield,
                "P_total" => row.P_total,
                "entropy" => row.entropy,
                "epsilon" => row.epsilon,
                "trace_anomaly" => row.trace_anomaly,
                "P_meson_over_P_total" => row.P_meson_over_P_total,
                "delta_P_vs_no_meson" => row.delta_P_vs_no_meson,
                "P_pi_qp" => row.P_pi_qp,
                "P_pi_ld" => row.P_pi_ld,
                "P_K_qp" => row.P_K_qp,
                "P_K_ld" => row.P_K_ld,
                "P_primary" => row.P_primary,
                "P_secondary" => row.P_secondary,
                "P_primary_qp" => row.P_primary_qp,
                "P_primary_ld" => row.P_primary_ld,
                "P_secondary_qp" => row.P_secondary_qp,
                "P_secondary_ld" => row.P_secondary_ld,
                "equilibrium_converged" => row.equilibrium_converged,
                "phase_structure" => row.phase_structure,
                "phase_shift_variant" => row.phase_shift_variant,
                "thermo_derivation_mode" => row.thermo_derivation_mode,
                "ld_cutoff" => row.ld_cutoff,
                "ld_cutoff_mode" => row.ld_cutoff_mode,
                "ld_threshold_mode" => row.ld_threshold_mode,
                "scheme" => String(opts.scheme),
                "qmax" => opts.qmax,
                "q_nodes" => opts.q_nodes,
                "omega_min" => opts.omega_min,
                "omega_max" => opts.omega_max,
                "omega_nodes" => opts.omega_nodes,
                "eta" => opts.eta,
            )

            println(io, join(_row_to_values(output_columns, output_row), ','))
            flush(io)

            continuation_state = hasproperty(point, :continuation_state) ? point.continuation_state : continuation_state
            T += opts.tstep_mev
        end
    end

    stats = ProvenanceMetadata.csv_stats(out_csv; converged_col="equilibrium_converged")
    _write_readme(readme_path, opts, out_csv, stats)
    ctx = ProvenanceMetadata.new_run_context("scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl", copy(ARGS))
    ProvenanceMetadata.write_run_sidecars(
        outdir;
        ctx=ctx,
        effective_config=_effective_config(opts, out_csv),
        artifacts=[out_csv, readme_path],
        summary=Dict{String,Any}(
            "points_total" => stats.points_total,
            "success_count" => stats.success_count,
            "error_count" => stats.error_count,
            "workflow_entry" => "Models.solve_gap_and_phase_shift_meson_thermo_point",
            "canonical_case" => "muB0",
        ),
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
