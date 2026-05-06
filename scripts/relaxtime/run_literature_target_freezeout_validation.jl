"""
对单个 literature target 执行 charged/freeze-out meson-density workflow 扫描，
并生成：

1. workflow 原始输出 CSV
2. 与 target 对齐的比较 CSV
3. 一页最小结果说明

当前用途：
- 首个样本：`Blaschke:2019col` Figure 4 right panel
- 目标曲线：`K^-/pi^-`, `mu_pi = 100 MeV`

治理定位：
- 当前是 workflow reproduction / staged comparison
- 不是 strict literature reconstruction
- 不承担正式 validation pass/fail 门槛
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Models: run_freezeout_meson_density_scan

const DEFAULT_TARGET_CSV = joinpath(
    "data", "outputs", "results", "relaxtime", "literature", "meson_density_targets",
    "blaschke2019col_kminus_piminus_mu_pi_100_fig4_right.csv",
)

const DEFAULT_OUTPUT_DIR = joinpath(
    "data", "outputs", "results", "relaxtime", "meson_density", "freezeout_validation",
    "blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default",
)

struct ValidationOptions
    target_csv::String
    output_dir::String
    freezeout_profile::String
    path_profile::String
    flavor_chemical_profile::String
    meson_chemical_profile::String
    regime::Symbol
    overwrite::Bool
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
    println("Usage: julia --project=. scripts/relaxtime/run_literature_target_freezeout_validation.jl [options]\n")
    println("Options:")
    println("  --target-csv <path>                 literature target CSV")
    println("  --output-dir <path>                 输出目录")
    println("  --freezeout-profile <name>          freeze-out profile (default default)")
    println("  --path-profile <name>               path profile (default baseline_freezeout)")
    println("  --flavor-chemical-profile <name>    flavor chemical profile (default default)")
    println("  --meson-chemical-profile <name>     meson chemical profile (default blaschke2019_mu_pi_100)")
    println("  --regime <phase_shift_gbu|phase_shift_current|stable|strict_bw_stage1|strict_bw_stage2>")
    println("  --overwrite                         覆盖 scan/comparison/summary")
    println("  --p-num <int>                       平衡态动量节点数 (default 24)")
    println("  --t-num <int>                       平衡态角度节点数 (default 8)")
    println("  --max-iter <int>                    solver/mass iterations 上限 (default 40)")
    println("  --qmax <float>                      q 上限 (default 12)")
    println("  --q-nodes <int>                     q 节点数 (default 48)")
    println("  --omega-min <float>                 phase-shift omega 下限 (default 0.05)")
    println("  --omega-max <float>                 omega 上限 (default 10)")
    println("  --omega-nodes <int>                 omega 节点数 (default 48)")
    println("  --eta <float>                       phase-shift eta (default 1e-6)")
    println("  -h, --help                          显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :target_csv => DEFAULT_TARGET_CSV,
        :output_dir => DEFAULT_OUTPUT_DIR,
        :freezeout_profile => "default",
        :path_profile => "baseline_freezeout",
        :flavor_chemical_profile => "default",
        :meson_chemical_profile => "blaschke2019_mu_pi_100",
        :regime => :phase_shift_gbu,
        :overwrite => false,
        :p_num => 24,
        :t_num => 8,
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

        if arg == "--target-csv"
            opts[:target_csv] = require_value()
        elseif arg == "--output-dir"
            opts[:output_dir] = require_value()
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
        elseif arg == "--overwrite"
            opts[:overwrite] = true
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
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return ValidationOptions(
        String(opts[:target_csv]),
        String(opts[:output_dir]),
        String(opts[:freezeout_profile]),
        String(opts[:path_profile]),
        String(opts[:flavor_chemical_profile]),
        String(opts[:meson_chemical_profile]),
        Symbol(opts[:regime]),
        Bool(opts[:overwrite]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Float64(opts[:qmax]),
        Int(opts[:q_nodes]),
        Float64(opts[:omega_min]),
        Float64(opts[:omega_max]),
        Int(opts[:omega_nodes]),
        Float64(opts[:eta]),
    )
end

function read_xy_csv(path::String)
    lines = readlines(path)
    isempty(lines) && error("target CSV is empty: $path")
    strip(lines[1]) == "x_value,y_value" || error("unexpected target CSV header: $(lines[1])")

    xs = Float64[]
    ys = Float64[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        length(cols) == 2 || error("invalid target row: $line")
        push!(xs, parse(Float64, strip(cols[1])))
        push!(ys, parse(Float64, strip(cols[2])))
    end
    isempty(xs) && error("target CSV has no data rows: $path")
    return xs, ys
end

function read_scan_csv(path::String)
    lines = readlines(path)
    isempty(lines) && error("scan CSV is empty: $path")
    header = split(strip(lines[1]), ',')
    index = Dict(name => i for (i, name) in enumerate(header))

    required = ("sqrt_s_NN_GeV", "T_MeV", "muB_MeV", "kpi_ratio", "equilibrium_converged", "message")
    for key in required
        haskey(index, key) || error("scan CSV missing column: $key")
    end

    rows = Dict{Float64, NamedTuple}()
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        x = key_x(parse(Float64, cols[index["sqrt_s_NN_GeV"]]))
        rows[x] = (
            T_MeV=parse(Float64, cols[index["T_MeV"]]),
            muB_MeV=parse(Float64, cols[index["muB_MeV"]]),
            kpi_ratio=parse(Float64, cols[index["kpi_ratio"]]),
            converged=strip(cols[index["equilibrium_converged"]]) == "true",
            message=strip(cols[index["message"]]),
        )
    end
    return rows
end

fmt(x::Real) = string(x)
key_x(x::Real) = round(Float64(x); digits=6)

function write_comparison_csv(path::String, xs::Vector{Float64}, ys::Vector{Float64}, scan_rows::Dict{Float64, NamedTuple})
    open(path, "w") do io
        println(io, "x_value,target_y,model_y,abs_diff,rel_diff,T_MeV,muB_MeV,converged,message")
        for (x, y) in zip(xs, ys)
            row = get(scan_rows, key_x(x), nothing)
            if row === nothing
                println(io, join((fmt(x), fmt(y), "NaN", "NaN", "NaN", "NaN", "NaN", "false", "\"missing scan row\""), ','))
                continue
            end
            model_y = row.kpi_ratio
            abs_diff = abs(model_y - y)
            rel_diff = iszero(y) ? NaN : abs_diff / abs(y)
            message = replace(row.message, '"' => '\'')
            println(io, join((
                fmt(x),
                fmt(y),
                fmt(model_y),
                fmt(abs_diff),
                fmt(rel_diff),
                fmt(row.T_MeV),
                fmt(row.muB_MeV),
                string(row.converged),
                "\"$(message)\"",
            ), ','))
        end
    end
end

function summarize_comparison(xs::Vector{Float64}, ys::Vector{Float64}, scan_rows::Dict{Float64, NamedTuple})
    model_values = Float64[]
    abs_diffs = Float64[]
    rel_diffs = Float64[]
    converged_count = 0
    missing_count = 0
    for (x, y) in zip(xs, ys)
        row = get(scan_rows, key_x(x), nothing)
        if row === nothing
            missing_count += 1
            continue
        end
        row.converged && (converged_count += 1)
        push!(model_values, row.kpi_ratio)
        ad = abs(row.kpi_ratio - y)
        push!(abs_diffs, ad)
        if !iszero(y)
            push!(rel_diffs, ad / abs(y))
        end
    end

    return (
        points=length(xs),
        missing=missing_count,
        converged=converged_count,
        target_min=minimum(ys),
        target_max=maximum(ys),
        model_min=isempty(model_values) ? NaN : minimum(model_values),
        model_max=isempty(model_values) ? NaN : maximum(model_values),
        mean_abs_diff=isempty(abs_diffs) ? NaN : sum(abs_diffs) / length(abs_diffs),
        max_abs_diff=isempty(abs_diffs) ? NaN : maximum(abs_diffs),
        mean_rel_diff=isempty(rel_diffs) ? NaN : sum(rel_diffs) / length(rel_diffs),
        max_rel_diff=isempty(rel_diffs) ? NaN : maximum(rel_diffs),
    )
end

function write_summary_md(path::String, opts::ValidationOptions, target_csv::String, scan_csv::String, comparison_csv::String, summary)
    open(path, "w") do io
        println(io, "# charged/freeze-out literature workflow reproduction 最小结果说明")
        println(io)
        println(io, "更新日期：$(Dates.today())")
        println(io)
        println(io, "## 1. 运行口径")
        println(io)
        println(io, "- target: `$(target_csv)`")
        println(io, "- freezeout profile: `$(opts.freezeout_profile)`")
        println(io, "- path profile: `$(opts.path_profile)`")
        println(io, "- flavor chemical profile: `$(opts.flavor_chemical_profile)`")
        println(io, "- meson chemical profile: `$(opts.meson_chemical_profile)`")
        println(io, "- regime: `$(opts.regime)`")
        println(io, "- workflow output: `$(scan_csv)`")
        println(io, "- comparison output: `$(comparison_csv)`")
        println(io)
        println(io, "## 2. 当前结果定位")
        println(io)
        println(io, "这是一条 **workflow reproduction / staged comparison**，不是 strict literature reconstruction，也不是正式 validation gate。")
        println(io)
        println(io, "当前固定的是：")
        println(io)
        println(io, "- `sqrt(s_NN)` 直接取自 literature target 的横轴采样点")
        println(io, "- `sqrt(s_NN) -> (T, mu_B)` 使用本仓库 freeze-out baseline profile")
        println(io, "- path strategy 通过显式 path profile 作用在 baseline 之上")
        println(io, "- flavor-level `mu_u, mu_d, mu_s` 使用显式 flavor chemical profile")
        println(io, "- `mu_pi` 使用显式 meson chemical profile")
        println(io)
        println(io, "## 3. 数值摘要")
        println(io)
        println(io, "- points: `$(summary.points)`")
        println(io, "- missing rows: `$(summary.missing)`")
        println(io, "- converged rows: `$(summary.converged)`")
        println(io, "- target range: `[ $(summary.target_min), $(summary.target_max) ]`")
        println(io, "- model range: `[ $(summary.model_min), $(summary.model_max) ]`")
        println(io, "- mean abs diff: `$(summary.mean_abs_diff)`")
        println(io, "- max abs diff: `$(summary.max_abs_diff)`")
        println(io, "- mean rel diff: `$(summary.mean_rel_diff)`")
        println(io, "- max rel diff: `$(summary.max_rel_diff)`")
        println(io)
        println(io, "## 4. 解读约束")
        println(io)
        println(io, "- 该比较当前主要用于确认 charged/freeze-out workflow 与生产输出链已闭环可运行。")
        println(io, "- 当前不对相对误差设置 pass/fail 门槛，也不把该结果直接纳入 regression gate。")
        println(io, "- 若与文献量级仍有系统差异，应优先把它解释为 path/regime 语义差异，而不是直接回退 workflow 入口契约。")
    end
end

function main()
    opts = parse_args(ARGS)
    target_csv = normpath(joinpath(PROJECT_ROOT, opts.target_csv))
    isfile(target_csv) || error("target CSV not found: $target_csv")

    xs, ys = read_xy_csv(target_csv)
    outdir = normpath(joinpath(PROJECT_ROOT, opts.output_dir))
    mkpath(outdir)

    scan_csv = joinpath(outdir, "workflow_scan.csv")
    comparison_csv = joinpath(outdir, "comparison_vs_target.csv")
    summary_md = joinpath(outdir, "README.md")

    scan_result = run_freezeout_meson_density_scan(
        sqrt_s_NN_values=xs,
        xi_values=[0.0],
        freezeout_profile_name=opts.freezeout_profile,
        path_profile_name=opts.path_profile,
        flavor_chemical_profile_name=opts.flavor_chemical_profile,
        meson_chemical_profile_name=opts.meson_chemical_profile,
        regime=opts.regime,
        output_path=scan_csv,
        traversal=:as_given,
        overwrite=opts.overwrite,
        resume=!opts.overwrite,
        p_num=opts.p_num,
        t_num=opts.t_num,
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

    scan_rows = read_scan_csv(scan_csv)
    write_comparison_csv(comparison_csv, xs, ys, scan_rows)
    summary = summarize_comparison(xs, ys, scan_rows)
    write_summary_md(summary_md, opts, opts.target_csv, relpath(scan_csv, PROJECT_ROOT), relpath(comparison_csv, PROJECT_ROOT), summary)

    println(scan_result)
    println("comparison_csv=$(relpath(comparison_csv, PROJECT_ROOT))")
    println("summary_md=$(relpath(summary_md, PROJECT_ROOT))")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    using Dates
    main()
end
