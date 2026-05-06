"""
对单个 crossover-line literature target 执行 charged meson-density workflow 扫描，
并生成：

1. workflow 原始输出 CSV
2. 与 target 对齐的 comparison CSV（按 T/mu_B 线性插值）
3. 一页最小结果说明

当前首要用途：
- Friesen et al. 2019 Fig.5 top panel
- K+ / pi+ 与 K- / pi- 两条 charged ratio 曲线
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Models: run_crossover_meson_density_scan
using Dates

const DEFAULT_TARGET_CSV = raw"D:\Desktop\paper\dev\outputs\formalized\friesen2019_kpi_ratio_curves\friesen2019_fig5_top_charged_ratio_curves_mu055.csv"

struct ValidationOptions
    target_csv::String
    output_dir::String
    curve_id::String
    flavor_chemical_profile::String
    meson_chemical_profile::String
    regime::Symbol
    mu_min_MeV::Float64
    mu_max_MeV::Float64
    T_min_MeV::Float64
    T_max_MeV::Float64
    n_mu::Int
    xi::Float64
    method::Symbol
    variable::Symbol
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
    println("Usage: julia --project=. scripts/relaxtime/run_literature_target_crossover_validation.jl [options]\n")
    println("Options:")
    println("  --target-csv <path>                 target CSV")
    println("  --output-dir <path>                 output directory")
    println("  --curve-id <kplus_over_piplus|kminus_over_piminus>")
    println("  --flavor-chemical-profile <name>    flavor chemical profile")
    println("  --meson-chemical-profile <name>     meson chemical profile")
    println("  --regime <stable|phase_shift_current|phase_shift_gbu|strict_bw_stage1|strict_bw_stage2>")
    println("  --mu-min <MeV>                      lower mu_q bound")
    println("  --mu-max <MeV>                      upper mu_q bound")
    println("  --T-min <MeV>                       lower T bound for crossover search")
    println("  --T-max <MeV>                       upper T bound for crossover search")
    println("  --n-mu <int>                        mu_q scan points")
    println("  --xi <float>                        anisotropy xi")
    println("  --method <peak|inflection>          crossover locator")
    println("  --variable <symbol>                 crossover variable")
    println("  --overwrite                         overwrite output directory contents")
    println("  --p-num <int>                       equilibrium momentum nodes")
    println("  --t-num <int>                       equilibrium theta nodes")
    println("  --max-iter <int>                    solver/mass iterations")
    println("  --qmax <float>                      q upper bound")
    println("  --q-nodes <int>                     q nodes")
    println("  --omega-min <float>                 phase-shift omega lower bound")
    println("  --omega-max <float>                 omega upper bound")
    println("  --omega-nodes <int>                 omega nodes")
    println("  --eta <float>                       phase-shift eta")
    println("  -h, --help                          show help")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :target_csv => DEFAULT_TARGET_CSV,
        :output_dir => joinpath(
            "data", "outputs", "results", "relaxtime", "meson_density", "crossover_validation",
            "friesen2019_fig5_top_kplus_over_piplus_mu055_stable",
        ),
        :curve_id => "kplus_over_piplus",
        :flavor_chemical_profile => "friesen2019_mu_s_0p55",
        :meson_chemical_profile => "friesen2019_kplus_over_piplus_neutral",
        :regime => :stable,
        :mu_min_MeV => 1.0,
        :mu_max_MeV => 1200.0,
        :T_min_MeV => 1.0,
        :T_max_MeV => 220.0,
        :n_mu => 48,
        :xi => 0.0,
        :method => :peak,
        :variable => :phi_u,
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
        elseif arg == "--curve-id"
            opts[:curve_id] = require_value()
        elseif arg == "--flavor-chemical-profile"
            opts[:flavor_chemical_profile] = require_value()
        elseif arg == "--meson-chemical-profile"
            opts[:meson_chemical_profile] = require_value()
        elseif arg == "--regime"
            opts[:regime] = Symbol(require_value())
        elseif arg == "--mu-min"
            opts[:mu_min_MeV] = parse(Float64, require_value())
        elseif arg == "--mu-max"
            opts[:mu_max_MeV] = parse(Float64, require_value())
        elseif arg == "--T-min"
            opts[:T_min_MeV] = parse(Float64, require_value())
        elseif arg == "--T-max"
            opts[:T_max_MeV] = parse(Float64, require_value())
        elseif arg == "--n-mu"
            opts[:n_mu] = parse(Int, require_value())
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--method"
            opts[:method] = Symbol(require_value())
        elseif arg == "--variable"
            opts[:variable] = Symbol(require_value())
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
        String(opts[:curve_id]),
        String(opts[:flavor_chemical_profile]),
        String(opts[:meson_chemical_profile]),
        Symbol(opts[:regime]),
        Float64(opts[:mu_min_MeV]),
        Float64(opts[:mu_max_MeV]),
        Float64(opts[:T_min_MeV]),
        Float64(opts[:T_max_MeV]),
        Int(opts[:n_mu]),
        Float64(opts[:xi]),
        Symbol(opts[:method]),
        Symbol(opts[:variable]),
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

@inline function _target_column(curve_id::String)
    curve_id == "kplus_over_piplus" && return :Kplus_over_piplus
    curve_id == "kminus_over_piminus" && return :Kminus_over_piminus
    error("unsupported curve_id: $curve_id")
end

function read_target_curve(path::String, curve_id::String)
    lines = readlines(path)
    isempty(lines) && error("target CSV is empty: $path")
    header = split(strip(lines[1]), ',')
    index = Dict(name => i for (i, name) in enumerate(header))
    haskey(index, "curve_id") || error("target CSV missing curve_id column")
    haskey(index, "T_over_muB") || error("target CSV missing T_over_muB column")
    target_col = String(_target_column(curve_id))
    haskey(index, target_col) || error("target CSV missing $(target_col) column")

    xs = Float64[]
    ys = Float64[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        strip(cols[index["curve_id"]]) == curve_id || continue
        y_raw = strip(cols[index[target_col]])
        isempty(y_raw) && continue
        push!(xs, parse(Float64, strip(cols[index["T_over_muB"]])))
        push!(ys, parse(Float64, y_raw))
    end
    isempty(xs) && error("no rows found for curve_id=$(curve_id)")
    return xs, ys
end

function read_scan_rows(path::String)
    lines = readlines(path)
    isempty(lines) && error("scan CSV is empty: $path")
    header = split(strip(lines[1]), ',')
    index = Dict(name => i for (i, name) in enumerate(header))
    required = ("T_over_muB", "kpi_ratio", "muq_MeV", "muB_MeV", "T_MeV", "equilibrium_converged", "message")
    for key in required
        haskey(index, key) || error("scan CSV missing column: $key")
    end

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        cols = split(s, ',')
        x = try parse(Float64, cols[index["T_over_muB"]]) catch; NaN end
        y = try parse(Float64, cols[index["kpi_ratio"]]) catch; NaN end
        isfinite(x) && isfinite(y) || continue
        push!(rows, (
            x=x,
            y=y,
            muq_MeV=parse(Float64, cols[index["muq_MeV"]]),
            muB_MeV=parse(Float64, cols[index["muB_MeV"]]),
            T_MeV=parse(Float64, cols[index["T_MeV"]]),
            converged=strip(cols[index["equilibrium_converged"]]) == "true",
            message=strip(cols[index["message"]]),
        ))
    end
    isempty(rows) && error("scan CSV contains no finite rows")
    sort!(rows; by=r -> r.x)
    return rows
end

@inline fmt(x::Real) = string(x)

function interpolate_scan(rows, x::Float64)
    n = length(rows)
    x < rows[1].x && return nothing
    x > rows[end].x && return nothing
    for i in 1:n
        if rows[i].x == x
            return (y=rows[i].y, left=rows[i], right=rows[i], weight=0.0, exact=true)
        end
    end
    for i in 1:(n - 1)
        left = rows[i]
        right = rows[i + 1]
        if left.x <= x <= right.x
            dx = right.x - left.x
            dx > 0.0 || return (y=left.y, left=left, right=right, weight=0.0, exact=false)
            w = (x - left.x) / dx
            y = left.y + w * (right.y - left.y)
            return (y=y, left=left, right=right, weight=w, exact=false)
        end
    end
    return nothing
end

function write_comparison_csv(path::String, xs::Vector{Float64}, ys::Vector{Float64}, rows)
    open(path, "w") do io
        println(io, "x_value,target_y,model_y,abs_diff,rel_diff,interpolation,exact_match,left_x,right_x,left_y,right_y,left_T_MeV,right_T_MeV,left_muB_MeV,right_muB_MeV")
        for (x, y) in zip(xs, ys)
            interp = interpolate_scan(rows, x)
            if interp === nothing
                println(io, join((fmt(x), fmt(y), "NaN", "NaN", "NaN", "out_of_range", "false", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN"), ','))
                continue
            end
            model_y = interp.y
            abs_diff = abs(model_y - y)
            rel_diff = iszero(y) ? NaN : abs_diff / abs(y)
            println(io, join((
                fmt(x),
                fmt(y),
                fmt(model_y),
                fmt(abs_diff),
                fmt(rel_diff),
                interp.exact ? "exact" : "linear",
                string(interp.exact),
                fmt(interp.left.x),
                fmt(interp.right.x),
                fmt(interp.left.y),
                fmt(interp.right.y),
                fmt(interp.left.T_MeV),
                fmt(interp.right.T_MeV),
                fmt(interp.left.muB_MeV),
                fmt(interp.right.muB_MeV),
            ), ','))
        end
    end
end

function summarize_comparison(xs::Vector{Float64}, ys::Vector{Float64}, rows)
    model_values = Float64[]
    abs_diffs = Float64[]
    rel_diffs = Float64[]
    in_range = 0
    for (x, y) in zip(xs, ys)
        interp = interpolate_scan(rows, x)
        interp === nothing && continue
        in_range += 1
        push!(model_values, interp.y)
        ad = abs(interp.y - y)
        push!(abs_diffs, ad)
        !iszero(y) && push!(rel_diffs, ad / abs(y))
    end
    return (
        target_points=length(xs),
        in_range_points=in_range,
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

function write_summary_md(path::String, opts::ValidationOptions, scan_csv::String, comparison_csv::String, summary)
    open(path, "w") do io
        println(io, "# crossover literature workflow reproduction 最小结果说明")
        println(io)
        println(io, "更新日期：$(Dates.today())")
        println(io)
        println(io, "## 1. 运行口径")
        println(io)
        println(io, "- target: `$(opts.target_csv)`")
        println(io, "- curve_id: `$(opts.curve_id)`")
        println(io, "- flavor chemical profile: `$(opts.flavor_chemical_profile)`")
        println(io, "- meson chemical profile: `$(opts.meson_chemical_profile)`")
        println(io, "- regime: `$(opts.regime)`")
        println(io, "- crossover method / variable: `$(opts.method)` / `$(opts.variable)`")
        println(io, "- mu_q window [MeV]: `$(opts.mu_min_MeV)` -> `$(opts.mu_max_MeV)`")
        println(io, "- T window [MeV]: `$(opts.T_min_MeV)` -> `$(opts.T_max_MeV)`")
        println(io, "- n_mu: `$(opts.n_mu)`")
        println(io, "- workflow scan: `$(scan_csv)`")
        println(io, "- comparison: `$(comparison_csv)`")
        println(io)
        println(io, "## 2. 当前结果定位")
        println(io)
        println(io, "这是 `Friesen 2019 Fig.5 top` 的 crossover-line workflow reproduction。")
        println(io)
        println(io, "- 已固定 charged channel、`mu_s = 0.55 mu_u` 与 `T/mu_B` 对齐输出；")
        println(io, "- 但当前仍不直接升格为 formal validation gate；")
        println(io, "- 若量级差异明显，优先解释为 path / meson-mass scale / workflow semantics 差异，而不是先把 target 当作回归真值。")
        println(io)
        println(io, "## 3. 数值摘要")
        println(io)
        println(io, "- target points: `$(summary.target_points)`")
        println(io, "- in-range points: `$(summary.in_range_points)`")
        println(io, "- target range: `$(summary.target_min)` -> `$(summary.target_max)`")
        println(io, "- model range: `$(summary.model_min)` -> `$(summary.model_max)`")
        println(io, "- mean abs diff: `$(summary.mean_abs_diff)`")
        println(io, "- max abs diff: `$(summary.max_abs_diff)`")
        println(io, "- mean rel diff: `$(summary.mean_rel_diff)`")
        println(io, "- max rel diff: `$(summary.max_rel_diff)`")
    end
end

function main(args::Vector{String})
    opts = parse_args(args)

    opts.curve_id in ("kplus_over_piplus", "kminus_over_piminus") || error("unsupported curve_id $(opts.curve_id)")

    if !opts.overwrite && isdir(opts.output_dir)
        error("output_dir already exists: $(opts.output_dir)")
    end
    mkpath(opts.output_dir)

    scan_csv = joinpath(opts.output_dir, "workflow_scan.csv")
    comparison_csv = joinpath(opts.output_dir, "comparison_vs_target.csv")
    summary_md = joinpath(opts.output_dir, "README.md")

    run_crossover_meson_density_scan(
        mu_min_MeV=opts.mu_min_MeV,
        mu_max_MeV=opts.mu_max_MeV,
        T_min_MeV=opts.T_min_MeV,
        T_max_MeV=opts.T_max_MeV,
        n_mu=opts.n_mu,
        xi=opts.xi,
        method=opts.method,
        variable=opts.variable,
        flavor_chemical_profile_name=opts.flavor_chemical_profile,
        meson_chemical_profile_name=opts.meson_chemical_profile,
        regime=opts.regime,
        output_path=scan_csv,
        overwrite=true,
        p_num=opts.p_num,
        t_num=opts.t_num,
        max_iter=opts.max_iter,
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
    )

    xs, ys = read_target_curve(opts.target_csv, opts.curve_id)
    rows = read_scan_rows(scan_csv)
    write_comparison_csv(comparison_csv, xs, ys, rows)
    summary = summarize_comparison(xs, ys, rows)
    write_summary_md(summary_md, opts, scan_csv, comparison_csv, summary)

    println("crossover literature validation completed")
    println("scan_csv=$(scan_csv)")
    println("comparison_csv=$(comparison_csv)")
    println("summary_md=$(summary_md)")
    println("curve_id=$(opts.curve_id)")
    println("regime=$(opts.regime)")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
