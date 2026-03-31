#!/usr/bin/env julia

"""
Status: deprecate-candidate
Purpose: legacy one-off output bundle script for historical 2026-02-05 task context.
Replacement: prefer `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl` or `scripts/relaxtime/run_relaxtime_orchestrator.jl`.

为 docs/dev/active/2026_02_05_调用扫描脚本and画图脚本输出数据.md 的需求生成数据与图片。

覆盖的输出：
1) 散射截面：固定 (T, μ_B, ξ=0)，按过程输出 σ(√s)
2) 固定 μ_B=0：手征相变温度上下各一点（默认 190/210 MeV），扫描 ξ
3) 固定 μ_B=0、ξ=0：扫描 T

默认输出到：
- data/outputs/results/... （按子文件夹归档）
- data/outputs/figures/...（按子文件夹归档）

说明：
- 体粘滞 ζ 需要额外计算，默认开启（可能很慢）。可用 --no-bulk 关闭，但会缺少 ζ 相关列。
- 绘图通过 python scripts/plot_scan_csv.py 完成。

示例：
  julia --project=. scripts/relaxtime/run_outputs_2026_02_05.jl --overwrite
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

struct Options
    overwrite::Bool
    make_plots::Bool
    python_cmd::String
    integration_mode::String
    compute_bulk::Bool

    # cross section
    xs_T_mev::Float64
    xs_muB_mev::Float64
    xs_xi::Float64

    # xi scan at two T points
    xi_min::Float64
    xi_max::Float64
    xi_step::Float64
    T_below_mev::Float64
    T_above_mev::Float64

    # T scan at xi=0
    T_scan_min_mev::Float64
    T_scan_max_mev::Float64
    T_scan_step_mev::Float64
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_outputs_2026_02_05.jl [options]\n")
    println("Options:")
    println("  --overwrite          覆盖已有 CSV/图片")
    println("  --no-plots           只产出 CSV，不作图")
    println("  --python <cmd>       Python 命令（default: python）")
    println("  --mode <mode>        run_gap_transport_scan.jl 积分模式: semi_infinite|finite_15|finite_lambda (default finite_15)")
    println("  --no-bulk            关闭体粘滞计算（更快，但缺 ζ 相关输出）")
    println("\nCross section (σ vs √s):")
    println("  --xs-T <MeV>         default 200")
    println("  --xs-muB <MeV>       default 0")
    println("  --xs-xi <value>      default 0")
    println("\nXi scan points:")
    println("  --xi-min/--xi-max/--xi-step  default -0.6..0.6 step 0.1")
    println("  --T-below <MeV>      default 190")
    println("  --T-above <MeV>      default 210")
    println("\nT scan (μB=0, ξ=0):")
    println("  --Tmin/--Tmax/--Tstep default 120..350 step 10")
    println("  -h, --help           显示帮助")
end

function parse_args(args::Vector{String})::Options
    opts = Dict{Symbol,Any}(
        :overwrite => false,
        :make_plots => true,
        :python_cmd => "python",
        :integration_mode => "finite_15",
        :compute_bulk => true,
        :xs_T_mev => 200.0,
        :xs_muB_mev => 0.0,
        :xs_xi => 0.0,
        :xi_min => -0.6,
        :xi_max => 0.6,
        :xi_step => 0.1,
        :T_below_mev => 190.0,
        :T_above_mev => 210.0,
        :T_scan_min_mev => 120.0,
        :T_scan_max_mev => 350.0,
        :T_scan_step_mev => 10.0,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            v = args[i + 1]
            i += 1
            return v
        end

        if arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--no-plots"
            opts[:make_plots] = false
        elseif arg == "--python"
            opts[:python_cmd] = require_value()
        elseif arg == "--mode"
            opts[:integration_mode] = require_value()
        elseif arg == "--no-bulk"
            opts[:compute_bulk] = false
        elseif arg == "--xs-T"
            opts[:xs_T_mev] = parse(Float64, require_value())
        elseif arg == "--xs-muB"
            opts[:xs_muB_mev] = parse(Float64, require_value())
        elseif arg == "--xs-xi"
            opts[:xs_xi] = parse(Float64, require_value())
        elseif arg == "--xi-min"
            opts[:xi_min] = parse(Float64, require_value())
        elseif arg == "--xi-max"
            opts[:xi_max] = parse(Float64, require_value())
        elseif arg == "--xi-step"
            opts[:xi_step] = parse(Float64, require_value())
        elseif arg == "--T-below"
            opts[:T_below_mev] = parse(Float64, require_value())
        elseif arg == "--T-above"
            opts[:T_above_mev] = parse(Float64, require_value())
        elseif arg == "--Tmin"
            opts[:T_scan_min_mev] = parse(Float64, require_value())
        elseif arg == "--Tmax"
            opts[:T_scan_max_mev] = parse(Float64, require_value())
        elseif arg == "--Tstep"
            opts[:T_scan_step_mev] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage(); exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return Options(
        Bool(opts[:overwrite]),
        Bool(opts[:make_plots]),
        String(opts[:python_cmd]),
        String(opts[:integration_mode]),
        Bool(opts[:compute_bulk]),
        Float64(opts[:xs_T_mev]),
        Float64(opts[:xs_muB_mev]),
        Float64(opts[:xs_xi]),
        Float64(opts[:xi_min]),
        Float64(opts[:xi_max]),
        Float64(opts[:xi_step]),
        Float64(opts[:T_below_mev]),
        Float64(opts[:T_above_mev]),
        Float64(opts[:T_scan_min_mev]),
        Float64(opts[:T_scan_max_mev]),
        Float64(opts[:T_scan_step_mev]),
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function xi_list_string(xmin::Float64, xmax::Float64, xstep::Float64)::String
    xstep > 0 || error("xi-step must be positive")
    xmax >= xmin || error("xi-max must be >= xi-min")
    vals = collect(range(xmin; stop=xmax, step=xstep))
    # 保持稳定字符串格式，避免 -0.0
    fmt(x) = (abs(x) < 1e-12) ? "0.0" : string(x)
    return join(fmt.(vals), ",")
end

function run_cmd(cmd::Cmd)
    println("\n> ", cmd)
    run(cmd)
end

function julia_cmd(argv::Vector{String})::Cmd
    jc = Base.julia_cmd()
    cmd = Cmd(vcat(jc.exec, argv))
    # 兼容不同 Julia 版本：Cmd(Vector) 不支持 env/dir kwargs。
    # 工作目录由外层 cd(PROJECT_ROOT) 负责；环境变量如有需要用 setenv 注入。
    return jc.env === nothing ? cmd : setenv(cmd, jc.env)
end

function simple_cmd(exe::AbstractString, argv::Vector{String})::Cmd
    return Cmd(vcat([String(exe)], String.(argv)))
end

function maybe_rm(path::AbstractString; overwrite::Bool)
    if overwrite && isfile(path)
        rm(path; force=true)
    end
end

function maybe_rmdir(path::AbstractString; overwrite::Bool)
    if overwrite && isdir(path)
        rm(path; recursive=true, force=true)
    end
end

function read_csv_column_floats(path::AbstractString, col::AbstractString)::Vector{Float64}
    isfile(path) || error("CSV not found: $path")
    header_cols = nothing
    col_idx = nothing
    vals = Float64[]
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            (isempty(s) || startswith(s, "#")) && continue
            if header_cols === nothing
                header_cols = split(s, ',')
                col_idx = findfirst(==(col), header_cols)
                col_idx === nothing && error("CSV missing column '$col': $path")
                continue
            end
            parts = split(s, ',')
            (col_idx > length(parts)) && continue
            v = try
                parse(Float64, parts[col_idx])
            catch
                NaN
            end
            isfinite(v) && push!(vals, v)
        end
    end
    return vals
end

function three_ranges(vals::Vector{Float64})::Vector{Tuple{Float64,Float64}}
    isempty(vals) && error("cannot build ranges from empty values")
    v = sort(vals)
    n = length(v)
    a = v[1]
    b = v[clamp(Int(ceil(n / 3)), 1, n)]
    c = v[clamp(Int(ceil(2n / 3)), 1, n)]
    d = v[end]
    # De-duplicate boundaries if data is concentrated.
    r = Tuple{Float64,Float64}[]
    push!(r, (a, b))
    if c > b
        push!(r, (b, c))
    end
    if d > c
        push!(r, (c, d))
    end
    return r
end

function split_cross_section_csv_by_process_groups(xs_csv::AbstractString; overwrite::Bool)::Vector{Tuple{String,String}}
    # Returns [(group_tag, group_csv_path), ...]
    # Grouping rule: group processes by their effective sqrt_s ranges (based on per-process min sqrt_s),
    # so a single process is never split across multiple plots.

    isfile(xs_csv) || error("CSV not found: $xs_csv")
    meta_lines = String[]
    header_line = nothing
    data_lines = String[]
    open(xs_csv, "r") do io
        for line in eachline(io)
            s = chomp(line)
            if startswith(strip(s), "#")
                push!(meta_lines, s)
                continue
            end
            if header_line === nothing
                header_line = s
                continue
            end
            isempty(strip(s)) && continue
            push!(data_lines, s)
        end
    end
    header_line === nothing && error("CSV missing header: $xs_csv")
    cols = split(header_line, ',')
    idx_process = findfirst(==("process"), cols)
    idx_sqrt = findfirst(==("sqrt_s_MeV"), cols)
    idx_process === nothing && error("CSV missing column 'process': $xs_csv")
    idx_sqrt === nothing && error("CSV missing column 'sqrt_s_MeV': $xs_csv")

    proc_min = Dict{String,Float64}()
    proc_max = Dict{String,Float64}()
    proc_lines = Dict{String,Vector{String}}()
    for ln in data_lines
        parts = split(ln, ',')
        (idx_process > length(parts) || idx_sqrt > length(parts)) && continue
        p = parts[idx_process]
        v = try
            parse(Float64, parts[idx_sqrt])
        catch
            NaN
        end
        get!(proc_lines, p, String[])
        push!(proc_lines[p], ln)
        if isfinite(v)
            if !haskey(proc_min, p)
                proc_min[p] = v
                proc_max[p] = v
            else
                proc_min[p] = min(proc_min[p], v)
                proc_max[p] = max(proc_max[p], v)
            end
        end
    end

    procs = sort(collect(keys(proc_lines)); by=p->get(proc_min, p, Inf))
    isempty(procs) && error("No process data in CSV: $xs_csv")

    # Partition into up to 3 groups by locating the two largest gaps in sorted min(sqrt_s).
    # This keeps low-threshold processes together and avoids splitting a cluster by mere count.
    mins = [get(proc_min, p, Inf) for p in procs]
    n = length(procs)
    if n <= 3
        groups = [[p] for p in procs]
    else
        gaps = [(i, mins[i + 1] - mins[i]) for i in 1:(n - 1)]
        # pick up to 2 largest positive gaps
        gaps2 = sort(gaps; by=x->x[2], rev=true)
        cut_idxs = Int[]
        for (i, g) in gaps2
            g > 0 || continue
            push!(cut_idxs, i)
            length(cut_idxs) == 2 && break
        end
        cut_idxs = sort(unique(cut_idxs))
        if length(cut_idxs) == 0
            # fallback: all in one group
            groups = [procs]
        elseif length(cut_idxs) == 1
            c1 = cut_idxs[1]
            groups = [procs[1:c1], procs[(c1 + 1):end]]
        else
            c1, c2 = cut_idxs[1], cut_idxs[2]
            groups = [procs[1:c1], procs[(c1 + 1):c2], procs[(c2 + 1):end]]
        end
    end

    # Ensure at most 3 groups
    if length(groups) > 3
        groups = groups[1:3]
    end

    out = Tuple{String,String}[]
    for (gi, plist) in enumerate(groups)
        tag = "group$(gi)"
        group_csv = replace(xs_csv, ".csv" => "_$(tag).csv")
        if overwrite
            isfile(group_csv) && rm(group_csv; force=true)
        end
        open(group_csv, "w") do io
            for ml in meta_lines
                println(io, ml)
            end
            println(io, header_line)
            for p in plist
                for ln in proc_lines[p]
                    println(io, ln)
                end
            end
        end
        push!(out, (tag, group_csv))
    end
    return out
end

function main()
    opts = parse_args(copy(ARGS))

    cd(PROJECT_ROOT) do
        # ---- 1) cross section by process ----
        xs_csv = joinpath("data", "outputs", "results", "relaxtime", "cross_section",
            "xs_vs_s_by_process_T$(Int(round(opts.xs_T_mev)))_muB$(Int(round(opts.xs_muB_mev)))_xi$(opts.xs_xi).csv")
        xs_fig_dir = joinpath("data", "outputs", "figures", "relaxtime", "cross_section",
            "T$(Int(round(opts.xs_T_mev)))_muB$(Int(round(opts.xs_muB_mev)))_xi$(opts.xs_xi)")
        ensure_parent_dir(xs_csv)
        if opts.overwrite
            maybe_rm(xs_csv; overwrite=true)
            for gi in 1:3
                maybe_rm(replace(xs_csv, ".csv" => "_group$(gi).csv"); overwrite=true)
            end
            maybe_rmdir(xs_fig_dir; overwrite=true)
        end

        xs_args = String[
            "--project=.",
            "scripts/relaxtime/scan_cross_section_vs_s_by_process.jl",
            "--T-MeV", string(opts.xs_T_mev),
            "--muB-MeV", string(opts.xs_muB_mev),
            "--xi", string(opts.xs_xi),
            "--out", xs_csv,
        ]
        opts.overwrite && push!(xs_args, "--overwrite")
        run_cmd(julia_cmd(xs_args))

        if opts.make_plots
            ensure_parent_dir(joinpath(xs_fig_dir, "dummy.txt"))
            py = opts.python_cmd

            group_csvs = split_cross_section_csv_by_process_groups(xs_csv; overwrite=opts.overwrite)
            for (tag, group_csv) in group_csvs
                seg_dir = joinpath(xs_fig_dir, "sqrt_s_$(tag)")
                ensure_parent_dir(joinpath(seg_dir, "dummy.txt"))
                run_cmd(simple_cmd(py, [
                    "scripts/plot_scan_csv.py",
                    "--mode", "lines",
                    "--csv", group_csv,
                    "--x", "sqrt_s_MeV",
                    "--ys", "sigma",
                    "--group", "process",
                    "--out-dir", seg_dir,
                ]))

                # 额外一份：限制 y 轴，避免极端过程拉伸其它曲线
                seg_dir2 = joinpath(seg_dir, "ylim_0_10")
                ensure_parent_dir(joinpath(seg_dir2, "dummy.txt"))
                run_cmd(simple_cmd(py, [
                    "scripts/plot_scan_csv.py",
                    "--mode", "lines",
                    "--csv", group_csv,
                    "--x", "sqrt_s_MeV",
                    "--ys", "sigma",
                    "--group", "process",
                    "--ylim", "0", "10",
                    "--out-dir", seg_dir2,
                ]))
            end

            # 单独输出 ssbar_to_uubar，避免其影响整体 y 轴可读性
            ss_only_dir = joinpath(xs_fig_dir, "process_ssbar_to_uubar")
            ensure_parent_dir(joinpath(ss_only_dir, "dummy.txt"))
            run_cmd(simple_cmd(py, [
                "scripts/plot_scan_csv.py",
                "--mode", "lines",
                "--csv", xs_csv,
                "--x", "sqrt_s_MeV",
                "--ys", "sigma",
                "--where", "process=ssbar_to_uubar",
                "--out-dir", ss_only_dir,
            ]))
        end

        # ---- 2) muB=0: xi scan at two T points ----
        xi_list = xi_list_string(opts.xi_min, opts.xi_max, opts.xi_step)
        for (tag, T_mev) in (("below", opts.T_below_mev), ("above", opts.T_above_mev))
            out_csv = joinpath("data", "outputs", "results", "relaxtime", "transport_vs_xi", "muB0",
                "T$(Int(round(T_mev)))_$(tag)", "gap_transport_vs_xi.csv")
            out_fig = joinpath("data", "outputs", "figures", "relaxtime", "transport_vs_xi", "muB0",
                "T$(Int(round(T_mev)))_$(tag)")
            ensure_parent_dir(out_csv)
            if opts.overwrite
                maybe_rmdir(dirname(out_csv); overwrite=true)
                maybe_rmdir(out_fig; overwrite=true)
            end

            scan_args = String[
                "--project=.",
                "scripts/relaxtime/run_gap_transport_scan.jl",
                "--tmin", string(T_mev),
                "--tmax", string(T_mev),
                "--tstep", "1",
                "--mubmin", "0",
                "--mubmax", "0",
                "--mubstep", "1",
                "--xi-list", xi_list,
                "--mode", opts.integration_mode,
                "--output", out_csv,
            ]
            opts.compute_bulk && push!(scan_args, "--compute-bulk")
            opts.overwrite && push!(scan_args, "--overwrite")
            run_cmd(julia_cmd(scan_args))

            if opts.make_plots
                ensure_parent_dir(joinpath(out_fig, "dummy.txt"))
                py = opts.python_cmd
                run_cmd(simple_cmd(py, ["scripts/plot_scan_csv.py", "--mode", "lines", "--csv", out_csv, "--x", "xi", "--ys", "tau_u,tau_d,tau_s,tau_ubar,tau_dbar,tau_sbar", "--multi-y", "--out-dir", out_fig]))
                run_cmd(simple_cmd(py, ["scripts/plot_scan_csv.py", "--mode", "lines", "--csv", out_csv, "--x", "xi", "--ys", "sigma_over_T,eta_over_s,zeta_over_s", "--out-dir", out_fig]))
                run_cmd(simple_cmd(py, ["scripts/plot_scan_csv.py", "--mode", "lines", "--csv", out_csv, "--x", "xi", "--ys", "sigma_over_T_over_eta_over_s,zeta_over_s_over_eta_over_s", "--out-dir", out_fig]))
            end
        end

        # ---- 3) muB=0, xi=0: scan vs T ----
        out_csv_T = joinpath("data", "outputs", "results", "relaxtime", "transport_vs_T", "muB0", "xi0", "gap_transport_vs_T.csv")
        out_fig_T = joinpath("data", "outputs", "figures", "relaxtime", "transport_vs_T", "muB0", "xi0")
        ensure_parent_dir(out_csv_T)
        if opts.overwrite
            maybe_rmdir(dirname(out_csv_T); overwrite=true)
            maybe_rmdir(out_fig_T; overwrite=true)
        end

        scanT_args = String[
            "--project=.",
            "scripts/relaxtime/run_gap_transport_scan.jl",
            "--tmin", string(opts.T_scan_min_mev),
            "--tmax", string(opts.T_scan_max_mev),
            "--tstep", string(opts.T_scan_step_mev),
            "--mubmin", "0",
            "--mubmax", "0",
            "--mubstep", "1",
            "--xi-list", "0.0",
            "--mode", opts.integration_mode,
            "--output", out_csv_T,
        ]
        opts.compute_bulk && push!(scanT_args, "--compute-bulk")
        opts.overwrite && push!(scanT_args, "--overwrite")
        run_cmd(julia_cmd(scanT_args))

        if opts.make_plots
            ensure_parent_dir(joinpath(out_fig_T, "dummy.txt"))
            py = opts.python_cmd
            run_cmd(simple_cmd(py, ["scripts/plot_scan_csv.py", "--mode", "lines", "--csv", out_csv_T, "--x", "T_MeV", "--ys", "tau_u,tau_d,tau_s,tau_ubar,tau_dbar,tau_sbar", "--multi-y", "--yscale", "log", "--out-dir", out_fig_T]))
            run_cmd(simple_cmd(py, ["scripts/plot_scan_csv.py", "--mode", "lines", "--csv", out_csv_T, "--x", "T_MeV", "--ys", "sigma_over_T", "--out-dir", out_fig_T]))
            run_cmd(simple_cmd(py, ["scripts/plot_scan_csv.py", "--mode", "lines", "--csv", out_csv_T, "--x", "T_MeV", "--ys", "eta_over_s,zeta_over_s", "--multi-y", "--yscale", "log", "--out-dir", out_fig_T]))
            run_cmd(simple_cmd(py, ["scripts/plot_scan_csv.py", "--mode", "lines", "--csv", out_csv_T, "--x", "T_MeV", "--ys", "sigma_over_T_over_eta_over_s,zeta_over_s_over_eta_over_s", "--out-dir", out_fig_T]))
        end

        println("\nAll done. Results under data/outputs/results and figures under data/outputs/figures")
    end
end

main()
