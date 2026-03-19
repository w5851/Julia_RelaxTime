#!/usr/bin/env julia

"""
Manual relaxation-time workflow entrypoint.

This script is intended for explicit, user-triggered generation of three output
families under data/outputs/figures/relaxtime:
1. cross_section
2. plan_a
3. plan_b

It reuses the existing scan scripts and the generic plot_scan_csv.py plotting
utility rather than introducing a second implementation path.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

struct Options
    sections::Set{Symbol}
    overwrite::Bool
    make_plots::Bool
    python_cmd::String
    integration_mode::String
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
    sigma_grid_n::Int
    compute_bulk::Bool
    xs_T_mev::Float64
    xs_muB_mev::Float64
    xs_xi::Float64
    plan_a_Tmin_mev::Float64
    plan_a_Tmax_mev::Float64
    plan_a_Tstep_mev::Float64
    plan_b_T_list_mev::Vector{Float64}
    xi_min::Float64
    xi_max::Float64
    xi_step::Float64
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_manual_relaxation_scan_workflow.jl [options]\n")
    println("Options:")
    println("  --sections <list>      all|cross_section|plan_a|plan_b (comma separated, default all)")
    println("  --overwrite            overwrite existing CSVs/figures")
    println("  --no-plots             only write CSVs")
    println("  --python <cmd>         python executable (default python)")
    println("  --mode <mode>          semi_infinite|finite_15|finite_lambda (default finite_15)")
    println("  --tau-p-nodes <int>    τ 平均散射率动量节点 (default 20)")
    println("  --tau-angle-nodes <int> τ 平均散射率 cosθ 节点 (default 4)")
    println("  --tau-phi-nodes <int>  τ 平均散射率 φ 节点 (default 8)")
    println("  --tau-n-sigma <int>    τ 使用的 σ(s) 积分点数 (default 6)")
    println("  --sigma-grid-n <int>   τ 使用的 σ(s) 预计算网格点数 (default 60)")
    println("  --no-bulk              disable bulk-viscosity calculation")
    println("  --xs-T <MeV>           cross section T (default 200)")
    println("  --xs-muB <MeV>         cross section mu_B (default 0)")
    println("  --xs-xi <value>        cross section xi (default 0)")
    println("  --plan-a-Tmin <MeV>    plan_a Tmin (default 120)")
    println("  --plan-a-Tmax <MeV>    plan_a Tmax (default 350)")
    println("  --plan-a-Tstep <MeV>   plan_a Tstep (default 10)")
    println("  --plan-b-T-list <csv>  fixed-T list for plan_b (default 150,190,200,250)")
    println("  --xi-min <value>       plan_b xi min (default -0.5)")
    println("  --xi-max <value>       plan_b xi max (default 0.5)")
    println("  --xi-step <value>      plan_b xi step (default 0.02)")
    println("  -h, --help             show help")
end

function parse_sections(raw::AbstractString)::Set{Symbol}
    items = Set{Symbol}()
    for part in split(raw, ',')
        token = Symbol(strip(part))
        token == :all && return Set([:cross_section, :plan_a, :plan_b])
        token in (:cross_section, :plan_a, :plan_b) || error("unknown section: $(part)")
        push!(items, token)
    end
    isempty(items) && error("--sections cannot be empty")
    return items
end

function parse_float_list(raw::AbstractString)::Vector{Float64}
    vals = Float64[]
    for part in split(raw, ',')
        s = strip(part)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && error("empty numeric list")
    return vals
end

function parse_args(args::Vector{String})::Options
    opts = Dict{Symbol,Any}(
        :sections => Set([:cross_section, :plan_a, :plan_b]),
        :overwrite => false,
        :make_plots => true,
        :python_cmd => "python",
        :integration_mode => "finite_15",
        :tau_p_nodes => 20,
        :tau_angle_nodes => 4,
        :tau_phi_nodes => 8,
        :tau_n_sigma_points => 6,
        :sigma_grid_n => 60,
        :compute_bulk => true,
        :xs_T_mev => 200.0,
        :xs_muB_mev => 0.0,
        :xs_xi => 0.0,
        :plan_a_Tmin_mev => 120.0,
        :plan_a_Tmax_mev => 350.0,
        :plan_a_Tstep_mev => 10.0,
        :plan_b_T_list_mev => Float64[150.0, 190.0, 200.0, 250.0],
        :xi_min => -0.5,
        :xi_max => 0.5,
        :xi_step => 0.02,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            value = args[i + 1]
            i += 1
            return value
        end

        if arg == "--sections"
            opts[:sections] = parse_sections(require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--no-plots"
            opts[:make_plots] = false
        elseif arg == "--python"
            opts[:python_cmd] = require_value()
        elseif arg == "--mode"
            mode = require_value()
            mode in ("semi_infinite", "finite_15", "finite_lambda") || error("unknown mode: $mode")
            opts[:integration_mode] = mode
        elseif arg == "--tau-p-nodes"
            opts[:tau_p_nodes] = parse(Int, require_value())
        elseif arg == "--tau-angle-nodes"
            opts[:tau_angle_nodes] = parse(Int, require_value())
        elseif arg == "--tau-phi-nodes"
            opts[:tau_phi_nodes] = parse(Int, require_value())
        elseif arg == "--tau-n-sigma"
            opts[:tau_n_sigma_points] = parse(Int, require_value())
        elseif arg == "--sigma-grid-n"
            opts[:sigma_grid_n] = parse(Int, require_value())
        elseif arg == "--no-bulk"
            opts[:compute_bulk] = false
        elseif arg == "--xs-T"
            opts[:xs_T_mev] = parse(Float64, require_value())
        elseif arg == "--xs-muB"
            opts[:xs_muB_mev] = parse(Float64, require_value())
        elseif arg == "--xs-xi"
            opts[:xs_xi] = parse(Float64, require_value())
        elseif arg == "--plan-a-Tmin"
            opts[:plan_a_Tmin_mev] = parse(Float64, require_value())
        elseif arg == "--plan-a-Tmax"
            opts[:plan_a_Tmax_mev] = parse(Float64, require_value())
        elseif arg == "--plan-a-Tstep"
            opts[:plan_a_Tstep_mev] = parse(Float64, require_value())
        elseif arg == "--plan-b-T-list"
            opts[:plan_b_T_list_mev] = parse_float_list(require_value())
        elseif arg == "--xi-min"
            opts[:xi_min] = parse(Float64, require_value())
        elseif arg == "--xi-max"
            opts[:xi_max] = parse(Float64, require_value())
        elseif arg == "--xi-step"
            opts[:xi_step] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    opts[:xi_step] > 0 || error("xi-step must be positive")
    opts[:plan_a_Tstep_mev] > 0 || error("plan-a-Tstep must be positive")

    return Options(
        Set{Symbol}(opts[:sections]),
        Bool(opts[:overwrite]),
        Bool(opts[:make_plots]),
        String(opts[:python_cmd]),
        String(opts[:integration_mode]),
        Int(opts[:tau_p_nodes]),
        Int(opts[:tau_angle_nodes]),
        Int(opts[:tau_phi_nodes]),
        Int(opts[:tau_n_sigma_points]),
        Int(opts[:sigma_grid_n]),
        Bool(opts[:compute_bulk]),
        Float64(opts[:xs_T_mev]),
        Float64(opts[:xs_muB_mev]),
        Float64(opts[:xs_xi]),
        Float64(opts[:plan_a_Tmin_mev]),
        Float64(opts[:plan_a_Tmax_mev]),
        Float64(opts[:plan_a_Tstep_mev]),
        Float64.(opts[:plan_b_T_list_mev]),
        Float64(opts[:xi_min]),
        Float64(opts[:xi_max]),
        Float64(opts[:xi_step]),
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function maybe_rm(path::AbstractString; overwrite::Bool)
    overwrite && isfile(path) && rm(path; force=true)
end

function maybe_rmdir(path::AbstractString; overwrite::Bool)
    overwrite && isdir(path) && rm(path; recursive=true, force=true)
end

function run_cmd(cmd::Cmd)
    println("\n> ", cmd)
    run(cmd)
end

function julia_cmd(argv::Vector{String})::Cmd
    jc = Base.julia_cmd()
    cmd = Cmd(vcat(jc.exec, argv))
    return jc.env === nothing ? cmd : setenv(cmd, jc.env)
end

function simple_cmd(exe::AbstractString, argv::Vector{String})::Cmd
    return Cmd(vcat([String(exe)], String.(argv)))
end

function xi_list_string(xmin::Float64, xmax::Float64, xstep::Float64)::String
    vals = collect(range(xmin; stop=xmax, step=xstep))
    fmt(x) = abs(x) < 1e-12 ? "0.0" : string(x)
    return join(fmt.(vals), ",")
end

function split_cross_section_csv_by_process_groups(xs_csv::AbstractString; overwrite::Bool)::Vector{Tuple{String,String}}
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
    idx_process === nothing && error("missing process column: $xs_csv")
    idx_sqrt === nothing && error("missing sqrt_s_MeV column: $xs_csv")

    proc_lines = Dict{String,Vector{String}}()
    proc_min = Dict{String,Float64}()
    for line in data_lines
        parts = split(line, ',')
        p = parts[idx_process]
        value = parse(Float64, parts[idx_sqrt])
        push!(get!(proc_lines, p, String[]), line)
        proc_min[p] = min(get(proc_min, p, value), value)
    end

    procs = sort(collect(keys(proc_lines)); by=p -> proc_min[p])
    isempty(procs) && error("no process data in $xs_csv")

    groups = if length(procs) <= 3
        [[p] for p in procs]
    else
        mins = [proc_min[p] for p in procs]
        gaps = [(i, mins[i + 1] - mins[i]) for i in 1:(length(mins) - 1)]
        cuts = sort(first.(sort(gaps; by=last, rev=true)[1:min(2, length(gaps))]))
        if isempty(cuts)
            [procs]
        elseif length(cuts) == 1
            c1 = cuts[1]
            [procs[1:c1], procs[(c1 + 1):end]]
        else
            c1, c2 = cuts
            [procs[1:c1], procs[(c1 + 1):c2], procs[(c2 + 1):end]]
        end
    end

    outputs = Tuple{String,String}[]
    for (index, plist) in enumerate(groups)
        tag = "group$(index)"
        group_csv = replace(xs_csv, ".csv" => "_$(tag).csv")
        maybe_rm(group_csv; overwrite=overwrite)
        open(group_csv, "w") do io
            for line in meta_lines
                println(io, line)
            end
            println(io, header_line)
            for process in plist
                for line in proc_lines[process]
                    println(io, line)
                end
            end
        end
        push!(outputs, (tag, group_csv))
    end

    return outputs
end

function merge_csvs(inputs::Vector{String}, output::String; overwrite::Bool)
    isempty(inputs) && error("no input CSVs to merge")
    ensure_parent_dir(output)
    maybe_rm(output; overwrite=overwrite)

    meta_lines = String[]
    header_line = nothing
    rows = String[]

    for (index, path) in enumerate(inputs)
        open(path, "r") do io
            for line in eachline(io)
                s = chomp(line)
                isempty(strip(s)) && continue
                if startswith(strip(s), "#")
                    index == 1 && push!(meta_lines, s)
                    continue
                end
                if header_line === nothing
                    header_line = s
                    continue
                end
                if index == 1 || s != header_line
                    push!(rows, s)
                end
            end
        end
    end

    open(output, "w") do io
        for line in meta_lines
            println(io, line)
        end
        println(io, header_line)
        for row in rows
            println(io, row)
        end
    end
end

function run_cross_section(opts::Options)
    xs_csv = joinpath("data", "outputs", "results", "relaxtime", "cross_section",
        "xs_vs_s_by_process_T$(Int(round(opts.xs_T_mev)))_muB$(Int(round(opts.xs_muB_mev)))_xi$(opts.xs_xi).csv")
    xs_fig_dir = joinpath("data", "outputs", "figures", "relaxtime", "cross_section",
        "T$(Int(round(opts.xs_T_mev)))_muB$(Int(round(opts.xs_muB_mev)))_xi$(opts.xs_xi)")

    ensure_parent_dir(xs_csv)
    if opts.overwrite
        maybe_rm(xs_csv; overwrite=true)
        for index in 1:3
            maybe_rm(replace(xs_csv, ".csv" => "_group$(index).csv"); overwrite=true)
        end
        maybe_rmdir(xs_fig_dir; overwrite=true)
    end

    args = String[
        "--project=.",
        "scripts/relaxtime/scan_cross_section_vs_s_by_process.jl",
        "--T-MeV", string(opts.xs_T_mev),
        "--muB-MeV", string(opts.xs_muB_mev),
        "--xi", string(opts.xs_xi),
        "--out", xs_csv,
    ]
    opts.overwrite && push!(args, "--overwrite")
    run_cmd(julia_cmd(args))

    if opts.make_plots
        group_csvs = split_cross_section_csv_by_process_groups(xs_csv; overwrite=opts.overwrite)
        for (tag, group_csv) in group_csvs
            out_dir = joinpath(xs_fig_dir, "sqrt_s_$(tag)")
            ensure_parent_dir(joinpath(out_dir, "dummy.txt"))
            run_cmd(simple_cmd(opts.python_cmd, [
                "scripts/plot_scan_csv.py",
                "--mode", "lines",
                "--csv", group_csv,
                "--x", "sqrt_s_MeV",
                "--ys", "sigma",
                "--group", "process",
                "--out-dir", out_dir,
            ]))

            yzoom_dir = joinpath(out_dir, "ylim_0_10")
            ensure_parent_dir(joinpath(yzoom_dir, "dummy.txt"))
            run_cmd(simple_cmd(opts.python_cmd, [
                "scripts/plot_scan_csv.py",
                "--mode", "lines",
                "--csv", group_csv,
                "--x", "sqrt_s_MeV",
                "--ys", "sigma",
                "--group", "process",
                "--ylim", "0", "10",
                "--out-dir", yzoom_dir,
            ]))
        end

        ss_only_dir = joinpath(xs_fig_dir, "process_ssbar_to_uubar")
        ensure_parent_dir(joinpath(ss_only_dir, "dummy.txt"))
        run_cmd(simple_cmd(opts.python_cmd, [
            "scripts/plot_scan_csv.py",
            "--mode", "lines",
            "--csv", xs_csv,
            "--x", "sqrt_s_MeV",
            "--ys", "sigma",
            "--where", "process=ssbar_to_uubar",
            "--out-dir", ss_only_dir,
        ]))
    end
end

function run_plan_a(opts::Options)
    out_csv = joinpath("data", "outputs", "results", "relaxtime", "plan_a", "gap_transport_vs_T_muB0_xi0.csv")
    out_fig = joinpath("data", "outputs", "figures", "relaxtime", "plan_a")

    ensure_parent_dir(out_csv)
    if opts.overwrite
        maybe_rmdir(dirname(out_csv); overwrite=true)
        maybe_rmdir(out_fig; overwrite=true)
    end

    args = String[
        "--project=.",
        "scripts/relaxtime/run_gap_transport_scan.jl",
        "--tmin", string(opts.plan_a_Tmin_mev),
        "--tmax", string(opts.plan_a_Tmax_mev),
        "--tstep", string(opts.plan_a_Tstep_mev),
        "--mubmin", "0",
        "--mubmax", "0",
        "--mubstep", "1",
        "--xi-list", "0.0",
        "--mode", opts.integration_mode,
        "--tau-p-nodes", string(opts.tau_p_nodes),
        "--tau-angle-nodes", string(opts.tau_angle_nodes),
        "--tau-phi-nodes", string(opts.tau_phi_nodes),
        "--tau-n-sigma", string(opts.tau_n_sigma_points),
        "--sigma-grid-n", string(opts.sigma_grid_n),
        "--output", out_csv,
    ]
    opts.compute_bulk && push!(args, "--compute-bulk")
    opts.overwrite && push!(args, "--overwrite")
    run_cmd(julia_cmd(args))

    if opts.make_plots
        ensure_parent_dir(joinpath(out_fig, "dummy.txt"))
        run_cmd(simple_cmd(opts.python_cmd, [
            "scripts/plot_scan_csv.py",
            "--mode", "lines",
            "--csv", out_csv,
            "--x", "T_MeV",
            "--ys", "tau_u,tau_d,tau_s,tau_ubar,tau_dbar,tau_sbar",
            "--multi-y",
            "--yscale", "log",
            "--out-dir", out_fig,
        ]))
        run_cmd(simple_cmd(opts.python_cmd, [
            "scripts/plot_scan_csv.py",
            "--mode", "lines",
            "--csv", out_csv,
            "--x", "T_MeV",
            "--ys", "sigma_over_T",
            "--out-dir", out_fig,
        ]))
        run_cmd(simple_cmd(opts.python_cmd, [
            "scripts/plot_scan_csv.py",
            "--mode", "lines",
            "--csv", out_csv,
            "--x", "T_MeV",
            "--ys", "eta_over_s,zeta_over_s",
            "--multi-y",
            "--yscale", "log",
            "--out-dir", out_fig,
        ]))
        run_cmd(simple_cmd(opts.python_cmd, [
            "scripts/plot_scan_csv.py",
            "--mode", "lines",
            "--csv", out_csv,
            "--x", "T_MeV",
            "--ys", "sigma_over_T_over_eta_over_s,zeta_over_s_over_eta_over_s",
            "--out-dir", out_fig,
        ]))
    end
end

function run_plan_b(opts::Options)
    xi_list = xi_list_string(opts.xi_min, opts.xi_max, opts.xi_step)
    result_dir = joinpath("data", "outputs", "results", "relaxtime", "plan_b")
    figure_dir = joinpath("data", "outputs", "figures", "relaxtime", "plan_b")
    merged_csv = joinpath(result_dir, "plan_b_merged.csv")
    csv_paths = String[]

    if opts.overwrite
        maybe_rmdir(result_dir; overwrite=true)
        maybe_rmdir(figure_dir; overwrite=true)
    end

    for T_mev in opts.plan_b_T_list_mev
        csv_path = joinpath(result_dir, "transport_vs_xi_T$(Int(round(T_mev)))_muB0.csv")
        fig_path = joinpath(figure_dir, "T$(Int(round(T_mev)))")
        ensure_parent_dir(csv_path)

        args = String[
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
            "--tau-p-nodes", string(opts.tau_p_nodes),
            "--tau-angle-nodes", string(opts.tau_angle_nodes),
            "--tau-phi-nodes", string(opts.tau_phi_nodes),
            "--tau-n-sigma", string(opts.tau_n_sigma_points),
            "--sigma-grid-n", string(opts.sigma_grid_n),
            "--output", csv_path,
        ]
        opts.compute_bulk && push!(args, "--compute-bulk")
        opts.overwrite && push!(args, "--overwrite")
        run_cmd(julia_cmd(args))
        push!(csv_paths, csv_path)

        if opts.make_plots
            ensure_parent_dir(joinpath(fig_path, "dummy.txt"))
            run_cmd(simple_cmd(opts.python_cmd, [
                "scripts/plot_scan_csv.py",
                "--mode", "lines",
                "--csv", csv_path,
                "--x", "xi",
                "--ys", "tau_u,tau_ubar,tau_s,tau_sbar",
                "--multi-y",
                "--out-dir", fig_path,
            ]))
            for y in ("sigma_over_T", "eta_over_s", "zeta_over_s", "sigma_over_T_over_eta_over_s", "zeta_over_s_over_eta_over_s")
                run_cmd(simple_cmd(opts.python_cmd, [
                    "scripts/plot_scan_csv.py",
                    "--mode", "lines",
                    "--csv", csv_path,
                    "--x", "xi",
                    "--ys", y,
                    "--out-dir", fig_path,
                ]))
            end
        end
    end

    merge_csvs(csv_paths, merged_csv; overwrite=opts.overwrite)

    if opts.make_plots
        combined_dir = joinpath(figure_dir, "combined")
        ensure_parent_dir(joinpath(combined_dir, "dummy.txt"))
        for y in (
            "tau_u", "tau_ubar", "tau_s", "tau_sbar",
            "sigma_over_T", "eta_over_s", "zeta_over_s",
            "sigma_over_T_over_eta_over_s", "zeta_over_s_over_eta_over_s",
        )
            plot_args = String[
                "scripts/plot_scan_csv.py",
                "--mode", "lines",
                "--csv", merged_csv,
                "--x", "xi",
                "--ys", y,
                "--group", "T_MeV",
                "--out-dir", combined_dir,
            ]
            if y == "eta_over_s"
                append!(plot_args, ["--yscale", "log"])
            end
            run_cmd(simple_cmd(opts.python_cmd, plot_args))
        end
    end
end

function main()
    opts = parse_args(copy(ARGS))
    cd(PROJECT_ROOT) do
        :cross_section in opts.sections && run_cross_section(opts)
        :plan_a in opts.sections && run_plan_a(opts)
        :plan_b in opts.sections && run_plan_b(opts)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
