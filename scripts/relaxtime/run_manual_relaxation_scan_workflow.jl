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

@inline function _norm_slash(path::AbstractString)
    return replace(String(path), '\\' => '/')
end

function _repo_relpath(path::AbstractString)
    abs_path = normpath(abspath(String(path)))
    rel = try
        relpath(abs_path, PROJECT_ROOT)
    catch
        nothing
    end
    if rel !== nothing
        rel_s = String(rel)
        if rel_s == "." || !(startswith(rel_s, "..") || startswith(rel_s, string("..", Base.Filesystem.path_separator)))
            return _norm_slash(rel_s)
        end
        return _norm_slash(rel_s)
    end
    return _norm_slash(abs_path)
end

using Dates
using SHA

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
    base_output_dir::String
end

struct RunContext
    run_id::String
    timestamp_utc::String
    argv::Vector{String}
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
    println("  --base-output-dir <path>  base directory for outputs (default data/outputs)")
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
        :base_output_dir => joinpath("data", "outputs"),
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
        elseif arg == "--base-output-dir"
            opts[:base_output_dir] = require_value()
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
        String(opts[:base_output_dir]),
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

@inline function _json_escape(s::AbstractString)
    out = IOBuffer()
    for c in s
        if c == '"'
            print(out, "\\\"")
        elseif c == '\\'
            print(out, "\\\\")
        elseif c == '\n'
            print(out, "\\n")
        elseif c == '\r'
            print(out, "\\r")
        elseif c == '\t'
            print(out, "\\t")
        else
            print(out, c)
        end
    end
    return String(take!(out))
end

function _to_json(x)
    if x === nothing
        return "null"
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x isa Integer || x isa AbstractFloat
        return string(x)
    elseif x isa AbstractString
        return "\"$(_json_escape(x))\""
    elseif x isa AbstractVector
        return "[" * join((_to_json(v) for v in x), ",") * "]"
    elseif x isa AbstractDict
        pairs_sorted = sort(collect(pairs(x)); by=kv -> String(kv.first))
        parts = String[]
        for (k, v) in pairs_sorted
            push!(parts, "\"$(_json_escape(String(k)))\":" * _to_json(v))
        end
        return "{" * join(parts, ",") * "}"
    else
        return "\"$(_json_escape(string(x)))\""
    end
end

function _write_json(path::String, x)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, _to_json(x))
    end
end

function _current_git_commit()
    try
        return readchomp(`git -C $(PROJECT_ROOT) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _current_git_branch()
    try
        return readchomp(`git -C $(PROJECT_ROOT) rev-parse --abbrev-ref HEAD`)
    catch
        return "unknown"
    end
end

function _git_dirty()
    try
        return !isempty(readchomp(`git -C $(PROJECT_ROOT) status --porcelain`))
    catch
        return false
    end
end

function _sha256_file(path::String)
    open(path, "r") do io
        return bytes2hex(sha256(read(io)))
    end
end

function _count_csv_rows(path::String)
    rows = 0
    seen_header = false
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if !seen_header
                seen_header = true
                continue
            end
            rows += 1
        end
    end
    return rows
end

function _first_csv_header(path::String)
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            return [strip(x) for x in split(s, ',')]
        end
    end
    return String[]
end

function _csv_stats(path::String)
    header = String[]
    idx_seed = 0
    idx_conv = 0
    total = 0
    success = 0
    fallback = 0
    errors = 0

    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            if isempty(header)
                header = [strip(x) for x in split(s, ',')]
                idx_seed = findfirst(==("seed_source"), header)
                idx_conv = findfirst(==("converged"), header)
                continue
            end

            parts = split(s, ',')
            total += 1

            if idx_conv > 0 && idx_conv <= length(parts)
                v = lowercase(strip(parts[idx_conv]))
                if v == "true"
                    success += 1
                else
                    errors += 1
                end
            end

            if idx_seed > 0 && idx_seed <= length(parts)
                occursin("legacy_fallback", parts[idx_seed]) && (fallback += 1)
            end
        end
    end

    return (points_total=total, success_count=success, fallback_count=fallback, error_count=errors)
end

function _annotate_csv_with_run_id(path::String, run_id::String)
    tmp = path * ".tmp"
    open(path, "r") do src
        open(tmp, "w") do dst
            seen_header = false
            for line in eachline(src)
                s = chomp(line)
                if startswith(strip(s), "#") || isempty(strip(s))
                    println(dst, s)
                    continue
                end
                if !seen_header
                    println(dst, s, ",run_id")
                    seen_header = true
                else
                    println(dst, s, ",", run_id)
                end
            end
        end
    end
    mv(tmp, path; force=true)
end

function _source_T_from_path(path::String)
    m = match(r"transport_vs_xi_T([0-9]+)_muB0", basename(path))
    return m === nothing ? "" : m.captures[1]
end

function _artifact_entry(path::String)
    abs_path = abspath(path)
    return Dict{String,Any}(
        "path" => _repo_relpath(abs_path),
        "sha256" => _sha256_file(abs_path),
        "rows" => _count_csv_rows(abs_path),
        "schema_version" => "scan_csv_v1",
    )
end

function _write_provenance_sidecars(opts::Options, ctx::RunContext, section::String, out_dir::String, config::Dict{String,Any}, artifacts::Vector{String}, summary::Dict{String,Any})
    effective = Dict{String,Any}(
        "schema_version" => "v1",
        "workflow" => "manual_relaxation_scan",
        "section" => section,
        "options" => config,
    )
    effective_json = _to_json(effective)
    config_hash = bytes2hex(sha256(effective_json))

    project_toml = joinpath(PROJECT_ROOT, "Project.toml")
    manifest_toml = joinpath(PROJECT_ROOT, "Manifest.toml")

    manifest = Dict{String,Any}(
        "run_id" => ctx.run_id,
        "timestamp_utc" => ctx.timestamp_utc,
        "script" => "scripts/relaxtime/run_manual_relaxation_scan_workflow.jl",
        "argv" => copy(ctx.argv),
        "cwd" => _repo_relpath(pwd()),
        "project_path" => ".",
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
        "git_commit" => _current_git_commit(),
        "git_branch" => _current_git_branch(),
        "git_dirty" => _git_dirty(),
        "project_toml_hash" => isfile(project_toml) ? _sha256_file(project_toml) : "",
        "manifest_toml_hash" => isfile(manifest_toml) ? _sha256_file(manifest_toml) : "",
        "config_hash" => config_hash,
        "schema_version" => "v1",
        "artifacts" => [_artifact_entry(p) for p in artifacts],
        "summary" => summary,
    )

    _write_json(joinpath(out_dir, "effective_config.json"), effective)
    _write_json(joinpath(out_dir, "run_manifest.json"), manifest)
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

function merge_csvs(inputs::Vector{String}, output::String, run_id::String; overwrite::Bool)
    isempty(inputs) && error("no input CSVs to merge")
    ensure_parent_dir(output)
    maybe_rm(output; overwrite=overwrite)

    meta_lines = String[]
    header_line = nothing
    merged_header = nothing
    rows = String[]

    for (index, path) in enumerate(inputs)
        source_file = basename(path)
        source_T = _source_T_from_path(path)
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
                    merged_header = string(s, ",source_file,source_T_MeV")
                    continue
                end
                if s == header_line
                    continue
                end
                push!(rows, string(s, ",", source_file, ",", source_T))
            end
        end
    end

    header_line === nothing && error("merged CSV missing header")
    merged_header === nothing && error("failed to compose merged header")

    open(output, "w") do io
        for line in meta_lines
            println(io, line)
        end
        println(io, merged_header)
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

function run_plan_a(opts::Options, ctx::RunContext)
    out_dir = joinpath(opts.base_output_dir, "results", "relaxtime", "plan_a")
    out_csv = joinpath(out_dir, "gap_transport_vs_T_muB0_xi0.csv")
    out_fig = joinpath(opts.base_output_dir, "figures", "relaxtime", "plan_a")

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
    _annotate_csv_with_run_id(out_csv, ctx.run_id)

    stats = _csv_stats(out_csv)
    summary = Dict{String,Any}(
        "points_total" => stats.points_total,
        "success_count" => stats.success_count,
        "fallback_count" => stats.fallback_count,
        "error_count" => stats.error_count,
    )
    config = Dict{String,Any}(
        "plan_a_Tmin_mev" => opts.plan_a_Tmin_mev,
        "plan_a_Tmax_mev" => opts.plan_a_Tmax_mev,
        "plan_a_Tstep_mev" => opts.plan_a_Tstep_mev,
        "integration_mode" => opts.integration_mode,
        "tau_p_nodes" => opts.tau_p_nodes,
        "tau_angle_nodes" => opts.tau_angle_nodes,
        "tau_phi_nodes" => opts.tau_phi_nodes,
        "tau_n_sigma_points" => opts.tau_n_sigma_points,
        "sigma_grid_n" => opts.sigma_grid_n,
        "compute_bulk" => opts.compute_bulk,
        "seed_policy" => "phase_aware_xi_T_continuity",
    )
    _write_provenance_sidecars(opts, ctx, "plan_a", out_dir, config, [out_csv], summary)

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

function run_plan_b(opts::Options, ctx::RunContext)
    xi_list = xi_list_string(opts.xi_min, opts.xi_max, opts.xi_step)
    result_dir = joinpath(opts.base_output_dir, "results", "relaxtime", "plan_b")
    figure_dir = joinpath(opts.base_output_dir, "figures", "relaxtime", "plan_b")
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
        _annotate_csv_with_run_id(csv_path, ctx.run_id)
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

    merge_csvs(csv_paths, merged_csv, ctx.run_id; overwrite=opts.overwrite)

    merged_stats = _csv_stats(merged_csv)
    summary = Dict{String,Any}(
        "points_total" => merged_stats.points_total,
        "success_count" => merged_stats.success_count,
        "fallback_count" => merged_stats.fallback_count,
        "error_count" => merged_stats.error_count,
    )
    config = Dict{String,Any}(
        "plan_b_T_list_mev" => opts.plan_b_T_list_mev,
        "xi_min" => opts.xi_min,
        "xi_max" => opts.xi_max,
        "xi_step" => opts.xi_step,
        "integration_mode" => opts.integration_mode,
        "tau_p_nodes" => opts.tau_p_nodes,
        "tau_angle_nodes" => opts.tau_angle_nodes,
        "tau_phi_nodes" => opts.tau_phi_nodes,
        "tau_n_sigma_points" => opts.tau_n_sigma_points,
        "sigma_grid_n" => opts.sigma_grid_n,
        "compute_bulk" => opts.compute_bulk,
        "seed_policy" => "phase_aware_xi_T_continuity",
    )
    _write_provenance_sidecars(opts, ctx, "plan_b", result_dir, config, vcat(copy(csv_paths), [merged_csv]), summary)

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
    run_id = string(Dates.format(now(UTC), dateformat"yyyymmddTHHMMSS"), "_", bytes2hex(rand(UInt8, 4)))
    timestamp_utc = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ")
    ctx = RunContext(run_id, timestamp_utc, copy(ARGS))
    cd(PROJECT_ROOT) do
        :cross_section in opts.sections && run_cross_section(opts)
        :plan_a in opts.sections && run_plan_a(opts, ctx)
        :plan_b in opts.sections && run_plan_b(opts, ctx)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
