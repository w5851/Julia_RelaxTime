"""
离线补点：读取 scan CSV 中被 quality_flag 标记的点，
逐点调用 run_gap_transport_scan.jl 做高精重算，输出 patch CSV（不自动回写原文件）。
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "offline_patch_utils.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "scan_quality.jl"))
include(joinpath(PROJECT_ROOT, "scripts", "utils", "scan_csv.jl"))

using .OfflinePatchUtils: read_flagged_points, select_patch_points
using .ScanQuality: assess_point_quality
using .ScanCSV: ScanCSV

struct PatchOptions
    input::String
    output::String
    max_points::Int
    reason_filter::Union{Nothing, String}
    overwrite::Bool
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
    tr_p_nodes::Int
    tr_p_max_fm::Float64
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_offline_transport_patch.jl [options]")
    println("Options:")
    println("  --input <path>              输入 scan CSV（需含 quality_flag）")
    println("  --output <path>             输出 patch CSV")
    println("  --max-points <int>          最多重算点数（按 |quality_metric| 取 Top-K；0=全部）")
    println("  --reason-filter <csv>       仅重算指定 reason（逗号分隔）")
    println("  --overwrite                 覆盖输出文件")
    println("  --compute-bulk              重算时包含 bulk")
    println("  --p-num/--t-num/--max-iter  平衡态求解参数")
    println("  --tau-p-nodes <int>")
    println("  --tau-angle-nodes <int>")
    println("  --tau-phi-nodes <int>")
    println("  --tau-n-sigma <int>")
    println("  --tau-threshold-subtraction / --tau-no-threshold-subtraction")
    println("  --tau-asym-window <float>")
    println("  --tau-asym-fit-min-points <int>")
    println("  --tau-asym-extra-points <int>")
    println("  --tau-interpolation-mode <pchip|linear|direct|hybrid_threshold>")
    println("  --tr-p-nodes <int>")
    println("  --tr-p-max <fm^-1>")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :input => joinpath("data", "outputs", "results", "relaxtime", "scan", "gap_transport_scan.csv"),
        :output => joinpath("data", "outputs", "results", "relaxtime", "scan", "gap_transport_scan_patch.csv"),
        :max_points => 0,
        :reason_filter => nothing,
        :overwrite => false,
        :compute_bulk => false,
        :p_num => 12,
        :t_num => 6,
        :max_iter => 40,
        :tau_p_nodes => 36,
        :tau_angle_nodes => 8,
        :tau_phi_nodes => 16,
        :tau_n_sigma_points => 128,
        :tau_threshold_subtraction => true,
        :tau_asym_window => 0.6,
        :tau_asym_fit_min_points => 12,
        :tau_asym_extra_points => 16,
        :tau_interpolation_mode => :hybrid_threshold,
        :tr_p_nodes => 36,
        :tr_p_max => 8.0,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function need_value()
            i == length(args) && error("missing value for $arg")
            v = args[i + 1]
            i += 1
            return v
        end

        if arg == "--input"
            opts[:input] = need_value()
        elseif arg == "--output"
            opts[:output] = need_value()
        elseif arg == "--max-points"
            opts[:max_points] = parse(Int, need_value())
        elseif arg == "--reason-filter"
            opts[:reason_filter] = need_value()
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--compute-bulk"
            opts[:compute_bulk] = true
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, need_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, need_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, need_value())
        elseif arg == "--tau-p-nodes"
            opts[:tau_p_nodes] = parse(Int, need_value())
        elseif arg == "--tau-angle-nodes"
            opts[:tau_angle_nodes] = parse(Int, need_value())
        elseif arg == "--tau-phi-nodes"
            opts[:tau_phi_nodes] = parse(Int, need_value())
        elseif arg == "--tau-n-sigma"
            opts[:tau_n_sigma_points] = parse(Int, need_value())
        elseif arg == "--tau-threshold-subtraction"
            opts[:tau_threshold_subtraction] = true
        elseif arg == "--tau-no-threshold-subtraction"
            opts[:tau_threshold_subtraction] = false
        elseif arg == "--tau-asym-window"
            opts[:tau_asym_window] = parse(Float64, need_value())
        elseif arg == "--tau-asym-fit-min-points"
            opts[:tau_asym_fit_min_points] = parse(Int, need_value())
        elseif arg == "--tau-asym-extra-points"
            opts[:tau_asym_extra_points] = parse(Int, need_value())
        elseif arg == "--tau-interpolation-mode"
            mode = Symbol(need_value())
            mode in (:pchip, :linear, :direct, :hybrid_threshold) || error("unknown tau interpolation mode: $mode")
            opts[:tau_interpolation_mode] = mode
        elseif arg == "--tr-p-nodes"
            opts[:tr_p_nodes] = parse(Int, need_value())
        elseif arg == "--tr-p-max"
            opts[:tr_p_max] = parse(Float64, need_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return PatchOptions(
        String(opts[:input]),
        String(opts[:output]),
        Int(opts[:max_points]),
        isnothing(opts[:reason_filter]) ? nothing : String(opts[:reason_filter]),
        Bool(opts[:overwrite]),
        Bool(opts[:compute_bulk]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        Int(opts[:max_iter]),
        Int(opts[:tau_p_nodes]),
        Int(opts[:tau_angle_nodes]),
        Int(opts[:tau_phi_nodes]),
        Int(opts[:tau_n_sigma_points]),
        Bool(opts[:tau_threshold_subtraction]),
        Float64(opts[:tau_asym_window]),
        Int(opts[:tau_asym_fit_min_points]),
        Int(opts[:tau_asym_extra_points]),
        Symbol(opts[:tau_interpolation_mode]),
        Int(opts[:tr_p_nodes]),
        Float64(opts[:tr_p_max]),
    )
end

function ensure_parent_dir(path::AbstractString)
    parent = dirname(path)
    isempty(parent) && return
    isdir(parent) || mkpath(parent)
end

function current_git_commit()
    try
        return readchomp(`git rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _is_comment(line::AbstractString)
    s = strip(line)
    return isempty(s) || startswith(s, "#")
end

function _read_single_data_row(path::AbstractString)
    open(path, "r") do io
        header = nothing
        for line in eachline(io)
            _is_comment(line) && continue
            if header === nothing
                header = split(strip(line), ',')
                continue
            end
            vals = split(strip(line), ',')
            length(vals) == length(header) || continue
            return Dict{String, String}(String(strip(header[i])) => String(strip(vals[i])) for i in eachindex(header))
        end
    end
    return nothing
end

function write_patch_header_if_needed(io)
    cols = [
        "T_MeV", "muB_MeV", "xi",
        "quality_reason", "quality_metric",
        "tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
        "tauinv_u", "tauinv_d", "tauinv_s", "tauinv_ubar", "tauinv_dbar", "tauinv_sbar",
        "eta", "sigma", "zeta", "eta_over_s", "zeta_over_s", "sigma_over_T",
        "patch_quality_flag", "patch_quality_reason", "patch_quality_metric",
    ]
    ScanCSV.write_header(io, cols)
end

function _require_float(row::Dict{String, String}, key::String)
    haskey(row, key) || error("missing key in point csv: $key")
    v = tryparse(Float64, row[key])
    v === nothing && error("invalid float for key $key: $(row[key])")
    return Float64(v)
end

function _build_scan_cmd(opts::PatchOptions, p)
    scan_script = joinpath(PROJECT_ROOT, "scripts", "relaxtime", "run_gap_transport_scan.jl")
    cmd_parts = String[
        "julia", "--project=.", scan_script,
        "--tmin", string(p.T_MeV), "--tmax", string(p.T_MeV), "--tstep", "1",
        "--mubmin", string(p.muB_MeV), "--mubmax", string(p.muB_MeV), "--mubstep", "1",
        "--xi", string(p.xi),
        "--no-resume",
        "--overwrite",
        "--p-num", string(opts.p_num),
        "--t-num", string(opts.t_num),
        "--max-iter", string(opts.max_iter),
        "--tau-p-nodes", string(opts.tau_p_nodes),
        "--tau-angle-nodes", string(opts.tau_angle_nodes),
        "--tau-phi-nodes", string(opts.tau_phi_nodes),
        "--tau-n-sigma", string(opts.tau_n_sigma_points),
        "--tau-asym-window", string(opts.tau_asym_window),
        "--tau-asym-fit-min-points", string(opts.tau_asym_fit_min_points),
        "--tau-asym-extra-points", string(opts.tau_asym_extra_points),
        "--tau-interpolation-mode", string(opts.tau_interpolation_mode),
        "--tr-p-nodes", string(opts.tr_p_nodes),
        "--tr-p-max", string(opts.tr_p_max_fm),
    ]
    if opts.compute_bulk
        push!(cmd_parts, "--compute-bulk")
    end
    if opts.tau_threshold_subtraction
        push!(cmd_parts, "--tau-threshold-subtraction")
    else
        push!(cmd_parts, "--tau-no-threshold-subtraction")
    end
    return cmd_parts
end

function run_patch(opts::PatchOptions)
    flagged = read_flagged_points(opts.input)
    selected = select_patch_points(flagged; max_points=opts.max_points, reason_filter=opts.reason_filter)

    isempty(selected) && begin
        println("No flagged points selected. Nothing to do.")
        return
    end

    if isfile(opts.output)
        opts.overwrite || error("output exists: $(opts.output). Use --overwrite to replace it.")
        rm(opts.output)
    end
    ensure_parent_dir(opts.output)

    io = open(opts.output, "a")
    try
        ScanCSV.write_metadata(io, Dict(
            "schema" => "scan_csv_v1",
            "title" => "offline_transport_patch",
            "script" => "scripts/relaxtime/run_offline_transport_patch.jl",
            "git_commit" => current_git_commit(),
            "source_csv" => opts.input,
            "selection_count" => string(length(selected)),
        ))
        write_patch_header_if_needed(io)

        for (idx, p) in enumerate(selected)
            tmp_out = tempname() * ".csv"
            cmd_parts = _build_scan_cmd(opts, p)
            push!(cmd_parts, "--output", tmp_out)

            run(Cmd(cmd_parts; dir=PROJECT_ROOT))
            row = _read_single_data_row(tmp_out)
            rm(tmp_out; force=true)
            row === nothing && error("failed to read computed point row for T=$(p.T_MeV), muB=$(p.muB_MeV), xi=$(p.xi)")

            tau = (
                u=_require_float(row, "tau_u"),
                d=_require_float(row, "tau_d"),
                s=_require_float(row, "tau_s"),
                ubar=_require_float(row, "tau_ubar"),
                dbar=_require_float(row, "tau_dbar"),
                sbar=_require_float(row, "tau_sbar"),
            )
            eta_over_s = _require_float(row, "eta_over_s")
            sigma_over_T = _require_float(row, "sigma_over_T")
            patch_flag, patch_reason, patch_metric = assess_point_quality(tau, eta_over_s, sigma_over_T)

            out_row = join([
                string(p.T_MeV), string(p.muB_MeV), string(p.xi),
                p.quality_reason, string(p.quality_metric),
                row["tau_u"], row["tau_d"], row["tau_s"], row["tau_ubar"], row["tau_dbar"], row["tau_sbar"],
                row["tauinv_u"], row["tauinv_d"], row["tauinv_s"], row["tauinv_ubar"], row["tauinv_dbar"], row["tauinv_sbar"],
                row["eta"], row["sigma"], row["zeta"], row["eta_over_s"], row["zeta_over_s"], row["sigma_over_T"],
                patch_flag ? "true" : "false", patch_reason, string(patch_metric),
            ], ',')
            println(io, out_row)
            flush(io)

            if idx % 5 == 0
                println("patch progress: $(idx)/$(length(selected))")
            end
        end
    finally
        close(io)
    end

    println("Offline patch finished. Output: $(opts.output)")
end

function main()
    opts = parse_args(copy(ARGS))
    run_patch(opts)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
