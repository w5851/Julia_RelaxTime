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
    reason_filter::Union{Nothing,String}
    overwrite::Bool
    compute_bulk::Bool
    sampling_profile::String
    max_iter::Int
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_offline_transport_patch.jl [options]")
    println("  --input / --output <path>      Source scan and separate patch CSV")
    println("  --max-points <int>             Top K flagged points; 0 selects all")
    println("  --reason-filter <csv>          Restrict quality reasons")
    println("  --sampling-profile <name|path> Default thermal_weighted_v1; use a new profile for refinement")
    println("  --max-iter <int>               Solver iteration limit")
    println("  --no-compute-bulk              Explicitly omit bulk viscosity")
    println("  --overwrite                   Replace this patch output")
end

if !isdefined(Main,:TransportSamplingScanConfig)
    Base.include(Main,joinpath(@__DIR__,"transport_sampling_scan_config.jl"))
end

function parse_args(args::Vector{String})
    Main.TransportSamplingScanConfig.reject_retired_sampling_flags(args)
    opts=Dict{Symbol,Any}(
        :input=>joinpath("data","outputs","results","relaxtime","scan","gap_transport_scan.csv"),
        :output=>joinpath("data","outputs","results","relaxtime","scan","gap_transport_scan_patch.csv"),
        :max_points=>0,:reason_filter=>nothing,:overwrite=>false,:compute_bulk=>true,
        :sampling_profile=>"thermal_weighted_v1",:max_iter=>40)
    i=1
    while i<=length(args)
        arg=args[i]
        if arg in ("--overwrite","--compute-bulk","--no-compute-bulk")
            key=arg=="--overwrite" ? :overwrite : :compute_bulk
            opts[key]=arg!="--no-compute-bulk"
        elseif arg in ("-h","--help")
            print_usage();exit(0)
        elseif arg in ("--input","--output","--max-points","--reason-filter","--sampling-profile","--max-iter")
            i<length(args) || throw(ArgumentError("missing value for $arg"))
            i+=1;key=Symbol(replace(arg[3:end],'-'=>'_'))
            opts[key]=key in (:max_points,:max_iter) ? parse(Int,args[i]) : args[i]
        else
            throw(ArgumentError("unsupported option $arg; quadrature settings belong in --sampling-profile"))
        end
        i+=1
    end
    opts[:max_points]>=0 && opts[:max_iter]>0 || throw(ArgumentError("invalid patch limits"))
    abspath(opts[:input])==abspath(opts[:output]) && throw(ArgumentError("patch output must differ from its source"))
    return PatchOptions(opts[:input],opts[:output],opts[:max_points],opts[:reason_filter],
        opts[:overwrite],opts[:compute_bulk],opts[:sampling_profile],opts[:max_iter])
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

function _build_scan_cmd(opts::PatchOptions,p)
    scan_script=joinpath(PROJECT_ROOT,"scripts","relaxtime","run_gap_transport_scan.jl")
    cmd_parts=String[Base.julia_cmd().exec...,"--project=$PROJECT_ROOT",scan_script,
        "--tmin",string(p.T_MeV),"--tmax",string(p.T_MeV),"--tstep","1",
        "--mubmin",string(p.muB_MeV),"--mubmax",string(p.muB_MeV),"--mubstep","1",
        "--xi",string(p.xi),"--no-resume","--overwrite",
        "--sampling-profile",opts.sampling_profile,"--max-iter",string(opts.max_iter)]
    push!(cmd_parts,opts.compute_bulk ? "--compute-bulk" : "--no-compute-bulk")
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
            "sampling_profile" => opts.sampling_profile,
        ))
        write_patch_header_if_needed(io)

        for (idx, p) in enumerate(selected)
            point_dir = joinpath(splitext(opts.output)[1]*"_points","point_"*lpad(string(idx),4,'0'))
            mkpath(point_dir)
            tmp_out = joinpath(point_dir,"result.csv")
            cmd_parts = _build_scan_cmd(opts, p)
            push!(cmd_parts, "--output", tmp_out)

            run(Cmd(cmd_parts; dir=PROJECT_ROOT))
            row = _read_single_data_row(tmp_out)
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
