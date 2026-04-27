#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using JSON3

struct Options
    params_csv::String
    out_root::String
    xi_min::Float64
    xi_max::Float64
    xi_step::Float64
    dry_run::Bool
    resume::Bool
    overwrite::Bool
end

struct Sample
    sample_id::String
    T_MeV::Float64
    muB_MeV::Float64
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_xi_smoothness_batch.jl [options]\n")
    println("Options:")
    println("  --params <csv>      输入采样 CSV (required)")
    println("  --out-root <dir>    输出根目录 (required)")
    println("  --xi-min <float>    xi 扫描下界 (default -0.5)")
    println("  --xi-max <float>    xi 扫描上界 (default 0.5)")
    println("  --xi-step <float>   xi 扫描步长 (default 0.02)")
    println("  --dry-run           仅写 manifest，不执行计算")
    println("  --resume            已有非空 result_csv 时跳过")
    println("  --overwrite         覆盖已有结果文件")
    println("  -h, --help          显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :params => nothing,
        :out_root => nothing,
        :xi_min => -0.5,
        :xi_max => 0.5,
        :xi_step => 0.02,
        :dry_run => false,
        :resume => false,
        :overwrite => false,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $(arg)")
            value = args[i + 1]
            i += 1
            return value
        end

        if arg == "--params"
            opts[:params] = require_value()
        elseif arg == "--out-root"
            opts[:out_root] = require_value()
        elseif arg == "--xi-min"
            opts[:xi_min] = parse(Float64, require_value())
        elseif arg == "--xi-max"
            opts[:xi_max] = parse(Float64, require_value())
        elseif arg == "--xi-step"
            opts[:xi_step] = parse(Float64, require_value())
        elseif arg == "--dry-run"
            opts[:dry_run] = true
        elseif arg == "--resume"
            opts[:resume] = true
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $(arg)")
        end
        i += 1
    end

    isnothing(opts[:params]) && error("--params is required")
    isnothing(opts[:out_root]) && error("--out-root is required")

    params_csv = String(opts[:params])
    out_root = String(opts[:out_root])
    isfile(params_csv) || error("params csv not found: $(params_csv)")
    opts[:xi_step] > 0 || error("xi-step must be positive")
    opts[:xi_min] <= opts[:xi_max] || error("xi-min must be <= xi-max")

    return Options(
        params_csv,
        out_root,
        Float64(opts[:xi_min]),
        Float64(opts[:xi_max]),
        Float64(opts[:xi_step]),
        Bool(opts[:dry_run]),
        Bool(opts[:resume]),
        Bool(opts[:overwrite]),
    )
end

@inline _is_comment_or_empty(line::AbstractString) = isempty(strip(line)) || startswith(strip(line), "#")

function _parse_header_map(line::String)
    cols = [strip(c) for c in split(strip(line), ',')]
    return Dict{String, Int}(c => i for (i, c) in enumerate(cols))
end

function _parse_required_float(parts, idx::Int, label::String, line_no::Int)
    idx <= length(parts) || error("params csv line $(line_no): missing $(label)")
    raw = strip(parts[idx])
    val = tryparse(Float64, raw)
    val === nothing && error("params csv line $(line_no): invalid $(label): $(repr(raw))")
    return Float64(val)
end

function read_samples(path::String)
    samples = Sample[]
    seen_sample_lines = Dict{String, Int}()
    open(path, "r") do io
        header_map = Dict{String, Int}()
        header_seen = false
        line_no = 0
        for raw in eachline(io)
            line_no += 1
            _is_comment_or_empty(raw) && continue

            if !header_seen
                header_map = _parse_header_map(raw)
                haskey(header_map, "sample_id") || error("params csv missing required column: sample_id")
                haskey(header_map, "T_MeV") || error("params csv missing required column: T_MeV")
                has_muB = haskey(header_map, "muB_MeV")
                has_muq = haskey(header_map, "muq_MeV")
                (has_muB || has_muq) || error("params csv must include muB_MeV or muq_MeV")
                header_seen = true
                continue
            end

            parts = split(raw, ',')
            sample_id_idx = header_map["sample_id"]
            T_idx = header_map["T_MeV"]
            sample_id_idx <= length(parts) || error("params csv line $(line_no): missing sample_id")
            sid = strip(parts[sample_id_idx])
            isempty(sid) && error("params csv line $(line_no): empty sample_id")
            if haskey(seen_sample_lines, sid)
                first_line = seen_sample_lines[sid]
                error("params csv line $(line_no): duplicate sample_id '$(sid)' (first seen at line $(first_line))")
            end
            seen_sample_lines[sid] = line_no
            T_mev = _parse_required_float(parts, T_idx, "T_MeV", line_no)

            muB_mev = if haskey(header_map, "muB_MeV")
                _parse_required_float(parts, header_map["muB_MeV"], "muB_MeV", line_no)
            else
                3.0 * _parse_required_float(parts, header_map["muq_MeV"], "muq_MeV", line_no)
            end

            push!(samples, Sample(sid, T_mev, muB_mev))
        end
    end

    isempty(samples) && error("no samples found in params csv: $(path)")
    return samples
end

@inline _norm_slash(path::AbstractString) = replace(String(path), '\\' => '/')

function _repo_relpath(path::AbstractString)
    abs_path = normpath(abspath(String(path)))
    rel = try
        relpath(abs_path, PROJECT_ROOT)
    catch
        abs_path
    end
    return _norm_slash(String(rel))
end

function _fmt_xi(x::Float64)
    return abs(x) < 1e-12 ? "0.0" : string(x)
end

function _xi_list_string(opts::Options)
    vals = collect(range(opts.xi_min; stop=opts.xi_max, step=opts.xi_step))
    return join(_fmt_xi.(vals), ",")
end

function _julia_cmd(argv::Vector{String})
    jc = Base.julia_cmd()
    cmd = Cmd(vcat(jc.exec, argv))
    return jc.env === nothing ? cmd : setenv(cmd, jc.env)
end

function _run_sample(opts::Options, sample::Sample, xi_list::String)
    sample_dir = joinpath(opts.out_root, sample.sample_id)
    mkpath(sample_dir)
    result_csv = joinpath(sample_dir, "gap_transport_scan.csv")
    failed_points_csv = joinpath(sample_dir, "failed_points.csv")

    if opts.dry_run
        return Dict{String, Any}(
            "sample_id" => sample.sample_id,
            "T_MeV" => sample.T_MeV,
            "muB_MeV" => sample.muB_MeV,
            "status" => "skipped",
            "reason" => "dry_run",
            "result_csv" => _repo_relpath(result_csv),
            "failed_points_path" => _repo_relpath(failed_points_csv),
            "command" => "",
        )
    end

    if opts.resume && isfile(result_csv) && filesize(result_csv) > 0 && !opts.overwrite
        return Dict{String, Any}(
            "sample_id" => sample.sample_id,
            "T_MeV" => sample.T_MeV,
            "muB_MeV" => sample.muB_MeV,
            "status" => "skipped",
            "reason" => "resume_hit",
            "result_csv" => _repo_relpath(result_csv),
            "failed_points_path" => _repo_relpath(failed_points_csv),
            "command" => "",
        )
    end

    opts.overwrite && isfile(result_csv) && rm(result_csv; force=true)
    opts.overwrite && isfile(failed_points_csv) && rm(failed_points_csv; force=true)

    argv = String[
        "--project=.",
        "scripts/relaxtime/run_gap_transport_scan.jl",
        "--tmin", string(sample.T_MeV),
        "--tmax", string(sample.T_MeV),
        "--tstep", "1",
        "--mubmin", string(sample.muB_MeV),
        "--mubmax", string(sample.muB_MeV),
        "--mubstep", "1",
        "--xi-list", xi_list,
        "--output", result_csv,
        "--failed-points-output", failed_points_csv,
    ]
    opts.overwrite && push!(argv, "--overwrite")

    cmd = _julia_cmd(argv)
    cmd_text = sprint(show, cmd)

    status = "success"
    reason = ""
    try
        run(cmd)
    catch err
        status = "failed"
        reason = sprint(showerror, err)
    end

    return Dict{String, Any}(
        "sample_id" => sample.sample_id,
        "T_MeV" => sample.T_MeV,
        "muB_MeV" => sample.muB_MeV,
        "status" => status,
        "reason" => reason,
        "result_csv" => _repo_relpath(result_csv),
        "failed_points_path" => _repo_relpath(failed_points_csv),
        "command" => cmd_text,
    )
end

function _write_manifest(path::String, payload::Dict{String, Any})
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, JSON3.write(payload))
    end
end

function main(args::Vector{String}=copy(ARGS))
    opts = parse_args(args)
    samples = read_samples(opts.params_csv)
    xi_list = _xi_list_string(opts)

    mkpath(opts.out_root)
    entries = Dict{String, Any}[]
    success_count = 0
    failed_count = 0
    skipped_count = 0

    cd(PROJECT_ROOT) do
        for sample in samples
            entry = _run_sample(opts, sample, xi_list)
            push!(entries, entry)

            status = String(entry["status"])
            if status == "success"
                success_count += 1
            elseif status == "failed"
                failed_count += 1
            else
                skipped_count += 1
            end
        end
    end

    manifest = Dict{String, Any}(
        "schema_version" => "v1",
        "script" => "scripts/relaxtime/run_xi_smoothness_batch.jl",
        "params_csv" => _repo_relpath(opts.params_csv),
        "out_root" => _repo_relpath(opts.out_root),
        "dry_run" => opts.dry_run,
        "resume" => opts.resume,
        "overwrite" => opts.overwrite,
        "xi_min" => opts.xi_min,
        "xi_max" => opts.xi_max,
        "xi_step" => opts.xi_step,
        "summary" => Dict(
            "total" => length(entries),
            "success" => success_count,
            "failed" => failed_count,
            "skipped" => skipped_count,
        ),
        "samples" => entries,
    )

    manifest_path = joinpath(opts.out_root, "run_manifest.json")
    _write_manifest(manifest_path, manifest)
    println("Wrote batch manifest: " * manifest_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
