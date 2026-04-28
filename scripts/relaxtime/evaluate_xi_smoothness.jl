#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_OUT_ROOT = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "plan_c", "smoothness")
const DEFAULT_FIELDS = ("tau_u", "tau_s", "eta_over_s", "zeta_over_s", "sigma_over_T")

using JSON3
using Printf

struct Options
    manifest::String
    out_root::String
    smooth_s2::Float64
    suspect_s2::Float64
    smooth_s1jump::Float64
    suspect_s1jump::Float64
    suspect_max_spikes::Int
    spike_threshold::Float64
end

@inline _is_comment_or_empty(line::AbstractString) = isempty(strip(line)) || startswith(strip(line), "#")

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/evaluate_xi_smoothness.jl --manifest <run_manifest.json> [options]")
    println("Options:")
    println("  --manifest <path>            输入 run_manifest.json (required)")
    println("  --out-root <dir>             输出目录 (default data/outputs/results/relaxtime/plan_c/smoothness)")
    println("  --smooth-s2 <float>          smooth 的 S2 阈值 (default 0.05)")
    println("  --suspect-s2 <float>         suspect 的 S2 阈值 (default 0.15)")
    println("  --smooth-s1jump <float>      smooth 的 S1jump 阈值 (default 0.15)")
    println("  --suspect-s1jump <float>     suspect 的 S1jump 阈值 (default 0.40)")
    println("  --suspect-max-spikes <int>   suspect 的 N_spike 上限 (default 2)")
    println("  --spike-threshold <float>    N_spike 计数阈值 (default 1.0)")
    println("  -h, --help                   显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :manifest => nothing,
        :out_root => DEFAULT_OUT_ROOT,
        :smooth_s2 => 0.05,
        :suspect_s2 => 0.15,
        :smooth_s1jump => 0.15,
        :suspect_s1jump => 0.40,
        :suspect_max_spikes => 2,
        :spike_threshold => 1.0,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $(arg)"))
            i += 1
            return args[i]
        end

        if arg == "--manifest"
            opts[:manifest] = require_value()
        elseif arg == "--out-root"
            opts[:out_root] = require_value()
        elseif arg == "--smooth-s2"
            opts[:smooth_s2] = parse(Float64, require_value())
        elseif arg == "--suspect-s2"
            opts[:suspect_s2] = parse(Float64, require_value())
        elseif arg == "--smooth-s1jump"
            opts[:smooth_s1jump] = parse(Float64, require_value())
        elseif arg == "--suspect-s1jump"
            opts[:suspect_s1jump] = parse(Float64, require_value())
        elseif arg == "--suspect-max-spikes"
            opts[:suspect_max_spikes] = parse(Int, require_value())
        elseif arg == "--spike-threshold"
            opts[:spike_threshold] = parse(Float64, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $(arg)"))
        end
        i += 1
    end

    opts[:manifest] === nothing && throw(ArgumentError("--manifest is required"))
    opts[:smooth_s2] <= opts[:suspect_s2] || throw(ArgumentError("require smooth_s2 <= suspect_s2"))
    opts[:smooth_s1jump] <= opts[:suspect_s1jump] || throw(ArgumentError("require smooth_s1jump <= suspect_s1jump"))
    opts[:suspect_max_spikes] >= 0 || throw(ArgumentError("suspect_max_spikes must be non-negative"))
    opts[:spike_threshold] > 0 || throw(ArgumentError("spike_threshold must be positive"))

    return Options(
        normpath(abspath(String(opts[:manifest]))),
        normpath(abspath(String(opts[:out_root]))),
        Float64(opts[:smooth_s2]),
        Float64(opts[:suspect_s2]),
        Float64(opts[:smooth_s1jump]),
        Float64(opts[:suspect_s1jump]),
        Int(opts[:suspect_max_spikes]),
        Float64(opts[:spike_threshold]),
    )
end

@inline _norm_slash(path::AbstractString) = replace(String(path), '\\' => '/')

@inline function _series_scale(values::Vector{Float64})
    isempty(values) && return 1.0
    return max(sum(abs, values) / length(values), 1e-12)
end

function _normalized_second_differences(values::Vector{Float64})
    n = length(values)
    n < 3 && return Float64[]
    scale = _series_scale(values)
    out = Vector{Float64}(undef, n - 2)
    @inbounds for i in 2:(n - 1)
        out[i - 1] = (values[i + 1] - 2.0 * values[i] + values[i - 1]) / scale
    end
    return out
end

function compute_s2(values::Vector{Float64})
    d2 = _normalized_second_differences(values)
    isempty(d2) && return 0.0
    return sqrt(sum(abs2, d2) / length(d2))
end

function compute_s1jump(values::Vector{Float64})
    n = length(values)
    n < 2 && return 0.0
    scale = _series_scale(values)
    max_jump = 0.0
    @inbounds for i in 2:n
        jump = abs(values[i] - values[i - 1]) / scale
        max_jump = max(max_jump, jump)
    end
    return max_jump
end

function count_spikes(values::Vector{Float64}; spike_threshold::Float64=1.0)
    d2 = _normalized_second_differences(values)
    return count(v -> abs(v) > spike_threshold, d2)
end

function classify_smoothness(
    s2::Float64,
    s1jump::Float64,
    n_spike::Int;
    smooth_s2::Float64,
    suspect_s2::Float64,
    smooth_s1jump::Float64,
    suspect_s1jump::Float64,
    suspect_max_spikes::Int,
)
    if s2 <= smooth_s2 && s1jump <= smooth_s1jump && n_spike == 0
        return "smooth", "all_metrics_within_smooth_thresholds"
    end

    if s2 <= suspect_s2 && s1jump <= suspect_s1jump && n_spike <= suspect_max_spikes
        reasons = String[]
        s2 > smooth_s2 && push!(reasons, "S2_above_smooth")
        s1jump > smooth_s1jump && push!(reasons, "S1jump_above_smooth")
        n_spike > 0 && push!(reasons, "N_spike_nonzero")
        isempty(reasons) && push!(reasons, "within_suspect_thresholds")
        return "suspect", join(reasons, ";")
    end

    reasons = String[]
    s2 > suspect_s2 && push!(reasons, "S2_above_suspect")
    s1jump > suspect_s1jump && push!(reasons, "S1jump_above_suspect")
    n_spike > suspect_max_spikes && push!(reasons, "N_spike_above_suspect")
    isempty(reasons) && push!(reasons, "exceeds_dual_threshold_policy")
    return "not_smooth", join(reasons, ";")
end

@inline function _is_subpath(child::AbstractString, parent::AbstractString)
    abs_child = normpath(abspath(String(child)))
    abs_parent = normpath(abspath(String(parent)))
    rel = relpath(abs_child, abs_parent)
    return rel == "." || rel == abs_parent || !startswith(rel, "..")
end

function _resolve_result_csv(path_like::AbstractString, manifest_dir::String)
    allowed_roots = (normpath(abspath(manifest_dir)), normpath(abspath(PROJECT_ROOT)))

    function _validate_path(path::AbstractString)
        any(root -> _is_subpath(path, root), allowed_roots) ||
            throw(ArgumentError("result_csv path escapes allowed roots: $(path)"))
        return path
    end

    p = String(path_like)
    if isabspath(p)
        return _validate_path(normpath(abspath(p)))
    end
    cand_manifest = normpath(abspath(joinpath(manifest_dir, p)))
    isfile(cand_manifest) && return _validate_path(cand_manifest)
    return _validate_path(normpath(abspath(joinpath(PROJECT_ROOT, p))))
end

function _read_scan_field_series(path::String, fields::Tuple{Vararg{String}})
    isfile(path) || throw(ArgumentError("result_csv not found: $(path)"))

    header_map = Dict{String, Int}()
    header_seen = false
    xi_vals = Float64[]
    field_data = Dict{String, Vector{Float64}}(f => Float64[] for f in fields)

    open(path, "r") do io
        line_no = 0
        for raw in eachline(io)
            line_no += 1
            _is_comment_or_empty(raw) && continue

            if !header_seen
                cols = [strip(c) for c in split(strip(raw), ',')]
                header_map = Dict{String, Int}(c => i for (i, c) in enumerate(cols))
                haskey(header_map, "xi") || throw(ArgumentError("$(path): missing required column: xi"))
                for f in fields
                    haskey(header_map, f) || throw(ArgumentError("$(path): missing required column: $(f)"))
                end
                header_seen = true
                continue
            end

            parts = split(raw, ',')
            idx_xi = header_map["xi"]
            idx_xi <= length(parts) || throw(ArgumentError("$(path):$(line_no): missing xi value"))
            xi = tryparse(Float64, strip(parts[idx_xi]))
            xi === nothing && throw(ArgumentError("$(path):$(line_no): invalid xi value"))
            push!(xi_vals, Float64(xi))

            for f in fields
                idx = header_map[f]
                idx <= length(parts) || throw(ArgumentError("$(path):$(line_no): missing $(f) value"))
                v = tryparse(Float64, strip(parts[idx]))
                v === nothing && throw(ArgumentError("$(path):$(line_no): invalid $(f) value"))
                push!(field_data[f], Float64(v))
            end
        end
    end

    isempty(xi_vals) && throw(ArgumentError("$(path): no data rows"))
    p = sortperm(xi_vals)
    sorted = Dict{String, Vector{Float64}}()
    for f in fields
        sorted[f] = field_data[f][p]
    end
    return sorted
end

function _write_scores_csv(path::String, rows::Vector{NamedTuple})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "sample_id,field,S2,S1jump,N_spike,label,reason")
        for r in rows
            println(io, string(
                r.sample_id, ',', r.field, ',',
                @sprintf("%.8g", r.S2), ',',
                @sprintf("%.8g", r.S1jump), ',',
                r.N_spike, ',',
                r.label, ',',
                r.reason,
            ))
        end
    end
end

function _write_flags_csv(path::String, rows::Vector{NamedTuple})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "sample_id,field,label,reason,S2,S1jump,N_spike")
        for r in rows
            println(io, string(
                r.sample_id, ',', r.field, ',', r.label, ',', r.reason, ',',
                @sprintf("%.8g", r.S2), ',',
                @sprintf("%.8g", r.S1jump), ',',
                r.N_spike,
            ))
        end
    end
end

function _write_review_queue_csv(path::String, rows::Vector{NamedTuple})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "sample_id,field,label,reason,S2,S1jump,N_spike")
        for r in rows
            r.label == "smooth" && continue
            println(io, string(
                r.sample_id, ',', r.field, ',', r.label, ',', r.reason, ',',
                @sprintf("%.8g", r.S2), ',',
                @sprintf("%.8g", r.S1jump), ',',
                r.N_spike,
            ))
        end
    end
end

function _evaluate_manifest(opts::Options)
    isfile(opts.manifest) || throw(ArgumentError("manifest not found: $(opts.manifest)"))
    manifest_dir = dirname(opts.manifest)
    payload = JSON3.read(read(opts.manifest, String))
    haskey(payload, "samples") || throw(ArgumentError("manifest missing samples"))

    rows = NamedTuple[]
    skipped_samples = 0
    total_samples = 0
    for (idx, sample) in enumerate(payload["samples"])
        total_samples += 1
        sample_id = haskey(sample, "sample_id") ? strip(String(sample["sample_id"])) : "sample_$(idx)"
        isempty(sample_id) && (sample_id = "sample_$(idx)")

        if haskey(sample, "status")
            status = lowercase(strip(String(sample["status"])))
            if status != "success"
                skipped_samples += 1
                @warn "skip sample with non-success status" sample_id status
                continue
            end
        end

        if !haskey(sample, "result_csv")
            skipped_samples += 1
            @warn "skip sample: missing result_csv" sample_id
            continue
        end
        result_csv = _resolve_result_csv(String(sample["result_csv"]), manifest_dir)
        if !isfile(result_csv)
            skipped_samples += 1
            @warn "skip sample: result_csv not found" sample_id result_csv
            continue
        end

        series_map = try
            _read_scan_field_series(result_csv, DEFAULT_FIELDS)
        catch err
            skipped_samples += 1
            @warn "skip sample: failed to parse result_csv" sample_id result_csv exception=(err, catch_backtrace())
            continue
        end

        for field in DEFAULT_FIELDS
            ys = series_map[field]
            s2 = compute_s2(ys)
            s1 = compute_s1jump(ys)
            ns = count_spikes(ys; spike_threshold=opts.spike_threshold)
            label, reason = classify_smoothness(
                s2,
                s1,
                ns;
                smooth_s2=opts.smooth_s2,
                suspect_s2=opts.suspect_s2,
                smooth_s1jump=opts.smooth_s1jump,
                suspect_s1jump=opts.suspect_s1jump,
                suspect_max_spikes=opts.suspect_max_spikes,
            )
            push!(rows, (
                sample_id=sample_id,
                field=field,
                S2=s2,
                S1jump=s1,
                N_spike=ns,
                label=label,
                reason=reason,
            ))
        end
    end

    length(rows) == 0 && throw(ArgumentError(
        "all samples were skipped (skipped_samples=$(skipped_samples), total_samples=$(total_samples)); ensure manifest has at least one success sample with a readable result_csv",
    ))

    return rows, skipped_samples
end

function main(args::Vector{String}=copy(ARGS))
    opts = parse_args(args)
    rows, skipped_samples = _evaluate_manifest(opts)
    mkpath(opts.out_root)

    scores_csv = joinpath(opts.out_root, "smoothness_scores.csv")
    flags_csv = joinpath(opts.out_root, "smoothness_flags.csv")
    review_csv = joinpath(opts.out_root, "review_queue.csv")

    _write_scores_csv(scores_csv, rows)
    _write_flags_csv(flags_csv, rows)
    _write_review_queue_csv(review_csv, rows)

    smooth_count = count(r -> r.label == "smooth", rows)
    suspect_count = count(r -> r.label == "suspect", rows)
    not_smooth_count = count(r -> r.label == "not_smooth", rows)
    println("xi smoothness evaluation done: total=$(length(rows)), smooth=$(smooth_count), suspect=$(suspect_count), not_smooth=$(not_smooth_count), skipped_samples=$(skipped_samples)")
    println("scores: " * _norm_slash(scores_csv))
    println("flags: " * _norm_slash(flags_csv))
    println("review_queue: " * _norm_slash(review_csv))
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
