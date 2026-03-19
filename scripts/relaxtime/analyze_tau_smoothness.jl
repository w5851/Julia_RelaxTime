"""
分析 τ(T) 曲线中的局部不光滑点，并结合过程级贡献诊断给出可疑通道。

输入：
- scan CSV: scripts/relaxtime/scan_relaxation_times_vs_T.jl 输出的主扫描文件
- diagnostics CSV: 同一扫描可选输出的过程级贡献文件

输出：
- summary CSV: 每个 flagged kink 点一行
- channels CSV: 每个 flagged kink 点对应的 channel 变化明细
"""

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

struct Options
    scan_csv::String
    diagnostics_csv::String
    summary_out::String
    channels_out::String
    kink_threshold::Float64
    top_channels::Int
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/analyze_tau_smoothness.jl [options]\n")
    println("Options:")
    println("  --scan-csv <path>         主扫描 CSV")
    println("  --diagnostics-csv <path>  过程级诊断 CSV")
    println("  --summary-out <path>      kink 摘要输出 CSV")
    println("  --channels-out <path>     channel 明细输出 CSV")
    println("  --kink-threshold <float>  kink 阈值 (default 0.03)")
    println("  --top-channels <int>      摘要中保留的 top channel 数 (default 3)")
    println("  -h, --help                Show this help")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol,Any}(
        :scan_csv => joinpath("data", "outputs", "results", "relaxtime", "scan", "relaxation_times_vs_T_literature_compare.csv"),
        :diagnostics_csv => joinpath("data", "outputs", "results", "relaxtime", "scan", "relaxation_times_vs_T_process_diagnostics.csv"),
        :summary_out => joinpath("data", "outputs", "results", "relaxtime", "scan", "tau_smoothness_kinks_summary.csv"),
        :channels_out => joinpath("data", "outputs", "results", "relaxtime", "scan", "tau_smoothness_kinks_channels.csv"),
        :kink_threshold => 0.03,
        :top_channels => 3,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $(arg)")
            val = args[i + 1]
            i += 1
            return val
        end

        if arg == "--scan-csv"
            opts[:scan_csv] = require_value()
        elseif arg == "--diagnostics-csv"
            opts[:diagnostics_csv] = require_value()
        elseif arg == "--summary-out"
            opts[:summary_out] = require_value()
        elseif arg == "--channels-out"
            opts[:channels_out] = require_value()
        elseif arg == "--kink-threshold"
            opts[:kink_threshold] = parse(Float64, require_value())
        elseif arg == "--top-channels"
            opts[:top_channels] = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            print_usage(); exit(0)
        else
            error("unknown option: $(arg)")
        end
        i += 1
    end

    return Options(
        String(opts[:scan_csv]),
        String(opts[:diagnostics_csv]),
        String(opts[:summary_out]),
        String(opts[:channels_out]),
        Float64(opts[:kink_threshold]),
        Int(opts[:top_channels]),
    )
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function parse_csv_line(line::String)
    return split(chomp(line), ',')
end

function read_comment_csv(path::AbstractString)
    rows = Vector{Dict{String,String}}()
    open(path, "r") do io
        header = String[]
        for raw in eachline(io)
            line = strip(raw)
            isempty(line) && continue
            startswith(line, "#") && continue
            if isempty(header)
                header = parse_csv_line(raw)
                continue
            end
            values = parse_csv_line(raw)
            push!(rows, Dict(header[i] => values[i] for i in eachindex(header)))
        end
    end
    return rows
end

@inline function safe_float(row::Dict{String,String}, key::String)
    return parse(Float64, row[key])
end

@inline function kink_metric(prev::Float64, curr::Float64, next::Float64)
    denom = max(abs(curr), 1e-12)
    return abs(curr - 0.5 * (prev + next)) / denom
end

function main()
    opts = parse_args(ARGS)
    ensure_parent_dir(opts.summary_out)
    ensure_parent_dir(opts.channels_out)

    scan_rows = read_comment_csv(opts.scan_csv)
    diag_rows = read_comment_csv(opts.diagnostics_csv)

    scan_groups = Dict{Tuple{String,String}, Vector{Dict{String,String}}}()
    for row in scan_rows
        key = (row["muB_MeV"], row["xi"])
        push!(get!(scan_groups, key, Vector{Dict{String,String}}()), row)
    end
    for rows in values(scan_groups)
        sort!(rows, by=row -> safe_float(row, "T_MeV"))
    end

    diag_map = Dict{Tuple{String,String,String,String}, Float64}()
    for row in diag_rows
        key = (row["muB_MeV"], row["T_MeV"], row["species"], row["channel"])
        diag_map[key] = safe_float(row, "contribution")
    end

    species_fields = (
        (species="u", field="tau_u"),
        (species="s", field="tau_s"),
        (species="ubar", field="tau_ubar"),
        (species="sbar", field="tau_sbar"),
    )

    summary_io = open(opts.summary_out, "w")
    channels_io = open(opts.channels_out, "w")
    try
        println(summary_io, join([
            "muB_MeV", "xi", "species",
            "T_prev_MeV", "T_MeV", "T_next_MeV",
            "tau_prev", "tau_curr", "tau_next",
            "kink_metric", "top_channels",
        ], ','))
        println(channels_io, join([
            "muB_MeV", "xi", "species", "T_MeV", "channel",
            "current_contribution", "midpoint_contribution",
            "delta", "relative_delta",
        ], ','))

        flagged_count = 0
        for ((muB, xi), rows) in sort(collect(scan_groups); by=item -> parse(Float64, item[1][1]))
            length(rows) < 3 && continue
            for spec in species_fields
                for idx in 2:(length(rows) - 1)
                    prev = rows[idx - 1]
                    curr = rows[idx]
                    next = rows[idx + 1]
                    prev_tau = safe_float(prev, spec.field)
                    curr_tau = safe_float(curr, spec.field)
                    next_tau = safe_float(next, spec.field)
                    metric = kink_metric(prev_tau, curr_tau, next_tau)
                    metric < opts.kink_threshold && continue

                    channels = String[]
                    for row in diag_rows
                        row["species"] == spec.species || continue
                        row["muB_MeV"] == muB || continue
                        push!(channels, row["channel"])
                    end
                    unique!(channels)

                    channel_deltas = NamedTuple[]
                    for channel in channels
                        prev_key = (muB, prev["T_MeV"], spec.species, channel)
                        curr_key = (muB, curr["T_MeV"], spec.species, channel)
                        next_key = (muB, next["T_MeV"], spec.species, channel)
                        haskey(diag_map, curr_key) || continue
                        prev_val = get(diag_map, prev_key, NaN)
                        curr_val = diag_map[curr_key]
                        next_val = get(diag_map, next_key, NaN)
                        midpoint = (isfinite(prev_val) && isfinite(next_val)) ? 0.5 * (prev_val + next_val) : NaN
                        delta = isfinite(midpoint) ? curr_val - midpoint : NaN
                        rel = isfinite(delta) ? delta / max(abs(curr_val), 1e-12) : NaN
                        push!(channel_deltas, (channel=channel, current=curr_val, midpoint=midpoint, delta=delta, relative=rel))
                    end

                    sort!(channel_deltas, by=row -> isfinite(row.delta) ? -abs(row.delta) : Inf)
                    top = first(channel_deltas, min(opts.top_channels, length(channel_deltas)))
                    top_summary = join([@sprintf("%s:%.4e", row.channel, row.delta) for row in top if isfinite(row.delta)], ';')

                    println(summary_io, join([
                        muB, xi, spec.species,
                        prev["T_MeV"], curr["T_MeV"], next["T_MeV"],
                        string(prev_tau), string(curr_tau), string(next_tau),
                        string(metric), top_summary,
                    ], ','))

                    for row in channel_deltas
                        println(channels_io, join([
                            muB, xi, spec.species, curr["T_MeV"], row.channel,
                            string(row.current), string(row.midpoint),
                            string(row.delta), string(row.relative),
                        ], ','))
                    end

                    flagged_count += 1
                end
            end
        end

        @printf("Flagged kink points: %d\n", flagged_count)
        @printf("Saved summary: %s\n", opts.summary_out)
        @printf("Saved channel details: %s\n", opts.channels_out)
    finally
        close(summary_io)
        close(channels_io)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end