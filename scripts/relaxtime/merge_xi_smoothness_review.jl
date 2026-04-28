#!/usr/bin/env julia

const ALLOWED_MANUAL_LABELS = Set(["confirm_smooth", "confirm_not_smooth", "undecided"])

struct Options
    flags::String
    manual_review::String
    out::String
end

@inline _is_comment_or_empty(line::AbstractString) = isempty(strip(line)) || startswith(strip(line), "#")

function _parse_csv_line(raw::AbstractString)
    vals = String[]
    buf = IOBuffer()
    in_quotes = false
    i = firstindex(raw)
    while i <= lastindex(raw)
        c = raw[i]
        if c == '"'
            if in_quotes
                ni = nextind(raw, i)
                if ni <= lastindex(raw) && raw[ni] == '"'
                    write(buf, '"')
                    i = ni
                else
                    in_quotes = false
                end
            else
                in_quotes = true
            end
        elseif c == ',' && !in_quotes
            push!(vals, String(take!(buf)))
        else
            write(buf, c)
        end
        i = nextind(raw, i)
    end
    in_quotes && throw(ArgumentError("unterminated quoted field in csv line"))
    push!(vals, String(take!(buf)))
    return vals
end

function _csv_escape(x::AbstractString)
    s = String(x)
    if occursin('"', s) || occursin(',', s) || occursin('\n', s) || occursin('\r', s)
        return '"' * replace(s, '"' => "\"\"") * '"'
    end
    return s
end

function print_usage()
    println("Usage: julia --project=. scripts/relaxtime/merge_xi_smoothness_review.jl --flags <smoothness_flags.csv> --manual-review <manual_review.csv> --out <smoothness_final.csv>")
    println("Options:")
    println("  --flags <path>          自动标注 flags csv (required)")
    println("  --manual-review <path>  人工复核 csv (required)")
    println("  --out <path>            合并输出 csv (required)")
    println("  -h, --help              显示帮助")
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Union{Nothing, String}}(
        :flags => nothing,
        :manual_review => nothing,
        :out => nothing,
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $(arg)"))
            i += 1
            return args[i]
        end

        if arg == "--flags"
            opts[:flags] = require_value()
        elseif arg == "--manual-review"
            opts[:manual_review] = require_value()
        elseif arg == "--out"
            opts[:out] = require_value()
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $(arg)"))
        end
        i += 1
    end

    opts[:flags] === nothing && throw(ArgumentError("--flags is required"))
    opts[:manual_review] === nothing && throw(ArgumentError("--manual-review is required"))
    opts[:out] === nothing && throw(ArgumentError("--out is required"))

    flags = normpath(abspath(opts[:flags]))
    manual_review = normpath(abspath(opts[:manual_review]))
    out = normpath(abspath(opts[:out]))

    isfile(flags) || throw(ArgumentError("flags csv not found: $(flags)"))
    isfile(manual_review) || throw(ArgumentError("manual review csv not found: $(manual_review)"))

    return Options(flags, manual_review, out)
end

function _parse_csv(path::String)
    header = String[]
    rows = Vector{Dict{String, String}}()
    open(path, "r") do io
        header_seen = false
        for (line_no, raw) in enumerate(eachline(io))
            _is_comment_or_empty(raw) && continue
            if !header_seen
                header = [strip(x) for x in _parse_csv_line(raw)]
                header_seen = true
                continue
            end
            vals = [strip(x) for x in _parse_csv_line(raw)]
            row = Dict{String, String}()
            for (i, key) in enumerate(header)
                row[key] = i <= length(vals) ? vals[i] : ""
            end
            row["__line__"] = string(line_no)
            push!(rows, row)
        end
    end
    isempty(header) && throw(ArgumentError("$(path): missing csv header"))
    return header, rows
end

function _require_columns(path::String, header::AbstractVector{<:AbstractString}, required::AbstractVector{<:AbstractString})
    for col in required
        col in header || throw(ArgumentError("$(path): missing required column: $(col)"))
    end
end

function _row_key(path::String, row::Dict{String, String})
    sample_id = get(row, "sample_id", "")
    field = get(row, "field", "")
    isempty(sample_id) && throw(ArgumentError("$(path):$(row["__line__"]): sample_id must be non-empty"))
    isempty(field) && throw(ArgumentError("$(path):$(row["__line__"]): field must be non-empty"))
    return sample_id, field
end

function _read_flags(path::String)
    header, rows = _parse_csv(path)
    _require_columns(path, header, ["sample_id", "field", "label", "reason"])

    out = Dict{Tuple{String, String}, NamedTuple{(:sample_id, :field, :auto_label, :reason), Tuple{String, String, String, String}}}()
    for row in rows
        key = _row_key(path, row)
        haskey(out, key) && throw(ArgumentError("$(path):$(row["__line__"]): duplicate sample_id+field key: $(key[1]),$(key[2])"))
        auto_label = get(row, "label", "")
        isempty(auto_label) && throw(ArgumentError("$(path):$(row["__line__"]): label must be non-empty"))
        reason = get(row, "reason", "")
        out[key] = (
            sample_id=key[1],
            field=key[2],
            auto_label=auto_label,
            reason=reason,
        )
    end
    return out
end

function _read_manual_review(path::String)
    header, rows = _parse_csv(path)
    _require_columns(path, header, ["sample_id", "field", "manual_label"])

    out = Dict{Tuple{String, String}, NamedTuple{(:manual_label, :reason), Tuple{String, String}}}()
    for row in rows
        key = _row_key(path, row)
        haskey(out, key) && throw(ArgumentError("$(path):$(row["__line__"]): duplicate sample_id+field key: $(key[1]),$(key[2])"))
        manual_label = get(row, "manual_label", "")
        manual_label in ALLOWED_MANUAL_LABELS || throw(ArgumentError("$(path):$(row["__line__"]): invalid manual_label=$(repr(manual_label)); allowed: confirm_smooth, confirm_not_smooth, undecided"))
        reason = get(row, "reason", "")
        out[key] = (manual_label=manual_label, reason=reason)
    end
    return out
end

function _merge_rows(flags_rows, manual_rows)
    final_rows = NamedTuple[]
    ordered_keys = sort!(collect(keys(flags_rows)); by=x -> (x[1], x[2]))

    for key in ordered_keys
        flag_row = flags_rows[key]
        if haskey(manual_rows, key)
            manual = manual_rows[key]
            final_label = manual.manual_label
            reason = isempty(manual.reason) ? flag_row.reason : manual.reason
            push!(final_rows, (
                sample_id=flag_row.sample_id,
                field=flag_row.field,
                auto_label=flag_row.auto_label,
                manual_label=manual.manual_label,
                final_label=final_label,
                reason=reason,
            ))
        else
            push!(final_rows, (
                sample_id=flag_row.sample_id,
                field=flag_row.field,
                auto_label=flag_row.auto_label,
                manual_label="",
                final_label=flag_row.auto_label,
                reason=flag_row.reason,
            ))
        end
    end

    for key in keys(manual_rows)
        haskey(flags_rows, key) || throw(ArgumentError("manual review row not found in flags csv: $(key[1]),$(key[2])"))
    end

    return final_rows
end

function _write_final_csv(path::String, rows::Vector{<:NamedTuple})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "sample_id,field,auto_label,manual_label,final_label,reason")
        for r in rows
            println(io, join((_csv_escape(r.sample_id), _csv_escape(r.field), _csv_escape(r.auto_label), _csv_escape(r.manual_label), _csv_escape(r.final_label), _csv_escape(r.reason)), ','))
        end
    end
end

function main(args::Vector{String}=copy(ARGS))
    opts = parse_args(args)
    flags_rows = _read_flags(opts.flags)
    manual_rows = _read_manual_review(opts.manual_review)
    merged = _merge_rows(flags_rows, manual_rows)
    _write_final_csv(opts.out, merged)
    println("merged xi smoothness review rows: $(length(merged))")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
