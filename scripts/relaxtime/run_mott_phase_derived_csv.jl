#!/usr/bin/env julia

function _print_usage()
    println("Usage: julia --project=. scripts/relaxtime/run_mott_phase_derived_csv.jl --in <scan_csv> --out <derived_csv>")
end

function _parse_args(args::Vector{String})
    input = nothing
    output = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $arg"))
            i += 1
            return args[i]
        end

        if arg == "--in"
            input = require_value()
        elseif arg == "--out"
            output = require_value()
        elseif arg in ("-h", "--help")
            _print_usage()
            exit(0)
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end

    input === nothing && throw(ArgumentError("--in is required"))
    output === nothing && throw(ArgumentError("--out is required"))
    return String(input), String(output)
end

@inline function _safe_parse_float(s::AbstractString)
    x = tryparse(Float64, strip(s))
    return x === nothing ? NaN : x
end

function _derived_vals(row::Dict{String,String})
    mu = _safe_parse_float(get(row, "m_u", "NaN"))
    md = _safe_parse_float(get(row, "m_d", "NaN"))
    ms = _safe_parse_float(get(row, "m_s", "NaN"))
    return mu + md, mu + ms
end

function _col_exists(header::AbstractVector{<:AbstractString}, name::String)
    return any(==(name), header)
end

function _process(input_csv::String, output_csv::String)
    lines = readlines(input_csv)

    metadata = String[]
    header = String[]
    data_rows = Vector{Dict{String,String}}()
    header_seen = false

    for line in lines
        s = strip(line)
        if isempty(s)
            continue
        elseif startswith(s, "#") && !header_seen
            push!(metadata, line)
            continue
        elseif !header_seen
            header = [strip(x) for x in split(line, ',')]
            header_seen = true
            continue
        end

        vals = split(line, ',')
        row = Dict{String,String}()
        for (idx, key) in enumerate(header)
            row[key] = idx <= length(vals) ? strip(vals[idx]) : ""
        end
        push!(data_rows, row)
    end

    isempty(header) && throw(ArgumentError("input csv has no header"))

    out_header = copy(header)
    if !("M_u_plus_M_d" in out_header)
        push!(out_header, "M_u_plus_M_d")
    end
    if !("M_u_plus_M_s" in out_header)
        push!(out_header, "M_u_plus_M_s")
    end

    if _col_exists(header, "M_eta")
        if !("M_eta_plus_M_eta_prime" in out_header)
            push!(out_header, "M_eta_plus_M_eta_prime")
        end
    end

    mkpath(dirname(output_csv))
    open(output_csv, "w") do io
        for m in metadata
            println(io, m)
        end
        println(io, "# y_unit.M_u_plus_M_d: fm^-1")
        println(io, "# y_unit.M_u_plus_M_s: fm^-1")
        println(io, join(out_header, ','))

        for row in data_rows
            mud, mus = _derived_vals(row)
            row["M_u_plus_M_d"] = string(mud)
            row["M_u_plus_M_s"] = string(mus)

            if haskey(row, "M_eta") && haskey(row, "M_eta_prime")
                meta = _safe_parse_float(get(row, "M_eta", "NaN"))
                metap = _safe_parse_float(get(row, "M_eta_prime", "NaN"))
                row["M_eta_plus_M_eta_prime"] = string(meta + metap)
            end

            vals = [get(row, c, "") for c in out_header]
            println(io, join(vals, ','))
        end
    end
end

function main()
    input_csv, output_csv = _parse_args(ARGS)
    _process(input_csv, output_csv)
    println("Wrote derived CSV: ", output_csv)
end

main()
