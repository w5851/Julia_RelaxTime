"""
Utility helpers for a lightweight scan CSV format with metadata header.

Format (scan_csv_v1):
- Any line starting with `#` is treated as metadata/comment and ignored by CSV readers.
- Metadata lines should use `# key: value` (recommended) or `# key=value`.
- The first non-comment line is the CSV header.

This intentionally avoids external deps (no CSV.jl) and assumes numeric/identifier
fields without embedded commas.
"""

module ScanCSV

export write_metadata, write_header, read_existing_keys

"""Write metadata lines as `# key: value` (one per line)."""
function write_metadata(io, metadata::AbstractDict{<:AbstractString,<:AbstractString})
    for (k, v) in metadata
        println(io, "# ", k, ": ", v)
    end
end

"""Write a single CSV header line from a vector of column names."""
function write_header(io, cols::AbstractVector{<:AbstractString})
    println(io, join(cols, ','))
end

function _is_comment(line::AbstractString)
    s = strip(line)
    return isempty(s) || startswith(s, "#")
end

function _first_non_comment_line(io)
    for line in eachline(io)
        _is_comment(line) && continue
        return line
    end
    return nothing
end

"""Parse header into a mapping colname=>index (1-based)."""
function _parse_header_map(header_line::AbstractString)
    cols = split(strip(header_line), ',')
    return Dict{String,Int}(String(strip(c)) => i for (i, c) in enumerate(cols))
end

"""Read existing keys from a scan CSV file.

Keys are tuples of Float64 values extracted from the specified `key_cols`.
Skips metadata/comment lines.
"""
function read_existing_keys(path::AbstractString, key_cols::AbstractVector{<:AbstractString})
    keys = Set{Tuple{Vararg{Float64}}}()
    isfile(path) || return keys

    open(path, "r") do io
        header = _first_non_comment_line(io)
        header === nothing && return keys
        colmap = _parse_header_map(header)

        idxs = Int[]
        for c in key_cols
            haskey(colmap, String(c)) || return keys
            push!(idxs, colmap[String(c)])
        end

        for line in eachline(io)
            _is_comment(line) && continue
            parts = split(line, ',')
            length(parts) < maximum(idxs) && continue

            vals = Float64[]
            ok = true
            for idx in idxs
                v = tryparse(Float64, strip(parts[idx]))
                if v === nothing
                    ok = false
                    break
                end
                push!(vals, v)
            end
            ok || continue

            push!(keys, Tuple(vals))
        end
    end

    return keys
end

end # module
