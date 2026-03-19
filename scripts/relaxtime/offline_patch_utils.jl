module OfflinePatchUtils

export PatchPoint, read_flagged_points, select_patch_points

struct PatchPoint
    T_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    quality_reason::String
    quality_metric::Float64
end

@inline _is_comment(line::AbstractString) = isempty(strip(line)) || startswith(strip(line), "#")

function _first_non_comment_line(io)
    for line in eachline(io)
        _is_comment(line) && continue
        return line
    end
    return nothing
end

function _parse_header_map(header_line::AbstractString)
    cols = split(strip(header_line), ',')
    return Dict{String, Int}(String(strip(c)) => i for (i, c) in enumerate(cols))
end

function _parse_bool_or_nothing(s::AbstractString)
    t = lowercase(strip(s))
    t in ("true", "1", "yes", "y") && return true
    t in ("false", "0", "no", "n") && return false
    return nothing
end

function _parse_float_or_nan(s::AbstractString)
    v = tryparse(Float64, strip(s))
    return v === nothing ? NaN : Float64(v)
end

function read_flagged_points(path::AbstractString; require_quality_flag::Bool=true)
    points = PatchPoint[]
    isfile(path) || error("input csv does not exist: $path")

    open(path, "r") do io
        header = _first_non_comment_line(io)
        header === nothing && return points

        colmap = _parse_header_map(header)
        for col in ("T_MeV", "muB_MeV", "xi")
            haskey(colmap, col) || error("missing required column: $col")
        end
        if require_quality_flag && !haskey(colmap, "quality_flag")
            error("missing required column: quality_flag")
        end

        idx_T = colmap["T_MeV"]
        idx_muB = colmap["muB_MeV"]
        idx_xi = colmap["xi"]
        idx_flag = get(colmap, "quality_flag", 0)
        idx_reason = get(colmap, "quality_reason", 0)
        idx_metric = get(colmap, "quality_metric", 0)

        for line in eachline(io)
            _is_comment(line) && continue
            parts = split(line, ',')
            length(parts) < max(idx_xi, idx_metric) && continue

            T_MeV = tryparse(Float64, strip(parts[idx_T]))
            muB_MeV = tryparse(Float64, strip(parts[idx_muB]))
            xi = tryparse(Float64, strip(parts[idx_xi]))
            (T_MeV === nothing || muB_MeV === nothing || xi === nothing) && continue

            if idx_flag > 0
                flag_raw = strip(parts[idx_flag])
                flag = _parse_bool_or_nothing(flag_raw)
                flag === nothing && continue
                flag || continue
            elseif require_quality_flag
                continue
            end

            reason = idx_reason > 0 && idx_reason <= length(parts) ? strip(parts[idx_reason]) : ""
            metric = idx_metric > 0 && idx_metric <= length(parts) ? _parse_float_or_nan(parts[idx_metric]) : NaN
            push!(points, PatchPoint(Float64(T_MeV), Float64(muB_MeV), Float64(xi), reason, metric))
        end
    end

    # Deduplicate by (T, muB, xi), keeping the row with larger |metric|.
    by_key = Dict{Tuple{Float64, Float64, Float64}, PatchPoint}()
    for p in points
        key = (p.T_MeV, p.muB_MeV, p.xi)
        if !haskey(by_key, key)
            by_key[key] = p
            continue
        end
        prev = by_key[key]
        prev_abs = isfinite(prev.quality_metric) ? abs(prev.quality_metric) : -Inf
        now_abs = isfinite(p.quality_metric) ? abs(p.quality_metric) : -Inf
        now_abs > prev_abs && (by_key[key] = p)
    end

    dedup = collect(values(by_key))
    sort!(dedup, by=p -> (p.T_MeV, p.muB_MeV, p.xi))
    return dedup
end

function select_patch_points(points::Vector{PatchPoint}; max_points::Int=0, reason_filter::Union{Nothing, String}=nothing)
    selected = points
    if reason_filter !== nothing
        wanted = Set{String}(filter(!isempty, strip.(split(reason_filter, ','))))
        selected = [p for p in selected if p.quality_reason in wanted]
    end

    if max_points > 0 && length(selected) > max_points
        sort!(selected, by=p -> (isfinite(p.quality_metric) ? abs(p.quality_metric) : -Inf), rev=true)
        return selected[1:max_points]
    end

    return selected
end

end
