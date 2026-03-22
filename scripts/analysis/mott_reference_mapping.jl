module MottReferenceMapping

export load_reference_table, validate_reference_schema
export estimate_mott_temperature, refine_mott_temperature_bisection
export supported_theoretical_mesons, collect_mesons_from_rows

const REQUIRED_COLUMNS = (
    :record_id,
    :source_impl,
    :T_MeV,
    :muB_MeV,
    :xi,
    :meson,
    :mass_MeV,
    :threshold_MeV,
    :gap_MeV,
    :mott_flag,
    :solver_status,
)

const THEORETICAL_MESONS = (
    :pi,
    :K,
    :eta,
    :eta_prime,
    :sigma_pi,
    :sigma_K,
    :sigma,
    :sigma_prime,
)

supported_theoretical_mesons() = collect(THEORETICAL_MESONS)

function collect_mesons_from_rows(rows::Vector{Dict{Symbol,Any}})
    out = Symbol[]
    seen = Set{Symbol}()
    for row in rows
        haskey(row, :meson) || continue
        m = Symbol(row[:meson])
        m in seen && continue
        push!(out, m)
        push!(seen, m)
    end
    sort!(out; by=String)
    return out
end

function _parse_header(header_line::String)
    cols = split(strip(header_line), ',')
    return Symbol.(strip.(cols))
end

function _parse_value(key::Symbol, raw::AbstractString)
    v = strip(String(raw))
    if key in (:record_id, :source_impl, :meson, :solver_status)
        return v
    elseif key == :mott_flag
        return parse(Int, v)
    else
        return parse(Float64, v)
    end
end

function load_reference_table(path::String)
    isfile(path) || throw(ArgumentError("reference CSV not found: $path"))
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("reference CSV is empty: $path"))

    header = _parse_header(lines[1])
    rows = Dict{Symbol,Any}[]

    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, '#') && continue
        cols = split(s, ',')
        length(cols) == length(header) || throw(ArgumentError("invalid row column count: $line"))
        row = Dict{Symbol,Any}()
        for (idx, key) in enumerate(header)
            row[key] = _parse_value(key, cols[idx])
        end
        push!(rows, row)
    end

    return rows
end

function validate_reference_schema(rows::Vector{Dict{Symbol,Any}})
    isempty(rows) && throw(ArgumentError("reference rows are empty"))
    for (i, row) in enumerate(rows)
        for key in REQUIRED_COLUMNS
            haskey(row, key) || throw(ArgumentError("row $i missing required key: $key"))
        end
    end
    return rows
end

function _sorted_pairs(Ts::Vector{Float64}, gaps::Vector{Float64})
    length(Ts) == length(gaps) || throw(ArgumentError("Ts and gaps length mismatch"))
    length(Ts) >= 2 || throw(ArgumentError("need at least 2 points to estimate mott temperature"))
    pairs = sort(collect(zip(Ts, gaps)); by=x -> x[1])
    return pairs
end

function estimate_mott_temperature(Ts::Vector{Float64}, gaps::Vector{Float64})
    pairs = _sorted_pairs(Ts, gaps)

    for i in 1:(length(pairs)-1)
        t1, g1 = pairs[i]
        t2, g2 = pairs[i + 1]
        if g1 == 0.0
            return (T_mott_MeV=t1, method=:linear, approx=false)
        end
        if g2 == 0.0
            return (T_mott_MeV=t2, method=:linear, approx=false)
        end
        if signbit(g1) != signbit(g2)
            t = t1 - g1 * (t2 - t1) / (g2 - g1)
            return (T_mott_MeV=t, method=:linear, approx=false)
        end
    end

    min_idx = argmin(abs.([p[2] for p in pairs]))
    return (T_mott_MeV=pairs[min_idx][1], method=:minabs_approx, approx=true)
end

function refine_mott_temperature_bisection(gap_func, T_lo::Float64, T_hi::Float64;
                                           tol::Float64=1e-4, max_iter::Int=50)
    T_lo < T_hi || throw(ArgumentError("invalid bracket: T_lo must be < T_hi"))
    g_lo = gap_func(T_lo)
    g_hi = gap_func(T_hi)
    if g_lo == 0.0
        return (T_mott_MeV=T_lo, converged=true, iterations=0, method=:bisection)
    end
    if g_hi == 0.0
        return (T_mott_MeV=T_hi, converged=true, iterations=0, method=:bisection)
    end
    signbit(g_lo) != signbit(g_hi) || throw(ArgumentError("bisection requires sign change bracket"))

    lo = T_lo
    hi = T_hi
    gleft = g_lo

    for iter in 1:max_iter
        mid = 0.5 * (lo + hi)
        gmid = gap_func(mid)
        if abs(gmid) <= tol
            return (T_mott_MeV=mid, converged=true, iterations=iter, method=:bisection)
        end
        if signbit(gleft) != signbit(gmid)
            hi = mid
        else
            lo = mid
            gleft = gmid
        end
    end

    mid = 0.5 * (lo + hi)
    return (T_mott_MeV=mid, converged=false, iterations=max_iter, method=:bisection)
end

end
