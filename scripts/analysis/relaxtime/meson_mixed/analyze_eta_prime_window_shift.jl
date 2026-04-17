#!/usr/bin/env julia

function _parse_args(args::Vector{String})
    csv = nothing
    jump_threshold = 0.25
    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && throw(ArgumentError("missing value for $arg"))
            i += 1
            return args[i]
        end
        if arg == "--csv"
            csv = require_value()
        elseif arg == "--jump-threshold"
            jump_threshold = parse(Float64, require_value())
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end
    csv === nothing && throw(ArgumentError("--csv is required"))
    return String(csv), jump_threshold
end

function _load_rows(path::String)
    lines = readlines(path)
    header = String[]
    rows = Vector{Dict{String,Any}}()
    need = [
        "T_MeV", "xi",
        "M_eta_prime", "Gamma_eta_prime",
        "root_quality_eta_prime", "governance_selection_reason_eta_prime",
    ]
    for ln in lines
        s = strip(ln)
        isempty(s) && continue
        startswith(s, "#") && continue
        if isempty(header)
            header = [strip(x) for x in split(s, ',')]
            continue
        end
        vals = split(s, ',')
        d = Dict{String,Any}()
        ok = true
        for k in need
            idx = findfirst(==(k), header)
            idx === nothing && (ok = false; break)
            idx > length(vals) && (ok = false; break)
            v = strip(vals[idx])
            x = tryparse(Float64, v)
            d[k] = (x === nothing ? v : x)
        end
        ok || continue
        push!(rows, d)
    end
    return rows
end

function _group_by_xi(rows)
    g = Dict{Float64,Vector{Dict{String,Any}}}()
    for r in rows
        xi = Float64(r["xi"])
        push!(get!(g, xi, Dict{String,Any}[]), r)
    end
    for rr in values(g)
        sort!(rr; by=r -> Float64(r["T_MeV"]))
    end
    return g
end

function _jump_events(rr::Vector{Dict{String,Any}}, obs::String, thr::Float64)
    ev = NamedTuple[]
    for i in 2:length(rr)
        prev = Float64(rr[i - 1][obs])
        curr = Float64(rr[i][obs])
        dj = abs(curr - prev)
        if dj >= thr
            push!(ev, (
                T=Float64(rr[i]["T_MeV"]),
                jump=dj,
                prev=prev,
                curr=curr,
                quality=String(rr[i]["root_quality_eta_prime"]),
                gov=String(rr[i]["governance_selection_reason_eta_prime"]),
            ))
        end
    end
    return ev
end

function _estimate_highT_window_center(rr::Vector{Dict{String,Any}}, thr::Float64)
    ev = _jump_events(rr, "Gamma_eta_prime", thr)
    high = filter(e -> e.T >= 220.0, ev)
    isempty(high) && return NaN
    # weighted center by jump amplitude
    wsum = sum(e.jump for e in high)
    wsum <= 0 && return NaN
    return sum(e.T * e.jump for e in high) / wsum
end

function main()
    csv, thr = _parse_args(ARGS)
    rows = _load_rows(csv)
    grouped = _group_by_xi(rows)

    println("row_count=", length(rows))
    println("jump_threshold=", thr)
    println("xi,event_count_M,event_count_G,max_jump_M,max_jump_G,highT_center_G")

    xis = sort(collect(keys(grouped)))
    centers = NamedTuple[]
    for xi in xis
        rr = grouped[xi]
        em = _jump_events(rr, "M_eta_prime", thr)
        eg = _jump_events(rr, "Gamma_eta_prime", thr)
        max_m = isempty(em) ? 0.0 : maximum(e.jump for e in em)
        max_g = isempty(eg) ? 0.0 : maximum(e.jump for e in eg)
        center = _estimate_highT_window_center(rr, thr)
        push!(centers, (xi=xi, center=center))
        println(string(xi, ",", length(em), ",", length(eg), ",", max_m, ",", max_g, ",", center))
    end

    valid = filter(c -> isfinite(c.center), centers)
    println("highT_window_centers=")
    for c in valid
        println(string(c.xi, ",", c.center))
    end

    if length(valid) >= 2
        # rough monotonic trend check: center increases with xi?
        sorted = sort(valid; by=x -> x.xi)
        nondecreasing = all(sorted[i].center <= sorted[i + 1].center + 1e-12 for i in 1:length(sorted)-1)
        nonincreasing = all(sorted[i].center >= sorted[i + 1].center - 1e-12 for i in 1:length(sorted)-1)
        trend = nondecreasing ? "center_increases_with_xi" : (nonincreasing ? "center_decreases_with_xi" : "non_monotonic")
        println("trend=", trend)
    else
        println("trend=insufficient_data")
    end
end

main()
