#!/usr/bin/env julia

function _parse_args(args::Vector{String})
    csv_path = nothing
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
            csv_path = require_value()
        elseif arg == "--jump-threshold"
            jump_threshold = parse(Float64, require_value())
        else
            throw(ArgumentError("unknown option: $arg"))
        end
        i += 1
    end
    csv_path === nothing && throw(ArgumentError("--csv is required"))
    return String(csv_path), jump_threshold
end

function _load_rows(csv_path::String)
    lines = readlines(csv_path)
    header = String[]
    rows = Vector{Dict{String,Float64}}()
    needed = Set([
        "T_MeV", "xi",
        "M_eta", "M_eta_prime", "M_sigma", "M_sigma_prime",
        "Gamma_eta", "Gamma_eta_prime", "Gamma_sigma", "Gamma_sigma_prime",
    ])
    for ln in lines
        s = strip(ln)
        isempty(s) && continue
        startswith(s, "#") && continue
        if isempty(header)
            header = split(s, ',')
            continue
        end
        cols = split(s, ',')
        length(cols) == length(header) || continue
        r = Dict{String,Float64}()
        ok = true
        for (k, v) in zip(header, cols)
            k in needed || continue
            x = tryparse(Float64, strip(v))
            if x === nothing
                ok = false
                break
            end
            r[k] = x
        end
        ok || continue
        all(haskey(r, k) for k in needed) || continue
        push!(rows, r)
    end
    return rows
end

function _by_xi(rows)
    out = Dict{Float64,Vector{Dict{String,Float64}}}()
    for r in rows
        xi = r["xi"]
        get!(out, xi, Dict{String,Float64}[])
        push!(out[xi], r)
    end
    for v in values(out)
        sort!(v; by=r -> r["T_MeV"])
    end
    return out
end

function _analyze(rows, jump_threshold)
    obs = (
        "M_eta",
        "M_eta_prime",
        "M_sigma",
        "M_sigma_prime",
        "Gamma_eta",
        "Gamma_eta_prime",
        "Gamma_sigma",
        "Gamma_sigma_prime",
    )

    grouped = _by_xi(rows)
    println("xi,observable,T_MeV,abs_jump")
    for (xi, rr) in sort!(collect(grouped); by=first)
        for o in obs
            for i in 2:length(rr)
                dj = abs(rr[i][o] - rr[i - 1][o])
                if dj >= jump_threshold
                    println(string(xi, ",", o, ",", rr[i]["T_MeV"], ",", dj))
                end
            end
        end
    end
end

function main()
    csv_path, jump_threshold = _parse_args(ARGS)
    rows = _load_rows(csv_path)
    println("row_count=", length(rows))
    _analyze(rows, jump_threshold)
end

main()
