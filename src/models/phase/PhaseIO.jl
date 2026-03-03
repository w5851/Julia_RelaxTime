function group_curves_by_temperature(rows; xi::Real=0.0, tol::Real=1e-6)
    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    for row in rows
        T = try parse(Float64, row["T_MeV"]) catch; continue end
        xi_val = haskey(row, "xi") ? (try parse(Float64, row["xi"]) catch; NaN end) : NaN
        (isnan(xi_val) || abs(xi_val - xi) > tol) && continue
        rho = try parse(Float64, row["rho"]) catch; continue end
        mu = haskey(row, "mu_avg_MeV") ? (try parse(Float64, row["mu_avg_MeV"]) catch; NaN end) : NaN
        isnan(mu) && (mu = haskey(row, "mu_MeV") ? (try parse(Float64, row["mu_MeV"]) catch; NaN end) : NaN)
        (isnan(mu) || !isfinite(rho)) && continue
        bucket = get!(grouped, T) do
            Vector{Tuple{Float64, Float64}}()
        end
        push!(bucket, (mu, rho))
    end
    return grouped
end

function load_curves_from_trho_csv(path::String; xi::Real=0.0, tol::Real=1e-6, min_points::Int=3)
    curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    isfile(path) || return curves

    lines = readlines(path)
    isempty(lines) && return curves

    header = split(lines[1], ',')
    idx_T = findfirst(==("T_MeV"), header)
    idx_rho = findfirst(==("rho"), header)
    idx_xi = findfirst(==("xi"), header)
    idx_mu_avg = findfirst(==("mu_avg_MeV"), header)
    idx_mu = idx_mu_avg === nothing ? findfirst(==("mu_MeV"), header) : idx_mu_avg
    idx_conv = findfirst(==("converged"), header)

    if idx_T === nothing || idx_rho === nothing || idx_mu === nothing
        return curves
    end

    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()

    for line in lines[2:end]
        isempty(strip(line)) && continue
        cols = split(line, ',')
        length(cols) < max(idx_T, idx_rho, idx_mu) && continue

        T = tryparse(Float64, cols[idx_T])
        rho = tryparse(Float64, cols[idx_rho])
        mu = tryparse(Float64, cols[idx_mu])
        (T === nothing || rho === nothing || mu === nothing) && continue

        if idx_xi !== nothing && idx_xi <= length(cols)
            xi_val = tryparse(Float64, cols[idx_xi])
            (xi_val === nothing || abs(xi_val - xi) > tol) && continue
        end

        if idx_conv !== nothing && idx_conv <= length(cols)
            conv_raw = lowercase(strip(cols[idx_conv]))
            conv = conv_raw in ("true", "1", "yes")
            conv || continue
        end

        bucket = get!(grouped, Float64(T)) do
            Vector{Tuple{Float64, Float64}}()
        end
        push!(bucket, (Float64(mu), Float64(rho)))
    end

    for (T, pts) in grouped
        length(pts) < min_points && continue
        mu_vals = Float64[first(p) for p in pts]
        rho_vals = Float64[last(p) for p in pts]
        curves[T] = (mu_vals, rho_vals)
    end

    return curves
end
