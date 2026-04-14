#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))

using .Constants_PNJL: ħc_MeV_fm
using Main.GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using Main.OneLoopIntegrals: A
using Main.PolarizationAniso: polarization_with_width

function _parse_args(args::Vector{String})
    opts = Dict{String,String}()
    i = 1
    while i <= length(args)
        a = args[i]
        if startswith(a, "--")
            i == length(args) && throw(ArgumentError("missing value for $a"))
            opts[a] = args[i + 1]
            i += 2
        else
            throw(ArgumentError("unknown arg: $a"))
        end
    end
    for k in ("--derived-csv", "--scan-csv", "--gap-csv", "--out-csv")
        haskey(opts, k) || throw(ArgumentError("missing required arg: $k"))
    end
    return (
        derived_csv=opts["--derived-csv"],
        scan_csv=opts["--scan-csv"],
        gap_csv=opts["--gap-csv"],
        out_csv=opts["--out-csv"],
    )
end

function _read_csv(path::String)
    lines = String[]
    open(path, "r") do io
        for raw in eachline(io)
            s = strip(raw)
            isempty(s) && continue
            startswith(s, "#") && continue
            push!(lines, s)
        end
    end
    isempty(lines) && throw(ArgumentError("empty csv: $path"))
    header = [strip(x) for x in split(lines[1], ',')]
    rows = Vector{Dict{String,String}}()
    for ln in lines[2:end]
        parts = split(ln, ',')
        length(parts) < length(header) && continue
        row = Dict{String,String}()
        for (i, c) in enumerate(header)
            row[c] = strip(parts[i])
        end
        push!(rows, row)
    end
    return rows
end

@inline _f(row::Dict{String,String}, key::String) = parse(Float64, row[key])

function _group_by_xi(rows::Vector{Dict{String,String}})
    d = Dict{Float64,Vector{Dict{String,String}}}()
    for r in rows
        xi = _f(r, "xi")
        push!(get!(d, xi, Dict{String,String}[]), r)
    end
    for (_, rr) in d
        sort!(rr, by=r -> _f(r, "T_MeV"))
    end
    return d
end

function _interp_t(rows::Vector{Dict{String,String}}, col::String, t::Float64)
    ts = [_f(r, "T_MeV") for r in rows]
    ys = [_f(r, col) for r in rows]
    if t <= ts[1]
        return ys[1]
    elseif t >= ts[end]
        return ys[end]
    end
    for i in 1:(length(ts)-1)
        t0, t1 = ts[i], ts[i+1]
        if t0 <= t <= t1
            y0, y1 = ys[i], ys[i+1]
            return y0 + (t - t0) * (y1 - y0) / (t1 - t0)
        end
    end
    return ys[end]
end

function _interp_gap(gap_by_xi::Dict{Float64,Vector{Dict{String,String}}}, xi::Float64, col::String, t::Float64)
    if haskey(gap_by_xi, xi)
        return _interp_t(gap_by_xi[xi], col, t)
    end
    xis = sort(collect(keys(gap_by_xi)))
    xi < xis[1] && return _interp_t(gap_by_xi[xis[1]], col, t)
    xi > xis[end] && return _interp_t(gap_by_xi[xis[end]], col, t)
    for i in 1:(length(xis)-1)
        x0, x1 = xis[i], xis[i+1]
        if x0 <= xi <= x1
            y0 = _interp_t(gap_by_xi[x0], col, t)
            y1 = _interp_t(gap_by_xi[x1], col, t)
            return y0 + (xi - x0) * (y1 - y0) / (x1 - x0)
        end
    end
    return _interp_t(gap_by_xi[xis[end]], col, t)
end

function _estimate_crossing(rows::Vector{Dict{String,String}}, mes_col::String, thr_col::String)
    ts = [_f(r, "T_MeV") for r in rows]
    ds = [_f(r, mes_col) - _f(r, thr_col) for r in rows]
    for i in 1:(length(ts)-1)
        d0, d1 = ds[i], ds[i+1]
        if d0 == 0.0
            return ts[i]
        end
        if d0 * d1 < 0.0
            t0, t1 = ts[i], ts[i+1]
            return t0 + (0.0 - d0) * (t1 - t0) / (d1 - d0)
        end
    end
    return ts[end]
end

function _compute_point(channel::Symbol, xi::Float64, t_mev::Float64,
        scan_rows::Vector{Dict{String,String}}, gap_by_xi::Dict{Float64,Vector{Dict{String,String}}})
    k0 = channel == :pi ? _interp_t(scan_rows, "M_pi", t_mev) : _interp_t(scan_rows, "M_K", t_mev)
    gamma = channel == :pi ? _interp_t(scan_rows, "Gamma_pi", t_mev) : _interp_t(scan_rows, "Gamma_K", t_mev)
    m_u = _interp_t(scan_rows, "m_u", t_mev)
    m_s = _interp_t(scan_rows, "m_s", t_mev)
    phi = _interp_gap(gap_by_xi, xi, "Phi", t_mev)
    phibar = _interp_gap(gap_by_xi, xi, "Phibar", t_mev)
    T_fm = t_mev / ħc_MeV_fm

    A_u = A(m_u, 0.0, T_fm, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(m_s, 0.0, T_fm, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)

    if channel == :pi
        full_re, full_im = polarization_with_width(:P, k0, gamma, 0.0, m_u, m_u, 0.0, 0.0, T_fm, phi, phibar, xi, A_u, A_u, 0)
        ind_re, ind_im = polarization_with_width(:P, k0, gamma, 0.0, m_u, m_u, 0.0, 0.0, T_fm, phi, phibar, 0.0, A_u, A_u, 0)
    else
        full_re, full_im = polarization_with_width(:P, k0, gamma, 0.0, m_u, m_s, 0.0, 0.0, T_fm, phi, phibar, xi, A_u, A_s, 1)
        ind_re, ind_im = polarization_with_width(:P, k0, gamma, 0.0, m_u, m_s, 0.0, 0.0, T_fm, phi, phibar, 0.0, A_u, A_s, 1)
    end
    return (
        re_full=full_re,
        im_full=full_im,
        re_indirect=ind_re,
        im_indirect=ind_im,
        re_direct=full_re - ind_re,
        im_direct=full_im - ind_im,
        t_probe=t_mev,
        k0=k0,
        gamma=gamma,
        phi=phi,
        phibar=phibar,
    )
end

@inline function _finite_decomp(pt)
    vals = (
        pt.re_full, pt.im_full,
        pt.re_indirect, pt.im_indirect,
        pt.re_direct, pt.im_direct,
        pt.k0, pt.gamma,
    )
    return all(isfinite, vals)
end

function _compute_point_best(channel::Symbol, xi::Float64, t_target::Float64,
        scan_rows::Vector{Dict{String,String}}, gap_by_xi::Dict{Float64,Vector{Dict{String,String}}})
    ts = unique([_f(r, "T_MeV") for r in scan_rows])
    sort!(ts, by=t -> abs(t - t_target))
    # 优先在目标附近寻找有限值，避免正好踩在阈值奇点。
    for t in ts
        pt = _compute_point(channel, xi, t, scan_rows, gap_by_xi)
        if _finite_decomp(pt)
            return pt
        end
    end
    # 若全部非有限，退回目标点结果（用于显式暴露问题，而不是静默吞掉）
    return _compute_point(channel, xi, t_target, scan_rows, gap_by_xi)
end

function main()
    opts = _parse_args(ARGS)
    derived = _read_csv(opts.derived_csv)
    scan = _read_csv(opts.scan_csv)
    gap = _read_csv(opts.gap_csv)

    derived_by_xi = _group_by_xi(derived)
    scan_by_xi = _group_by_xi(scan)
    gap_by_xi = _group_by_xi(gap)
    xis = sort(collect(intersect(keys(derived_by_xi), keys(scan_by_xi))))

    mkpath(dirname(opts.out_csv))
    open(opts.out_csv, "w") do io
        println(io, "channel,xi,t_probe_MeV,re_full,im_full,re_indirect,im_indirect,re_direct,im_direct,k0,gamma,Phi,Phibar")
        for xi in xis
            drows = derived_by_xi[xi]
            t_pi = _estimate_crossing(drows, "M_pi", "M_u_plus_M_d")
            t_k = _estimate_crossing(drows, "M_K", "M_u_plus_M_s")

            pi_pt = _compute_point_best(:pi, xi, t_pi, scan_by_xi[xi], gap_by_xi)
            k_pt = _compute_point_best(:K, xi, t_k, scan_by_xi[xi], gap_by_xi)

            println(io, join([
                "pi", string(xi), string(pi_pt.t_probe), string(pi_pt.re_full), string(pi_pt.im_full),
                string(pi_pt.re_indirect), string(pi_pt.im_indirect), string(pi_pt.re_direct), string(pi_pt.im_direct),
                string(pi_pt.k0), string(pi_pt.gamma), string(pi_pt.phi), string(pi_pt.phibar)
            ], ','))
            println(io, join([
                "K", string(xi), string(k_pt.t_probe), string(k_pt.re_full), string(k_pt.im_full),
                string(k_pt.re_indirect), string(k_pt.im_indirect), string(k_pt.re_direct), string(k_pt.im_direct),
                string(k_pt.k0), string(k_pt.gamma), string(k_pt.phi), string(k_pt.phibar)
            ], ','))
        end
    end

    println("Wrote polarization decomposition csv: " * opts.out_csv)
end

main()
