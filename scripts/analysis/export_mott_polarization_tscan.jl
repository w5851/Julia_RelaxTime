#!/usr/bin/env julia
# Export polarization function ReΠ and ImΠ across a temperature scan,
# with direct/indirect decomposition, for auxiliary figure generation.
#
# Usage:
#   julia --project=. scripts/analysis/export_mott_polarization_tscan.jl \
#       --derived-csv data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/mott_phase_derived.csv \
#       --gap-csv data/outputs/results/relaxtime/scan/gap_transport_scan_step5_muB0_xi-0p3to0p3.csv \
#       --out-csv data/outputs/results/relaxtime/mott_phase/reference_100_300_fine/analysis_exports/polarization_tscan.csv

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))

using .Constants_PNJL: ħc_MeV_fm
using Main.GaussLegendre: gauleg
using Main.OneLoopIntegralsCorrection: A_aniso
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
    for k in ("--derived-csv", "--gap-csv", "--out-csv")
        haskey(opts, k) || throw(ArgumentError("missing required arg: $k"))
    end
    return opts
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
        for (ic, c) in enumerate(header)
            row[c] = strip(parts[ic])
        end
        push!(rows, row)
    end
    return rows
end

_f(row, key) = parse(Float64, row[key])

function _group_by_xi(rows)
    d = Dict{Float64,Vector{Dict{String,String}}}()
    for r in rows
        xi = _f(r, "xi")
        push!(get!(d, xi, Dict{String,String}[]), r)
    end
    return d
end

function _interp_t(rows, col, t)
    ts = [_f(r, "T_MeV") for r in rows]
    ys = [_f(r, col) for r in rows]
    if t <= ts[1]; return ys[1]; end
    if t >= ts[end]; return ys[end]; end
    for i in 1:(length(ts)-1)
        if ts[i] <= t <= ts[i+1]
            return ys[i] + (t - ts[i]) * (ys[i+1] - ys[i]) / (ts[i+1] - ts[i])
        end
    end
    return ys[end]
end

function _interp_gap(gap_by_xi, xi, col, t)
    if haskey(gap_by_xi, xi)
        return _interp_t(gap_by_xi[xi], col, t)
    end
    xis = sort(collect(keys(gap_by_xi)))
    xi < xis[1] && return _interp_t(gap_by_xi[xis[1]], col, t)
    xi > xis[end] && return _interp_t(gap_by_xi[xis[end]], col, t)
    for i in 1:(length(xis)-1)
        if xis[i] <= xi <= xis[i+1]
            y0 = _interp_t(gap_by_xi[xis[i]], col, t)
            y1 = _interp_t(gap_by_xi[xis[i+1]], col, t)
            return y0 + (xi - xis[i]) * (y1 - y0) / (xis[i+1] - xis[i])
        end
    end
    return _interp_t(gap_by_xi[xis[end]], col, t)
end

# Grid parameters matching solver's build_A_triplet defaults:
# p_nodes=16, p_max=20.0, cos_nodes=8
const A_P_NODES = 16
const A_P_MAX   = 20.0
const A_COS_N   = 8
const A_nodes_p, A_weights_p = gauleg(0.0, A_P_MAX, A_P_NODES)
const A_nodes_cos, A_weights_cos = gauleg(-1.0, 1.0, A_COS_N)

function _compute_point(channel::Symbol, xi, t_mev, scan_rows, gap_by_xi)
    k0 = channel == :pi ? _interp_t(scan_rows, "M_pi", t_mev) : _interp_t(scan_rows, "M_K", t_mev)
    gamma = channel == :pi ? _interp_t(scan_rows, "Gamma_pi", t_mev) : _interp_t(scan_rows, "Gamma_K", t_mev)
    m_u = _interp_t(scan_rows, "m_u", t_mev)
    m_s = _interp_t(scan_rows, "m_s", t_mev)
    phi = _interp_gap(gap_by_xi, xi, "Phi", t_mev)
    phibar = _interp_gap(gap_by_xi, xi, "Phibar", t_mev)
    T_fm = t_mev / ħc_MeV_fm

    # Use A_aniso (matching solver's build_A_triplet) with consistent grid
    A_u = A_aniso(m_u, 0.0, T_fm, phi, phibar, xi, A_nodes_p, A_weights_p, A_nodes_cos, A_weights_cos)
    A_s = A_aniso(m_s, 0.0, T_fm, phi, phibar, xi, A_nodes_p, A_weights_p, A_nodes_cos, A_weights_cos)

    # Indirect A: isotropic distribution at the ξ-affected (m, phi) state
    A_u_ind = A_aniso(m_u, 0.0, T_fm, phi, phibar, 0.0, A_nodes_p, A_weights_p, A_nodes_cos, A_weights_cos)
    A_s_ind = A_aniso(m_s, 0.0, T_fm, phi, phibar, 0.0, A_nodes_p, A_weights_p, A_nodes_cos, A_weights_cos)

    if channel == :pi
        full_re, full_im = polarization_with_width(:P, k0, gamma, 0.0, m_u, m_u, 0.0, 0.0, T_fm, phi, phibar, xi, A_u, A_u, 0)
        ind_re, ind_im  = polarization_with_width(:P, k0, gamma, 0.0, m_u, m_u, 0.0, 0.0, T_fm, phi, phibar, 0.0, A_u_ind, A_u_ind, 0)
    else
        full_re, full_im = polarization_with_width(:P, k0, gamma, 0.0, m_u, m_s, 0.0, 0.0, T_fm, phi, phibar, xi, A_u, A_s, 1)
        ind_re, ind_im  = polarization_with_width(:P, k0, gamma, 0.0, m_u, m_s, 0.0, 0.0, T_fm, phi, phibar, 0.0, A_u_ind, A_s_ind, 1)
    end
    return (
        re_full=full_re, im_full=full_im,
        re_indirect=ind_re, im_indirect=ind_im,
        re_direct=full_re - ind_re,
        im_direct=full_im - ind_im,
        k0=k0, gamma=gamma, phi=phi, phibar=phibar,
    )
end

function main()
    opts = _parse_args(ARGS)
    derived = _read_csv(opts["--derived-csv"])
    gap = _read_csv(opts["--gap-csv"])
    gap_by_xi = _group_by_xi(gap)

    derived_by_xi = _group_by_xi(derived)
    # Target xi values (matching core figures)
    target_xis = [-0.3, 0.0, 0.3]
    # Downsample T: every 5 MeV is enough for smooth curves
    T_step = 5.0

    mkpath(dirname(opts["--out-csv"]))
    open(opts["--out-csv"], "w") do io
        println(io, "channel,xi,T_MeV,re_full,im_full,re_indirect,im_indirect,re_direct,im_direct,k0,gamma,Phi,Phibar")
        for xi in target_xis
            haskey(derived_by_xi, xi) || continue
            drows = derived_by_xi[xi]
            all_T = sort(unique([_f(r, "T_MeV") for r in drows]))
            # Take every T_step'th point
            probe_Ts = all_T[1]:T_step:all_T[end]
            for t in probe_Ts
                pt = _compute_point(:pi, xi, t, drows, gap_by_xi)
                println(io, join([
                    "pi", string(xi), string(Int(t)),
                    string(pt.re_full), string(pt.im_full),
                    string(pt.re_indirect), string(pt.im_indirect),
                    string(pt.re_direct), string(pt.im_direct),
                    string(pt.k0), string(pt.gamma), string(pt.phi), string(pt.phibar)
                ], ','))
            end
            for t in probe_Ts
                pt = _compute_point(:K, xi, t, drows, gap_by_xi)
                println(io, join([
                    "K", string(xi), string(Int(t)),
                    string(pt.re_full), string(pt.im_full),
                    string(pt.re_indirect), string(pt.im_indirect),
                    string(pt.re_direct), string(pt.im_direct),
                    string(pt.k0), string(pt.gamma), string(pt.phi), string(pt.phibar)
                ], ','))
            end
        end
    end

    println("Wrote T-scan polarization csv: " * opts["--out-csv"])
end

main()
