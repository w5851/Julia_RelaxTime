#!/usr/bin/env julia

using DelimitedFiles

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FORTRAN_ROOT = raw"D:\Desktop\fortran_mott\PNJL-mott-mu_T\PNJL-mu-T"

const CONSTANTS_SCRIPT = joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl")
const MODELS_SCRIPT = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")

isfile(CONSTANTS_SCRIPT) || error("missing constants script: $CONSTANTS_SCRIPT")
isfile(MODELS_SCRIPT) || error("missing models script: $MODELS_SCRIPT")

if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, CONSTANTS_SCRIPT)
end
if !isdefined(Main, :Models)
    Base.include(Main, MODELS_SCRIPT)
end

using Main.OneLoopIntegrals: B0
using Main.PolarizationAniso: polarization_with_width
using Main.MesonMass: ensure_quark_params_has_A
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings, mixing_matrix_elements
using Main.Constants_PNJL: G_fm2, K_fm5, ħc_MeV_fm

const HBARC_MEV_FM = ħc_MeV_fm

function _parse_env_float_list(name::String, default::Vector{Float64})
    raw = get(ENV, name, "")
    isempty(strip(raw)) && return default
    vals = Float64[]
    for part in split(raw, ',')
        s = strip(part)
        isempty(s) && continue
        v = tryparse(Float64, s)
        v === nothing && error("invalid numeric value '$s' in ENV[$name]")
        push!(vals, v)
    end
    isempty(vals) && error("ENV[$name] did not provide any numeric values")
    return vals
end

function _muB_output_tag(muBs_MeV::Vector{Float64})
    if length(muBs_MeV) == 1
        v = muBs_MeV[1]
        iv = round(Int, v)
        if isapprox(v, iv; atol=1e-9)
            return "muB$(iv)"
        end
    end
    return "muBmulti"
end

function _read_numeric_rows(path::String)
    rows = Vector{Vector{Float64}}()
    for ln in eachline(path)
        s = strip(ln)
        isempty(s) && continue
        cols = split(s)
        vals = Float64[]
        ok = true
        for c in cols
            v = tryparse(Float64, c)
            if v === nothing
                ok = false
                break
            end
            push!(vals, v)
        end
        ok || continue
        push!(rows, vals)
    end
    return rows
end

function _parameter_row_by_T(path::String)
    m = Dict{Float64,NamedTuple}()
    for r in _read_numeric_rows(path)
        length(r) >= 8 || continue
        T = r[1]
        m[T] = (
            phi_u=r[2],
            phi_d=r[3],
            phi_s=r[4],
            Phi=r[5],
            Phibar=r[6],
            m_u_MeV=r[7],
            m_s_MeV=r[8],
        )
    end
    return m
end

function _root_diag_rows(path::String; muBs_MeV::Vector{Float64}=[600.0], Ts::Vector{Float64}=[220.0, 240.0])
    out = Vector{NamedTuple}()
    for r in _read_numeric_rows(path)
        length(r) >= 14 || continue
        T = r[1]
        idx = Int(round(r[2]))
        muB = r[14]
        (idx == 3 || idx == 4) || continue
        any(t -> isapprox(T, t; atol=1e-9), Ts) || continue
        any(x -> isapprox(muB, x; atol=1e-9), muBs_MeV) || continue
        meson = idx == 3 ? :eta : :eta_prime
        is_sign = idx == 3 ? 1 : -1
        push!(out, (
            T_MeV=T,
            muB_MeV=muB,
            meson=meson,
            is_sign=is_sign,
            mass_MeV=r[3],
            gamma_MeV=r[4],
            fortran_resid=r[5],
            fortran_iters=Int(round(r[6])),
            fortran_converged=Int(round(r[7])) == 1,
        ))
    end
    return out
end

function _write_b0_probe_input(path::String, roots, params)
    open(path, "w") do io
        write(io, "lam,k,m1,m2,mu1,mu2,T,Phi,Phibar\n")
        for rr in roots
            p = params[rr.T_MeV]
            λ = rr.mass_MeV / HBARC_MEV_FM
            μq = (rr.muB_MeV / 3.0) / HBARC_MEV_FM
            m_u = p.m_u_MeV / HBARC_MEV_FM
            m_s = p.m_s_MeV / HBARC_MEV_FM
            T = rr.T_MeV / HBARC_MEV_FM
            write(io, string(λ, ',', 0.0, ',', m_u, ',', m_u, ',', μq, ',', μq, ',', T, ',', p.Phi, ',', p.Phibar, '\n'))
            write(io, string(λ, ',', 0.0, ',', m_s, ',', m_s, ',', μq, ',', μq, ',', T, ',', p.Phi, ',', p.Phibar, '\n'))
        end
    end
end

function _read_b0_probe_output(path::String)
    lines = readlines(path)
    out = Vector{NamedTuple}()
    for ln in lines[2:end]
        s = strip(ln)
        isempty(s) && continue
        cols = split(s, ',')
        length(cols) >= 11 || continue
        vals = map(x -> parse(Float64, strip(x)), cols)
        push!(out, (
            lam=vals[1],
            k=vals[2],
            m1=vals[3],
            m2=vals[4],
            mu1=vals[5],
            mu2=vals[6],
            T=vals[7],
            Phi=vals[8],
            Phibar=vals[9],
            b0_re=vals[10],
            b0_im=vals[11],
        ))
    end
    return out
end

function _read_pi_diag(path::String)
    rows = Dict{Tuple{Float64,Symbol,Float64},NamedTuple}()
    for r in _read_numeric_rows(path)
        length(r) >= 7 || continue
        T = r[1]
        idx = Int(round(r[2]))
        meson = idx == 3 ? :eta : idx == 4 ? :eta_prime : nothing
        meson === nothing && continue
        muB = r[7]
        rows[(T, meson, muB)] = (
            pi_uu_re=r[3],
            pi_uu_im=r[4],
            pi_ss_re=r[5],
            pi_ss_im=r[6],
        )
    end
    return rows
end

function _pick_b0(rows, λ, m)
    best = nothing
    bestd = Inf
    for r in rows
        d = abs(r.lam - λ) + abs(r.m1 - m) + abs(r.m2 - m) + abs(r.k)
        if d < bestd
            bestd = d
            best = r
        end
    end
    return best
end

function _julia_layers(rr, p, b0_fortran_rows, pi_diag_rows)
    T = rr.T_MeV / HBARC_MEV_FM
    μq = (rr.muB_MeV / 3.0) / HBARC_MEV_FM
    m_u = p.m_u_MeV / HBARC_MEV_FM
    m_s = p.m_s_MeV / HBARC_MEV_FM
    k0 = rr.mass_MeV / HBARC_MEV_FM
    γ = rr.gamma_MeV / HBARC_MEV_FM

    qp = (
        m=(u=m_u, d=m_u, s=m_s),
        μ=(u=μq, d=μq, s=μq),
    )
    tp = (T=T, Φ=p.Phi, Φbar=p.Phibar, ξ=0.0)

    qpA = ensure_quark_params_has_A(qp, tp)
    G_u = calculate_G_from_A(qpA.A.u, qpA.m.u)
    G_s = calculate_G_from_A(qpA.A.s, qpA.m.s)
    Kc = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    λ = k0
    b0_uu_re, b0_uu_im = B0(λ, 0.0, m_u, μq, m_u, μq, T; Φ=p.Phi, Φbar=p.Phibar)
    b0_ss_re, b0_ss_im = B0(λ, 0.0, m_s, μq, m_s, μq, T; Φ=p.Phi, Φbar=p.Phibar)

    puu_re, puu_im = polarization_with_width(:P, k0, γ, 0.0, m_u, m_u, μq, μq, T, p.Phi, p.Phibar, 0.0, qpA.A.u, qpA.A.u, 0)
    pss_re, pss_im = polarization_with_width(:P, k0, γ, 0.0, m_s, m_s, μq, μq, T, p.Phi, p.Phibar, 0.0, qpA.A.s, qpA.A.s, 2)
    Πuu = ComplexF64(puu_re, puu_im)
    Πss = ComplexF64(pss_re, pss_im)

    elems = mixing_matrix_elements(Πuu, Πss, Kc, :P)
    M00 = elems.M00
    M08 = elems.M08
    M88 = elems.M88
    sqrt_term = sqrt((M00 - M88)^2 + 4.0 * M08^2)
    branch_val = (M00 + M88) + rr.is_sign * sqrt_term
    branch_score = hypot(real(branch_val), imag(branch_val))

    b0f_uu = _pick_b0(b0_fortran_rows, λ, m_u)
    b0f_ss = _pick_b0(b0_fortran_rows, λ, m_s)

    return (
        b0_uu_jl_re=b0_uu_re,
        b0_uu_jl_im=b0_uu_im,
        b0_ss_jl_re=b0_ss_re,
        b0_ss_jl_im=b0_ss_im,
        b0_uu_ft_re=b0f_uu === nothing ? NaN : b0f_uu.b0_re,
        b0_uu_ft_im=b0f_uu === nothing ? NaN : b0f_uu.b0_im,
        b0_ss_ft_re=b0f_ss === nothing ? NaN : b0f_ss.b0_re,
        b0_ss_ft_im=b0f_ss === nothing ? NaN : b0f_ss.b0_im,
        puu_re=puu_re,
        puu_im=puu_im,
        pss_re=pss_re,
        pss_im=pss_im,
        puu_ft_re=get(pi_diag_rows, (rr.T_MeV, rr.meson, rr.muB_MeV), (pi_uu_re=NaN,)).pi_uu_re,
        puu_ft_im=get(pi_diag_rows, (rr.T_MeV, rr.meson, rr.muB_MeV), (pi_uu_im=NaN,)).pi_uu_im,
        pss_ft_re=get(pi_diag_rows, (rr.T_MeV, rr.meson, rr.muB_MeV), (pi_ss_re=NaN,)).pi_ss_re,
        pss_ft_im=get(pi_diag_rows, (rr.T_MeV, rr.meson, rr.muB_MeV), (pi_ss_im=NaN,)).pi_ss_im,
        M00_re=real(M00),
        M00_im=imag(M00),
        M08_re=real(M08),
        M08_im=imag(M08),
        M88_re=real(M88),
        M88_im=imag(M88),
        branch_re=real(branch_val),
        branch_im=imag(branch_val),
        branch_score=branch_score,
    )
end

function main()
    root_diag_path = joinpath(FORTRAN_ROOT, "quark_phase", "root_diag.dat")
    param_path = joinpath(FORTRAN_ROOT, "quark_phase", "parameter.dat")
    b0_driver = joinpath(FORTRAN_ROOT, "b0_driver.exe")
    pi_diag_path = joinpath(FORTRAN_ROOT, "quark_phase", "pi_diag.dat")

    isfile(root_diag_path) || error("missing fortran root diagnostics: $root_diag_path")
    isfile(param_path) || error("missing fortran parameter file: $param_path")
    isfile(b0_driver) || error("missing fortran b0 driver: $b0_driver")
    isfile(pi_diag_path) || error("missing fortran pi diagnostics: $pi_diag_path")

    muBs = _parse_env_float_list("RECON_MUBS", [600.0])
    Ts = _parse_env_float_list("RECON_TS", [220.0, 240.0])
    out_tag = get(ENV, "RECON_OUT_TAG", _muB_output_tag(muBs))

    params = _parameter_row_by_T(param_path)
    roots = _root_diag_rows(root_diag_path; muBs_MeV=muBs, Ts=Ts)
    isempty(roots) && error("no eta/eta_prime diagnostics matched in $root_diag_path")

    out_dir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics")
    mkpath(out_dir)
    b0_in = joinpath(out_dir, "layer_reconcile_b0_input.csv")
    b0_out = joinpath(out_dir, "layer_reconcile_b0_fortran_output.csv")

    _write_b0_probe_input(b0_in, roots, params)
    run(`$b0_driver $b0_in $b0_out`)
    b0_fortran_rows = _read_b0_probe_output(b0_out)
    pi_diag_rows = _read_pi_diag(pi_diag_path)

    out_csv = joinpath(out_dir, "fortran_julia_layer_reconcile_$(out_tag)_eta_etap_pi.csv")
    open(out_csv, "w") do io
        write(io, "muB_MeV,T_MeV,meson,mass_MeV,gamma_MeV,fortran_resid,fortran_iters,fortran_converged,")
        write(io, "b0_uu_ft_re,b0_uu_ft_im,b0_uu_jl_re,b0_uu_jl_im,b0_uu_abs_diff,")
        write(io, "b0_ss_ft_re,b0_ss_ft_im,b0_ss_jl_re,b0_ss_jl_im,b0_ss_abs_diff,")
        write(io, "puu_ft_re,puu_ft_im,puu_re,puu_im,puu_abs_diff,")
        write(io, "pss_ft_re,pss_ft_im,pss_re,pss_im,pss_abs_diff,")
        write(io, "M00_re,M00_im,M08_re,M08_im,M88_re,M88_im,branch_re,branch_im,branch_score\n")

        for rr in sort(roots; by=r -> (r.T_MeV, String(r.meson)))
            p = params[rr.T_MeV]
            jl = _julia_layers(rr, p, b0_fortran_rows, pi_diag_rows)
            b0_uu_diff = hypot(jl.b0_uu_jl_re - jl.b0_uu_ft_re, jl.b0_uu_jl_im - jl.b0_uu_ft_im)
            b0_ss_diff = hypot(jl.b0_ss_jl_re - jl.b0_ss_ft_re, jl.b0_ss_jl_im - jl.b0_ss_ft_im)
            puu_diff = hypot(jl.puu_re - jl.puu_ft_re, jl.puu_im - jl.puu_ft_im)
            pss_diff = hypot(jl.pss_re - jl.pss_ft_re, jl.pss_im - jl.pss_ft_im)
            write(io, string(
                rr.muB_MeV, ',', rr.T_MeV, ',', rr.meson, ',', rr.mass_MeV, ',', rr.gamma_MeV, ',', rr.fortran_resid, ',', rr.fortran_iters, ',', rr.fortran_converged, ',',
                jl.b0_uu_ft_re, ',', jl.b0_uu_ft_im, ',', jl.b0_uu_jl_re, ',', jl.b0_uu_jl_im, ',', b0_uu_diff, ',',
                jl.b0_ss_ft_re, ',', jl.b0_ss_ft_im, ',', jl.b0_ss_jl_re, ',', jl.b0_ss_jl_im, ',', b0_ss_diff, ',',
                jl.puu_ft_re, ',', jl.puu_ft_im, ',', jl.puu_re, ',', jl.puu_im, ',', puu_diff, ',',
                jl.pss_ft_re, ',', jl.pss_ft_im, ',', jl.pss_re, ',', jl.pss_im, ',', pss_diff, ',',
                jl.M00_re, ',', jl.M00_im, ',', jl.M08_re, ',', jl.M08_im, ',', jl.M88_re, ',', jl.M88_im, ',', jl.branch_re, ',', jl.branch_im, ',', jl.branch_score,
                '\n',
            ))
        end
    end

    println("Wrote layer reconcile table: $out_csv")
end

main()
