#!/usr/bin/env julia

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

using Main.MesonMass: ensure_quark_params_has_A
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings, mixing_matrix_elements
using Main.PolarizationAniso: polarization_with_width
using Main.Constants_PNJL: G_fm2, K_fm5, ħc_MeV_fm

const HBARC_MEV_FM = ħc_MeV_fm

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

function _load_param_by_T(path::String)
    d = Dict{Float64,NamedTuple}()
    for r in _read_numeric_rows(path)
        length(r) >= 8 || continue
        T = r[1]
        d[T] = (
            Phi=r[5],
            Phibar=r[6],
            m_u_MeV=r[7],
            m_s_MeV=r[8],
        )
    end
    return d
end

function _load_roots(path::String; muB_MeV::Float64=600.0)
    out = Vector{NamedTuple}()
    for r in _read_numeric_rows(path)
        length(r) >= 14 || continue
        idx = Int(round(r[2]))
        idx in (3, 4) || continue
        isapprox(r[14], muB_MeV; atol=1e-9) || continue
        meson = idx == 3 ? :eta : :eta_prime
        push!(out, (
            T_MeV=r[1],
            meson=meson,
            idx=idx,
            mass_MeV=r[3],
            gamma_MeV=r[4],
            fortran_resid=r[5],
            fortran_iters=Int(round(r[6])),
            fortran_converged=Int(round(r[7])) == 1,
            fortran_sign=Int(round(r[10])),
            muB_MeV=r[14],
        ))
    end
    return out
end

function _branch_scores(row, param)
    T = row.T_MeV / HBARC_MEV_FM
    muq = (row.muB_MeV / 3.0) / HBARC_MEV_FM
    m_u = param.m_u_MeV / HBARC_MEV_FM
    m_s = param.m_s_MeV / HBARC_MEV_FM
    k0 = row.mass_MeV / HBARC_MEV_FM
    g = row.gamma_MeV / HBARC_MEV_FM

    qp = (
        m=(u=m_u, d=m_u, s=m_s),
        μ=(u=muq, d=muq, s=muq),
    )
    tp = (T=T, Φ=param.Phi, Φbar=param.Phibar, ξ=0.0)

    qpA = ensure_quark_params_has_A(qp, tp)
    Gu = calculate_G_from_A(qpA.A.u, qpA.m.u)
    Gs = calculate_G_from_A(qpA.A.s, qpA.m.s)
    Kc = calculate_effective_couplings(G_fm2, K_fm5, Gu, Gs)

    pu = polarization_with_width(:P, k0, g, 0.0, m_u, m_u, muq, muq, T, param.Phi, param.Phibar, 0.0, qpA.A.u, qpA.A.u, 0)
    ps = polarization_with_width(:P, k0, g, 0.0, m_s, m_s, muq, muq, T, param.Phi, param.Phibar, 0.0, qpA.A.s, qpA.A.s, 2)
    Πuu = ComplexF64(pu[1], pu[2])
    Πss = ComplexF64(ps[1], ps[2])

    elems = mixing_matrix_elements(Πuu, Πss, Kc, :P)
    rt = sqrt((elems.M00 - elems.M88)^2 + 4.0 * elems.M08^2)
    plus_val = (elems.M00 + elems.M88) + rt
    minus_val = (elems.M00 + elems.M88) - rt
    score_plus = hypot(real(plus_val), imag(plus_val))
    score_minus = hypot(real(minus_val), imag(minus_val))
    return score_plus, score_minus
end

function main()
    param_path = joinpath(FORTRAN_ROOT, "quark_phase", "parameter.dat")
    root_path = joinpath(FORTRAN_ROOT, "quark_phase", "root_diag.dat")
    isfile(param_path) || error("missing parameter.dat")
    isfile(root_path) || error("missing root_diag.dat")

    params = _load_param_by_T(param_path)
    rows = _load_roots(root_path; muB_MeV=600.0)
    isempty(rows) && error("no mixed-root rows found")

    out_dir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics")
    mkpath(out_dir)
    out_csv = joinpath(out_dir, "experimental_mixed_branch_identity_muB600.csv")

    open(out_csv, "w") do io
        write(io, "T_MeV,meson,fortran_sign,fortran_resid,fortran_conv,score_plus,score_minus,score_label_align,score_fortran_sign,score_identity_align,preferred_sign\n")
        for r in sort(rows; by=x -> (x.T_MeV, String(x.meson)))
            p = params[r.T_MeV]
            splus, sminus = _branch_scores(r, p)
            score_label = r.meson == :eta ? splus : sminus
            score_sign = r.fortran_sign == 1 ? splus : sminus
            score_id = min(splus, sminus)
            pref_sign = splus <= sminus ? 1 : -1
            write(io, string(
                r.T_MeV, ',', r.meson, ',', r.fortran_sign, ',', r.fortran_resid, ',', r.fortran_converged, ',',
                splus, ',', sminus, ',', score_label, ',', score_sign, ',', score_id, ',', pref_sign, '\n',
            ))
        end
    end

    println("Wrote experimental branch-identity table: $out_csv")
end

main()
