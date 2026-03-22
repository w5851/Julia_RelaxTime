#!/usr/bin/env julia

using DelimitedFiles

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

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

const HBARC_MEV_FM = Main.Constants_PNJL.ħc_MeV_fm
const MESON_WORKFLOW_MOD = Main.Models.meson_workflow_module()

using Main.MesonMass: meson_mass_equation, ensure_quark_params_has_A
using Main.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using Main.Constants_PNJL: G_fm2, K_fm5

function _as_namedtuple_params(quark, thermo)
    qp = (
        m=(u=Float64(quark.m.u), d=Float64(quark.m.d), s=Float64(quark.m.s)),
        μ=(u=Float64(quark.μ.u), d=Float64(quark.μ.d), s=Float64(quark.μ.s)),
    )
    tp = (
        T=Float64(thermo.T),
        Φ=Float64(thermo.Φ),
        Φbar=Float64(thermo.Φbar),
        ξ=Float64(thermo.ξ),
    )
    return qp, tp
end

function _prepare_k_coeffs(qp::NamedTuple, tp::NamedTuple)
    qpA = ensure_quark_params_has_A(qp, tp)
    G_u = calculate_G_from_A(qpA.A.u, qpA.m.u)
    G_s = calculate_G_from_A(qpA.A.s, qpA.m.s)
    Kc = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
    return qpA, Kc
end

function _scan_landscape(; muB_MeV::Float64, T_MeV::Float64,
    m_min_MeV::Float64, m_max_MeV::Float64, m_steps::Int,
    g_min_MeV::Float64, g_max_MeV::Float64, g_steps::Int)

    T_fm = T_MeV / HBARC_MEV_FM
    muq_fm = (muB_MeV / 3.0) / HBARC_MEV_FM

    base = MESON_WORKFLOW_MOD.solve_gap_and_meson_point(
        T_fm,
        muq_fm;
        xi=0.0,
        mesons=(:eta, :eta_prime),
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=25,),
        mass_kwargs=(iterations=40,),
    )

    qp, tp = _as_namedtuple_params(base.quark_params, base.thermo_params)
    qpA, Kc = _prepare_k_coeffs(qp, tp)

    m_grid = collect(range(m_min_MeV / HBARC_MEV_FM, m_max_MeV / HBARC_MEV_FM, length=m_steps))
    g_grid = collect(range(g_min_MeV / HBARC_MEV_FM, g_max_MeV / HBARC_MEV_FM, length=g_steps))

    rows = Vector{NTuple{7,Float64}}()
    best = (score=Inf, mass=NaN, gamma=NaN, re=NaN, im=NaN)

    for m in m_grid
        for g in g_grid
            f = meson_mass_equation(:eta_prime, m, g, 0.0, qpA, tp, Kc)
            re = real(f)
            im = imag(f)
            score = hypot(re, im)
            push!(rows, (
                muB_MeV,
                T_MeV,
                m * HBARC_MEV_FM,
                g * HBARC_MEV_FM,
                re,
                im,
                score,
            ))
            if score < best.score
                best = (score=score, mass=m * HBARC_MEV_FM, gamma=g * HBARC_MEV_FM, re=re, im=im)
            end
        end
    end

    return rows, best, base
end

function main()
    out_dir = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "diagnostics")
    mkpath(out_dir)

    cases = [
        (muB_MeV=0.0, T_MeV=220.0, m_min_MeV=380.0, m_max_MeV=820.0, m_steps=301, g_min_MeV=0.0, g_max_MeV=120.0, g_steps=241),
        (muB_MeV=0.0, T_MeV=240.0, m_min_MeV=380.0, m_max_MeV=820.0, m_steps=301, g_min_MeV=0.0, g_max_MeV=120.0, g_steps=241),
    ]

    summary_lines = String[
        "muB_MeV,T_MeV,base_mass_MeV,base_gamma_MeV,base_residual,base_converged,best_mass_MeV,best_gamma_MeV,best_re,best_im,best_score",
    ]

    for c in cases
        rows, best, base = _scan_landscape(; c...)

        out_csv = joinpath(out_dir, "eta_prime_residual_landscape_muB$(Int(round(c.muB_MeV)))_T$(Int(round(c.T_MeV))).csv")
        open(out_csv, "w") do io
            write(io, "muB_MeV,T_MeV,mass_MeV,gamma_MeV,residual_re,residual_im,residual_score\n")
            for r in rows
                write(io, string(r[1], ',', r[2], ',', r[3], ',', r[4], ',', r[5], ',', r[6], ',', r[7], '\n'))
            end
        end

        b = base.meson_results[:eta_prime]
        push!(summary_lines,
            string(
                c.muB_MeV, ',', c.T_MeV, ',',
                b.mass * HBARC_MEV_FM, ',', b.gamma * HBARC_MEV_FM, ',', b.residual, ',', b.converged, ',',
                best.mass, ',', best.gamma, ',', best.re, ',', best.im, ',', best.score,
            ),
        )
    end

    summary_csv = joinpath(out_dir, "eta_prime_residual_landscape_summary.csv")
    write(summary_csv, join(summary_lines, "\n") * "\n")
    println("Wrote diagnostics summary: $summary_csv")
end

main()
