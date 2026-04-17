if !isdefined(Main, :Constants_PNJL)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "constants", "Constants_PNJL.jl"))
end

if !isdefined(Main, :Models)
    include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
end

const HBARC = Main.Constants_PNJL.ħc_MeV_fm
const WF = Main.Models.meson_workflow_module()

const _MESONS = (:pi, :K, :eta, :eta_prime, :sigma_pi, :sigma_K, :sigma, :sigma_prime)
const _RTOL_BY_MESON = Dict(
    :pi => 1e-2,
    :K => 1e-2,
    :eta => 4e-2,
    :eta_prime => 3e-1,
    :sigma_pi => 1e-2,
    :sigma_K => 1e-2,
    :sigma => 1e-2,
    :sigma_prime => 1e-2,
)

@inline _julia_symbol_for_legacy(meson::Symbol) =
    meson == :eta ? :eta_prime :
    meson == :eta_prime ? :eta :
    meson == :sigma ? :sigma_prime :
    meson == :sigma_prime ? :sigma : meson

@inline _fortran_ps(meson::Symbol) = (meson in (:pi, :K, :eta, :eta_prime)) ? -1 : 1
@inline _fortran_is(meson::Symbol) = (meson in (:eta, :sigma)) ? -1 : 1

function _fortran_flavor_triplet(meson::Symbol)
    if meson == :pi
        return (1, 2, 3)
    elseif meson == :K
        # mu_T main uses (3,1,2)
        return (3, 1, 2)
    elseif meson == :sigma_pi
        return (1, 2, 3)
    elseif meson == :sigma_K
        # mu_T main uses (1,3,2)
        return (1, 3, 2)
    elseif meson in (:eta, :eta_prime, :sigma, :sigma_prime)
        return (1, 2, 3)
    end
    throw(ArgumentError("unsupported meson: $meson"))
end

function _read_csv_rows(path::String)
    lines = readlines(path)
    header = strip.(split(strip(lines[1]), ','))
    rows = Dict{String, String}[]
    for ln in lines[2:end]
        s = strip(ln)
        isempty(s) && continue
        c = strip.(split(s, ','))
        length(c) == length(header) || continue
        r = Dict{String, String}()
        for i in eachindex(header)
            r[header[i]] = c[i]
        end
        push!(rows, r)
    end
    return rows
end

function _collect_legacy_fail_rows()
    paths = (
        joinpath(@__DIR__, "..", "..", "..", "tests", "validation", "data", "targets", "relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB0_v1.csv"),
        joinpath(@__DIR__, "..", "..", "..", "tests", "validation", "data", "targets", "relaxtime", "legacy", "meson", "legacy_meson_scan_fortran_muB600_v1.csv"),
    )

    targets = Dict{Tuple{Float64, Float64, Symbol}, Float64}()
    for p in paths
        for row in _read_csv_rows(p)
            row["source_impl"] == "fortran" || continue
            parse(Float64, row["xi"]) == 0.0 || continue
            meson = Symbol(row["meson"])
            meson in _MESONS || continue
            key = (parse(Float64, row["muB_MeV"]), parse(Float64, row["T_MeV"]), meson)
            targets[key] = parse(Float64, row["mass_MeV"])
        end
    end

    failed = NamedTuple[]
    for muB in (0.0, 600.0), T in (120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0)
        res = WF.solve_gap_and_meson_point(
            T / HBARC,
            (muB / 3.0) / HBARC;
            xi=0.0,
            mesons=_MESONS,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=25,),
            mass_kwargs=(iterations=25,),
        )
        for meson in _MESONS
            expected = targets[(muB, T, meson)]
            julia_meson = _julia_symbol_for_legacy(meson)
            actual = Float64(res.meson_results[julia_meson].mass) * HBARC
            if !(isfinite(actual) && isapprox(actual, expected; rtol=_RTOL_BY_MESON[meson], atol=0.0))
                push!(failed, (
                    muB_MeV=muB,
                    T_MeV=T,
                    legacy_meson=meson,
                    julia_meson=julia_meson,
                    target_mass_fm=expected / HBARC,
                    qp=res.quark_params,
                    tp=res.thermo_params,
                ))
            end
        end
    end
    return failed
end

function _kcoeffs_from_qp(qp)
    Gu = Main.EffectiveCouplings.calculate_G_from_A(qp.A.u, qp.m.u)
    Gs = Main.EffectiveCouplings.calculate_G_from_A(qp.A.s, qp.m.s)
    return Main.EffectiveCouplings.calculate_effective_couplings(
        Main.Constants_PNJL.G_fm2,
        Main.Constants_PNJL.K_fm5,
        Gu,
        Gs,
    )
end

function _pi_nonmixed_from_julia_lowlevel(meson::Symbol, k0::Float64, gamma::Float64, qp, tp)
    pp = Main.MesonMass._meson_polarization_params(meson, qp)
    re, im = Main.PolarizationAniso.polarization_with_width(
        pp.channel,
        k0,
        gamma,
        0.0,
        pp.m1,
        pp.m2,
        pp.μ1,
        pp.μ2,
        tp.T,
        tp.Φ,
        tp.Φbar,
        tp.ξ,
        pp.A1,
        pp.A2,
        pp.num_s_quark,
    )
    return ComplexF64(re, im)
end

function _fortran_pi_ij(meson::Symbol, qp)
    ps = _fortran_ps(meson)
    (_, _, kk) = _fortran_flavor_triplet(meson)
    A = (Float64(qp.A.u), Float64(qp.A.d), Float64(qp.A.s))
    m = (Float64(qp.m.u), Float64(qp.m.d), Float64(qp.m.s))
    return Main.Constants_PNJL.G_fm2 - ps * 0.5 * Main.Constants_PNJL.K_fm5 * A[kk] * (-3.0 / (4.0 * pi^2) * m[kk])
end

function _hard_compare_nonmixed(meson::Symbol, k0::Float64, gamma::Float64, qp, tp, Kc)
    Π = _pi_nonmixed_from_julia_lowlevel(meson, k0, gamma, qp, tp)
    K_fortran = _fortran_pi_ij(meson, qp)
    K_julia = Main.MesonMass._meson_coupling(meson, Kc)

    f_fortran = ComplexF64(1.0 - 4.0 * K_fortran * real(Π), imag(Π))
    f_julia = Main.MesonMass.meson_mass_equation(meson, k0, gamma, 0.0, qp, tp, Kc)

    return (
        K_fortran=K_fortran,
        K_julia=K_julia,
        delta_K=K_fortran - K_julia,
        f_fortran=f_fortran,
        f_julia=f_julia,
        delta_f_re=real(f_fortran) - real(f_julia),
        imag_scale=(abs(imag(f_fortran)) > 0 ? imag(f_julia) / imag(f_fortran) : NaN),
    )
end

function _fortran_mixed_value(meson::Symbol, k0::Float64, gamma::Float64, qp, tp)
    ps = _fortran_ps(meson)
    sign_is = _fortran_is(meson)

    m_u = Float64(qp.m.u)
    m_s = Float64(qp.m.s)
    μ_u = Float64(qp.μ.u)
    μ_s = Float64(qp.μ.s)
    A_u = Float64(qp.A.u)
    A_s = Float64(qp.A.s)

    coeff = (-3.0 / (4.0 * pi^2))
    G = Main.Constants_PNJL.G_fm2
    K = Main.Constants_PNJL.K_fm5

    k00 = G + ps * K * (2.0 * A_u * m_u + A_s * m_s) / 3.0 * coeff
    k88 = G - ps * K * (4.0 * A_u * m_u - A_s * m_s) / 6.0 * coeff
    k08 = -ps * sqrt(2.0) * K * (A_u * m_u - A_s * m_s) / 6.0 * coeff
    detk = k00 * k88 - k08^2

    channel = (meson in (:eta, :eta_prime)) ? :P : :S
    Πuu_re, Πuu_im = Main.PolarizationAniso.polarization_with_width(
        channel, k0, gamma, 0.0,
        m_u, m_u,
        μ_u, μ_u,
        tp.T, tp.Φ, tp.Φbar, tp.ξ,
        A_u, A_u,
        0,
    )
    Πss_re, Πss_im = Main.PolarizationAniso.polarization_with_width(
        channel, k0, gamma, 0.0,
        m_s, m_s,
        μ_s, μ_s,
        tp.T, tp.Φ, tp.Φbar, tp.ξ,
        A_s, A_s,
        2,
    )
    Πuu = ComplexF64(Πuu_re, Πuu_im)
    Πss = ComplexF64(Πss_re, Πss_im)

    m00 = k00 - (4.0 / 3.0) * (Πuu + 2.0 * Πss) * detk
    m08 = k08 + (4.0 / 3.0) * sqrt(2.0) * (Πuu - Πss) * detk
    m88 = k88 - (4.0 / 3.0) * (2.0 * Πuu + Πss) * detk

    root = sqrt((m00 - m88)^2 + 4.0 * m08^2)
    f_fortran = m00 + m88 + sign_is * root
    return (f_fortran=f_fortran, k00=k00, k08=k08, k88=k88, detk=detk)
end

function _hard_compare_mixed(meson::Symbol, k0::Float64, gamma::Float64, qp, tp, Kc)
    f_fortran = _fortran_mixed_value(meson, k0, gamma, qp, tp)
    f_julia = Main.MesonMass.meson_mass_equation(meson, k0, gamma, 0.0, qp, tp, Kc)
    return (
        f_fortran=f_fortran.f_fortran,
        f_julia=f_julia,
        delta_f=f_fortran.f_fortran - f_julia,
        k00=f_fortran.k00,
        k08=f_fortran.k08,
        k88=f_fortran.k88,
        detk=f_fortran.detk,
    )
end

function main()
    failed = _collect_legacy_fail_rows()
    println("hard_source_check_rows=", length(failed))
    println("muB_MeV,T_MeV,meson,type,delta_core_re,delta_core_im,note")

    max_nonmixed_deltaK = 0.0
    max_mixed_delta = 0.0
    nonmixed_rows = 0
    mixed_rows = 0

    for row in failed
        meson = row.julia_meson
        k0 = Float64(row.target_mass_fm)
        gamma = 0.0
        qp = Main.MesonMass.ensure_quark_params_has_A((
            m=(
                u=Float64(row.qp.m.u),
                d=Float64(row.qp.m.d),
                s=Float64(row.qp.m.s),
            ),
            μ=(
                u=Float64(row.qp.μ.u),
                d=Float64(row.qp.μ.d),
                s=Float64(row.qp.μ.s),
            ),
        ), (
            T=Float64(row.tp.T),
            Φ=Float64(row.tp.Φ),
            Φbar=Float64(row.tp.Φbar),
            ξ=Float64(row.tp.ξ),
        ))
        tp = (
            T=Float64(row.tp.T),
            Φ=Float64(row.tp.Φ),
            Φbar=Float64(row.tp.Φbar),
            ξ=Float64(row.tp.ξ),
        )
        Kc = _kcoeffs_from_qp(qp)

        if meson in (:eta, :eta_prime, :sigma, :sigma_prime)
            mixed_rows += 1
            cmp = _hard_compare_mixed(meson, k0, gamma, qp, tp, Kc)
            d = cmp.delta_f
            max_mixed_delta = max(max_mixed_delta, abs(d))
            note = abs(d) > 1e-8 ? "HARD_MISMATCH" : "formula_equivalent"
            println(string(row.muB_MeV, ",", row.T_MeV, ",", meson, ",mixed,", real(d), ",", imag(d), ",", note))
        else
            nonmixed_rows += 1
            cmp = _hard_compare_nonmixed(meson, k0, gamma, qp, tp, Kc)
            max_nonmixed_deltaK = max(max_nonmixed_deltaK, abs(cmp.delta_K))
            note = abs(cmp.delta_K) > 1e-10 ? "HARD_COUPLING_MISMATCH" : "formula_equivalent_core"
            if isfinite(cmp.imag_scale)
                note = string(note, ";imag_scale_julia_over_fortranF2=", cmp.imag_scale)
            end
            println(string(row.muB_MeV, ",", row.T_MeV, ",", meson, ",nonmixed,", cmp.delta_f_re, ",", 0.0, ",", note))
        end
    end

    println("nonmixed_rows=", nonmixed_rows)
    println("mixed_rows=", mixed_rows)
    println("max_abs_deltaK_nonmixed=", max_nonmixed_deltaK)
    println("max_abs_delta_mixed_complex=", max_mixed_delta)
end

main()
