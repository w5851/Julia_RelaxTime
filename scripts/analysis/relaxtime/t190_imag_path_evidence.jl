const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

using .Constants_PNJL: N_color

const OUTDIR = raw"D:\Desktop\Temp\relaxtime_t190_window"
const OUT_DETAIL = joinpath(OUTDIR, "t190_imag_path_evidence_detail.csv")
const OUT_SUMMARY = joinpath(OUTDIR, "t190_imag_path_evidence_summary.csv")

function ensure_outdir(path::AbstractString)
    isdir(path) || mkpath(path)
end

function infer_b0_from_pi(pi::ComplexF64, channel::Symbol, k0::Float64, k_norm::Float64,
    m1::Float64, m2::Float64, μ1::Float64, μ2::Float64, A1::Float64, A2::Float64)
    factor = -N_color / (8π^2)
    λ = k0 + μ1 - μ2
    pref = k_norm^2 - λ^2
    if channel === :P
        pref += (m1 - m2)^2
    elseif channel === :S
        pref += (m1 + m2)^2
    else
        error("unsupported channel: $channel")
    end
    return ((pi / factor) - (A1 + A2)) / pref
end

function single_tilde_term(sign_flag::Symbol, λ::Float64, m::Float64, mprime::Float64, μ::Float64,
    T::Float64, Φ::Float64, Φbar::Float64)
    m_pos = max(m, 0.0)
    mprime_pos = max(mprime, 0.0)
    Emin = m_pos
    Emax = Main.OneLoopIntegrals.energy_cutoff(m_pos)
    denominator_term = (λ^2 + m_pos^2 - mprime_pos^2) / 2.0
    poles = Main.OneLoopIntegrals.singularity_k_zero(λ, Emin, Emax, denominator_term)
    t_re, t_im = Main.OneLoopIntegrals.tilde_B0_k_zero(sign_flag, λ, m_pos, mprime_pos, μ, T, Φ, Φbar)
    has_pole = !isempty(poles)
    E0 = has_pole ? poles[1] : NaN
    return (has_pole=has_pole, E0=E0, Emin=Emin, Emax=Emax, t_re=t_re, t_im=t_im)
end

function b0_term_decomp(λ::Float64, m1::Float64, μ1::Float64, m2::Float64, μ2::Float64,
    T::Float64, Φ::Float64, Φbar::Float64)
    term1 = single_tilde_term(:plus, -λ, m1, m2, μ1, T, Φ, Φbar)
    term2 = single_tilde_term(:minus, λ, m1, m2, μ1, T, Φ, Φbar)
    term3 = single_tilde_term(:plus, λ, m2, m1, μ2, T, Φ, Φbar)
    term4 = single_tilde_term(:minus, -λ, m2, m1, μ2, T, Φ, Φbar)
    b0_re = term1.t_re - term2.t_re + term3.t_re - term4.t_re
    b0_im = term1.t_im - term2.t_im + term3.t_im - term4.t_im
    return (term1=term1, term2=term2, term3=term3, term4=term4, b0_re=b0_re, b0_im=b0_im)
end

function tmid_for_process(process::Symbol, s::Float64, st)
    th = process_threshold_info(process, st.quark_params)
    tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
    return 0.5 * (tb.t_min + tb.t_max)
end

function write_simple_rows!(io, T_MeV::Float64, muB_MeV::Float64, xi::Float64, ds::Float64)
    process = :usbar_to_usbar
    st = build_state_point(T_MeV, muB_MeV, xi)
    th = process_threshold_info(process, st.quark_params)
    s = th.s_th + ds
    t = tmid_for_process(process, s, st)
    cms = Main.TotalPropagator.calculate_cms_momentum(process, s, t, :s, st.quark_params)

    m1 = st.quark_params.m.u
    m2 = st.quark_params.m.s
    μ1 = st.quark_params.μ.u
    μ2 = st.quark_params.μ.s
    A1 = st.quark_params.A.u
    A2 = st.quark_params.A.s

    Π_re, Π_im = Main.PolarizationCache.polarization_aniso_cached(
        :P, cms.k0, cms.k,
        m1, m2,
        μ1, μ2,
        st.thermo_params.T, st.thermo_params.Φ, st.thermo_params.Φbar, st.thermo_params.ξ,
        A1, A2,
        1,
    )
    Π = ComplexF64(Π_re, Π_im)
    den = 1.0 - 4.0 * st.K_coeffs.K4567_plus * Π

    λ = cms.k0 + μ1 - μ2
    b0_pos = b0_term_decomp(λ, m1, μ1, m2, μ2, st.thermo_params.T, st.thermo_params.Φ, st.thermo_params.Φbar)
    b0_neg = b0_term_decomp(-λ, m1, μ1, m2, μ2, st.thermo_params.T, st.thermo_params.Φ, st.thermo_params.Φbar)
    b0_im_used = 0.5 * (b0_pos.b0_im + b0_neg.b0_im)
    b0_from_pi = infer_b0_from_pi(Π, :P, cms.k0, cms.k, m1, m2, μ1, μ2, A1, A2)

    for (branch_tag, info) in (("lambda_pos", b0_pos), ("lambda_neg", b0_neg))
        for (term_tag, term) in (("term1", info.term1), ("term2", info.term2), ("term3", info.term3), ("term4", info.term4))
            println(io, join((
                "simple",
                string(process),
                string(xi),
                string(ds),
                string(cms.k0),
                string(cms.k),
                branch_tag,
                term_tag,
                string(term.has_pole),
                string(term.E0),
                string(term.Emin),
                string(term.Emax),
                string(term.t_im),
                string(real(den)),
                string(imag(den)),
                string(Π_im),
                string(real(b0_from_pi)),
                string(imag(b0_from_pi)),
                string(b0_im_used),
            ), ','))
        end
    end
end

function write_mixed_rows!(io, T_MeV::Float64, muB_MeV::Float64, xi::Float64, ds::Float64)
    process = :uubar_to_ddbar
    st = build_state_point(T_MeV, muB_MeV, xi)
    th = process_threshold_info(process, st.quark_params)
    s = th.s_th + ds
    t = tmid_for_process(process, s, st)

    decomp = decompose_mixed_p_propagator_chain(process, s, t, st.quark_params, st.thermo_params, st.K_coeffs)
    cms = Main.TotalPropagator.calculate_cms_momentum(process, s, t, :s, st.quark_params)

    λ_uu = cms.k0 + st.quark_params.μ.u - st.quark_params.μ.u
    λ_ss = cms.k0 + st.quark_params.μ.s - st.quark_params.μ.s

    b0_uu = b0_term_decomp(λ_uu,
        st.quark_params.m.u, st.quark_params.μ.u,
        st.quark_params.m.u, st.quark_params.μ.u,
        st.thermo_params.T, st.thermo_params.Φ, st.thermo_params.Φbar)
    b0_ss = b0_term_decomp(λ_ss,
        st.quark_params.m.s, st.quark_params.μ.s,
        st.quark_params.m.s, st.quark_params.μ.s,
        st.thermo_params.T, st.thermo_params.Φ, st.thermo_params.Φbar)

    for (flavor_tag, info, pi_val, m, μ, A) in (
        ("uu", b0_uu, decomp.Π_uu, st.quark_params.m.u, st.quark_params.μ.u, st.quark_params.A.u),
        ("ss", b0_ss, decomp.Π_ss, st.quark_params.m.s, st.quark_params.μ.s, st.quark_params.A.s),
    )
        b0_from_pi = infer_b0_from_pi(pi_val, :P, cms.k0, cms.k, m, m, μ, μ, A, A)
        for (term_tag, term) in (("term1", info.term1), ("term2", info.term2), ("term3", info.term3), ("term4", info.term4))
            println(io, join((
                "mixed_" * flavor_tag,
                string(process),
                string(xi),
                string(ds),
                string(cms.k0),
                string(cms.k),
                "lambda_pos",
                term_tag,
                string(term.has_pole),
                string(term.E0),
                string(term.Emin),
                string(term.Emax),
                string(term.t_im),
                string(real(decomp.detM)),
                string(imag(decomp.detM)),
                string(imag(pi_val)),
                string(real(b0_from_pi)),
                string(imag(b0_from_pi)),
                string(info.b0_im),
            ), ','))
        end
    end
end

function build_summary(detail_path::String, out_path::String)
    rows = readlines(detail_path)
    header = split(rows[1], ',')
    idx = Dict(name => i for (i, name) in enumerate(header))

    groups = Dict{Tuple{String, String}, Vector{Vector{SubString{String}}}}()
    for line in rows[2:end]
        cols = split(line, ',')
        key = (String(cols[idx["branch"]]), String(cols[idx["xi"]]))
        push!(get!(groups, key, Vector{Vector{SubString{String}}}()), cols)
    end

    open(out_path, "w") do io
        println(io, "branch,xi,k_min,k_max,has_pole_terms,all_terms,pole_fraction,den_re,den_im,pi_im_abs_max,b0_im_used_abs_max")
        for ((branch, xi), cols_vec) in sort(collect(groups); by=x -> (x[1][1], parse(Float64, x[1][2])))
            k_vals = [parse(Float64, String(c[idx["k_norm"]])) for c in cols_vec]
            has_pole = [lowercase(String(c[idx["has_pole"]])) == "true" for c in cols_vec]
            den_re = parse(Float64, String(cols_vec[1][idx["den_re"]]))
            den_im = parse(Float64, String(cols_vec[1][idx["den_im"]]))
            pi_im_vals = [abs(parse(Float64, String(c[idx["pi_im"]]))) for c in cols_vec]
            b0_used_vals = [abs(parse(Float64, String(c[idx["b0_im_used"]]))) for c in cols_vec]
            pole_n = count(identity, has_pole)
            total_n = length(has_pole)
            println(io, join((
                branch,
                xi,
                string(minimum(k_vals)),
                string(maximum(k_vals)),
                string(pole_n),
                string(total_n),
                string(pole_n / total_n),
                string(den_re),
                string(den_im),
                string(maximum(pi_im_vals)),
                string(maximum(b0_used_vals)),
            ), ','))
        end
    end
end

function main()
    ensure_outdir(OUTDIR)
    T_MeV = 190.0
    muB_MeV = 0.0
    ds = 1.0e-8

    xi_simple = [-0.46, -0.44, -0.42, -0.40]
    xi_mixed = [-0.12, -0.10, -0.08]

    open(OUT_DETAIL, "w") do io
        println(io, "branch,process,xi,delta_s,k0_s,k_norm,lambda_branch,term,has_pole,E0,Emin,Emax,term_im,den_re,den_im,pi_im,b0_re_from_pi,b0_im_from_pi,b0_im_used")
        for xi in xi_simple
            write_simple_rows!(io, T_MeV, muB_MeV, xi, ds)
        end
        for xi in xi_mixed
            write_mixed_rows!(io, T_MeV, muB_MeV, xi, ds)
        end
    end

    build_summary(OUT_DETAIL, OUT_SUMMARY)
    println("Wrote detail:  $OUT_DETAIL")
    println("Wrote summary: $OUT_SUMMARY")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
