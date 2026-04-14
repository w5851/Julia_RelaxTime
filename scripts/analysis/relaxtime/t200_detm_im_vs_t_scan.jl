#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

const OUT_CSV = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_vs_t_scan.csv"
const XI_LIST = [0.34, 0.36, 0.38]
const DS_POINTS = [2.0e-2, 2.0e-1, 2.0, 10.0]
const T_ALPHA_POINTS = collect(range(0.0, 1.0; length=41))

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function main()
    process = :uubar_to_uubar
    T_MeV = 200.0
    muB_MeV = 0.0

    ensure_parent_dir(OUT_CSV)
    io = open(OUT_CSV, "w")
    try
        println(
            io,
            "scenario,process,xi,ds,t_alpha,s_value,t_min,t_max,t_value,k0_s,k_s,pi_uu_im,pi_ss_im,im_M00M88,im_M08sq,detM_im,invabs_detM,abs_Dmixed_sq",
        )

        for xi in XI_LIST
            st = build_state_point(T_MeV, muB_MeV, xi)
            th = process_threshold_info(process, st.quark_params)

            for ds in DS_POINTS
                s = th.s_th + ds
                tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)

                for t_alpha in T_ALPHA_POINTS
                    t = tb.t_min + t_alpha * (tb.t_max - tb.t_min)
                    cms_s = Main.TotalPropagator.calculate_cms_momentum(process, s, t, :s, st.quark_params)

                    Π_uu_re, Π_uu_im = Main.PolarizationCache.polarization_aniso_cached(
                        :P,
                        cms_s.k0,
                        cms_s.k,
                        st.quark_params.m.u,
                        st.quark_params.m.u,
                        st.quark_params.μ.u,
                        st.quark_params.μ.u,
                        st.thermo_params.T,
                        st.thermo_params.Φ,
                        st.thermo_params.Φbar,
                        st.thermo_params.ξ,
                        st.quark_params.A.u,
                        st.quark_params.A.u,
                        0,
                    )
                    Π_ss_re, Π_ss_im = Main.PolarizationCache.polarization_aniso_cached(
                        :P,
                        cms_s.k0,
                        cms_s.k,
                        st.quark_params.m.s,
                        st.quark_params.m.s,
                        st.quark_params.μ.s,
                        st.quark_params.μ.s,
                        st.thermo_params.T,
                        st.thermo_params.Φ,
                        st.thermo_params.Φbar,
                        st.thermo_params.ξ,
                        st.quark_params.A.s,
                        st.quark_params.A.s,
                        2,
                    )

                    Π_uu = ComplexF64(Π_uu_re, Π_uu_im)
                    Π_ss = ComplexF64(Π_ss_re, Π_ss_im)
                    M00, M08, M88 = Main.MesonPropagator.calculate_coupling_elements(Π_uu, Π_ss, st.K_coeffs, :P)

                    term_a = M00 * M88
                    term_b = M08 * M08
                    detM = term_a - term_b

                    det_K = st.K_coeffs.det_K_plus
                    q1, q2, q3, q4 = Main.TotalPropagator.parse_scattering_process(process)
                    D_mixed = Main.MesonPropagator.meson_propagator_mixed(det_K, M00, M08, M88, q1, q2, q3, q4, :s)

                    println(
                        io,
                        join(
                            (
                                "tauu_pos_uubaruubar",
                                process,
                                xi,
                                ds,
                                t_alpha,
                                s,
                                tb.t_min,
                                tb.t_max,
                                t,
                                cms_s.k0,
                                cms_s.k,
                                Π_uu_im,
                                Π_ss_im,
                                imag(term_a),
                                imag(term_b),
                                imag(detM),
                                inv(abs2(detM)),
                                abs2(D_mixed),
                            ),
                            ',',
                        ),
                    )
                end
            end
        end
    finally
        close(io)
    end

    println("Wrote t-scan trace to: " * OUT_CSV)
end

main()
