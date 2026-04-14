#!/usr/bin/env julia

using Statistics: mean

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

const OUT_TRACE_CSV = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_term_switch_trace.csv"
const OUT_SUMMARY_CSV = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_detm_im_upstream_term_switch_summary.csv"
const XI_LIST = [0.34, 0.36, 0.38]
const DS_POINTS = [1.0e-3, 5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0]

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function phase_tag(ds::Float64)
    if ds <= 5.0
        return "decrease"
    elseif ds <= 10.0
        return "recover"
    end
    return "grow"
end

function dominant_label(c1::Float64, c2::Float64, c3::Float64)
    a1 = abs(c1)
    a2 = abs(c2)
    a3 = abs(c3)
    if a1 >= a2 && a1 >= a3
        return "ReM00*ImM88"
    elseif a2 >= a1 && a2 >= a3
        return "ImM00*ReM88"
    end
    return "-2ReM08*ImM08"
end

function main()
    process = :uubar_to_uubar
    T_MeV = 200.0
    muB_MeV = 0.0

    ensure_parent_dir(OUT_TRACE_CSV)
    ensure_parent_dir(OUT_SUMMARY_CSV)

    rows = NamedTuple[]

    io = open(OUT_TRACE_CSV, "w")
    try
        println(
            io,
            "scenario,process,xi,phase,ds,s_value,k0,k_norm,ReM00,ImM00,ReM08,ImM08,ReM88,ImM88,c1_ReM00_ImM88,c2_ImM00_ReM88,c3_minus2ReM08ImM08,im_M00M88,im_M08sq,detM_im,detM_im_rebuild,rebuild_error,dom_term,imM00_from_Piuu,imM00_from_Piss,imM08_from_Piuu,imM08_from_Piss,imM88_from_Piuu,imM88_from_Piss",
        )

        for xi in XI_LIST
            st = build_state_point(T_MeV, muB_MeV, xi)
            th = process_threshold_info(process, st.quark_params)

            for ds in DS_POINTS
                s = th.s_th + ds
                tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
                t = 0.5 * (tb.t_min + tb.t_max)

                decomp = decompose_mixed_p_propagator_chain(
                    process,
                    s,
                    t,
                    st.quark_params,
                    st.thermo_params,
                    st.K_coeffs,
                )

                M00 = decomp.M00
                M08 = decomp.M08
                M88 = decomp.M88

                c1 = real(M00) * imag(M88)
                c2 = imag(M00) * real(M88)
                c3 = -2.0 * real(M08) * imag(M08)
                rebuild = c1 + c2 + c3
                detM_im = imag(decomp.detM)
                err = rebuild - detM_im

                dom = dominant_label(c1, c2, c3)
                ptag = phase_tag(ds)

                push!(
                    rows,
                    (
                        xi = xi,
                        phase = ptag,
                        ds = ds,
                        detM_im = detM_im,
                        c1 = c1,
                        c2 = c2,
                        c3 = c3,
                        dom = dom,
                    ),
                )

                println(
                    io,
                    join(
                        (
                            "tauu_pos_uubaruubar",
                            process,
                            xi,
                            ptag,
                            ds,
                            s,
                            decomp.k0,
                            decomp.k_norm,
                            real(M00),
                            imag(M00),
                            real(M08),
                            imag(M08),
                            real(M88),
                            imag(M88),
                            c1,
                            c2,
                            c3,
                            imag(decomp.detM_term_M00M88),
                            imag(decomp.detM_term_M08sq),
                            detM_im,
                            rebuild,
                            err,
                            dom,
                            imag(decomp.M00_from_Piuu),
                            imag(decomp.M00_from_Piss),
                            imag(decomp.M08_from_Piuu),
                            imag(decomp.M08_from_Piss),
                            imag(decomp.M88_from_Piuu),
                            imag(decomp.M88_from_Piss),
                        ),
                        ',',
                    ),
                )
            end
        end
    finally
        close(io)
    end

    io_sum = open(OUT_SUMMARY_CSV, "w")
    try
        println(
            io_sum,
            "scenario,xi,phase,detM_im_mean,detM_im_min,detM_im_max,c1_mean,c2_mean,c3_mean,dom_term_meanabs,ds_at_min_detM_im,detM_im_at_min,c1_at_min,c2_at_min,c3_at_min,dom_term_at_min",
        )

        for xi in XI_LIST
            sub_xi = [r for r in rows if r.xi == xi]
            min_row = sub_xi[argmin(map(r -> r.detM_im, sub_xi))]

            for ph in ("decrease", "recover", "grow")
                sub = [r for r in sub_xi if r.phase == ph]
                isempty(sub) && continue
                det_vals = map(r -> r.detM_im, sub)
                c1_vals = map(r -> r.c1, sub)
                c2_vals = map(r -> r.c2, sub)
                c3_vals = map(r -> r.c3, sub)

                c1_abs = mean(abs.(c1_vals))
                c2_abs = mean(abs.(c2_vals))
                c3_abs = mean(abs.(c3_vals))
                dom_mean = dominant_label(c1_abs, c2_abs, c3_abs)

                println(
                    io_sum,
                    join(
                        (
                            "tauu_pos_uubaruubar",
                            xi,
                            ph,
                            mean(det_vals),
                            minimum(det_vals),
                            maximum(det_vals),
                            mean(c1_vals),
                            mean(c2_vals),
                            mean(c3_vals),
                            dom_mean,
                            min_row.ds,
                            min_row.detM_im,
                            min_row.c1,
                            min_row.c2,
                            min_row.c3,
                            min_row.dom,
                        ),
                        ',',
                    ),
                )
            end
        end
    finally
        close(io_sum)
    end

    println("Wrote upstream term-switch trace to: " * OUT_TRACE_CSV)
    println("Wrote upstream term-switch summary to: " * OUT_SUMMARY_CSV)
end

main()
