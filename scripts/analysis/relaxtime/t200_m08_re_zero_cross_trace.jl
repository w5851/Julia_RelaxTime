#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

const OUT_TRACE_CSV = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_m08_re_zero_cross_trace.csv"
const OUT_SUMMARY_CSV = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_m08_re_zero_cross_summary.csv"
const XI_LIST = [0.34, 0.36, 0.38]
const DS_SCAN = vcat(collect(range(1.0e-3, 0.2; length=36)), collect(range(0.25, 5.0; length=40)), collect(range(5.5, 20.0; length=60)))

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function locate_zero_cross(ds_vals::Vector{Float64}, y_vals::Vector{Float64})
    for i in 2:length(ds_vals)
        y1 = y_vals[i - 1]
        y2 = y_vals[i]
        if y1 == 0.0
            return ds_vals[i - 1]
        end
        if signbit(y1) != signbit(y2)
            t = y1 / (y1 - y2)
            return ds_vals[i - 1] + t * (ds_vals[i] - ds_vals[i - 1])
        end
    end
    return NaN
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
            "scenario,process,xi,ds,s_value,k0,k_norm,K08_plus,detK_plus,c_offdiag,re_Piuu,re_Piss,re_Piuu_minus_Piss,re_Piuu_minus_Piss_target,re_M08_from_K08,re_M08_from_Piuu,re_M08_from_Piss,re_M08_total,re_M08_formula,residual_formula",
        )

        for xi in XI_LIST
            st = build_state_point(T_MeV, muB_MeV, xi)
            th = process_threshold_info(process, st.quark_params)
            c_offdiag = (4.0 / 3.0) * sqrt(2.0) * st.K_coeffs.det_K_plus
            target_delta = -st.K_coeffs.K08_plus / c_offdiag

            for ds in DS_SCAN
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

                re_pi_uu = real(decomp.Π_uu)
                re_pi_ss = real(decomp.Π_ss)
                delta_re_pi = re_pi_uu - re_pi_ss

                re_m08_total = real(decomp.M08)
                re_m08_formula = st.K_coeffs.K08_plus + c_offdiag * delta_re_pi
                resid = re_m08_total - re_m08_formula

                push!(
                    rows,
                    (
                        xi = xi,
                        ds = ds,
                        re_m08_total = re_m08_total,
                        delta_re_pi = delta_re_pi,
                        target_delta = target_delta,
                    ),
                )

                println(
                    io,
                    join(
                        (
                            "tauu_pos_uubaruubar",
                            process,
                            xi,
                            ds,
                            s,
                            decomp.k0,
                            decomp.k_norm,
                            st.K_coeffs.K08_plus,
                            st.K_coeffs.det_K_plus,
                            c_offdiag,
                            re_pi_uu,
                            re_pi_ss,
                            delta_re_pi,
                            target_delta,
                            real(decomp.M08_from_K08),
                            real(decomp.M08_from_Piuu),
                            real(decomp.M08_from_Piss),
                            re_m08_total,
                            re_m08_formula,
                            resid,
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
            "scenario,xi,ds_zero_cross_linear,reM08_min,reM08_max,deltaRePi_at_minAbsM08,target_deltaRePi,min_abs_M08,min_abs_ds",
        )
        for xi in XI_LIST
            sub = sort([r for r in rows if r.xi == xi], by = r -> r.ds)
            ds_vals = map(r -> r.ds, sub)
            m08_vals = map(r -> r.re_m08_total, sub)
            zero_ds = locate_zero_cross(ds_vals, m08_vals)
            min_idx = argmin(abs.(m08_vals))
            min_row = sub[min_idx]
            println(
                io_sum,
                join(
                    (
                        "tauu_pos_uubaruubar",
                        xi,
                        zero_ds,
                        minimum(m08_vals),
                        maximum(m08_vals),
                        min_row.delta_re_pi,
                        min_row.target_delta,
                        abs(min_row.re_m08_total),
                        min_row.ds,
                    ),
                    ',',
                ),
            )
        end
    finally
        close(io_sum)
    end

    println("Wrote M08 zero-cross trace to: " * OUT_TRACE_CSV)
    println("Wrote M08 zero-cross summary to: " * OUT_SUMMARY_CSV)
end

main()
