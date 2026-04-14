#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

function trapz(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == length(y) || error("trapz: length mismatch")
    n <= 1 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        acc += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return acc
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function main()
    T_MeV = 190.0
    muB_MeV = 0.0
    xi_A = -0.10
    xi_B = -0.08
    processes = Symbol[:uubar_to_ddbar, :uubar_to_uubar]
    delta_s_max = 2.0
    n_s = 36

    out = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_mixed_p_propagator_chain_decomposition.csv"
    out_sum = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_mixed_p_propagator_chain_decomposition_summary.csv"
    ensure_parent_dir(out)
    ensure_parent_dir(out_sum)

    stA = build_state_point(T_MeV, muB_MeV, xi_A)
    stB = build_state_point(T_MeV, muB_MeV, xi_B)
    ds_vals = collect(range(1e-6, delta_s_max; length=n_s))

    io = open(out, "w")
    io_sum = open(out_sum, "w")
    try
        println(io, "process,xi_state,s,s_minus_sth,t_mid,k0_s,k_s,detK_plus,abs_detM_sq,abs_detM_term_M00M88_sq,abs_detM_term_M08sq_sq,detM_cross_term,real_detM_term_M00M88,imag_detM_term_M00M88,real_detM_term_M08sq,imag_detM_term_M08sq,abs_M00_from_K0_sq,abs_M00_from_Piuu_sq,abs_M00_from_Piss_sq,abs_M08_from_K08_sq,abs_M08_from_Piuu_sq,abs_M08_from_Piss_sq,abs_M88_from_K8_sq,abs_M88_from_Piuu_sq,abs_M88_from_Piss_sq,abs_JMJ_sq,abs_D_mixed_P_sq,abs_Piuu_sq,abs_Piss_sq,abs_M00_sq,abs_M08_sq,abs_M88_sq")
        println(io_sum, "process,delta_s_max,area_abs_detM_sq_A,area_abs_detM_sq_B,ratio_detM_B_over_A,area_abs_detM_term_M00M88_sq_A,area_abs_detM_term_M00M88_sq_B,ratio_detM_term_M00M88_B_over_A,area_abs_detM_term_M08sq_sq_A,area_abs_detM_term_M08sq_sq_B,ratio_detM_term_M08sq_B_over_A,area_detM_cross_term_A,area_detM_cross_term_B,ratio_detM_cross_term_B_over_A,area_abs_M00_from_K0_sq_A,area_abs_M00_from_K0_sq_B,ratio_M00_from_K0_B_over_A,area_abs_M00_from_Piuu_sq_A,area_abs_M00_from_Piuu_sq_B,ratio_M00_from_Piuu_B_over_A,area_abs_M00_from_Piss_sq_A,area_abs_M00_from_Piss_sq_B,ratio_M00_from_Piss_B_over_A,area_abs_M08_from_K08_sq_A,area_abs_M08_from_K08_sq_B,ratio_M08_from_K08_B_over_A,area_abs_M08_from_Piuu_sq_A,area_abs_M08_from_Piuu_sq_B,ratio_M08_from_Piuu_B_over_A,area_abs_M08_from_Piss_sq_A,area_abs_M08_from_Piss_sq_B,ratio_M08_from_Piss_B_over_A,area_abs_M88_from_K8_sq_A,area_abs_M88_from_K8_sq_B,ratio_M88_from_K8_B_over_A,area_abs_M88_from_Piuu_sq_A,area_abs_M88_from_Piuu_sq_B,ratio_M88_from_Piuu_B_over_A,area_abs_M88_from_Piss_sq_A,area_abs_M88_from_Piss_sq_B,ratio_M88_from_Piss_B_over_A,area_abs_JMJ_sq_A,area_abs_JMJ_sq_B,ratio_JMJ_B_over_A,area_abs_Dmixed_sq_A,area_abs_Dmixed_sq_B,ratio_Dmixed_B_over_A,area_detK_plus_A,area_detK_plus_B,ratio_detK_B_over_A,area_abs_Piuu_sq_A,area_abs_Piuu_sq_B,ratio_Piuu_B_over_A,area_abs_Piss_sq_A,area_abs_Piss_sq_B,ratio_Piss_B_over_A")

        for process in processes
            thA = process_threshold_info(process, stA.quark_params)
            thB = process_threshold_info(process, stB.quark_params)

            detM_A = Float64[]
            detM_B = Float64[]
            JMJ_A = Float64[]
            JMJ_B = Float64[]
            Dmix_A = Float64[]
            Dmix_B = Float64[]
            detK_A = Float64[]
            detK_B = Float64[]
            Piuu_A = Float64[]
            Piuu_B = Float64[]
            Piss_A = Float64[]
            Piss_B = Float64[]
            detM_term_M00M88_A = Float64[]
            detM_term_M00M88_B = Float64[]
            detM_term_M08sq_A = Float64[]
            detM_term_M08sq_B = Float64[]
            detM_cross_A = Float64[]
            detM_cross_B = Float64[]
            M00_K0_A = Float64[]
            M00_K0_B = Float64[]
            M00_Piuu_A = Float64[]
            M00_Piuu_B = Float64[]
            M00_Piss_A = Float64[]
            M00_Piss_B = Float64[]
            M08_K08_A = Float64[]
            M08_K08_B = Float64[]
            M08_Piuu_A = Float64[]
            M08_Piuu_B = Float64[]
            M08_Piss_A = Float64[]
            M08_Piss_B = Float64[]
            M88_K8_A = Float64[]
            M88_K8_B = Float64[]
            M88_Piuu_A = Float64[]
            M88_Piuu_B = Float64[]
            M88_Piss_A = Float64[]
            M88_Piss_B = Float64[]

            for ds in ds_vals
                sA = thA.s_th + ds
                sB = thB.s_th + ds
                tbA = Main.TotalCrossSection.calculate_t_bounds(sA, thA.mi, thA.mj, thA.mc, thA.md)
                tbB = Main.TotalCrossSection.calculate_t_bounds(sB, thB.mi, thB.mj, thB.mc, thB.md)
                tA = 0.5 * (tbA.t_min + tbA.t_max)
                tB = 0.5 * (tbB.t_min + tbB.t_max)

                pA = decompose_mixed_p_propagator_chain(process, sA, tA, stA.quark_params, stA.thermo_params, stA.K_coeffs)
                pB = decompose_mixed_p_propagator_chain(process, sB, tB, stB.quark_params, stB.thermo_params, stB.K_coeffs)

                push!(detM_A, pA.abs_detM_sq)
                push!(detM_B, pB.abs_detM_sq)
                push!(JMJ_A, pA.abs_JMJ_sq)
                push!(JMJ_B, pB.abs_JMJ_sq)
                push!(Dmix_A, pA.abs_D_mixed_P_sq)
                push!(Dmix_B, pB.abs_D_mixed_P_sq)
                push!(detK_A, pA.detK_plus)
                push!(detK_B, pB.detK_plus)
                push!(Piuu_A, abs2(pA.Π_uu))
                push!(Piuu_B, abs2(pB.Π_uu))
                push!(Piss_A, abs2(pA.Π_ss))
                push!(Piss_B, abs2(pB.Π_ss))
                push!(detM_term_M00M88_A, pA.abs_detM_term_M00M88_sq)
                push!(detM_term_M00M88_B, pB.abs_detM_term_M00M88_sq)
                push!(detM_term_M08sq_A, pA.abs_detM_term_M08sq_sq)
                push!(detM_term_M08sq_B, pB.abs_detM_term_M08sq_sq)
                push!(detM_cross_A, pA.detM_cross_term)
                push!(detM_cross_B, pB.detM_cross_term)
                push!(M00_K0_A, abs2(pA.M00_from_K0))
                push!(M00_K0_B, abs2(pB.M00_from_K0))
                push!(M00_Piuu_A, abs2(pA.M00_from_Piuu))
                push!(M00_Piuu_B, abs2(pB.M00_from_Piuu))
                push!(M00_Piss_A, abs2(pA.M00_from_Piss))
                push!(M00_Piss_B, abs2(pB.M00_from_Piss))
                push!(M08_K08_A, abs2(pA.M08_from_K08))
                push!(M08_K08_B, abs2(pB.M08_from_K08))
                push!(M08_Piuu_A, abs2(pA.M08_from_Piuu))
                push!(M08_Piuu_B, abs2(pB.M08_from_Piuu))
                push!(M08_Piss_A, abs2(pA.M08_from_Piss))
                push!(M08_Piss_B, abs2(pB.M08_from_Piss))
                push!(M88_K8_A, abs2(pA.M88_from_K8))
                push!(M88_K8_B, abs2(pB.M88_from_K8))
                push!(M88_Piuu_A, abs2(pA.M88_from_Piuu))
                push!(M88_Piuu_B, abs2(pB.M88_from_Piuu))
                push!(M88_Piss_A, abs2(pA.M88_from_Piss))
                push!(M88_Piss_B, abs2(pB.M88_from_Piss))

                println(io, join((string(process), string(xi_A), sA, ds, tA, pA.k0, pA.k_norm, pA.detK_plus, pA.abs_detM_sq, pA.abs_detM_term_M00M88_sq, pA.abs_detM_term_M08sq_sq, pA.detM_cross_term, real(pA.detM_term_M00M88), imag(pA.detM_term_M00M88), real(pA.detM_term_M08sq), imag(pA.detM_term_M08sq), abs2(pA.M00_from_K0), abs2(pA.M00_from_Piuu), abs2(pA.M00_from_Piss), abs2(pA.M08_from_K08), abs2(pA.M08_from_Piuu), abs2(pA.M08_from_Piss), abs2(pA.M88_from_K8), abs2(pA.M88_from_Piuu), abs2(pA.M88_from_Piss), pA.abs_JMJ_sq, pA.abs_D_mixed_P_sq, abs2(pA.Π_uu), abs2(pA.Π_ss), abs2(pA.M00), abs2(pA.M08), abs2(pA.M88)), ','))
                println(io, join((string(process), string(xi_B), sB, ds, tB, pB.k0, pB.k_norm, pB.detK_plus, pB.abs_detM_sq, pB.abs_detM_term_M00M88_sq, pB.abs_detM_term_M08sq_sq, pB.detM_cross_term, real(pB.detM_term_M00M88), imag(pB.detM_term_M00M88), real(pB.detM_term_M08sq), imag(pB.detM_term_M08sq), abs2(pB.M00_from_K0), abs2(pB.M00_from_Piuu), abs2(pB.M00_from_Piss), abs2(pB.M08_from_K08), abs2(pB.M08_from_Piuu), abs2(pB.M08_from_Piss), abs2(pB.M88_from_K8), abs2(pB.M88_from_Piuu), abs2(pB.M88_from_Piss), pB.abs_JMJ_sq, pB.abs_D_mixed_P_sq, abs2(pB.Π_uu), abs2(pB.Π_ss), abs2(pB.M00), abs2(pB.M08), abs2(pB.M88)), ','))
            end

            area_detM_A = trapz(ds_vals, detM_A)
            area_detM_B = trapz(ds_vals, detM_B)
            area_JMJ_A = trapz(ds_vals, JMJ_A)
            area_JMJ_B = trapz(ds_vals, JMJ_B)
            area_Dmix_A = trapz(ds_vals, Dmix_A)
            area_Dmix_B = trapz(ds_vals, Dmix_B)
            area_detK_A = trapz(ds_vals, detK_A)
            area_detK_B = trapz(ds_vals, detK_B)
            area_Piuu_A = trapz(ds_vals, Piuu_A)
            area_Piuu_B = trapz(ds_vals, Piuu_B)
            area_Piss_A = trapz(ds_vals, Piss_A)
            area_Piss_B = trapz(ds_vals, Piss_B)
            area_detM_term_M00M88_A = trapz(ds_vals, detM_term_M00M88_A)
            area_detM_term_M00M88_B = trapz(ds_vals, detM_term_M00M88_B)
            area_detM_term_M08sq_A = trapz(ds_vals, detM_term_M08sq_A)
            area_detM_term_M08sq_B = trapz(ds_vals, detM_term_M08sq_B)
            area_detM_cross_A = trapz(ds_vals, detM_cross_A)
            area_detM_cross_B = trapz(ds_vals, detM_cross_B)
            area_M00_K0_A = trapz(ds_vals, M00_K0_A)
            area_M00_K0_B = trapz(ds_vals, M00_K0_B)
            area_M00_Piuu_A = trapz(ds_vals, M00_Piuu_A)
            area_M00_Piuu_B = trapz(ds_vals, M00_Piuu_B)
            area_M00_Piss_A = trapz(ds_vals, M00_Piss_A)
            area_M00_Piss_B = trapz(ds_vals, M00_Piss_B)
            area_M08_K08_A = trapz(ds_vals, M08_K08_A)
            area_M08_K08_B = trapz(ds_vals, M08_K08_B)
            area_M08_Piuu_A = trapz(ds_vals, M08_Piuu_A)
            area_M08_Piuu_B = trapz(ds_vals, M08_Piuu_B)
            area_M08_Piss_A = trapz(ds_vals, M08_Piss_A)
            area_M08_Piss_B = trapz(ds_vals, M08_Piss_B)
            area_M88_K8_A = trapz(ds_vals, M88_K8_A)
            area_M88_K8_B = trapz(ds_vals, M88_K8_B)
            area_M88_Piuu_A = trapz(ds_vals, M88_Piuu_A)
            area_M88_Piuu_B = trapz(ds_vals, M88_Piuu_B)
            area_M88_Piss_A = trapz(ds_vals, M88_Piss_A)
            area_M88_Piss_B = trapz(ds_vals, M88_Piss_B)

            ratio_detM = area_detM_A > 0 ? area_detM_B / area_detM_A : NaN
            ratio_JMJ = area_JMJ_A > 0 ? area_JMJ_B / area_JMJ_A : NaN
            ratio_Dmix = area_Dmix_A > 0 ? area_Dmix_B / area_Dmix_A : NaN
            ratio_detK = area_detK_A != 0 ? area_detK_B / area_detK_A : NaN
            ratio_Piuu = area_Piuu_A > 0 ? area_Piuu_B / area_Piuu_A : NaN
            ratio_Piss = area_Piss_A > 0 ? area_Piss_B / area_Piss_A : NaN
            ratio_detM_term_M00M88 = area_detM_term_M00M88_A > 0 ? area_detM_term_M00M88_B / area_detM_term_M00M88_A : NaN
            ratio_detM_term_M08sq = area_detM_term_M08sq_A > 0 ? area_detM_term_M08sq_B / area_detM_term_M08sq_A : NaN
            ratio_detM_cross = area_detM_cross_A != 0 ? area_detM_cross_B / area_detM_cross_A : NaN
            ratio_M00_K0 = area_M00_K0_A > 0 ? area_M00_K0_B / area_M00_K0_A : NaN
            ratio_M00_Piuu = area_M00_Piuu_A > 0 ? area_M00_Piuu_B / area_M00_Piuu_A : NaN
            ratio_M00_Piss = area_M00_Piss_A > 0 ? area_M00_Piss_B / area_M00_Piss_A : NaN
            ratio_M08_K08 = area_M08_K08_A > 0 ? area_M08_K08_B / area_M08_K08_A : NaN
            ratio_M08_Piuu = area_M08_Piuu_A > 0 ? area_M08_Piuu_B / area_M08_Piuu_A : NaN
            ratio_M08_Piss = area_M08_Piss_A > 0 ? area_M08_Piss_B / area_M08_Piss_A : NaN
            ratio_M88_K8 = area_M88_K8_A > 0 ? area_M88_K8_B / area_M88_K8_A : NaN
            ratio_M88_Piuu = area_M88_Piuu_A > 0 ? area_M88_Piuu_B / area_M88_Piuu_A : NaN
            ratio_M88_Piss = area_M88_Piss_A > 0 ? area_M88_Piss_B / area_M88_Piss_A : NaN

            println(io_sum, join((
                string(process),
                delta_s_max,
                area_detM_A,
                area_detM_B,
                ratio_detM,
                area_detM_term_M00M88_A,
                area_detM_term_M00M88_B,
                ratio_detM_term_M00M88,
                area_detM_term_M08sq_A,
                area_detM_term_M08sq_B,
                ratio_detM_term_M08sq,
                area_detM_cross_A,
                area_detM_cross_B,
                ratio_detM_cross,
                area_M00_K0_A,
                area_M00_K0_B,
                ratio_M00_K0,
                area_M00_Piuu_A,
                area_M00_Piuu_B,
                ratio_M00_Piuu,
                area_M00_Piss_A,
                area_M00_Piss_B,
                ratio_M00_Piss,
                area_M08_K08_A,
                area_M08_K08_B,
                ratio_M08_K08,
                area_M08_Piuu_A,
                area_M08_Piuu_B,
                ratio_M08_Piuu,
                area_M08_Piss_A,
                area_M08_Piss_B,
                ratio_M08_Piss,
                area_M88_K8_A,
                area_M88_K8_B,
                ratio_M88_K8,
                area_M88_Piuu_A,
                area_M88_Piuu_B,
                ratio_M88_Piuu,
                area_M88_Piss_A,
                area_M88_Piss_B,
                ratio_M88_Piss,
                area_JMJ_A,
                area_JMJ_B,
                ratio_JMJ,
                area_Dmix_A,
                area_Dmix_B,
                ratio_Dmix,
                area_detK_A,
                area_detK_B,
                ratio_detK,
                area_Piuu_A,
                area_Piuu_B,
                ratio_Piuu,
                area_Piss_A,
                area_Piss_B,
                ratio_Piss,
            ), ','))
            flush(io)
            flush(io_sum)
        end
    finally
        close(io)
        close(io_sum)
    end
end

main()
