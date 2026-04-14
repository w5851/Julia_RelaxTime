#!/usr/bin/env julia

using Statistics: mean

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

    out = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_sigma_amplitude_decomposition.csv"
    out_summary = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_sigma_amplitude_decomposition_summary.csv"
    ensure_parent_dir(out)
    ensure_parent_dir(out_summary)

    stA = build_state_point(T_MeV, muB_MeV, xi_A)
    stB = build_state_point(T_MeV, muB_MeV, xi_B)

    io = open(out, "w")
    io_sum = open(out_summary, "w")
    try
        println(io, "process,xi_state,s,s_minus_sth,t_mid,M2_total,M_s_sq,M_t_sq,M_interf,M_s_S,M_s_P,M_t_S,M_t_P,M2_total_Kswap")
        println(io_sum, "process,delta_s_max,area_M2_A,area_M2_B,ratio_B_over_A,area_s_sq_A,area_s_sq_B,area_t_sq_A,area_t_sq_B,area_interf_A,area_interf_B,s_channel_share_A,s_channel_share_B,t_channel_share_A,t_channel_share_B,interf_share_A,interf_share_B,k_swap_rel_A,k_swap_rel_B,t_mid_mean_A,t_mid_mean_B")

        ds_vals = collect(range(1e-6, delta_s_max; length=n_s))

        for process in processes
            thA = process_threshold_info(process, stA.quark_params)
            thB = process_threshold_info(process, stB.quark_params)

            M2A = Float64[]; M2B = Float64[]
            SsA = Float64[]; SsB = Float64[]
            TsA = Float64[]; TsB = Float64[]
            IfA = Float64[]; IfB = Float64[]
            KswA = Float64[]; KswB = Float64[]
            tmA = Float64[]; tmB = Float64[]

            for ds in ds_vals
                sA = thA.s_th + ds
                sB = thB.s_th + ds
                tbA = Main.TotalCrossSection.calculate_t_bounds(sA, thA.mi, thA.mj, thA.mc, thA.md)
                tbB = Main.TotalCrossSection.calculate_t_bounds(sB, thB.mi, thB.mj, thB.mc, thB.md)
                tA = 0.5 * (tbA.t_min + tbA.t_max)
                tB = 0.5 * (tbB.t_min + tbB.t_max)

                pA = decompose_qqbar_amplitude_terms(process, sA, tA, stA.quark_params, stA.thermo_params, stA.K_coeffs)
                pB = decompose_qqbar_amplitude_terms(process, sB, tB, stB.quark_params, stB.thermo_params, stB.K_coeffs)
                pA_k = decompose_qqbar_amplitude_terms(process, sA, tA, stA.quark_params, stA.thermo_params, stB.K_coeffs)
                pB_k = decompose_qqbar_amplitude_terms(process, sB, tB, stB.quark_params, stB.thermo_params, stA.K_coeffs)

                push!(M2A, pA.M2_total); push!(M2B, pB.M2_total)
                push!(SsA, pA.M_s_sq); push!(SsB, pB.M_s_sq)
                push!(TsA, pA.M_t_sq); push!(TsB, pB.M_t_sq)
                push!(IfA, pA.M_interf); push!(IfB, pB.M_interf)
                push!(KswA, pA_k.M2_total); push!(KswB, pB_k.M2_total)
                push!(tmA, tA); push!(tmB, tB)

                println(io, join((string(process), string(xi_A), sA, ds, tA, pA.M2_total, pA.M_s_sq, pA.M_t_sq, pA.M_interf, pA.M_s_S, pA.M_s_P, pA.M_t_S, pA.M_t_P, pA_k.M2_total), ','))
                println(io, join((string(process), string(xi_B), sB, ds, tB, pB.M2_total, pB.M_s_sq, pB.M_t_sq, pB.M_interf, pB.M_s_S, pB.M_s_P, pB.M_t_S, pB.M_t_P, pB_k.M2_total), ','))
            end

            area_M2_A = trapz(ds_vals, M2A)
            area_M2_B = trapz(ds_vals, M2B)
            area_S_A = trapz(ds_vals, SsA)
            area_S_B = trapz(ds_vals, SsB)
            area_T_A = trapz(ds_vals, TsA)
            area_T_B = trapz(ds_vals, TsB)
            area_I_A = trapz(ds_vals, IfA)
            area_I_B = trapz(ds_vals, IfB)
            area_K_A = trapz(ds_vals, KswA)
            area_K_B = trapz(ds_vals, KswB)

            ratio = area_M2_A > 0 ? area_M2_B / area_M2_A : NaN
            s_share_A = area_M2_A != 0 ? area_S_A / area_M2_A : NaN
            s_share_B = area_M2_B != 0 ? area_S_B / area_M2_B : NaN
            t_share_A = area_M2_A != 0 ? area_T_A / area_M2_A : NaN
            t_share_B = area_M2_B != 0 ? area_T_B / area_M2_B : NaN
            i_share_A = area_M2_A != 0 ? area_I_A / area_M2_A : NaN
            i_share_B = area_M2_B != 0 ? area_I_B / area_M2_B : NaN
            k_rel_A = area_M2_A > 0 ? (area_K_A / area_M2_A - 1.0) : NaN
            k_rel_B = area_M2_B > 0 ? (area_K_B / area_M2_B - 1.0) : NaN

            println(io_sum, join((
                string(process),
                delta_s_max,
                area_M2_A,
                area_M2_B,
                ratio,
                area_S_A,
                area_S_B,
                area_T_A,
                area_T_B,
                area_I_A,
                area_I_B,
                s_share_A,
                s_share_B,
                t_share_A,
                t_share_B,
                i_share_A,
                i_share_B,
                k_rel_A,
                k_rel_B,
                mean(tmA),
                mean(tmB),
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
