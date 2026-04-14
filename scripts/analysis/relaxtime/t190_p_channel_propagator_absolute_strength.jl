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

    out = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_p_channel_propagator_absolute_strength.csv"
    out_sum = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_p_channel_propagator_absolute_strength_summary.csv"
    ensure_parent_dir(out)
    ensure_parent_dir(out_sum)

    stA = build_state_point(T_MeV, muB_MeV, xi_A)
    stB = build_state_point(T_MeV, muB_MeV, xi_B)
    ds_vals = collect(range(1e-6, delta_s_max; length=n_s))

    io = open(out, "w")
    io_sum = open(out_sum, "w")
    try
        println(io, "process,xi_state,s,s_minus_sth,t_mid,k0_s,k_s,abs_D_s_P_sq_simple,abs_D_s_P_sq_mixed,abs_D_s_P_sq_total,kin_s_P,M_s_P_total")
        println(io_sum, "process,delta_s_max,area_absDsP_simple_A,area_absDsP_simple_B,ratio_simple_B_over_A,area_absDsP_mixed_A,area_absDsP_mixed_B,ratio_mixed_B_over_A,area_absDsP_total_A,area_absDsP_total_B,ratio_total_B_over_A,mixed_share_A,mixed_share_B,area_kin_s_P_A,area_kin_s_P_B,ratio_kin_B_over_A,area_MsP_A,area_MsP_B,ratio_MsP_B_over_A")

        for process in processes
            thA = process_threshold_info(process, stA.quark_params)
            thB = process_threshold_info(process, stB.quark_params)

            dps_simple_A = Float64[]
            dps_simple_B = Float64[]
            dps_mixed_A = Float64[]
            dps_mixed_B = Float64[]
            dps_total_A = Float64[]
            dps_total_B = Float64[]
            kin_A = Float64[]
            kin_B = Float64[]
            msp_A = Float64[]
            msp_B = Float64[]

            for ds in ds_vals
                sA = thA.s_th + ds
                sB = thB.s_th + ds
                tbA = Main.TotalCrossSection.calculate_t_bounds(sA, thA.mi, thA.mj, thA.mc, thA.md)
                tbB = Main.TotalCrossSection.calculate_t_bounds(sB, thB.mi, thB.mj, thB.mc, thB.md)
                tA = 0.5 * (tbA.t_min + tbA.t_max)
                tB = 0.5 * (tbB.t_min + tbB.t_max)

                pA = decompose_p_channel_propagator_strength(process, sA, tA, stA.quark_params, stA.thermo_params, stA.K_coeffs)
                pB = decompose_p_channel_propagator_strength(process, sB, tB, stB.quark_params, stB.thermo_params, stB.K_coeffs)

                push!(dps_simple_A, pA.abs_D_s_P_sq_simple)
                push!(dps_simple_B, pB.abs_D_s_P_sq_simple)
                push!(dps_mixed_A, pA.abs_D_s_P_sq_mixed)
                push!(dps_mixed_B, pB.abs_D_s_P_sq_mixed)
                push!(dps_total_A, pA.abs_D_s_P_sq_total)
                push!(dps_total_B, pB.abs_D_s_P_sq_total)
                push!(kin_A, pA.kin_s_P)
                push!(kin_B, pB.kin_s_P)
                push!(msp_A, pA.M_s_P_total)
                push!(msp_B, pB.M_s_P_total)

                println(io, join((string(process), string(xi_A), sA, ds, tA, pA.k0, pA.k_norm, pA.abs_D_s_P_sq_simple, pA.abs_D_s_P_sq_mixed, pA.abs_D_s_P_sq_total, pA.kin_s_P, pA.M_s_P_total), ','))
                println(io, join((string(process), string(xi_B), sB, ds, tB, pB.k0, pB.k_norm, pB.abs_D_s_P_sq_simple, pB.abs_D_s_P_sq_mixed, pB.abs_D_s_P_sq_total, pB.kin_s_P, pB.M_s_P_total), ','))
            end

            area_simple_A = trapz(ds_vals, dps_simple_A)
            area_simple_B = trapz(ds_vals, dps_simple_B)
            area_mixed_A = trapz(ds_vals, dps_mixed_A)
            area_mixed_B = trapz(ds_vals, dps_mixed_B)
            area_total_A = trapz(ds_vals, dps_total_A)
            area_total_B = trapz(ds_vals, dps_total_B)
            area_kin_A = trapz(ds_vals, kin_A)
            area_kin_B = trapz(ds_vals, kin_B)
            area_msp_A = trapz(ds_vals, msp_A)
            area_msp_B = trapz(ds_vals, msp_B)

            ratio_simple = area_simple_A > 0 ? area_simple_B / area_simple_A : NaN
            ratio_mixed = area_mixed_A > 0 ? area_mixed_B / area_mixed_A : NaN
            ratio_total = area_total_A > 0 ? area_total_B / area_total_A : NaN
            mixed_share_A = area_total_A > 0 ? area_mixed_A / area_total_A : NaN
            mixed_share_B = area_total_B > 0 ? area_mixed_B / area_total_B : NaN
            ratio_kin = area_kin_A != 0 ? area_kin_B / area_kin_A : NaN
            ratio_msp = area_msp_A > 0 ? area_msp_B / area_msp_A : NaN

            println(io_sum, join((
                string(process),
                delta_s_max,
                area_simple_A,
                area_simple_B,
                ratio_simple,
                area_mixed_A,
                area_mixed_B,
                ratio_mixed,
                area_total_A,
                area_total_B,
                ratio_total,
                mixed_share_A,
                mixed_share_B,
                area_kin_A,
                area_kin_B,
                ratio_kin,
                area_msp_A,
                area_msp_B,
                ratio_msp,
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
