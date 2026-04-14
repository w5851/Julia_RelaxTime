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

    out = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_s_channel_rootcause_decomposition.csv"
    out_sum = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_s_channel_rootcause_decomposition_summary.csv"
    ensure_parent_dir(out)
    ensure_parent_dir(out_sum)

    stA = build_state_point(T_MeV, muB_MeV, xi_A)
    stB = build_state_point(T_MeV, muB_MeV, xi_B)
    ds_vals = collect(range(1e-6, delta_s_max; length=n_s))

    io = open(out, "w")
    io_sum = open(out_sum, "w")
    try
        println(io, "process,xi_state,s,s_minus_sth,t_mid,abs_D_s_S_sq,abs_D_s_P_sq,kin_s_S,kin_s_P,M_s_S,M_s_P,M_s_sq")
        println(io_sum, "process,delta_s_max,area_Ms_A,area_Ms_B,ratio_B_over_A,area_absDsS_A,area_absDsS_B,area_absDsP_A,area_absDsP_B,area_kinS_A,area_kinS_B,area_kinP_A,area_kinP_B,share_S_A,share_S_B,share_P_A,share_P_B")

        for process in processes
            thA = process_threshold_info(process, stA.quark_params)
            thB = process_threshold_info(process, stB.quark_params)

            MsA = Float64[]; MsB = Float64[]
            dssA = Float64[]; dssB = Float64[]
            dspA = Float64[]; dspB = Float64[]
            kinSA = Float64[]; kinSB = Float64[]
            kinPA = Float64[]; kinPB = Float64[]
            mssA = Float64[]; mssB = Float64[]
            mspA = Float64[]; mspB = Float64[]

            for ds in ds_vals
                sA = thA.s_th + ds
                sB = thB.s_th + ds
                tbA = Main.TotalCrossSection.calculate_t_bounds(sA, thA.mi, thA.mj, thA.mc, thA.md)
                tbB = Main.TotalCrossSection.calculate_t_bounds(sB, thB.mi, thB.mj, thB.mc, thB.md)
                tA = 0.5 * (tbA.t_min + tbA.t_max)
                tB = 0.5 * (tbB.t_min + tbB.t_max)

                fa = decompose_qqbar_s_channel_factors(process, sA, tA, stA.quark_params, stA.thermo_params, stA.K_coeffs)
                fb = decompose_qqbar_s_channel_factors(process, sB, tB, stB.quark_params, stB.thermo_params, stB.K_coeffs)

                MsqA = fa.M_s_S + fa.M_s_P
                MsqB = fb.M_s_S + fb.M_s_P

                push!(MsA, MsqA); push!(MsB, MsqB)
                push!(dssA, fa.abs_D_s_S_sq); push!(dssB, fb.abs_D_s_S_sq)
                push!(dspA, fa.abs_D_s_P_sq); push!(dspB, fb.abs_D_s_P_sq)
                push!(kinSA, fa.kin_s_S); push!(kinSB, fb.kin_s_S)
                push!(kinPA, fa.kin_s_P); push!(kinPB, fb.kin_s_P)
                push!(mssA, fa.M_s_S); push!(mssB, fb.M_s_S)
                push!(mspA, fa.M_s_P); push!(mspB, fb.M_s_P)

                println(io, join((string(process), string(xi_A), sA, ds, tA, fa.abs_D_s_S_sq, fa.abs_D_s_P_sq, fa.kin_s_S, fa.kin_s_P, fa.M_s_S, fa.M_s_P, MsqA), ','))
                println(io, join((string(process), string(xi_B), sB, ds, tB, fb.abs_D_s_S_sq, fb.abs_D_s_P_sq, fb.kin_s_S, fb.kin_s_P, fb.M_s_S, fb.M_s_P, MsqB), ','))
            end

            area_Ms_A = trapz(ds_vals, MsA)
            area_Ms_B = trapz(ds_vals, MsB)
            area_absDsS_A = trapz(ds_vals, dssA)
            area_absDsS_B = trapz(ds_vals, dssB)
            area_absDsP_A = trapz(ds_vals, dspA)
            area_absDsP_B = trapz(ds_vals, dspB)
            area_kinS_A = trapz(ds_vals, kinSA)
            area_kinS_B = trapz(ds_vals, kinSB)
            area_kinP_A = trapz(ds_vals, kinPA)
            area_kinP_B = trapz(ds_vals, kinPB)
            area_MsS_A = trapz(ds_vals, mssA)
            area_MsS_B = trapz(ds_vals, mssB)
            area_MsP_A = trapz(ds_vals, mspA)
            area_MsP_B = trapz(ds_vals, mspB)

            ratio = area_Ms_A > 0 ? area_Ms_B / area_Ms_A : NaN
            shareS_A = area_Ms_A != 0 ? area_MsS_A / area_Ms_A : NaN
            shareS_B = area_Ms_B != 0 ? area_MsS_B / area_Ms_B : NaN
            shareP_A = area_Ms_A != 0 ? area_MsP_A / area_Ms_A : NaN
            shareP_B = area_Ms_B != 0 ? area_MsP_B / area_Ms_B : NaN

            println(io_sum, join((
                string(process),
                delta_s_max,
                area_Ms_A,
                area_Ms_B,
                ratio,
                area_absDsS_A,
                area_absDsS_B,
                area_absDsP_A,
                area_absDsP_B,
                area_kinS_A,
                area_kinS_B,
                area_kinP_A,
                area_kinP_B,
                shareS_A,
                shareS_B,
                shareP_A,
                shareP_B,
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
