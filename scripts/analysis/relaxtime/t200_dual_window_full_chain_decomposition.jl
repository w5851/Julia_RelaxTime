#!/usr/bin/env julia

using CSV

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
    T_MeV = 200.0
    muB_MeV = 0.0
    xi_values = sort(vcat(collect(range(-0.4, -0.2; length=11)), collect(range(0.2, 0.4; length=11))))
    processes = Symbol[:uubar_to_ddbar, :uubar_to_uubar]
    ds_vals = collect(range(1e-6, 2.0; length=36))

    merged_csv = raw"D:\Desktop\Julia_RelaxTime\.worktrees\repro-main-oldparams\data\outputs\tmp\repro_main_oldparams\results\relaxtime\fixed_temperature_xi_scan_muB0\fixed_temperature_xi_scan_muB0_merged.csv"
    out_detail = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_full_chain_detail.csv"
    out_summary = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_full_chain_summary.csv"
    out_adjacent = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_adjacent_transition_summary.csv"
    ensure_parent_dir(out_detail)
    ensure_parent_dir(out_summary)
    ensure_parent_dir(out_adjacent)

    states = Dict{Float64, Any}()
    for xi in xi_values
        states[xi] = build_state_point(T_MeV, muB_MeV, xi)
    end

    io = open(out_detail, "w")
    io_sum = open(out_summary, "w")
    io_adj = open(out_adjacent, "w")
    try
        println(io, "process,xi,s_minus_sth,abs_detM_sq,abs_D_mixed_P_sq,abs_JMJ_sq,detK_plus,abs_DsP_total_sq,abs_DsP_simple_sq,abs_DsP_mixed_sq,MsP_total,M2_total,M_s_sq,M_t_sq,M_interf,sigma_total")
        println(io_sum, "process,xi,area_abs_detM_sq,area_abs_Dmixed_sq,area_abs_JMJ_sq,area_detK_plus,area_abs_DsP_total_sq,area_abs_DsP_simple_sq,area_abs_DsP_mixed_sq,area_MsP_total,area_M2_total,area_M_s_sq,area_M_t_sq,area_M_interf,area_sigma_total")
        println(io_adj, "process,xi_left,xi_right,ratio_area_sigma,right_over_left_ratio_sigma_over_T,right_over_left_ratio_tau_u,ratio_area_abs_Dmixed_sq,ratio_area_abs_detM_sq,ratio_area_abs_DsP_mixed_sq,ratio_area_abs_DsP_simple_sq,ratio_area_M_s_sq,ratio_area_M_t_sq,ratio_area_M_interf,ratio_area_abs_JMJ_sq,ratio_area_detK_plus")

        summary_by_process = Dict{Symbol, Dict{Float64, Dict{String, Float64}}}()
        for process in processes
            summary_by_process[process] = Dict{Float64, Dict{String, Float64}}()
            for xi in xi_values
                st = states[xi]
                th = process_threshold_info(process, st.quark_params)

                detM_v = Float64[]
                Dmix_v = Float64[]
                JMJ_v = Float64[]
                detK_v = Float64[]
                DsP_total_v = Float64[]
                DsP_simple_v = Float64[]
                DsP_mixed_v = Float64[]
                MsP_v = Float64[]
                M2_v = Float64[]
                Ms_v = Float64[]
                Mt_v = Float64[]
                Mi_v = Float64[]
                sigma_v = Float64[]

                for ds in ds_vals
                    s = th.s_th + ds
                    tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
                    t = 0.5 * (tb.t_min + tb.t_max)

                    p_chain = decompose_mixed_p_propagator_chain(process, s, t, st.quark_params, st.thermo_params, st.K_coeffs)
                    p_p = decompose_p_channel_propagator_strength(process, s, t, st.quark_params, st.thermo_params, st.K_coeffs)
                    p_m = decompose_qqbar_amplitude_terms(process, s, t, st.quark_params, st.thermo_params, st.K_coeffs)
                    sigma_val = Main.TotalCrossSection.total_cross_section(process, s, st.quark_params, st.thermo_params, st.K_coeffs; n_points=24)

                    push!(detM_v, p_chain.abs_detM_sq)
                    push!(Dmix_v, p_chain.abs_D_mixed_P_sq)
                    push!(JMJ_v, p_chain.abs_JMJ_sq)
                    push!(detK_v, p_chain.detK_plus)
                    push!(DsP_total_v, p_p.abs_D_s_P_sq_total)
                    push!(DsP_simple_v, p_p.abs_D_s_P_sq_simple)
                    push!(DsP_mixed_v, p_p.abs_D_s_P_sq_mixed)
                    push!(MsP_v, p_p.M_s_P_total)
                    push!(M2_v, p_m.M2_total)
                    push!(Ms_v, p_m.M_s_sq)
                    push!(Mt_v, p_m.M_t_sq)
                    push!(Mi_v, p_m.M_interf)
                    push!(sigma_v, sigma_val)

                    println(io, join((string(process), xi, ds, p_chain.abs_detM_sq, p_chain.abs_D_mixed_P_sq, p_chain.abs_JMJ_sq, p_chain.detK_plus, p_p.abs_D_s_P_sq_total, p_p.abs_D_s_P_sq_simple, p_p.abs_D_s_P_sq_mixed, p_p.M_s_P_total, p_m.M2_total, p_m.M_s_sq, p_m.M_t_sq, p_m.M_interf, sigma_val), ','))
                end

                metrics = Dict(
                    "area_abs_detM_sq" => trapz(ds_vals, detM_v),
                    "area_abs_Dmixed_sq" => trapz(ds_vals, Dmix_v),
                    "area_abs_JMJ_sq" => trapz(ds_vals, JMJ_v),
                    "area_detK_plus" => trapz(ds_vals, detK_v),
                    "area_abs_DsP_total_sq" => trapz(ds_vals, DsP_total_v),
                    "area_abs_DsP_simple_sq" => trapz(ds_vals, DsP_simple_v),
                    "area_abs_DsP_mixed_sq" => trapz(ds_vals, DsP_mixed_v),
                    "area_MsP_total" => trapz(ds_vals, MsP_v),
                    "area_M2_total" => trapz(ds_vals, M2_v),
                    "area_M_s_sq" => trapz(ds_vals, Ms_v),
                    "area_M_t_sq" => trapz(ds_vals, Mt_v),
                    "area_M_interf" => trapz(ds_vals, Mi_v),
                    "area_sigma_total" => trapz(ds_vals, sigma_v),
                )
                summary_by_process[process][xi] = metrics

                println(io_sum, join((
                    string(process),
                    xi,
                    metrics["area_abs_detM_sq"],
                    metrics["area_abs_Dmixed_sq"],
                    metrics["area_abs_JMJ_sq"],
                    metrics["area_detK_plus"],
                    metrics["area_abs_DsP_total_sq"],
                    metrics["area_abs_DsP_simple_sq"],
                    metrics["area_abs_DsP_mixed_sq"],
                    metrics["area_MsP_total"],
                    metrics["area_M2_total"],
                    metrics["area_M_s_sq"],
                    metrics["area_M_t_sq"],
                    metrics["area_M_interf"],
                    metrics["area_sigma_total"],
                ), ','))
            end
        end

        rows = collect(CSV.File(merged_csv; comment="#"))
        t200_rows = filter(r -> Float64(r.T_MeV) == 200.0 && Float64(r.muB_MeV) == 0.0, rows)
        sort!(t200_rows, by=r -> Float64(r.xi))
        sigma_over_T_by_xi = Dict{Float64, Float64}()
        tau_u_by_xi = Dict{Float64, Float64}()
        for r in t200_rows
            x = Float64(r.xi)
            sigma_over_T_by_xi[x] = Float64(r.sigma_over_T)
            tau_u_by_xi[x] = Float64(r.tau_u)
        end

        for process in processes
            for i in 1:(length(xi_values) - 1)
                xl = xi_values[i]
                xr = xi_values[i + 1]
                (xl < -0.4 || xr > -0.2) && (xl < 0.2 || xr > 0.4) && continue

                left = summary_by_process[process][xl]
                right = summary_by_process[process][xr]

                ratio_sigma = right["area_sigma_total"] / left["area_sigma_total"]
                ratio_sigma_over_T = sigma_over_T_by_xi[xr] / sigma_over_T_by_xi[xl]
                ratio_tau_u = tau_u_by_xi[xr] / tau_u_by_xi[xl]

                println(io_adj, join((
                    string(process),
                    xl,
                    xr,
                    ratio_sigma,
                    ratio_sigma_over_T,
                    ratio_tau_u,
                    right["area_abs_Dmixed_sq"] / left["area_abs_Dmixed_sq"],
                    right["area_abs_detM_sq"] / left["area_abs_detM_sq"],
                    right["area_abs_DsP_mixed_sq"] / left["area_abs_DsP_mixed_sq"],
                    right["area_abs_DsP_simple_sq"] / left["area_abs_DsP_simple_sq"],
                    right["area_M_s_sq"] / left["area_M_s_sq"],
                    right["area_M_t_sq"] / left["area_M_t_sq"],
                    right["area_M_interf"] / left["area_M_interf"],
                    right["area_abs_JMJ_sq"] / left["area_abs_JMJ_sq"],
                    right["area_detK_plus"] / left["area_detK_plus"],
                ), ','))
            end
        end
    finally
        close(io)
        close(io_sum)
        close(io_adj)
    end
end

main()
