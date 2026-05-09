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
    xi_values = collect(range(0.0, 0.5; length=26))
    ds_vals = collect(range(1e-6, 2.0; length=36))
    process = :uubar_to_ddbar

    merged_csv = raw"D:\Desktop\Julia_RelaxTime\.worktrees\repro-main-oldparams\data\outputs\tmp\repro_main_oldparams\results\relaxtime\fixed_temperature_xi_scan_muB0\fixed_temperature_xi_scan_muB0_merged.csv"
    out = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_xi_positive_chain_contrast.csv"
    ensure_parent_dir(out)

    states = Dict{Float64, Any}()
    for xi in xi_values
        states[xi] = build_state_point(190.0, 0.0, xi)
    end

    rows = collect(CSV.File(merged_csv; comment="#"))
    t190 = filter(r -> Float64(r.T_MeV) == 190.0 && Float64(r.muB_MeV) == 0.0 && Float64(r.xi) >= 0.0 && Float64(r.xi) <= 0.5, rows)
    sort!(t190, by=r -> Float64(r.xi))
    sigma_over_T = Dict{Float64, Float64}(Float64(r.xi) => Float64(r.sigma_over_T) for r in t190)

    summary = Dict{Float64, Dict{String, Float64}}()
    for xi in xi_values
        st = states[xi]
        th = process_threshold_info(process, st.quark_params)

        detM_v = Float64[]
        Dmix_v = Float64[]
        DsP_mixed_v = Float64[]
        DsP_simple_v = Float64[]
        Ms_v = Float64[]
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
            push!(DsP_mixed_v, p_p.abs_D_s_P_sq_mixed)
            push!(DsP_simple_v, p_p.abs_D_s_P_sq_simple)
            push!(Ms_v, p_m.M_s_sq)
            push!(sigma_v, sigma_val)
        end

        summary[xi] = Dict(
            "area_abs_detM_sq" => trapz(ds_vals, detM_v),
            "area_abs_Dmixed_sq" => trapz(ds_vals, Dmix_v),
            "area_abs_DsP_mixed_sq" => trapz(ds_vals, DsP_mixed_v),
            "area_abs_DsP_simple_sq" => trapz(ds_vals, DsP_simple_v),
            "area_M_s_sq" => trapz(ds_vals, Ms_v),
            "area_sigma_total" => trapz(ds_vals, sigma_v),
        )
    end

    io = open(out, "w")
    try
        println(io, "xi_left,xi_right,ratio_area_sigma,right_over_left_ratio_sigma_over_T,ratio_area_abs_Dmixed_sq,ratio_area_abs_detM_sq,ratio_area_abs_DsP_mixed_sq,ratio_area_abs_DsP_simple_sq,ratio_area_M_s_sq")
        for i in 1:(length(xi_values) - 1)
            xl = xi_values[i]
            xr = xi_values[i + 1]
            left = summary[xl]
            right = summary[xr]
            println(io, join((
                xl,
                xr,
                right["area_sigma_total"] / left["area_sigma_total"],
                sigma_over_T[xr] / sigma_over_T[xl],
                right["area_abs_Dmixed_sq"] / left["area_abs_Dmixed_sq"],
                right["area_abs_detM_sq"] / left["area_abs_detM_sq"],
                right["area_abs_DsP_mixed_sq"] / left["area_abs_DsP_mixed_sq"],
                right["area_abs_DsP_simple_sq"] / left["area_abs_DsP_simple_sq"],
                right["area_M_s_sq"] / left["area_M_s_sq"],
            ), ','))
        end
    finally
        close(io)
    end
end

main()
