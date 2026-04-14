#!/usr/bin/env julia

using Printf
using StaticArrays
using Statistics: mean

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "ScatteringAmplitude.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: standard_nodes_weights
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .ScatteringAmplitude: prepare_scattering_context, scattering_amplitude_squared_prepared
using .TotalCrossSection: total_cross_section, calculate_t_bounds, calculate_final_state_energies, combined_final_state_factor

const PNJL = Models.pnjl_module()
const solve = getproperty(PNJL, :solve)
const FixedMu = getproperty(PNJL, :FixedMu)
const Integrals = getproperty(PNJL, :Integrals)
const DEFAULT_MOMENTUM_NODES = getproperty(Integrals, :DEFAULT_MOMENTUM_NODES)
const DEFAULT_MOMENTUM_WEIGHTS = getproperty(Integrals, :DEFAULT_MOMENTUM_WEIGHTS)

struct StatePoint
    xi::Float64
    quark_params::NamedTuple
    thermo_params::NamedTuple
    K_coeffs::NamedTuple
end

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function trapz(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == length(y) || error("trapz: x/y length mismatch")
    n <= 1 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        acc += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return acc
end

@inline function get_mass(flavor::Symbol, quark_params::NamedTuple)
    if flavor === :u || flavor === :ubar
        return quark_params.m.u
    elseif flavor === :d || flavor === :dbar
        return quark_params.m.d
    elseif flavor === :s || flavor === :sbar
        return quark_params.m.s
    end
    error("unknown flavor: $flavor")
end

@inline function get_mu(flavor::Symbol, quark_params::NamedTuple)
    if flavor === :u || flavor === :ubar
        return quark_params.μ.u
    elseif flavor === :d || flavor === :dbar
        return quark_params.μ.d
    elseif flavor === :s || flavor === :sbar
        return quark_params.μ.s
    end
    error("unknown flavor: $flavor")
end

function build_state(T_MeV::Float64, muB_MeV::Float64, xi::Float64)
    T_fm = T_MeV / ħc_MeV_fm
    muq_fm = (muB_MeV / 3.0) / ħc_MeV_fm

    base = solve(FixedMu(), T_fm, muq_fm; xi=xi, p_num=12, t_num=6)
    Bool(base.converged) || error("equilibrium not converged for xi=$xi")

    Phi = Float64(base.x_state[4])
    Phibar = Float64(base.x_state[5])
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

    A_u = A(masses.u, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(masses.s, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    quark_params = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm), A=(u=A_u, d=A_u, s=A_s))
    thermo_params = (T=T_fm, Φ=Phi, Φbar=Phibar, ξ=xi)
    return StatePoint(xi, quark_params, thermo_params, K_coeffs)
end

function process_threshold(process::Symbol, quark_params::NamedTuple)
    i, j, c, d = Main.AverageScatteringRate.parse_particles_from_process(process)
    mi = get_mass(i, quark_params)
    mj = get_mass(j, quark_params)
    mc = get_mass(c, quark_params)
    md = get_mass(d, quark_params)
    s_th = max((mi + mj)^2, (mc + md)^2)
    return (i=i, j=j, c=c, d=d, mi=mi, mj=mj, mc=mc, md=md, s_th=s_th)
end

function sigma_no_blocking(process::Symbol, s::Float64, quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple; n_points::Int)
    ctx = prepare_scattering_context(process, quark_params)
    mi, mj, mc, md = ctx.m1, ctx.m2, ctx.m3, ctx.m4
    s <= max((mi + mj)^2, (mc + md)^2) && return 0.0

    tb = calculate_t_bounds(s, mi, mj, mc, md)
    t_min, t_max = tb.t_min, tb.t_max
    if ctx.particle_c == ctx.particle_d
        t_min /= 2.0
    end
    (t_max - t_min) <= 1e-14 && return 0.0

    nodes, weights = standard_nodes_weights(n_points)
    scale_t = (t_max - t_min) / 2.0
    shift_t = (t_max + t_min) / 2.0
    s_12_plus = s - (mi + mj)^2
    s_12_minus = s - (mi - mj)^2

    sigma = 0.0
    @inbounds for k in 1:n_points
        t = scale_t * nodes[k] + shift_t
        w = scale_t * weights[k]
        M2 = scattering_amplitude_squared_prepared(ctx, s, t, quark_params, thermo_params, K_coeffs)
        if !(isfinite(M2) && M2 > 0.0)
            continue
        end
        dsdt = Main.DifferentialCrossSection.differential_cross_section(s_12_plus, s_12_minus, M2)
        isfinite(dsdt) || continue
        sigma += w * max(dsdt, 0.0)
    end
    return sigma
end

function mean_blocking_factor(process::Symbol, s::Float64, quark_params::NamedTuple, thermo_params::NamedTuple; n_points::Int)
    i, j, c, d = Main.AverageScatteringRate.parse_particles_from_process(process)
    mi, mj = get_mass(i, quark_params), get_mass(j, quark_params)
    mc, md = get_mass(c, quark_params), get_mass(d, quark_params)
    mu_c, mu_d = get_mu(c, quark_params), get_mu(d, quark_params)
    s <= max((mi + mj)^2, (mc + md)^2) && return NaN

    tb = calculate_t_bounds(s, mi, mj, mc, md)
    t_min, t_max = tb.t_min, tb.t_max
    if c == d
        t_min /= 2.0
    end
    nodes, weights = standard_nodes_weights(n_points)
    scale_t = (t_max - t_min) / 2.0
    shift_t = (t_max + t_min) / 2.0
    E_c, E_d = calculate_final_state_energies(s, mc, md)

    z = 0.0
    bz = 0.0
    @inbounds for k in 1:n_points
        t = scale_t * nodes[k] + shift_t
        w = scale_t * weights[k]
        cos_star = Main.TotalCrossSection.cos_theta_star(s, t, mi, mj, mc, md)
        b = combined_final_state_factor(c, d, E_c, E_d, mc, md, mu_c, mu_d,
            thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ, cos_star)
        isfinite(b) || continue
        z += w
        bz += w * b
    end
    z > 0 ? (bz / z) : NaN
end

function main()
    T_MeV = 190.0
    muB_MeV = 0.0
    xi_A = -0.10
    xi_B = -0.08
    processes = Symbol[:uubar_to_ddbar, :uubar_to_uubar]

    delta_s_max = 2.0
    n_s = 36
    n_sigma_points = 24
    out = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_sigma_chain_decomposition.csv"
    out_summary = raw"D:\Desktop\Temp\relaxtime_t190_window\t190_sigma_chain_decomposition_summary.csv"

    ensure_parent_dir(out)
    ensure_parent_dir(out_summary)

    state_A = build_state(T_MeV, muB_MeV, xi_A)
    state_B = build_state(T_MeV, muB_MeV, xi_B)

    io = open(out, "w")
    io_sum = open(out_summary, "w")
    try
        println(io, "process,xi_state,s_th,s,s_minus_sth,t_width,sigma_base,sigma_no_block,sigma_Kswap,mean_blocking")
        println(io_sum, "process,delta_s_max,area_base_A,area_base_B,ratio_base_B_over_A,area_no_block_A,area_no_block_B,ratio_no_block_B_over_A,block_factor_area_A,block_factor_area_B,area_A_with_KB,area_B_with_KA,K_sensitivity_A_rel,K_sensitivity_B_rel,t_width_mean_A,t_width_mean_B")

        for process in processes
            thA = process_threshold(process, state_A.quark_params)
            thB = process_threshold(process, state_B.quark_params)

            ds_vals = collect(range(1e-6, delta_s_max; length=n_s))

            sA_vals = [thA.s_th + ds for ds in ds_vals]
            sB_vals = [thB.s_th + ds for ds in ds_vals]

            baseA = Float64[]
            baseB = Float64[]
            nobA = Float64[]
            nobB = Float64[]
            kswA = Float64[]
            kswB = Float64[]
            blkA = Float64[]
            blkB = Float64[]
            twA = Float64[]
            twB = Float64[]

            for (idx, ds) in enumerate(ds_vals)
                sA = sA_vals[idx]
                sB = sB_vals[idx]

                tbA = calculate_t_bounds(sA, thA.mi, thA.mj, thA.mc, thA.md)
                tbB = calculate_t_bounds(sB, thB.mi, thB.mj, thB.mc, thB.md)
                widthA = tbA.t_max - (thA.c == thA.d ? tbA.t_min / 2.0 : tbA.t_min)
                widthB = tbB.t_max - (thB.c == thB.d ? tbB.t_min / 2.0 : tbB.t_min)

                s_base_A = total_cross_section(process, sA, state_A.quark_params, state_A.thermo_params, state_A.K_coeffs; n_points=n_sigma_points)
                s_base_B = total_cross_section(process, sB, state_B.quark_params, state_B.thermo_params, state_B.K_coeffs; n_points=n_sigma_points)
                s_nob_A = sigma_no_blocking(process, sA, state_A.quark_params, state_A.thermo_params, state_A.K_coeffs; n_points=n_sigma_points)
                s_nob_B = sigma_no_blocking(process, sB, state_B.quark_params, state_B.thermo_params, state_B.K_coeffs; n_points=n_sigma_points)
                s_ksw_A = total_cross_section(process, sA, state_A.quark_params, state_A.thermo_params, state_B.K_coeffs; n_points=n_sigma_points)
                s_ksw_B = total_cross_section(process, sB, state_B.quark_params, state_B.thermo_params, state_A.K_coeffs; n_points=n_sigma_points)
                bA = mean_blocking_factor(process, sA, state_A.quark_params, state_A.thermo_params; n_points=n_sigma_points)
                bB = mean_blocking_factor(process, sB, state_B.quark_params, state_B.thermo_params; n_points=n_sigma_points)

                push!(baseA, s_base_A); push!(baseB, s_base_B)
                push!(nobA, s_nob_A); push!(nobB, s_nob_B)
                push!(kswA, s_ksw_A); push!(kswB, s_ksw_B)
                push!(blkA, bA); push!(blkB, bB)
                push!(twA, widthA); push!(twB, widthB)

                println(io, join((string(process), string(xi_A), thA.s_th, sA, ds, widthA, s_base_A, s_nob_A, s_ksw_A, bA), ','))
                println(io, join((string(process), string(xi_B), thB.s_th, sB, ds, widthB, s_base_B, s_nob_B, s_ksw_B, bB), ','))
            end

            area_base_A = trapz(ds_vals, baseA)
            area_base_B = trapz(ds_vals, baseB)
            area_nob_A = trapz(ds_vals, nobA)
            area_nob_B = trapz(ds_vals, nobB)
            area_ksw_A = trapz(ds_vals, kswA)
            area_ksw_B = trapz(ds_vals, kswB)

            ratio_base = area_base_A > 0 ? area_base_B / area_base_A : NaN
            ratio_nob = area_nob_A > 0 ? area_nob_B / area_nob_A : NaN
            block_factor_A = area_nob_A > 0 ? area_base_A / area_nob_A : NaN
            block_factor_B = area_nob_B > 0 ? area_base_B / area_nob_B : NaN
            sens_A = area_base_A > 0 ? (area_ksw_A / area_base_A - 1.0) : NaN
            sens_B = area_ksw_B > 0 ? (area_base_B / area_ksw_B - 1.0) : NaN

            println(io_sum, join((
                string(process),
                delta_s_max,
                area_base_A,
                area_base_B,
                ratio_base,
                area_nob_A,
                area_nob_B,
                ratio_nob,
                block_factor_A,
                block_factor_B,
                area_ksw_A,
                area_ksw_B,
                sens_A,
                sens_B,
                mean(twA),
                mean(twB),
            ), ','))
            flush(io)
            flush(io_sum)
            @printf("process=%s done\n", string(process))
        end
    finally
        close(io)
        close(io_sum)
    end

    @printf("Wrote detail: %s\n", out)
    @printf("Wrote summary: %s\n", out_summary)
end

main()
