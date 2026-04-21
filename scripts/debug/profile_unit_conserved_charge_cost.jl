#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using Main.Constants_PNJL: ħc_MeV_fm
const PNJL = Models.pnjl_module()

function _time_block(f::Function, name::String)
    t = @elapsed f()
    return (name=name, sec=t)
end

function _profile_once()
    rows = NamedTuple[]

    push!(rows, _time_block("chi1..chi4 finite (48x12)") do
        T_fm = 140.0 / ħc_MeV_fm
        muB_fm = 360.0 / ħc_MeV_fm
        PNJL.chi1_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)
        PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)
        PNJL.chi3_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)
        PNJL.chi4_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)
    end)

    push!(rows, _time_block("chi1_B vs rho/T^3 (48x12)") do
        T_fm = 130.0 / ħc_MeV_fm
        muB_fm = 300.0 / ħc_MeV_fm
        muq_fm = muB_fm / 3
        PNJL.thermo_derivatives(T_fm, muq_fm; xi=0.0, p_num=48, t_num=12)
        PNJL.chi1_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)
    end)

    push!(rows, _time_block("cumulant+ratios consistency (40x10)") do
        T_fm = 0.55
        muB_fm = 0.9
        V = 125.0
        PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)
        PNJL.chi3_B(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)
        PNJL.chi4_B(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)
        PNJL.cumulant_B(T_fm, muB_fm, V; order=4, xi=0.0, p_num=40, t_num=10)
        PNJL.baryon_Ssigma(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)
        PNJL.baryon_kappa_sigma2(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)
    end)

    push!(rows, _time_block("odd-order near zero (32x10)") do
        T_fm = 0.5
        PNJL.chi1_B(T_fm, 0.0; xi=0.0, p_num=32, t_num=10)
        PNJL.chi3_B(T_fm, 0.0; xi=0.0, p_num=32, t_num=10)
    end)

    push!(rows, _time_block("flavor pressure derivatives (order=2)") do
        T_fm = 0.55
        mu_vec = SVector(0.10, 0.06, 0.02)
        PNJL.flavor_pressure_derivatives(T_fm, mu_vec; order=2, xi=0.0, p_num=32, t_num=10)
    end)

    push!(rows, _time_block("2nd-order conserved set (32x10)") do
        T_fm = 0.55
        muB_fm = 0.30
        muQ_fm = 0.06
        muS_fm = 0.02
        PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=32, t_num=10)
        PNJL.chi2_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10)
        PNJL.chi2_S(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10)
        PNJL.chi11_BQ(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10)
        PNJL.chi11_BS(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10)
        PNJL.chi11_QS(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10)
    end)

    push!(rows, _time_block("chi2_B mapping check (32x10)") do
        T_fm = 0.60
        muB_fm = 0.24
        PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=32, t_num=10)
        PNJL.conserved_charge_susceptibility(T_fm, muB_fm, 0.0, 0.0; orders=(2, 0, 0), xi=0.0, p_num=32, t_num=10)
    end)

    push!(rows, _time_block("chi_BQS wrappers parity (32x10)") do
        T_fm = 0.58
        muB_fm = 0.21
        muQ_fm = 0.04
        muS_fm = 0.01
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 2, 0), xi=0.0, p_num=32, t_num=10)
        PNJL.chi2_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10)
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), xi=0.0, p_num=32, t_num=10)
        PNJL.chi11_BQ(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10)
        PNJL.chi_BQS(T_fm, muB_fm, 0.0, 0.0; orders=(2, 0, 0), xi=0.0, p_num=32, t_num=10)
        PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=32, t_num=10)
    end)

    push!(rows, _time_block("pure Q/S higher-order (24x8)") do
        T_fm = 0.57
        muB_fm = 0.18
        muQ_fm = 0.05
        muS_fm = 0.02
        PNJL.chi1_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8)
        PNJL.chi3_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8)
        PNJL.chi4_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8)
        PNJL.chi1_S(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8)
        PNJL.chi3_S(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8)
        PNJL.chi4_S(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8)
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 3, 0), xi=0.0, p_num=24, t_num=8)
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 0, 4), xi=0.0, p_num=24, t_num=8)
    end)

    push!(rows, _time_block("cumulant_BQS VT^3 checks (32x10)") do
        T_fm = 0.57
        muB_fm = 0.18
        muQ_fm = 0.05
        muS_fm = 0.02
        V = 64.0
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 2, 0), xi=0.0, p_num=32, t_num=10)
        PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(0, 2, 0), xi=0.0, p_num=32, t_num=10)
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), xi=0.0, p_num=32, t_num=10)
        PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(1, 1, 0), xi=0.0, p_num=32, t_num=10)
    end)

    push!(rows, _time_block("cumulant_B vs BQS baryon path (32x10)") do
        T_fm = 0.54
        muB_fm = 0.33
        V = 50.0
        PNJL.cumulant_B(T_fm, muB_fm, V; order=4, xi=0.0, p_num=32, t_num=10)
        PNJL.cumulant_BQS(T_fm, muB_fm, 0.0, 0.0, V; orders=(4, 0, 0), xi=0.0, p_num=32, t_num=10)
    end)

    return rows
end

function _print_report(title::String, rows)
    total = sum(r -> r.sec, rows)
    println("\n" * title)
    println(repeat("-", 86))
    println(rpad("block", 56) * lpad("sec", 12) * lpad("share", 12))
    for r in sort(rows; by=x -> x.sec, rev=true)
        share = total > 0 ? 100 * r.sec / total : 0.0
        println(rpad(r.name, 56) * lpad(@sprintf("%.2f", r.sec), 12) * lpad(@sprintf("%.1f%%", share), 12))
    end
    println(repeat("-", 86))
    println(rpad("total", 56) * lpad(@sprintf("%.2f", total), 12) * lpad("100.0%", 12))
    return total
end

println("Profiling high-order conserved-charge unit workflow")
println("Project root: " * PROJECT_ROOT)

cold_rows = _profile_once()
hot_rows = _profile_once()

cold_total = _print_report("Cold pass (compile + run)", cold_rows)
hot_total = _print_report("Hot pass (mostly run-time)", hot_rows)

compile_est = max(cold_total - hot_total, 0.0)
compile_share = cold_total > 0 ? 100 * compile_est / cold_total : 0.0
runtime_share = cold_total > 0 ? 100 * hot_total / cold_total : 0.0

println("\nEstimated split")
println(repeat("-", 86))
println(@sprintf("cold total      : %.2f s", cold_total))
println(@sprintf("hot total       : %.2f s", hot_total))
println(@sprintf("compile estimate: %.2f s (%.1f%%)", compile_est, compile_share))
println(@sprintf("runtime estimate: %.2f s (%.1f%%)", hot_total, runtime_share))
