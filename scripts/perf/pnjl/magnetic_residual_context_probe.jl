#!/usr/bin/env julia

"""
    magnetic_residual_context_probe.jl

Diagnostic-only A/B for the magnetic finite-difference residual.

The baseline computes the center Omega even when `n_max` is explicitly supplied;
the context variant skips that unused center evaluation. Both variants use the
same perturbed states, quadrature controls, and explicit Landau cutoff.
This script never calls `solve_magnetic_gap` and does not write production data.
"""

using LinearAlgebra: norm
using Printf
using Statistics: median
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .Models

const X_STATE = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)
const TEMPERATURES_MEV = (1.0, 50.0, 150.0, 240.0)
const CHEMICAL_POTENTIALS_MEV = (0.0, 60.0, 240.0)
const MAGNETIC_FIELDS_FM2 = (0.1, 0.5, 2.0)
const P_NUM = 8
const T_NUM = 4
const PZ_MAX = 5.0
const N_MAX = 8
const CUTOFF_N = 10
const REPEATS = 5

function residual_skip_explicit_nmax(
    model,
    x,
    T_fm,
    mu_vec;
    xi::Real=0.0,
    p_num::Int=P_NUM,
    t_num::Int=T_NUM,
    pz_max::Real=PZ_MAX,
    n_max::Int=N_MAX,
    cutoff_N::Int=CUTOFF_N,
    finite_difference_step::Real=1e-5,
)
    _ = xi
    _ = t_num
    μ = SVector{3, Float64}(Tuple(mu_vec))
    x0 = SVector{5, Float64}(Tuple(x))
    h0 = Float64(finite_difference_step)
    gradient = MVector{5, Float64}(0.0, 0.0, 0.0, 0.0, 0.0)

    # With explicit n_max, the center Omega is not needed to construct the
    # finite-difference residual; only the ten +/- perturbations are required.
    @inbounds for i in 1:5
        h = h0 * max(1.0, abs(x0[i]))
        xp = MVector{5, Float64}(x0)
        xm = MVector{5, Float64}(x0)
        xp[i] += h
        xm[i] -= h
        ωp = Models.calculate_magnetic_omega_components(
            SVector{5, Float64}(xp), μ, T_fm, model.magnetic;
            xi=xi, p_num=p_num, t_num=t_num, pz_max=pz_max, n_max=n_max, cutoff_N=cutoff_N,
        ).omega
        ωm = Models.calculate_magnetic_omega_components(
            SVector{5, Float64}(xm), μ, T_fm, model.magnetic;
            xi=xi, p_num=p_num, t_num=t_num, pz_max=pz_max, n_max=n_max, cutoff_N=cutoff_N,
        ).omega
        gradient[i] = (ωp - ωm) / (2 * h)
    end
    return SVector{5, Float64}(gradient)
end

function timed_median(f; repeats::Int=REPEATS)
    f()
    samples_ms = Float64[]
    for _ in 1:repeats
        t0 = time_ns()
        f()
        push!(samples_ms, (time_ns() - t0) / 1.0e6)
    end
    return median(samples_ms)
end

function main()
    println("A/B-M1 magnetic residual explicit-n_max center-elision probe")
    println("controls: p_num=$(P_NUM), t_num=$(T_NUM), pz_max=$(PZ_MAX), n_max=$(N_MAX), repeats=$(REPEATS)")

    baseline_total = Float64[]
    variant_total = Float64[]
    for T_MeV in TEMPERATURES_MEV
        for μ_MeV in CHEMICAL_POTENTIALS_MEV
            for eB_fm2 in MAGNETIC_FIELDS_FM2
                model = Models.PNJLMagneticModel(
                    eB_fm2=eB_fm2,
                    p_num=P_NUM,
                    pz_max=PZ_MAX,
                    n_max=N_MAX,
                    cutoff_N=CUTOFF_N,
                )
                T_fm = T_MeV / ħc_MeV_fm
                μ = SVector{3, Float64}(fill(μ_MeV / ħc_MeV_fm, 3))
                baseline = () -> Models.magnetic_gap_residual(
                    model,
                    X_STATE,
                    T_fm,
                    μ;
                    xi=0.0,
                    p_num=P_NUM,
                    t_num=T_NUM,
                    pz_max=PZ_MAX,
                    n_max=N_MAX,
                    cutoff_N=CUTOFF_N,
                )
                variant = () -> residual_skip_explicit_nmax(
                    model,
                    X_STATE,
                    T_fm,
                    μ;
                    xi=0.0,
                    p_num=P_NUM,
                    t_num=T_NUM,
                    pz_max=PZ_MAX,
                    n_max=N_MAX,
                    cutoff_N=CUTOFF_N,
                )

                r_base = baseline()
                r_variant = variant()
                parity = maximum(abs.(r_base .- r_variant))
                base_ms = timed_median(baseline)
                variant_ms = timed_median(variant)
                push!(baseline_total, base_ms)
                push!(variant_total, variant_ms)
                speedup = base_ms / variant_ms
                @printf(
                    "point T=%g MeV mu=%g MeV eB=%g fm^-2 | baseline=%.6f ms variant=%.6f ms speedup=%.3fx parity=%.3e\n",
                    T_MeV, μ_MeV, eB_fm2, base_ms, variant_ms, speedup, parity,
                )
            end
        end
    end

    base_med = median(baseline_total)
    variant_med = median(variant_total)
    @printf("summary baseline_median=%.6f ms variant_median=%.6f ms speedup=%.3fx\n", base_med, variant_med, base_med / variant_med)
    println("status=diagnostic_only solver_called=false production_written=false")
end

main()
