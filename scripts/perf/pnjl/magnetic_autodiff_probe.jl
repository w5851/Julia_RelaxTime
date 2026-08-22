#!/usr/bin/env julia

"""
    magnetic_autodiff_probe.jl

Diagnostic-only A/B probe for the magnetic five-field stationarity residual.

The baseline uses the current finite-difference residual and NLsolve finite
Jacobian. The variant differentiates a fixed-`n_max` Landau Omega with
ForwardDiff and lets NLsolve differentiate that residual. This script uses one
explicit seed, low quadrature controls, and never calls `solve_magnetic_gap`,
production scans, or writes numerical artifacts.

Set `MAGNETIC_AD_PROBE_MAX_POINTS=1` for the shortest smoke probe.
"""

using LinearAlgebra: norm
using NLsolve
using Printf
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
const P_NUM = 4
const T_NUM = 4
const PZ_MAX = 5.0
const N_MAX = 1
const CUTOFF_N = 10
const ITERATIONS = 12

function _max_points()
    raw = get(ENV, "MAGNETIC_AD_PROBE_MAX_POINTS", "0")
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("MAGNETIC_AD_PROBE_MAX_POINTS must be an integer, got $(raw)"))
    end
    value >= 0 || throw(ArgumentError("MAGNETIC_AD_PROBE_MAX_POINTS must be >= 0, got $(value)"))
    return value
end

function _residual_function(model, residual_method::Symbol, T_fm, mu_vec)
    residual_method === :finite && return (x -> Models.magnetic_gap_residual(
        model, x, T_fm, mu_vec;
        p_num=P_NUM, t_num=T_NUM, pz_max=PZ_MAX, n_max=N_MAX, cutoff_N=CUTOFF_N, xi=0.0,
    ))
    residual_method === :forward && return (x -> Models.magnetic_gap_residual_autodiff(
        model, x, T_fm, mu_vec;
        p_num=P_NUM, t_num=T_NUM, pz_max=PZ_MAX, n_max=N_MAX, cutoff_N=CUTOFF_N, xi=0.0,
    ))
    throw(ArgumentError("unsupported residual method $(residual_method)"))
end

function _solve_once(model, T_fm, mu_vec, residual_method::Symbol)
    residual_fn = _residual_function(model, residual_method, T_fm, mu_vec)
    f! = function (F, x)
        r = residual_fn(x)
        @inbounds for i in 1:5
            F[i] = r[i]
        end
        return nothing
    end
    autodiff = residual_method === :forward ? :forward : :finite
    t0 = time_ns()
    try
        result = nlsolve(
            f!,
            collect(X_STATE);
            autodiff=autodiff,
            method=:trust_region,
            xtol=1e-8,
            ftol=1e-8,
            iterations=ITERATIONS,
        )
        wall_ms = (time_ns() - t0) / 1.0e6
        final_residual = residual_fn(result.zero)
        return (
            status=result.f_converged ? "converged" : "not_converged",
            wall_ms=wall_ms,
            f_calls=Int(result.f_calls),
            g_calls=Int(result.g_calls),
            iterations=Int(result.iterations),
            residual_norm=norm(final_residual),
        )
    catch err
        wall_ms = (time_ns() - t0) / 1.0e6
        return (
            status="error:" * string(nameof(typeof(err))),
            wall_ms=wall_ms,
            f_calls=-1,
            g_calls=-1,
            iterations=-1,
            residual_norm=Inf,
        )
    end
end

function main()
    max_points = _max_points()
    point_count = 0
    println("magnetic AD residual A/B probe (diagnostic-only)")
    println("controls: p_num=$(P_NUM), t_num=$(T_NUM), pz_max=$(PZ_MAX), n_max=$(N_MAX), iterations=$(ITERATIONS)")

    for T_MeV in TEMPERATURES_MEV
        for μ_MeV in CHEMICAL_POTENTIALS_MEV
            for eB_fm2 in MAGNETIC_FIELDS_FM2
                max_points > 0 && point_count >= max_points && break
                point_count += 1
                model = Models.PNJLMagneticModel(
                    eB_fm2=eB_fm2,
                    p_num=P_NUM,
                    pz_max=PZ_MAX,
                    n_max=N_MAX,
                    cutoff_N=CUTOFF_N,
                )
                T_fm = T_MeV / ħc_MeV_fm
                mu_vec = SVector{3, Float64}(fill(μ_MeV / ħc_MeV_fm, 3))
                finite_fn = _residual_function(model, :finite, T_fm, mu_vec)
                forward_fn = _residual_function(model, :forward, T_fm, mu_vec)
                parity = maximum(abs.(finite_fn(X_STATE) .- forward_fn(X_STATE)))
                finite = _solve_once(model, T_fm, mu_vec, :finite)
                forward = _solve_once(model, T_fm, mu_vec, :forward)
                @printf(
                    "point T=%g MeV mu=%g MeV eB=%g fm^-2 | finite wall=%.3fms f=%d g=%d iter=%d residual=%.3e status=%s | forward wall=%.3fms f=%d g=%d iter=%d residual=%.3e status=%s | initial_residual_max_diff=%.3e\n",
                    T_MeV, μ_MeV, eB_fm2,
                    finite.wall_ms, finite.f_calls, finite.g_calls, finite.iterations,
                    finite.residual_norm, finite.status,
                    forward.wall_ms, forward.f_calls, forward.g_calls, forward.iterations,
                    forward.residual_norm, forward.status,
                    parity,
                )
            end
            max_points > 0 && point_count >= max_points && break
        end
        max_points > 0 && point_count >= max_points && break
    end
    println("status=diagnostic_only points=$(point_count) production_scan_called=false production_written=false")
end

main()
