#!/usr/bin/env julia

"""
    magnetic_autodiff_probe.jl

Diagnostic-only steady-state probe for the magnetic five-field AD production
path. The first call is reported separately because it includes compilation;
the repeated calls in the same Julia process are the algorithm timing.

This probe uses one explicit Landau cutoff and one seed. It never calls the
production scan or writes numerical artifacts.
"""

using LinearAlgebra: norm
using NLsolve
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

function _repeats()
    raw = get(ENV, "MAGNETIC_AD_PROBE_REPEATS", "5")
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("MAGNETIC_AD_PROBE_REPEATS must be an integer, got $(raw)"))
    end
    value >= 1 || throw(ArgumentError("MAGNETIC_AD_PROBE_REPEATS must be >= 1, got $(value)"))
    return value
end

function _residual_function(model, T_fm, mu_vec)
    return x -> Models.magnetic_gap_residual(
        model, x, T_fm, mu_vec;
        p_num=P_NUM, t_num=T_NUM, pz_max=PZ_MAX, n_max=N_MAX, cutoff_N=CUTOFF_N, xi=0.0,
    )
end

function _solve_once(model, T_fm, mu_vec)
    residual_fn = _residual_function(model, T_fm, mu_vec)
    f! = function (F, x)
        r = residual_fn(x)
        @inbounds for i in 1:5
            F[i] = r[i]
        end
        return nothing
    end
    t0 = time_ns()
    try
        result = nlsolve(
            f!,
            collect(X_STATE);
            autodiff=:forward,
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

function _measure(model, T_fm, mu_vec, repeats::Int)
    first = _solve_once(model, T_fm, mu_vec)
    steady = [_solve_once(model, T_fm, mu_vec) for _ in 1:repeats]
    return (
        first=first,
        steady_wall_median=median(getfield.(steady, :wall_ms)),
        steady_f_calls=median(getfield.(steady, :f_calls)),
        steady_g_calls=median(getfield.(steady, :g_calls)),
        steady_iterations=median(getfield.(steady, :iterations)),
        steady_residual_norm=median(getfield.(steady, :residual_norm)),
        steady_status=join(unique(getfield.(steady, :status)), ","),
    )
end

function main()
    max_points = _max_points()
    repeats = _repeats()
    point_count = 0
    println("magnetic AD production-path probe (diagnostic-only)")
    println("controls: p_num=$(P_NUM), t_num=$(T_NUM), pz_max=$(PZ_MAX), n_max=$(N_MAX), iterations=$(ITERATIONS)")
    println("process_startup_excluded=true first_call_includes_variant_jit=true repeats=$(repeats)")

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
                result = _measure(model, T_fm, mu_vec, repeats)
                @printf(
                    "point T=%g MeV mu=%g MeV eB=%g fm^-2 | AD first_ms=%.3f first_f=%d first_g=%d first_iter=%d first_status=%s steady_median_ms=%.3f steady_f=%.1f steady_g=%.1f steady_iter=%.1f steady_residual=%.3e steady_status=%s\n",
                    T_MeV, μ_MeV, eB_fm2,
                    result.first.wall_ms, result.first.f_calls, result.first.g_calls, result.first.iterations,
                    result.first.status, result.steady_wall_median, result.steady_f_calls,
                    result.steady_g_calls, result.steady_iterations, result.steady_residual_norm,
                    result.steady_status,
                )
            end
            max_points > 0 && point_count >= max_points && break
        end
        max_points > 0 && point_count >= max_points && break
    end
    println("status=diagnostic_only points=$(point_count) production_scan_called=false production_written=false")
end

main()
