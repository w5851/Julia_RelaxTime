#!/usr/bin/env julia

"""
    ordinary_pnjl_ad_fd_probe.jl

Diagnostic-only AD/FD probe for the ordinary five-field PNJL residual.

The four variants form a small factorial comparison:

* residual derivative: finite difference or ForwardDiff;
* NLsolve Jacobian: finite difference or ForwardDiff.

Each variant reports its first call separately from repeated calls in the
same Julia process.  The first call is a JIT-inclusive diagnostic; the steady
median is the algorithm comparison and excludes process startup/JIT effects
that are not repeated for every scan point.

This probe never calls a production scan, writes no numerical artifacts, and
does not change the production PNJL solver defaults.

Environment controls:

* `PNJL_AD_FD_MAX_POINTS=1` limits the number of points (default: all points).
* `PNJL_AD_FD_REPEATS=5` controls steady-state repetitions per variant.
"""

using LinearAlgebra: norm
using NLsolve
using Printf
using Statistics: median
using StaticArrays
using ForwardDiff

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .Models

const MODEL = Models.PNJLModel()
const X_STATE = SVector{5, Float64}(-1.84329, -1.84329, -2.22701, 0.5, 0.5)
const TEMPERATURES_MEV = (1.0, 50.0, 150.0, 240.0)
const CHEMICAL_POTENTIALS_MEV = (0.0, 120.0, 240.0)
const P_NUM = 4
const T_NUM = 4
const ITERATIONS = 12
const FD_STEP = 1.0e-5

struct Variant
    residual::Symbol
    jacobian::Symbol
end

const VARIANTS = (
    Variant(:finite, :finite),
    Variant(:finite, :forward),
    Variant(:forward, :finite),
    Variant(:forward, :forward),
)

function _env_int(name::String, default::Int; minimum::Int=0)
    raw = get(ENV, name, string(default))
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$(name) must be an integer, got $(raw)"))
    end
    value >= minimum || throw(ArgumentError("$(name) must be >= $(minimum), got $(value)"))
    return value
end

_primal_value(x::ForwardDiff.Dual) = ForwardDiff.value(x)
_primal_value(x::Real) = Float64(x)

@inline function _state_vector(x)
    T = eltype(x)
    return SVector{5, T}(x[1], x[2], x[3], x[4], x[5])
end

@inline function _omega(model, x, T_fm, mu_vec)
    return Models.omega(
        model,
        Models.MeanFieldState(_state_vector(x)),
        T_fm,
        mu_vec;
        p_num=P_NUM,
        t_num=T_NUM,
        xi=0.0,
    )
end

function _finite_difference_residual(model, x, T_fm, mu_vec)
    base = _state_vector(x)
    values = ntuple(5) do j
        h = FD_STEP * max(1.0, abs(_primal_value(base[j])))
        xp = setindex(base, base[j] + h, j)
        xm = setindex(base, base[j] - h, j)
        (_omega(model, xp, T_fm, mu_vec) - _omega(model, xm, T_fm, mu_vec)) / (2h)
    end
    return SVector{5, eltype(values)}(values)
end

function _residual_function(model, variant::Variant, T_fm, mu_vec)
    if variant.residual === :finite
        return x -> _finite_difference_residual(model, x, T_fm, mu_vec)
    elseif variant.residual === :forward
        return x -> Models.gap_residual(
            model,
            x,
            T_fm,
            mu_vec;
            p_num=P_NUM,
            t_num=T_NUM,
            xi=0.0,
        )
    end
    throw(ArgumentError("unsupported residual method $(variant.residual)"))
end

function _solve_once(model, variant::Variant, T_fm, mu_vec)
    residual_fn = _residual_function(model, variant, T_fm, mu_vec)
    f! = function (F, x)
        residual = residual_fn(x)
        @inbounds for i in 1:5
            F[i] = residual[i]
        end
        return nothing
    end

    t0 = time_ns()
    try
        result = nlsolve(
            f!,
            collect(X_STATE);
            autodiff=variant.jacobian,
            method=:trust_region,
            xtol=1e-8,
            ftol=1e-8,
            iterations=ITERATIONS,
        )
        wall_ms = (time_ns() - t0) / 1.0e6
        residual_norm = norm(residual_fn(result.zero))
        return (
            status=result.f_converged ? "converged" : "not_converged",
            wall_ms=wall_ms,
            f_calls=Int(result.f_calls),
            g_calls=Int(result.g_calls),
            iterations=Int(result.iterations),
            residual_norm=residual_norm,
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

function _measure_variant(model, variant::Variant, T_fm, mu_vec, repeats::Int)
    first = _solve_once(model, variant, T_fm, mu_vec)
    steady = [_solve_once(model, variant, T_fm, mu_vec) for _ in 1:repeats]
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
    max_points = _env_int("PNJL_AD_FD_MAX_POINTS", 0)
    repeats = _env_int("PNJL_AD_FD_REPEATS", 5; minimum=1)
    point_count = 0

    println("ordinary PNJL AD/FD residual/Jacobian A/B probe (diagnostic-only)")
    println("controls: p_num=$(P_NUM), t_num=$(T_NUM), iterations=$(ITERATIONS), fd_step=$(FD_STEP), repeats=$(repeats)")
    println("process_startup_excluded=true first_call_includes_variant_jit=true")

    for T_MeV in TEMPERATURES_MEV
        for μ_MeV in CHEMICAL_POTENTIALS_MEV
            max_points > 0 && point_count >= max_points && break
            point_count += 1
            T_fm = T_MeV / ħc_MeV_fm
            mu_vec = SVector{3, Float64}(fill(μ_MeV / ħc_MeV_fm, 3))
            println("point T=$(T_MeV) MeV mu=$(μ_MeV) MeV")
            finite_fn = _residual_function(MODEL, Variant(:finite, :finite), T_fm, mu_vec)
            forward_fn = _residual_function(MODEL, Variant(:forward, :forward), T_fm, mu_vec)
            parity = maximum(abs.(finite_fn(X_STATE) .- forward_fn(X_STATE)))

            for variant in VARIANTS
                measured = _measure_variant(MODEL, variant, T_fm, mu_vec, repeats)
                first = measured.first
                @printf(
                    "  residual=%s jacobian=%s | first_ms=%.3f first_f=%d first_g=%d first_iter=%d first_status=%s | steady_median_ms=%.3f steady_f=%.1f steady_g=%.1f steady_iter=%.1f steady_residual=%.3e steady_status=%s | initial_residual_max_diff=%.3e\n",
                    variant.residual,
                    variant.jacobian,
                    first.wall_ms,
                    first.f_calls,
                    first.g_calls,
                    first.iterations,
                    first.status,
                    measured.steady_wall_median,
                    measured.steady_f_calls,
                    measured.steady_g_calls,
                    measured.steady_iterations,
                    measured.steady_residual_norm,
                    measured.steady_status,
                    parity,
                )
            end
        end
        max_points > 0 && point_count >= max_points && break
    end

    println("status=diagnostic_only points=$(point_count) production_scan_called=false production_written=false")
end

main()
