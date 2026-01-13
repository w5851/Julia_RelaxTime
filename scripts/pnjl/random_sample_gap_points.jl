#!/usr/bin/env julia
# Random sampler to stress-test single-point PNJL equilibrium solve.
#
# Goal: find parameter points where single-point solve returns a problematic
# (non-finite / non-physical) result.
#
# Sampling ranges (MeV):
#   T ∈ [0, 350],  muB ∈ [0, 1800],  xi ∈ [-0.8, 0.8]
#
# Usage:
#   julia --project=. scripts/pnjl/random_sample_gap_points.jl
#
# Optional env knobs:
#   N_SAMPLES=200
#   SEED=1
#   P_NUM=24
#   T_NUM=8
#   T_MIN=0    T_MAX=350
#   MUB_MIN=0  MUB_MAX=1800
#   XI_MIN=-0.8 XI_MAX=0.8
#   STOP_ON_FIRST=1
#   OUT_CSV= (empty = no csv)

using Random
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .PNJL

@inline function env_int(key::String, default::Int)
    return parse(Int, get(ENV, key, string(default)))
end

@inline function env_float(key::String, default::Float64)
    return parse(Float64, get(ENV, key, string(default)))
end

@inline function env_bool(key::String, default::Bool=false)
    raw = get(ENV, key, default ? "1" : "0")
    return raw in ("1", "true", "TRUE", "yes", "YES")
end

"""Primary physicality rule for this sampler: masses must be positive and finite."""
@inline function mass_ok(masses)
    return all(isfinite, masses) && all(>(0.0), masses)
end

@inline function phi_ok(x_state; tol=1e-8)
    Φ = x_state[4]
    Φbar = x_state[5]
    return isfinite(Φ) && isfinite(Φbar) && (-tol <= Φ <= 1 + tol) && (-tol <= Φbar <= 1 + tol)
end

function write_csv_header(io)
    println(io, join([
        "idx",
        "T_MeV","muB_MeV","muq_MeV","xi",
        "p_num","t_num",
        "converged","residual_norm","iterations",
        "Phi","Phibar",
        "m_u_fm","m_d_fm","m_s_fm",
        "m_u_MeV","m_d_MeV","m_s_MeV",
        "omega","pressure","rho_norm",
        "mass_ok","phi_ok",
    ], ','))
end

function write_csv_row(io, idx, T_MeV, muB_MeV, xi, p_num, t_num, res)
    muq_MeV = muB_MeV / 3.0
    Φ = Float64(res.x_state[4])
    Φbar = Float64(res.x_state[5])
    mu = Float64(res.masses[1]); md = Float64(res.masses[2]); ms = Float64(res.masses[3])
    mu_MeV = mu / ħc_MeV_fm
    md_MeV = md / ħc_MeV_fm
    ms_MeV = ms / ħc_MeV_fm
    m_ok = mass_ok(res.masses)
    p_ok = phi_ok(res.x_state)

    row = [
        idx,
        T_MeV, muB_MeV, muq_MeV, xi,
        p_num, t_num,
        res.converged, res.residual_norm, res.iterations,
        Φ, Φbar,
        mu, md, ms,
        mu_MeV, md_MeV, ms_MeV,
        res.omega, res.pressure, res.rho_norm,
        m_ok, p_ok,
    ]
    println(io, join(string.(row), ','))
end

function main()
    n_samples = env_int("N_SAMPLES", 20000)
    seed = env_int("SEED", 1)
    p_num = env_int("P_NUM", 24)
    t_num = env_int("T_NUM", 8)

    T_min = env_float("T_MIN", 0.0)
    T_max = env_float("T_MAX", 350.0)
    muB_min = env_float("MUB_MIN", 0.0)
    muB_max = env_float("MUB_MAX", 1800.0)
    xi_min = env_float("XI_MIN", -0.8)
    xi_max = env_float("XI_MAX", 0.8)

    stop_on_first = env_bool("STOP_ON_FIRST", true)
    out_csv = get(ENV, "OUT_CSV", "")

    rng = MersenneTwister(seed)

    @printf("random_sample_gap_points: N=%d seed=%d\n", n_samples, seed)
    @printf("ranges: T=[%.1f, %.1f] MeV  muB=[%.1f, %.1f] MeV  xi=[%.3f, %.3f]\n", T_min, T_max, muB_min, muB_max, xi_min, xi_max)
    @printf("grid: p_num=%d t_num=%d\n", p_num, t_num)

    io = nothing
    if !isempty(out_csv)
        isdir(dirname(out_csv)) || mkpath(dirname(out_csv))
        io = open(out_csv, "w")
        write_csv_header(io)
        @printf("csv: %s\n", out_csv)
    end

    n_converged = 0
    n_mass_bad = 0
    n_phi_bad = 0
    n_nonfinite = 0
    n_failed = 0

    try
        for i in 1:n_samples
            T_MeV = T_min + rand(rng) * (T_max - T_min)
            muB_MeV = muB_min + rand(rng) * (muB_max - muB_min)
            xi = xi_min + rand(rng) * (xi_max - xi_min)

            # Convert to internal units (fm^-1)
            # Avoid exact T=0 which can cause derivatives/thermo to be ill-defined.
            T_fm = max(T_MeV, 1e-6) / ħc_MeV_fm
            muq_fm = (muB_MeV / 3.0) / ħc_MeV_fm

            local res
            try
                res = PNJL.solve(PNJL.FixedMu(), T_fm, muq_fm; xi=xi, p_num=p_num, t_num=t_num)
            catch err
                n_failed += 1
                @printf("[%d] ERROR: solve threw exception at T=%.3f muB=%.3f xi=%.3f : %s\n", i, T_MeV, muB_MeV, xi, string(err))
                if stop_on_first
                    @printf("repro: julia --project=. scripts/pnjl/diagnose_gap_point.jl --T=%.6f --muB=%.6f --xi=%.6f --p_num=%d --t_num=%d --verbose\n",
                        T_MeV, muB_MeV, xi, p_num, t_num)
                    return 1
                end
                continue
            end

            n_converged += res.converged ? 1 : 0

            m_ok = mass_ok(res.masses)
            p_ok = phi_ok(res.x_state)
            all_finite = isfinite(res.omega) && isfinite(res.pressure) && isfinite(res.rho_norm) && isfinite(res.residual_norm)

            if !all_finite
                n_nonfinite += 1
            end
            if !m_ok
                n_mass_bad += 1
            end
            if !p_ok
                n_phi_bad += 1
            end

            if io !== nothing
                write_csv_row(io, i, T_MeV, muB_MeV, xi, p_num, t_num, res)
            end

            problem = (res.converged && (!m_ok || !all_finite)) || (!res.converged)
            # We consider Phi/Phibar out-of-range as a separate signal but not primary (user rule focuses on masses).
            if problem
                @printf("[%d] PROBLEM: converged=%s |res|=%.3e  mass_ok=%s phi_ok=%s finite=%s\n",
                    i, string(res.converged), res.residual_norm, string(m_ok), string(p_ok), string(all_finite))
                @printf("     T=%.3f MeV  muB=%.3f MeV  xi=%.3f\n", T_MeV, muB_MeV, xi)
                @printf("     Phi=%.6f Phibar=%.6f  masses(fm^-1)=[%.6f %.6f %.6f]\n",
                    res.x_state[4], res.x_state[5], res.masses[1], res.masses[2], res.masses[3])
                @printf("     omega=%.6e pressure=%.6e rho_norm=%.6e\n", res.omega, res.pressure, res.rho_norm)
                @printf("repro: julia --project=. scripts/pnjl/diagnose_gap_point.jl --T=%.6f --muB=%.6f --xi=%.6f --p_num=%d --t_num=%d --verbose\n",
                    T_MeV, muB_MeV, xi, p_num, t_num)
                if stop_on_first
                    return 2
                end
            end
        end
    finally
        io === nothing || close(io)
    end

    println("\nsummary:")
    @printf("  converged: %d/%d\n", n_converged, n_samples)
    @printf("  failed(exceptions): %d\n", n_failed)
    @printf("  nonfinite(thermo/residual): %d\n", n_nonfinite)
    @printf("  mass_bad(m<=0 or NaN): %d\n", n_mass_bad)
    @printf("  phi_bad(out of [0,1] tol): %d\n", n_phi_bad)

    return 0
end

exit(main())
