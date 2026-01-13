#!/usr/bin/env julia
# Find parameter points where PNJL.solve_multi(FixedMu) fails to find *any* physical solution.
#
# This script is intentionally separate from unit tests:
# - It explores parameter space.
# - It may be slow.
# - It is used to collect counterexamples (repro points) for improving seed banks / solver robustness.
#
# Run:
#   julia --project=. scripts/pnjl/find_multiseed_failure_points.jl
#
# Optional ENV knobs (defaults are moderate):
#   N_SAMPLES=200
#   SEED=1
#   P_NUM=12
#   T_NUM=4
#   T_MIN=0    T_MAX=350
#   MUB_MIN=0  MUB_MAX=1800
#   XI_MIN=-0.8 XI_MAX=0.8
#   STOP_ON_FIRST=1
#   OUT_CSV= (empty = no csv)
#
# Output:
# - Prints any failure points with a one-line repro command.
# - If OUT_CSV is set, writes all sampled points and whether solve_multi succeeded.

using Random
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .PNJL
using .PNJL.ConstraintModes: FixedMu

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

function write_csv_header(io)
    println(io, join([
        "idx",
        "T_MeV","muB_MeV","muq_MeV","xi",
        "p_num","t_num",
        "ok",
        "err",
        "Phi","Phibar",
        "m_u_fm","m_d_fm","m_s_fm",
        "omega","pressure","rho_norm",
        "residual_norm","iterations",
    ], ','))
end

function write_csv_row(io, idx, T_MeV, muB_MeV, xi, p_num, t_num; ok::Bool, err::String,
    Phi::Float64=NaN, Phibar::Float64=NaN,
    m_u::Float64=NaN, m_d::Float64=NaN, m_s::Float64=NaN,
    omega::Float64=NaN, pressure::Float64=NaN, rho_norm::Float64=NaN,
    residual_norm::Float64=NaN, iterations::Int=0)

    muq_MeV = muB_MeV / 3.0

    row = Any[
        idx,
        T_MeV, muB_MeV, muq_MeV, xi,
        p_num, t_num,
        ok,
        replace(err, '\n' => ' '),
        Phi, Phibar,
        m_u, m_d, m_s,
        omega, pressure, rho_norm,
        residual_norm, iterations,
    ]
    println(io, join(string.(row), ','))
end

function main()
    n_samples = env_int("N_SAMPLES", 20000)
    seed = env_int("SEED", 1)
    p_num = env_int("P_NUM", 12)
    t_num = env_int("T_NUM", 4)

    T_min = env_float("T_MIN", 0.0)
    T_max = env_float("T_MAX", 350.0)
    muB_min = env_float("MUB_MIN", 0.0)
    muB_max = env_float("MUB_MAX", 1800.0)
    xi_min = env_float("XI_MIN", -0.8)
    xi_max = env_float("XI_MAX", 0.8)

    stop_on_first = env_bool("STOP_ON_FIRST", true)
    out_csv = get(ENV, "OUT_CSV", "")

    rng = MersenneTwister(seed)

    @printf("find_multiseed_failure_points: N=%d seed=%d\n", n_samples, seed)
    @printf("ranges: T=[%.1f, %.1f] MeV  muB=[%.1f, %.1f] MeV  xi=[%.3f, %.3f]\n", T_min, T_max, muB_min, muB_max, xi_min, xi_max)
    @printf("grid: p_num=%d t_num=%d\n", p_num, t_num)

    io = nothing
    if !isempty(out_csv)
        isdir(dirname(out_csv)) || mkpath(dirname(out_csv))
        io = open(out_csv, "w")
        write_csv_header(io)
        @printf("csv: %s\n", out_csv)
    end

    n_ok = 0
    n_fail = 0

    try
        for i in 1:n_samples
            T_MeV = T_min + rand(rng) * (T_max - T_min)
            muB_MeV = muB_min + rand(rng) * (muB_max - muB_min)
            xi = xi_min + rand(rng) * (xi_max - xi_min)

            # Convert to internal units (fm^-1)
            # Avoid exact T=0.
            T_fm = max(T_MeV, 1e-6) / ħc_MeV_fm
            muq_fm = (muB_MeV / 3.0) / ħc_MeV_fm

            ok = false
            err_str = ""

            local res
            try
                res = PNJL.solve_multi(FixedMu(), T_fm, muq_fm; xi=xi, p_num=p_num, t_num=t_num)
                ok = true
                n_ok += 1

                if io !== nothing
                    write_csv_row(io, i, T_MeV, muB_MeV, xi, p_num, t_num;
                        ok=true, err="",
                        Phi=Float64(res.x_state[4]), Phibar=Float64(res.x_state[5]),
                        m_u=Float64(res.masses[1]), m_d=Float64(res.masses[2]), m_s=Float64(res.masses[3]),
                        omega=Float64(res.omega), pressure=Float64(res.pressure), rho_norm=Float64(res.rho_norm),
                        residual_norm=Float64(res.residual_norm), iterations=Int(res.iterations))
                end
            catch err
                n_fail += 1
                err_str = string(err)

                @printf("[%d] FAIL: solve_multi error at T=%.6f MeV muB=%.6f MeV xi=%.6f : %s\n",
                    i, T_MeV, muB_MeV, xi, err_str)
                @printf("repro: julia --project=. scripts/pnjl/diagnose_gap_point.jl --T=%.6f --muB=%.6f --xi=%.6f --p_num=%d --t_num=%d --verbose\n",
                    T_MeV, muB_MeV, xi, p_num, t_num)

                if io !== nothing
                    write_csv_row(io, i, T_MeV, muB_MeV, xi, p_num, t_num; ok=false, err=err_str)
                end

                if stop_on_first
                    return 2
                end
            end
        end
    finally
        io === nothing || close(io)
    end

    println("\nsummary:")
    @printf("  ok:   %d/%d\n", n_ok, n_samples)
    @printf("  fail: %d\n", n_fail)

    return (n_fail == 0) ? 0 : 1
end

exit(main())
