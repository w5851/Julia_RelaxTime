#!/usr/bin/env julia

using FastGaussQuadrature
using Printf

println("Starting quadrature_sigma_singularity.jl ...")
flush(stdout)

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

"""
Sweep over Δ values to test sensitivity of the s-integral to the integration interval width.
For each Δ in `delta_values` compute the hybrid reference (with given n_t) and
print the value and elapsed time.
"""
function test_delta_sweep(; process=:ssbar_to_uubar, delta_values=[0.01,0.1,0.5,1.0,5.0,10.0],
    n_t=16, n_ref_jacobi=128, n_ref_leg=256, split_frac=0.25, progress_every::Int=0, max_seconds::Float64=Inf)
    @printf("Building params (A, couplings) ...\n")
    flush(stdout)
    params = build_params()
    m_i = params.quark_params.m.s
    m_j = params.quark_params.m.s

    println("Δ\tvalue\telapsed_s")
    for Δ in delta_values
        a = (m_i + m_j)^2
        b = a + Δ
        t = @elapsed val = reference_integral_hybrid(a, b; process=process, params=params, split_frac=split_frac,
            n_ref_jacobi=n_ref_jacobi, n_ref_leg=n_ref_leg, n_t=n_t,
            progress_every=progress_every, max_seconds=max_seconds)
        @printf("%.3f\t%.10e\t%.3f\n", Δ, val, t); flush(stdout)
    end
end
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src", "relaxtime"))

println("Including Constants_PNJL.jl ..."); flush(stdout)
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
println("Including GaussLegendre.jl ..."); flush(stdout)
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
println("Including OneLoopIntegrals.jl ..."); flush(stdout)
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
println("Including EffectiveCouplings.jl ..."); flush(stdout)
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
println("Including TotalCrossSection.jl ..."); flush(stdout)
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: total_cross_section, DEFAULT_T_INTEGRAL_POINTS

println("Includes and using statements complete.")
flush(stdout)



"""
Build quark/thermo parameters matching scan_total_cross_section defaults.
"""
function build_params(; T_MeV=150.0, mu_u_MeV=0.0, mu_d_MeV=0.0, mu_s_MeV=0.0,
    m_u_MeV=300.0, m_d_MeV=300.0, m_s_MeV=500.0, phi=0.5, phibar=0.5, xi=0.0)
    T = T_MeV / ħc_MeV_fm
    μ_u = mu_u_MeV / ħc_MeV_fm
    μ_d = mu_d_MeV / ħc_MeV_fm
    μ_s = mu_s_MeV / ħc_MeV_fm
    m_u = m_u_MeV / ħc_MeV_fm
    m_d = m_d_MeV / ħc_MeV_fm
    m_s = m_s_MeV / ħc_MeV_fm

    nodes_p = DEFAULT_MOMENTUM_NODES
    weights_p = DEFAULT_MOMENTUM_WEIGHTS

    @printf("Computing A_u (m_u=%.6f) ...\n", m_u); flush(stdout)
    tA_u = @elapsed A_u = A(m_u, μ_u, T, phi, phibar, nodes_p, weights_p)
    @printf("Done A_u: elapsed=%.3fs\n", tA_u); flush(stdout)

    @printf("Computing A_d (m_d=%.6f) ...\n", m_d); flush(stdout)
    tA_d = @elapsed A_d = A(m_d, μ_d, T, phi, phibar, nodes_p, weights_p)
    @printf("Done A_d: elapsed=%.3fs\n", tA_d); flush(stdout)

    @printf("Computing A_s (m_s=%.6f) ...\n", m_s); flush(stdout)
    tA_s = @elapsed A_s = A(m_s, μ_s, T, phi, phibar, nodes_p, weights_p)
    @printf("Done A_s: elapsed=%.3fs\n", tA_s); flush(stdout)

    @printf("Calculating G_u/G_s from A_u/A_s ...\n"); flush(stdout)
    tG = @elapsed begin
        G_u = calculate_G_from_A(A_u)
        G_s = calculate_G_from_A(A_s)
    end
    @printf("Done G calc: elapsed=%.3fs\n", tG); flush(stdout)

    quark_params = (m=(u=m_u, d=m_d, s=m_s), μ=(u=μ_u, d=μ_d, s=μ_s), A=(u=A_u, d=A_d, s=A_s))
    thermo_params = (T=T, Φ=phi, Φbar=phibar, ξ=xi)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

"""
Reference: dense Gauss–Legendre integral of σ(s) over [a,b].
"""
function reference_integral_legendre(a, b; process=:ssbar_to_uubar, params, n_ref=256, n_t=DEFAULT_T_INTEGRAL_POINTS,
    progress_every::Int=0, max_seconds::Float64=Inf)
    x, w = gausslegendre(n_ref)
    s = @. (b - a) / 2 * (x + 1) + a
    jac = (b - a) / 2
    acc = 0.0
    start_t = time()
    total = length(w)
    for (idx, (si, wi)) in enumerate(zip(s, w))
        acc += wi * total_cross_section(process, si, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_t)
        if progress_every > 0 && idx % progress_every == 0
            @printf("  [leg-ref] %d/%d elapsed=%.1fs\r", idx, total, time() - start_t); flush(stdout)
        end
        if time() - start_t > max_seconds
            @warn "Legendre reference hit time limit" idx total max_seconds elapsed=time() - start_t
            break
        end
    end
    if progress_every > 0
        println()
    end
    return acc * jac
end

function reference_integral_quadgk(a, b; process=:ssbar_to_uubar, params, n_t=DEFAULT_T_INTEGRAL_POINTS)
    error("quadgk-based reference removed; use hybrid/legendre/jacobi reference instead")
end

"""
High-order Gauss–Jacobi reference on [a,b] with weight (s-a)^(-1/2).
This is efficient when the only singularity is the left endpoint.
"""
function reference_integral_jacobi(a, b; process=:ssbar_to_uubar, params, n_ref=256, n_t=DEFAULT_T_INTEGRAL_POINTS,
    progress_every::Int=0, max_seconds::Float64=Inf)
    x, w = gaussjacobi(n_ref, 0.0, -0.5)
    s = @. (b - a) / 2 * (x + 1) + a
    pref = sqrt((b - a) / 2)
    acc = 0.0
    start_t = time()
    total = length(w)
    for (idx, (si, wi)) in enumerate(zip(s, w))
        acc += wi * total_cross_section(process, si, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_t) * sqrt(si - a)
        if progress_every > 0 && idx % progress_every == 0
            @printf("  [jac-ref] %d/%d elapsed=%.1fs\r", idx, total, time() - start_t); flush(stdout)
        end
        if time() - start_t > max_seconds
            @warn "Jacobi reference hit time limit" idx total max_seconds elapsed=time() - start_t
            break
        end
    end
    if progress_every > 0
        println()
    end
    return pref * acc
end

"""
Hybrid reference: a small left panel with Gauss–Jacobi to absorb the (s-a)^(-1/2)
singularity, remainder with Gauss–Legendre. Avoids deep quadgk recursion.
"""
function reference_integral_hybrid(a, b; process=:ssbar_to_uubar, params, split_frac=0.25,
    n_ref_jacobi=128, n_ref_leg=512, n_t=DEFAULT_T_INTEGRAL_POINTS,
    progress_every::Int=0, max_seconds::Float64=Inf)
    mid = a + split_frac * (b - a)
    part1 = reference_integral_jacobi(a, mid; process=process, params=params, n_ref=n_ref_jacobi, n_t=n_t,
        progress_every=progress_every, max_seconds=max_seconds)
    part2 = reference_integral_legendre(mid, b; process=process, params=params, n_ref=n_ref_leg, n_t=n_t,
        progress_every=progress_every, max_seconds=max_seconds)
    return part1 + part2
end

"""
Adaptive Gauss–Legendre reference: start at n0 and double until relative & absolute
changes are within tolerances or max_n reached. Avoids quadgk recursion overhead.
"""
function reference_integral_legendre_adaptive(a, b; process=:ssbar_to_uubar, params,
    n0=64, max_n=4096, rtol=1e-6, atol=1e-8, n_t=DEFAULT_T_INTEGRAL_POINTS,
    progress_every::Int=0, max_seconds::Float64=Inf)
    n = n0
    start_t = time()
    prev = reference_integral_legendre(a, b; process=process, params=params, n_ref=n, n_t=n_t,
        progress_every=progress_every, max_seconds=max_seconds)
    while true
        n = min(2n, max_n)
        cur = reference_integral_legendre(a, b; process=process, params=params, n_ref=n, n_t=n_t,
            progress_every=progress_every, max_seconds=max_seconds)
        diff = abs(cur - prev)
        if diff <= max(rtol * abs(cur), atol)
            return cur
        end
        if n == max_n
            @warn "Adaptive Legendre did not reach tolerance" n_max=max_n diff=diff rtol=rtol atol=atol
            return cur
        end
        if time() - start_t > max_seconds
            @warn "Adaptive Legendre hit time limit" elapsed=time() - start_t max_seconds
            return cur
        end
        prev = cur
    end
end

"""
Plain Gauss–Legendre on σ(s).
"""
function legendre_integral(a, b, n; process=:ssbar_to_uubar, params, n_t=DEFAULT_T_INTEGRAL_POINTS)
    x, w = gausslegendre(n)
    s = @. (b - a) / 2 * (x + 1) + a
    jac = (b - a) / 2
    acc = 0.0
    for (si, wi) in zip(s, w)
        acc += wi * total_cross_section(process, si, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_t)
    end
    return acc * jac
end

"""
Gauss–Jacobi with weight (s - a)^(-1/2); integrand g(s) = σ(s) * sqrt(s - a).
"""
function jacobi_integral(a, b, n; process=:ssbar_to_uubar, params, n_t=DEFAULT_T_INTEGRAL_POINTS)
    x, w = gaussjacobi(n, 0.0, -0.5)  # weight = (1 + x)^(-1/2)
    s = @. (b - a) / 2 * (x + 1) + a
    pref = sqrt((b - a) / 2)
    acc = 0.0
    for (si, wi) in zip(s, w)
        acc += wi * total_cross_section(process, si, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_t) * sqrt(si - a)
    end
    return pref * acc
end

function run_experiment(; process=:ssbar_to_uubar, Δ=0.5, orders=(1,2,4,8,16,32,64,128), n_ref=256, n_t=16,
    ref_method=:jacobi, rtol=1e-6, atol=1e-8, split_frac=0.25, n_ref_jacobi=64, n_ref_leg=128,
    progress_every=32, max_seconds=60.0)
    @printf("Building params (A, couplings) ...\n")
    flush(stdout)
    t_build = @elapsed params = build_params()
    @printf("Done building params (elapsed=%.3fs)\n", t_build)
    flush(stdout)
    m_i = params.quark_params.m.s
    m_j = params.quark_params.m.s
    a = (m_i + m_j)^2  # initial-state threshold dominates for ssbar -> uubar
    b = a + Δ

    @printf("Process %s, interval [%.6f, %.6f], Δ=%.3f\n", string(process), a, b, Δ)
    flush(stdout)
    @printf("Using n_t=%d (t integral nodes)\n", n_t)
    flush(stdout)
    timed(f) = @elapsed f()

    ref = nothing
    if ref_method == :legendre
        @printf("Computing reference with Gauss-Legendre n_ref=%d ...\n", n_ref)
        t = timed() do
            ref = reference_integral_legendre(a, b; process=process, params=params, n_ref=n_ref, n_t=n_t,
                progress_every=progress_every, max_seconds=max_seconds)
        end
        @printf("reference time (Legendre n=%d): %.2fs\n", n_ref, t)
    elseif ref_method == :jacobi
        @printf("Computing reference with Gauss-Jacobi n_ref=%d ...\n", n_ref_jacobi)
        t = timed() do
            ref = reference_integral_jacobi(a, b; process=process, params=params, n_ref=n_ref_jacobi, n_t=n_t,
                progress_every=progress_every, max_seconds=max_seconds)
        end
        @printf("reference time (Jacobi n=%d): %.2fs\n", n_ref_jacobi, t)
    elseif ref_method == :hybrid
        @printf("Computing reference with hybrid (Jacobi left %.0f%%, Legendre right) n_ref_jacobi=%d, n_ref_leg=%d ...\n", split_frac*100, n_ref_jacobi, n_ref_leg)
        t = timed() do
            ref = reference_integral_hybrid(a, b; process=process, params=params, split_frac=split_frac,
                n_ref_jacobi=n_ref_jacobi, n_ref_leg=n_ref_leg, n_t=n_t,
                progress_every=progress_every, max_seconds=max_seconds)
        end
        @printf("reference time (hybrid): %.2fs\n", t)
    else
        @printf("Computing reference with adaptive Gauss-Legendre (start n0=%d, max_n=4096, rtol=%.1e, atol=%.1e)...\n", 64, rtol, atol)
        t = timed() do
            ref = reference_integral_legendre_adaptive(a, b; process=process, params=params, n0=64, max_n=4096, rtol=rtol, atol=atol, n_t=n_t,
                progress_every=progress_every, max_seconds=max_seconds)
        end
        @printf("reference time (adaptive Legendre): %.2fs\n", t)
    end
    @printf("Reference integral: %.10e\n", ref)
    println("n\tmethod\t\tapprox\t\trel_error")

    for n in orders
        approx = legendre_integral(a, b, n; process=process, params=params, n_t=n_t)
        rel = abs(approx - ref) / abs(ref)
        @printf("%d\tLegendre\t%.10e\t%.3e\n", n, approx, rel)
    end

    for n in orders
        approx = jacobi_integral(a, b, n; process=process, params=params, n_t=n_t)
        rel = abs(approx - ref) / abs(ref)
        @printf("%d\tJacobi\t\t%.10e\t%.3e\n", n, approx, rel)
    end
end



"""
Evaluate sensitivity of the s-integral result to the number of t-integration nodes.
Computes a high-precision baseline using `ref_n_t` and then evaluates the hybrid
reference at each `n_t_values`, printing differences and timings.
"""
function test_t_nodes(; process=:ssbar_to_uubar, Δ=0.5, n_t_values=[4,8,16,32,64],
    ref_n_t=256, n_ref_jacobi=128, n_ref_leg=256, split_frac=0.25,
    progress_every::Int=0, max_seconds::Float64=Inf)
    @printf("Building params (A, couplings) ...\n")
    flush(stdout)
    params = build_params()
    m_i = params.quark_params.m.s
    m_j = params.quark_params.m.s
    a = (m_i + m_j)^2
    b = a + Δ

    @printf("Computing baseline reference with ref_n_t=%d ...\n", ref_n_t); flush(stdout)
    tbase = @elapsed base = reference_integral_hybrid(a, b; process=process, params=params,
        split_frac=split_frac, n_ref_jacobi=n_ref_jacobi, n_ref_leg=n_ref_leg, n_t=ref_n_t,
        progress_every=progress_every, max_seconds=max_seconds)
    @printf("Baseline (n_t=%d): %.10e (elapsed=%.3fs)\n", ref_n_t, base, tbase); flush(stdout)

    println("n_t\tvalue\tdiff\trel_diff\telapsed_s")
    for nt in n_t_values
        t = @elapsed val = reference_integral_hybrid(a, b; process=process, params=params,
            split_frac=split_frac, n_ref_jacobi=n_ref_jacobi, n_ref_leg=n_ref_leg, n_t=nt,
            progress_every=progress_every, max_seconds=max_seconds)
        diff = val - base
        rel = diff / base
        @printf("%d\t%.10e\t%.3e\t%.3e\t%.3f\n", nt, val, diff, rel, t); flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Quick test run for developer feedback: smaller n_t/n_ref and a short max_seconds.
    # Override these if you need full-precision reference results.
    run_experiment(ref_method=:hybrid, Δ=10, n_ref=256, n_ref_jacobi=256, n_ref_leg=256, n_t=64,
                   progress_every=8, max_seconds=15.0)
end

"""
Sweep over Δ values and t-node counts. For each Δ compute a high-precision
baseline using `ref_n_t`, then evaluate the hybrid reference at each `n_t` in
`n_t_values`, printing value, absolute and relative differences and timing.
"""
function test_delta_vs_t(; process=:ssbar_to_uubar, delta_values=[0.01,0.1,0.5,1.0],
    n_t_values=[4,8,16,32,64], ref_n_t=256, n_ref_jacobi=128, n_ref_leg=256,
    split_frac=0.25, progress_every::Int=0, max_seconds::Float64=Inf)
    @printf("Building params (A, couplings) ...\n")
    flush(stdout)
    params = build_params()
    m_i = params.quark_params.m.s
    m_j = params.quark_params.m.s

    for Δ in delta_values
        a = (m_i + m_j)^2
        b = a + Δ
        @printf("\nΔ = %.6f: computing baseline with ref_n_t=%d ...\n", Δ, ref_n_t); flush(stdout)
        tbase = @elapsed base = reference_integral_hybrid(a, b; process=process, params=params,
            split_frac=split_frac, n_ref_jacobi=n_ref_jacobi, n_ref_leg=n_ref_leg, n_t=ref_n_t,
            progress_every=progress_every, max_seconds=max_seconds)
        @printf("Baseline (n_t=%d): %.10e (elapsed=%.3fs)\n", ref_n_t, base, tbase); flush(stdout)

        println("n_t\tvalue\tabs_diff\trel_diff\telapsed_s")
        for nt in n_t_values
            t = @elapsed val = reference_integral_hybrid(a, b; process=process, params=params,
                split_frac=split_frac, n_ref_jacobi=n_ref_jacobi, n_ref_leg=n_ref_leg, n_t=nt,
                progress_every=progress_every, max_seconds=max_seconds)
            diff = val - base
            rel = diff / base
            @printf("%d\t%.10e\t%.3e\t%.3e\t%.3f\n", nt, val, diff, rel, t); flush(stdout)
        end
    end
end