#!/usr/bin/env julia

using Random
using Printf

function find_project_root(start_dir::AbstractString)
    dir = abspath(start_dir)
    while true
        if isfile(joinpath(dir, "Project.toml"))
            return dir
        end
        parent = dirname(dir)
        parent == dir && error("Could not find Project.toml from: $start_dir")
        dir = parent
    end
end

const PROJECT_ROOT = find_project_root(@__DIR__)

push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src", "relaxtime"))

include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
using .OneLoopIntegrals

"Reproduce the E1/E2 formula used in OneLoopIntegrals.singularity_k_positive for k>0."
function compute_E1E2(λ::Float64, k::Float64, m::Float64, mprime::Float64)
    d0 = k^2 - λ^2
    denom = 2.0 * d0
    A = λ^2 - k^2 + m^2 - mprime^2
    disc = A^2 + 4.0 * m^2 * d0
    disc < 0 && return (NaN, NaN)
    sqrt_disc = sqrt(disc)
    E1 = (λ * A - k * sqrt_disc) / denom
    E2 = (λ * A + k * sqrt_disc) / denom
    if E1 > E2
        E1, E2 = E2, E1
    end
    return (E1, E2)
end

"Emulate the contiguous (bo,up) interval selection in C++ BPM for the k>0 imaginary part."
function cpp_like_interval(E1::Float64, E2::Float64, Emin::Float64, Emax::Float64)
    in1 = (E1 - Emax) * (E1 - Emin) < 0.0
    in2 = (E2 - Emax) * (E2 - Emin) < 0.0
    if !(in1 || in2)
        return nothing
    end
    if in1 && in2
        bo = min(E1, E2)
        up = max(E1, E2)
    else
        up = Emax
        bo = in1 ? E1 : E2
    end
    return (bo, up)
end

function main()
    rng = MersenneTwister(20251216)

    function try_find(; ntry::Int, emax_mode::Symbol)
        for it in 1:ntry
            # Choose masses (fm^-1)
            m = rand(rng) * 3.0 + 0.05
            mprime = rand(rng) * 5.0 + 0.05

            # Choose k>0 (fm^-1)
            k = rand(rng) * 10.0 + 0.05

            # Force a = λ^2 - k^2 < 0 but keep |λ| close to k.
            ratio = 0.98 + rand(rng) * (0.999999 - 0.98)
            λ = (rand(rng, Bool) ? 1.0 : -1.0) * ratio * k

            Emin = m
            Emax = emax_mode === :physical ? OneLoopIntegrals.energy_cutoff(m) : 100.0

            intervals, sign_type = OneLoopIntegrals.singularity_k_positive(λ, k, m, mprime, Emin, Emax)
            if sign_type != :outside || length(intervals) != 2
                continue
            end

            E1, E2 = compute_E1E2(λ, k, m, mprime)
            if !isfinite(E1) || !isfinite(E2) || !(E1 < E2)
                continue
            end

            # Require both roots to lie strictly inside the integration window
            if !(E1 > Emin && E2 < Emax)
                continue
            end

            return (λ=λ, k=k, m=m, mprime=mprime, Emin=Emin, Emax=Emax, E1=E1, E2=E2, intervals=intervals)
        end
        return nothing
    end

    function stats(; ntry::Int, emax_mode::Symbol)
        n_outside = 0
        n_between = 0
        n_none = 0
        n_outside_len0 = 0
        n_outside_len1 = 0
        n_outside_len2 = 0

        max_E1_minus_Emin = -Inf
        min_E1_minus_Emin = +Inf
        max_E2_minus_Emax = -Inf
        min_E2_minus_Emax = +Inf

        for it in 1:ntry
            m = rand(rng) * 3.0 + 0.05
            mprime = rand(rng) * 5.0 + 0.05
            k = rand(rng) * 10.0 + 0.05
            ratio = 0.98 + rand(rng) * (0.999999 - 0.98)
            λ = (rand(rng, Bool) ? 1.0 : -1.0) * ratio * k

            Emin = m
            Emax = emax_mode === :physical ? OneLoopIntegrals.energy_cutoff(m) : 100.0

            intervals, sign_type = OneLoopIntegrals.singularity_k_positive(λ, k, m, mprime, Emin, Emax)
            if sign_type === :outside
                n_outside += 1
                if length(intervals) == 0
                    n_outside_len0 += 1
                elseif length(intervals) == 1
                    n_outside_len1 += 1
                elseif length(intervals) == 2
                    n_outside_len2 += 1
                end
            elseif sign_type === :between
                n_between += 1
            else
                n_none += 1
            end

            E1, E2 = compute_E1E2(λ, k, m, mprime)
            if isfinite(E1)
                max_E1_minus_Emin = max(max_E1_minus_Emin, E1 - Emin)
                min_E1_minus_Emin = min(min_E1_minus_Emin, E1 - Emin)
            end
            if isfinite(E2)
                max_E2_minus_Emax = max(max_E2_minus_Emax, E2 - Emax)
                min_E2_minus_Emax = min(min_E2_minus_Emax, E2 - Emax)
            end
        end

        return (ntry=ntry, emax_mode=emax_mode,
            n_outside=n_outside, n_between=n_between, n_none=n_none,
            n_outside_len0=n_outside_len0, n_outside_len1=n_outside_len1, n_outside_len2=n_outside_len2,
            max_E1_minus_Emin=max_E1_minus_Emin, min_E1_minus_Emin=min_E1_minus_Emin,
            max_E2_minus_Emax=max_E2_minus_Emax, min_E2_minus_Emax=min_E2_minus_Emax)
    end

    println("[1/2] Searching two-interval :outside case under physical cutoff Emax=sqrt(m^2+Λ^2)...")
    found = try_find(ntry=2_000_000, emax_mode=:physical)
    if found === nothing
        println("No two-interval case found under physical cutoff in 2e6 trials.")

        println("\nSTAT under physical cutoff (same sampling):")
        st1 = stats(ntry=200_000, emax_mode=:physical)
        @printf("ntry=%d  outside=%d  between=%d  none=%d\n", st1.ntry, st1.n_outside, st1.n_between, st1.n_none)
        @printf("outside interval counts: len0=%d len1=%d len2=%d\n", st1.n_outside_len0, st1.n_outside_len1, st1.n_outside_len2)
        @printf("E1-Emin range: [%.6f, %.6f]\n", st1.min_E1_minus_Emin, st1.max_E1_minus_Emin)
        @printf("E2-Emax range: [%.6f, %.6f]\n", st1.min_E2_minus_Emax, st1.max_E2_minus_Emax)

        println("\n[2/2] Now searching with an enlarged cutoff Emax=100...")
        found = try_find(ntry=2_000_000, emax_mode=:enlarged)
        if found === nothing
            println("No two-interval case found even with enlarged cutoff in 2e6 trials.")
            println("STAT under enlarged cutoff (same sampling):")
            st2 = stats(ntry=200_000, emax_mode=:enlarged)
            @printf("ntry=%d  outside=%d  between=%d  none=%d\n", st2.ntry, st2.n_outside, st2.n_between, st2.n_none)
            @printf("outside interval counts: len0=%d len1=%d len2=%d\n", st2.n_outside_len0, st2.n_outside_len1, st2.n_outside_len2)
            @printf("E1-Emin range: [%.6f, %.6f]\n", st2.min_E1_minus_Emin, st2.max_E1_minus_Emin)
            @printf("E2-Emax range: [%.6f, %.6f]\n", st2.min_E2_minus_Emax, st2.max_E2_minus_Emax)
            return
        end
    end

    a = found.λ^2 - found.k^2
    @printf("\nFOUND example (a<0 gives OUTSIDE intervals)\n")
    @printf("λ=%.6f  k=%.6f  a=λ^2-k^2=%.6f (<0)\n", found.λ, found.k, a)
    @printf("m=%.6f  m' =%.6f\n", found.m, found.mprime)
    @printf("Emin=m=%.6f  Emax=%.6f\n", found.Emin, found.Emax)
    @printf("E1=%.6f  E2=%.6f (E1<E2, both inside [Emin,Emax])\n", found.E1, found.E2)
    @printf("Julia intervals (outside): [%.6f, %.6f] U [%.6f, %.6f]\n",
        found.intervals[1][1], found.intervals[1][2], found.intervals[2][1], found.intervals[2][2])

    cpp = cpp_like_interval(found.E1, found.E2, found.Emin, found.Emax)
    if cpp === nothing
        @printf("C++-like: no interval (would set B_im=0)\n")
    else
        @printf("C++-like contiguous interval: (%.6f, %.6f)  (always one piece)\n", cpp[1], cpp[2])
        @printf("NOTE: In this example, C++ would integrate between E1 and E2,\n")
        @printf("      while Julia says the condition holds outside [E1,E2].\n")
    end
end

main()
