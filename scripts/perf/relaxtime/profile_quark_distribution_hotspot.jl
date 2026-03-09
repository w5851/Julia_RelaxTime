using BenchmarkTools
using Profile

const _QUARK_DISTRIBUTION_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "QuarkDistribution.jl"))
if !isdefined(Main, :PNJLQuarkDistributions)
    Base.include(Main, _QUARK_DISTRIBUTION_PATH)
end

const QD = Main.PNJLQuarkDistributions

@inline function _pnjl_quark_distribution_rescaled(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    x = (E_inv_fm - μ_inv_fm) / T_inv_fm
    if x >= 0.0
        z = exp(-x)
        z2 = z * z
        z3 = z2 * z
        return (Φ * z + 2.0 * Φbar * z2 + z3) / (1.0 + 3.0 * Φ * z + 3.0 * Φbar * z2 + z3)
    else
        y = exp(x)
        y2 = y * y
        y3 = y2 * y
        return (Φ * y2 + 2.0 * Φbar * y + 1.0) / (y3 + 3.0 * Φ * y2 + 3.0 * Φbar * y + 1.0)
    end
end

@inline function _pnjl_antiquark_distribution_rescaled(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    x = (E_inv_fm + μ_inv_fm) / T_inv_fm
    if x >= 0.0
        z = exp(-x)
        z2 = z * z
        z3 = z2 * z
        return (Φbar * z + 2.0 * Φ * z2 + z3) / (1.0 + 3.0 * Φbar * z + 3.0 * Φ * z2 + z3)
    else
        y = exp(x)
        y2 = y * y
        y3 = y2 * y
        return (Φbar * y2 + 2.0 * Φ * y + 1.0) / (y3 + 3.0 * Φbar * y2 + 3.0 * Φ * y + 1.0)
    end
end

function _build_cases()
    Φ = 0.45
    Φbar = 0.43
    moderate = [(E=E, μ=μ, T=T, Φ=Φ, Φbar=Φbar) for E in range(0.4, 4.0; length=96), μ in (0.0, 0.15, 0.3), T in (0.12, 0.18)]
    extreme = [
        (E=0.2, μ=2.0, T=0.005, Φ=0.2, Φbar=0.25),
        (E=8.0, μ=0.0, T=0.01, Φ=0.3, Φbar=0.35),
        (E=1.0, μ=1.8, T=0.02, Φ=0.6, Φbar=0.55),
        (E=10.0, μ=0.3, T=0.02, Φ=0.8, Φbar=0.8),
    ]
    return moderate, extreme
end

const MODERATE_CASES, EXTREME_CASES = _build_cases()

function _sum_current_quark(cases)
    acc = 0.0
    @inbounds for c in cases
        acc += QD.quark_distribution(c.E, c.μ, c.T, c.Φ, c.Φbar)
    end
    return acc
end

function _sum_rescaled_quark(cases)
    acc = 0.0
    @inbounds for c in cases
        acc += _pnjl_quark_distribution_rescaled(c.E, c.μ, c.T, c.Φ, c.Φbar)
    end
    return acc
end

function _sum_current_antiquark(cases)
    acc = 0.0
    @inbounds for c in cases
        acc += QD.antiquark_distribution(c.E, c.μ, c.T, c.Φ, c.Φbar)
    end
    return acc
end

function _sum_rescaled_antiquark(cases)
    acc = 0.0
    @inbounds for c in cases
        acc += _pnjl_antiquark_distribution_rescaled(c.E, c.μ, c.T, c.Φ, c.Φbar)
    end
    return acc
end

function _max_abs_diff_quark(cases)
    maxdiff = 0.0
    @inbounds for c in cases
        cur = QD.quark_distribution(c.E, c.μ, c.T, c.Φ, c.Φbar)
        alt = _pnjl_quark_distribution_rescaled(c.E, c.μ, c.T, c.Φ, c.Φbar)
        maxdiff = max(maxdiff, abs(cur - alt))
    end
    return maxdiff
end

function _max_abs_diff_antiquark(cases)
    maxdiff = 0.0
    @inbounds for c in cases
        cur = QD.antiquark_distribution(c.E, c.μ, c.T, c.Φ, c.Φbar)
        alt = _pnjl_antiquark_distribution_rescaled(c.E, c.μ, c.T, c.Φ, c.Φbar)
        maxdiff = max(maxdiff, abs(cur - alt))
    end
    return maxdiff
end

function _profile_current_quark(cases, n_outer::Int)
    Profile.clear()
    @profile begin
        acc = 0.0
        for _ in 1:n_outer
            acc += _sum_current_quark(cases)
        end
        println("profile checksum(current quark) = ", acc)
    end
    Profile.print(maxdepth=12)
end

function _profile_rescaled_quark(cases, n_outer::Int)
    Profile.clear()
    @profile begin
        acc = 0.0
        for _ in 1:n_outer
            acc += _sum_rescaled_quark(cases)
        end
        println("profile checksum(rescaled quark) = ", acc)
    end
    Profile.print(maxdepth=12)
end

function main()
    println("QuarkDistribution hotspot profiling")
    println("  moderate cases = ", length(MODERATE_CASES))
    println("  extreme cases  = ", length(EXTREME_CASES))
    println()

    println("[correctness / stability]")
    println("  max |Δ quark| (moderate)      = ", _max_abs_diff_quark(MODERATE_CASES))
    println("  max |Δ antiquark| (moderate)  = ", _max_abs_diff_antiquark(MODERATE_CASES))
    println("  max |Δ quark| (extreme)       = ", _max_abs_diff_quark(EXTREME_CASES))
    println("  max |Δ antiquark| (extreme)   = ", _max_abs_diff_antiquark(EXTREME_CASES))
    println()

    println("[benchmark: quark moderate loop]")
    show(stdout, MIME("text/plain"), @benchmark _sum_current_quark($MODERATE_CASES))
    println()
    show(stdout, MIME("text/plain"), @benchmark _sum_rescaled_quark($MODERATE_CASES))
    println("\n")

    println("[benchmark: antiquark moderate loop]")
    show(stdout, MIME("text/plain"), @benchmark _sum_current_antiquark($MODERATE_CASES))
    println()
    show(stdout, MIME("text/plain"), @benchmark _sum_rescaled_antiquark($MODERATE_CASES))
    println("\n")

    println("[benchmark: quark extreme loop]")
    show(stdout, MIME("text/plain"), @benchmark _sum_current_quark($EXTREME_CASES))
    println()
    show(stdout, MIME("text/plain"), @benchmark _sum_rescaled_quark($EXTREME_CASES))
    println("\n")

    println("[benchmark: antiquark extreme loop]")
    show(stdout, MIME("text/plain"), @benchmark _sum_current_antiquark($EXTREME_CASES))
    println()
    show(stdout, MIME("text/plain"), @benchmark _sum_rescaled_antiquark($EXTREME_CASES))
    println("\n")

    println("[profile: current quark moderate loop]")
    _profile_current_quark(MODERATE_CASES, 20_000)
    println("\n[profile: rescaled quark moderate loop]")
    _profile_rescaled_quark(MODERATE_CASES, 20_000)
end

main()