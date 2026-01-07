using Pkg
Pkg.instantiate()

include(joinpath(@__DIR__, "..", "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(@__DIR__, "..", "src", "relaxtime", "OneLoopIntegralsAniso.jl"))

using .OneLoopIntegrals
using .OneLoopIntegralsCorrection

# parameters (roughly typical)
const λ = 0.35
const k = 0.45
const m1 = 0.30
const m2 = 0.45
const μ1 = 0.10
const μ2 = 0.05
const T = 0.16
const Φ = 0.2
const Φbar = 0.2
const ξ = 0.10

# warmup
r1 = OneLoopIntegrals.B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
r2 = OneLoopIntegralsCorrection.B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
println("warmup B0=", r1)
println("warmup B0corr=", r2)

function bench_B0(n::Int)
    s = 0.0
    t = @elapsed begin
        for _ in 1:n
            r = OneLoopIntegrals.B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
            s += r[1] + r[2]
        end
    end
    return t, s
end

function bench_B0corr(n::Int)
    s = 0.0
    t = @elapsed begin
        for _ in 1:n
            r = OneLoopIntegralsCorrection.B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
            s += r[1] + r[2]
        end
    end
    return t, s
end

n = 300
(t1, s1) = bench_B0(n)
(t2, s2) = bench_B0corr(n)
println("BENCH n=", n)
println("B0 time_s=", t1, " checksum=", s1)
println("B0corr time_s=", t2, " checksum=", s2)
