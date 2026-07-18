# 当前 B0 hybrid 路径与独立高节点固定 GL oracle 的性能/精度探针。
#
# 本脚本只报告测量结果，不设置 correctness gate。高节点 oracle 必须同时查看
# 1024 -> 2048 节点漂移，不能把单次固定节点值当成自带误差估计。
#
# 运行：
#   julia --project=. scripts/perf/relaxtime/performance_comparison.jl

using Printf
using Statistics
using FastGaussQuadrature: gausslegendre

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))

const OLI = Main.OneLoopIntegrals
const OLIC = Main.OneLoopIntegralsCorrection

function fixed_gl_reference(f, a::Float64, b::Float64, n::Int)
    nodes, weights = gausslegendre(n)
    half = (b - a) / 2
    center = (a + b) / 2
    total = 0.0
    @inbounds for i in eachindex(nodes)
        total += weights[i] * f(center + half * nodes[i])
    end
    return half * total
end

const TEST_CASES = [
    (name="有根(标准)", λ=-1.0, k=0.01, m=0.3, m_prime=0.3),
    (name="有根(大k)", λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
    (name="无根", λ=-0.5, k=0.05, m=0.3, m_prime=0.3),
    (name="无根(正λ)", λ=0.5, k=0.1, m=0.3, m_prime=0.35),
]

function main(; repeats::Int=100)
    repeats > 0 || throw(ArgumentError("repeats must be positive, got $(repeats)"))
    ξ = -0.2
    T = 0.15
    μ = 0.0
    Φ = 0.0
    Φbar = 0.0

    println("B0 fixed-hybrid performance vs fixed-GL oracle")
    println(@sprintf("%-14s %12s %12s %12s %12s", "case", "hybrid(ms)", "relerr", "oracle_delta", "roots"))

    for tc in TEST_CASES
        Emin = tc.m
        Emax = OLI.energy_cutoff(tc.m)
        integrand(E) = OLIC.real_integrand_k_positive(
            :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar,
        )
        ref_1024 = fixed_gl_reference(integrand, Emin, Emax, 1024)
        ref_2048 = fixed_gl_reference(integrand, Emin, Emax, 2048)
        oracle_delta = abs(ref_2048 - ref_1024) / max(abs(ref_2048), eps(Float64))

        OLIC.tilde_B0_correction_k_positive(
            :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ,
        )
        times_ms = Float64[]
        result = (0.0, 0.0)
        for _ in 1:repeats
            start = time_ns()
            result = OLIC.tilde_B0_correction_k_positive(
                :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ,
            )
            push!(times_ms, (time_ns() - start) / 1.0e6)
        end

        relative_error = abs(result[1] - ref_2048) / max(abs(ref_2048), eps(Float64))
        roots = length(OLIC.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax))
        println(@sprintf(
            "%-14s %12.5f %12.3e %12.3e %12d",
            tc.name,
            mean(times_ms),
            relative_error,
            oracle_delta,
            roots,
        ))
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    repeats = isempty(ARGS) ? 100 : parse(Int, ARGS[1])
    main(; repeats=repeats)
end
