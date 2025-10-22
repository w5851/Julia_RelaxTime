
using Test
include("../src/QuarkDistribution.jl")

using .PNJLQuarkDistributions

# 中心差分数值导数
function numerical_derivative(f, E, h)
    return (f(E + h) - f(E - h)) / (2h)
end

@testset "quark_distribution_antiderivative correctness" begin
    μ = 0.1
    T = 0.15
    Φ = 0.5
    Φbar = 0.5

    # 选取几个能量点测试
    for E in (0.2, 0.5, 1.0, 2.0)
        h = 1e-6
        anti = E -> quark_distribution_antiderivative(E, μ, T, Φ, Φbar)
        numeric = numerical_derivative(anti, E, h)
        analytic = quark_distribution(E, μ, T, Φ, Φbar)
        @test isapprox(numeric, analytic; atol=1e-6, rtol=1e-6)
    end
end
