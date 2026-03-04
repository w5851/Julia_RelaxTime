"""
QuarkDistribution (PNJLQuarkDistributions) 模块单元测试

合并原 test_pnjl_quark_distribution_sum.jl, test_quark_distribution_antiderivative.jl
"""

using Test

const _QUARK_DISTRIBUTION_PATH_QD = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "models", "pnjl_physics", "QuarkDistribution.jl"))
if !isdefined(Main, :PNJLQuarkDistributions)
    Base.include(Main, _QUARK_DISTRIBUTION_PATH_QD)
end

using Main.PNJLQuarkDistributions

# 中心差分数值导数
function _numerical_derivative(f, E, h)
    return (f(E + h) - f(E - h)) / (2h)
end

@testset "QuarkDistribution" begin
    @testset "quark + antiquark normalization" begin
        μ = 1.0
        T = 0.15
        Φ = 0.5
        Φbar = 0.5
        energies = (0.0, 0.2, 0.5, 1.0, 2.0)
        for E in energies
            quark = quark_distribution(E, μ, T, Φ, Φbar)
            antiquark = antiquark_distribution(-E, μ, T, Φ, Φbar)
            @test isapprox(quark + antiquark, 1.0; atol=1e-12, rtol=1e-12)
        end
    end

    @testset "antiderivative correctness" begin
        μ = 0.1
        T = 0.15
        Φ = 0.5
        Φbar = 0.5

        for E in (0.01, 0.2, 0.5, 1.0, 2.0, 5.0)
            h = 1e-5
            anti = x -> quark_distribution_antiderivative(x, μ, T, Φ, Φbar)
            numeric = _numerical_derivative(anti, E, h)
            analytic = quark_distribution(E, μ, T, Φ, Φbar)
            @test isapprox(numeric, analytic; atol=1e-5, rtol=1e-5)
        end
    end
end
