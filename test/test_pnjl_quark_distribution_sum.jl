"""
验证PNJL模型中相同参数下夸克与反夸克分布函数之和为1。
"""

using Test
include("../src/QuarkDistribution.jl")

using .NJLQuarkDistributions

@testset "NJL quark + antiquark normalization" begin
    μ = 0.0
    T = 0.15
    energies = (0.0, 0.2, 0.5, 1.0, 2.0)
    for E in energies
        quark = quark_distribution(E, μ, T)
        antiquark = antiquark_distribution(E, μ, T)
        @test isapprox(quark + antiquark, 1.0; atol=1e-12, rtol=1e-12)
    end
end
