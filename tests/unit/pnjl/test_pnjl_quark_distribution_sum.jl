"""
验证PNJL模型中相同参数下对E和μ取反后夸克与反夸克分布函数之和为1。
"""

using Test
include("../../../src/QuarkDistribution.jl")

using .PNJLQuarkDistributions

@testset "PNJL quark + antiquark normalization" begin
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


