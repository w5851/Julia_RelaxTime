
using Test
using Logging
include("../../../src/QuarkDistribution.jl")

# 避免跨文件导出冲突：不把符号 `using` 进 Main。
QD = PNJLQuarkDistributions

# 中心差分数值导数
function numerical_derivative(f, E, h)
    return (f(E + h) - f(E - h)) / (2h)
end

@testset "quark_distribution_antiderivative correctness" begin
    μ = 0.1
    T = 0.15
    Φ = 0.5
    Φbar = 0.5

    # 物理上 E 应为非负；同时避免在数值上不够光滑的区域做导数一致性检验。
    for E in (0.01, 0.2, 0.5, 1.0, 2.0, 5.0)
        h = 1e-5
        anti = E -> QD.quark_distribution_antiderivative(E, μ, T, Φ, Φbar)
        @info "quark_distribution_antiderivative" E=E value=anti(E)
        numeric = numerical_derivative(anti, E, h)
        analytic = QD.quark_distribution(E, μ, T, Φ, Φbar)
        @test isapprox(numeric, analytic; atol=1e-5, rtol=1e-5)
    end
end


