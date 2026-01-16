# average_scattering_rate 性能烟囱测试（smoke）
#
# 对应系统流程步骤：
# - `src/relaxtime/AverageScatteringRate.jl`
#
# 测试内容：
# - 构造常量截面 cache（避免把“截面计算”混入该基准）
# - `BenchmarkTools.@benchmark` 单样本计时
#
# 运行方式：
# - `julia --project=. tests/perf/relaxtime/test_average_scattering_rate_performance.jl`
using Test
using BenchmarkTools

include("../../../src/relaxtime/AverageScatteringRate.jl")

using .AverageScatteringRate

const QUARK_PARAMS = (
    m = (u = 1.52, d = 1.52, s = 2.50),
    μ = (u = 0.3, d = 0.3, s = 0.0)
)
const THERMO = (T = 0.15, Φ = 0.5, Φbar = 0.5, ξ = 0.0)
const K_COEFFS = (K_σπ=2.0, K_σK=2.0, K_σ=3.0, K_δπ=1.5, K_δK=1.5)

function constant_sigma_cache(process::Symbol; sigma::Float64=1.0)
    cache = CrossSectionCache(process)
    AverageScatteringRate.insert_sigma!(cache, 0.0, sigma)
    AverageScatteringRate.insert_sigma!(cache, 500.0, sigma)
    return cache
end

@testset "average_scattering_rate performance (smoke)" begin
    cache = constant_sigma_cache(:uu_to_uu; sigma=1.0)
    bench = @benchmark average_scattering_rate(
        :uu_to_uu,
        QUARK_PARAMS,
        THERMO,
        K_COEFFS;
        p_nodes=8,
        angle_nodes=4,
        phi_nodes=4,
        cs_cache=$cache,
        n_sigma_points=4,
    ) samples=1 evals=1

    @info "average_scattering_rate runtime (p=8,angle=4,phi=4)" time_ns = bench.times[1]
    @test bench.times[1] > 0
end
