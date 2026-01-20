# A 函数收敛性测试（p 上限 10，p_nodes=32）
# 运行：julia --project=. tests/analysis/convergence/test_A_convergence.jl

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .GaussLegendre: gauleg
using .OneLoopIntegrals: A

function rel_diff(a::Float64, b::Float64)
    denom = max(abs(b), 1e-12)
    return abs(a - b) / denom
end

const TOL = parse(Float64, get(ENV, "A_CONV_TOL", "5e-3"))

# 代表性参数（单位：fm）
T_fm = 200.0 / ħc_MeV_fm
μ_fm = 0.0
Φ = 0.5
Φbar = 0.5
m_u = 300.0 / ħc_MeV_fm
m_s = 500.0 / ħc_MeV_fm

# 参考值：较大的上限/节点数
nodes_ref, weights_ref = gauleg(0.0, 20.0, 256)
A_u_ref = A(m_u, μ_fm, T_fm, Φ, Φbar, nodes_ref, weights_ref)
A_s_ref = A(m_s, μ_fm, T_fm, Φ, Φbar, nodes_ref, weights_ref)

# 目标配置：上限 10，p_nodes=32
nodes_10_32, weights_10_32 = gauleg(0.0, 10.0, 32)
A_u_10_32 = A(m_u, μ_fm, T_fm, Φ, Φbar, nodes_10_32, weights_10_32)
A_s_10_32 = A(m_s, μ_fm, T_fm, Φ, Φbar, nodes_10_32, weights_10_32)

# 稍高配置：上限 10，p_nodes=64
nodes_10_64, weights_10_64 = gauleg(0.0, 10.0, 64)
A_u_10_64 = A(m_u, μ_fm, T_fm, Φ, Φbar, nodes_10_64, weights_10_64)
A_s_10_64 = A(m_s, μ_fm, T_fm, Φ, Φbar, nodes_10_64, weights_10_64)

@testset "A convergence (p<=10, p_nodes=32)" begin
    @test rel_diff(A_u_10_32, A_u_ref) < TOL
    @test rel_diff(A_s_10_32, A_s_ref) < TOL
    @test rel_diff(A_u_10_64, A_u_ref) < TOL
    @test rel_diff(A_s_10_64, A_s_ref) < TOL
end
