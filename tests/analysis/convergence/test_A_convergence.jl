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

function A_pmax(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64; pmax::Float64, n::Int)
    nodes, weights = gauleg(0.0, pmax, n)
    return A(m, μ, T, Φ, Φbar, nodes, weights)
end

function min_nodes_to_match_ref(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64;
    pmax::Float64, ref::Float64, tol::Float64, candidates::Vector{Int})
    for n in candidates
        v = A_pmax(m, μ, T, Φ, Φbar; pmax=pmax, n=n)
        if rel_diff(v, ref) < tol
            return n, v
        end
    end
    # 都未满足：返回最大候选
    n = candidates[end]
    return n, A_pmax(m, μ, T, Φ, Φbar; pmax=pmax, n=n)
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

# 覆盖 ak_compare 中更“敏感”的参数区：高温、轻质量、Φ 接近 1
# 这里更容易出现“pmax=10 不够、需要更大有效截断”的现象。
T2_fm = 1.5202959509451173
μ2_fm = 0.0
Φ2 = 0.840579181
Φ2bar = 0.840579181
m2_u = 0.0404709972
m2_s = 1.03140506

nodes2_ref, weights2_ref = gauleg(0.0, 20.0, 256)
A2_u_ref = A(m2_u, μ2_fm, T2_fm, Φ2, Φ2bar, nodes2_ref, weights2_ref)
A2_s_ref = A(m2_s, μ2_fm, T2_fm, Φ2, Φ2bar, nodes2_ref, weights2_ref)

nodes2_10_32, weights2_10_32 = gauleg(0.0, 10.0, 32)
A2_u_10_32 = A(m2_u, μ2_fm, T2_fm, Φ2, Φ2bar, nodes2_10_32, weights2_10_32)
A2_s_10_32 = A(m2_s, μ2_fm, T2_fm, Φ2, Φ2bar, nodes2_10_32, weights2_10_32)

nodes2_15_128, weights2_15_128 = gauleg(0.0, 15.0, 128)
A2_u_15_128 = A(m2_u, μ2_fm, T2_fm, Φ2, Φ2bar, nodes2_15_128, weights2_15_128)
A2_s_15_128 = A(m2_s, μ2_fm, T2_fm, Φ2, Φ2bar, nodes2_15_128, weights2_15_128)

@testset "A convergence (ak_compare-like, p<=15 improves)" begin
    @test rel_diff(A2_u_15_128, A2_u_ref) < TOL
    @test rel_diff(A2_s_15_128, A2_s_ref) < TOL
    @test rel_diff(A2_u_15_128, A2_u_ref) <= rel_diff(A2_u_10_32, A2_u_ref)
    @test rel_diff(A2_s_15_128, A2_s_ref) <= rel_diff(A2_s_10_32, A2_s_ref)
end

# 可选诊断输出：打印 pmax=15 vs 20 的差异，以及在候选节点数中达到 TOL 所需的最小节点数。
# 默认关闭，避免 CI/批量测试输出噪声。
const REPORT = get(ENV, "A_CONV_REPORT", "0") != "0"
if REPORT
    candidates = [16, 32, 64, 128, 256]
    ref_n = 512

    function report_case(name::String, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64, m_u::Float64, m_s::Float64)
        println("\n=== A convergence report: ", name, " ===")
        println("T=", T, " μ=", μ, " Φ=", Φ, " Φbar=", Φbar)
        println("candidates=", candidates, " ref(pmax=20,n=", ref_n, ")")

        Au_ref = A_pmax(m_u, μ, T, Φ, Φbar; pmax=20.0, n=ref_n)
        As_ref = A_pmax(m_s, μ, T, Φ, Φbar; pmax=20.0, n=ref_n)

        n15_u, Au15 = min_nodes_to_match_ref(m_u, μ, T, Φ, Φbar; pmax=15.0, ref=Au_ref, tol=TOL, candidates=candidates)
        n20_u, Au20 = min_nodes_to_match_ref(m_u, μ, T, Φ, Φbar; pmax=20.0, ref=Au_ref, tol=TOL, candidates=candidates)
        n15_s, As15 = min_nodes_to_match_ref(m_s, μ, T, Φ, Φbar; pmax=15.0, ref=As_ref, tol=TOL, candidates=candidates)
        n20_s, As20 = min_nodes_to_match_ref(m_s, μ, T, Φ, Φbar; pmax=20.0, ref=As_ref, tol=TOL, candidates=candidates)

        println("u: n_min(pmax=15)=", n15_u, "  rel(15 vs ref20)=", rel_diff(Au15, Au_ref))
        println("u: n_min(pmax=20)=", n20_u, "  rel(20 vs ref20)=", rel_diff(Au20, Au_ref))
        println("u: trunc(rel 15 vs 20_at_nmin)=", rel_diff(Au15, Au20), "  dA(15-20)=", (Au15 - Au20))

        println("s: n_min(pmax=15)=", n15_s, "  rel(15 vs ref20)=", rel_diff(As15, As_ref))
        println("s: n_min(pmax=20)=", n20_s, "  rel(20 vs ref20)=", rel_diff(As20, As_ref))
        println("s: trunc(rel 15 vs 20_at_nmin)=", rel_diff(As15, As20), "  dA(15-20)=", (As15 - As20))
    end

    report_case("baseline", T_fm, μ_fm, Φ, Φbar, m_u, m_s)
    report_case("ak_compare-like", T2_fm, μ2_fm, Φ2, Φ2bar, m2_u, m2_s)
end
