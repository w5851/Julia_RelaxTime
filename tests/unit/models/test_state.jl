# Models state.jl 单元测试
#
# 测试内容：
# 1. MeanFieldState 各构造方式
# 2. state_vector 往返一致性
# 3. normalize_mu_vec 各输入形式

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================
# MeanFieldState 构造
# ============================================================================

@testset "MeanFieldState" begin

    @testset "SVector + kwargs 构造" begin
        φ = SVector{3}(0.1, 0.2, 0.3)
        st = Models.MeanFieldState(φ; Phi=0.5, PhiBar=0.6)
        @test st.phi == φ
        @test st.Phi == 0.5
        @test st.PhiBar == 0.6
    end

    @testset "SVector 默认 Phi=PhiBar=1" begin
        φ = SVector{3}(0.1, 0.2, 0.3)
        st = Models.MeanFieldState(φ)
        @test st.Phi == 1.0
        @test st.PhiBar == 1.0
    end

    @testset "从 5 维向量构造" begin
        x = [0.1, 0.2, 0.3, 0.4, 0.5]
        st = Models.MeanFieldState(x)
        @test st.phi ≈ SVector{3}(0.1, 0.2, 0.3)
        @test st.Phi ≈ 0.4
        @test st.PhiBar ≈ 0.5
    end

    @testset "从 3 维向量构造 (NJL mode)" begin
        x = [0.1, 0.2, 0.3]
        st = Models.MeanFieldState(x)
        @test st.phi ≈ SVector{3}(0.1, 0.2, 0.3)
        @test st.Phi == 1.0
        @test st.PhiBar == 1.0
    end

    @testset "无效维度抛异常" begin
        x = [0.1, 0.2]
        @test_throws ArgumentError Models.MeanFieldState(x)
    end

    @testset "NamedTuple 构造 (φ key)" begin
        nt = (φ=SVector{3}(0.1, 0.2, 0.3), Φ=0.4, Φbar=0.5)
        st = Models.MeanFieldState(nt)
        @test st.phi ≈ SVector{3}(0.1, 0.2, 0.3)
        @test st.Phi ≈ 0.4
        @test st.PhiBar ≈ 0.5
    end
end

# ============================================================================
# state_vector 往返
# ============================================================================

@testset "state_vector" begin
    φ = SVector{3}(0.1, 0.2, 0.3)
    st = Models.MeanFieldState(φ; Phi=0.4, PhiBar=0.5)
    v = Models.state_vector(st)
    @test v ≈ SVector{5}(0.1, 0.2, 0.3, 0.4, 0.5)

    # Round-trip
    st2 = Models.MeanFieldState(v)
    @test st2.phi ≈ st.phi
    @test st2.Phi ≈ st.Phi
    @test st2.PhiBar ≈ st.PhiBar
end

# ============================================================================
# normalize_mu_vec
# ============================================================================

@testset "normalize_mu_vec" begin
    @testset "Real → SVector{3}" begin
        μ = Models.normalize_mu_vec(1.5)
        @test μ ≈ SVector{3}(1.5, 1.5, 1.5)
    end

    @testset "SVector{3} 直接传透" begin
        μ = SVector{3}(1.0, 2.0, 3.0)
        @test Models.normalize_mu_vec(μ) === μ
    end

    @testset "AbstractVector 长度 3" begin
        μ = Models.normalize_mu_vec([1.0, 2.0, 3.0])
        @test μ ≈ SVector{3}(1.0, 2.0, 3.0)
    end

    @testset "长度不为 3 抛异常" begin
        @test_throws ArgumentError Models.normalize_mu_vec([1.0, 2.0])
    end
end
