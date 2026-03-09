# PNJLDistributions 单元测试
#
# 测试内容：
# 1. pnjl_quark_distribution 基本值域
# 2. pnjl_antiquark_distribution 基本值域
# 3. 各向异性分布函数
# 4. 物理极限行为

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# 分布函数在 Models 模块作用域（不在 PNJLCore 中）
const _PD = Models

# ============================================================================
# pnjl_quark_distribution 测试
# ============================================================================

@testset "PNJLDistributions" begin

    @testset "pnjl_quark_distribution 基本" begin
        E = 2.0   # fm⁻¹
        μ = 0.5   # fm⁻¹
        T = 0.5   # fm⁻¹
        Φ = 0.3
        Φbar = 0.3

        f_q = _PD.pnjl_quark_distribution(E, μ, T, Φ, Φbar)
        @test isfinite(f_q)
        @test 0.0 <= f_q <= 1.0  # occupation should be ∈ [0,1]
    end

    @testset "pnjl_antiquark_distribution 基本" begin
        E = 2.0
        μ = 0.5
        T = 0.5
        Φ = 0.3
        Φbar = 0.3

        f_qbar = _PD.pnjl_antiquark_distribution(E, μ, T, Φ, Φbar)
        @test isfinite(f_qbar)
        @test 0.0 <= f_qbar <= 1.0
    end

    @testset "T→0 极限：分布趋零 (E>μ)" begin
        E = 5.0
        μ = 1.0
        T = 0.001  # 极低温
        Φ = 0.3
        Φbar = 0.3

        f_q = _PD.pnjl_quark_distribution(E, μ, T, Φ, Φbar)
        @test f_q < 0.01  # 高能态在低温下几乎无占据
    end

    @testset "T→∞ 极限：分布趋 1/3 (confined→free)" begin
        E = 1.0
        μ = 0.0
        T = 100.0  # 极高温
        Φ = 1.0    # 解禁闭
        Φbar = 1.0

        f_q = _PD.pnjl_quark_distribution(E, μ, T, Φ, Φbar)
        @test f_q > 0.3  # 接近自由气体极限
    end

    @testset "Φ=Φbar=1 退化为标准 Fermi-Dirac" begin
        E = 2.0
        μ = 1.0
        T = 0.5
        Φ = 1.0
        Φbar = 1.0

        f_q = _PD.pnjl_quark_distribution(E, μ, T, Φ, Φbar)
        f_fd = 1.0 / (1.0 + exp((E - μ) / T))  # standard Fermi-Dirac
        @test f_q ≈ f_fd rtol=1e-10
    end

    @testset "Φ=Φbar=0 禁闭极限" begin
        E = 2.0
        μ = 1.0
        T = 0.5
        Φ = 0.0
        Φbar = 0.0

        f_q = _PD.pnjl_quark_distribution(E, μ, T, Φ, Φbar)
        f_qbar = _PD.pnjl_antiquark_distribution(E, μ, T, Φ, Φbar)
        # 在 Φ=Φbar=0 时分布函数退化但仍有限
        @test isfinite(f_q)
        @test isfinite(f_qbar)
    end

    @testset "极端低温高化学势稳定性" begin
        E = 0.2
        μ = 2.0
        T = 0.005
        Φ = 0.2
        Φbar = 0.25

        f_q = _PD.pnjl_quark_distribution(E, μ, T, Φ, Φbar)
        f_qbar = _PD.pnjl_antiquark_distribution(E, μ, T, Φ, Φbar)

        @test isfinite(f_q)
        @test isfinite(f_qbar)
        @test 0.0 <= f_q <= 1.0
        @test 0.0 <= f_qbar <= 1.0
        @test f_q > 0.999999
    end

    @testset "各向异性分布 aniso" begin
        p = 2.0
        m = 1.0
        μ = 0.5
        T = 0.5
        Φ = 0.3
        Φbar = 0.3
        ξ = 0.0
        cosθ = 1.0

        # ξ=0 应退化为各向同性
        f_aniso = _PD.pnjl_quark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)
        E_iso = sqrt(p^2 + m^2)
        f_iso = _PD.pnjl_quark_distribution(E_iso, μ, T, Φ, Φbar)
        @test f_aniso ≈ f_iso rtol=1e-12

        # ξ > 0 应压低分布（动量空间拉伸）
        f_aniso_xi = _PD.pnjl_quark_distribution_aniso(p, m, μ, T, Φ, Φbar, 1.0, cosθ)
        @test f_aniso_xi < f_iso  # larger effective energy → smaller occupation
    end

    @testset "反夸克各向异性" begin
        p = 2.0
        m = 1.0
        μ = 0.5
        T = 0.5
        Φ = 0.3
        Φbar = 0.3

        f_aq0 = _PD.pnjl_antiquark_distribution_aniso(p, m, μ, T, Φ, Φbar, 0.0, 1.0)
        E_iso = sqrt(p^2 + m^2)
        f_aq_ref = _PD.pnjl_antiquark_distribution(E_iso, μ, T, Φ, Φbar)
        @test f_aq0 ≈ f_aq_ref rtol=1e-12
    end
end
