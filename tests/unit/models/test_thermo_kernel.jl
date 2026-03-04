# Models thermo_kernel.jl 单元测试
#
# 测试内容：
# 1. model_pressure 基本计算
# 2. model_rho 数密度向量 (ForwardDiff)
# 3. model_thermo 热力学量返回

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================
# model_pressure
# ============================================================================

@testset "thermo_kernel" begin

    @testset "model_pressure NJL" begin
        m = Models.create_model(:NJL)
        φ = SVector{3}(0.01, 0.02, 0.5)
        T = 0.8
        μ = SVector{3}(0.0, 0.0, 0.0)

        P = Models.model_pressure(m, φ, μ, T; p_num=24, t_num=6, xi=0.0)
        @test isfinite(P)
        # P = -Ω, 应与 omega 符号相反
        Ω = Models.omega(m, φ, T, μ; p_num=24, t_num=6, xi=0.0)
        @test P ≈ -Ω rtol=1e-12
    end

    # ============================================================================
    # model_rho
    # ============================================================================

    @testset "model_rho NJL" begin
        m = Models.create_model(:NJL)
        φ = SVector{3}(0.01, 0.02, 0.5)
        T = 0.8
        μ = SVector{3}(1.0, 1.0, 1.0)

        ρ = Models.model_rho(m, φ, μ, T; p_num=24, t_num=6, xi=0.0)
        @test length(ρ) == 3
        @test all(isfinite.(ρ))
    end

    @testset "model_rho μ=0 时密度为零" begin
        m = Models.create_model(:NJL)
        φ = SVector{3}(0.01, 0.02, 0.5)
        T = 0.8
        μ = SVector{3}(0.0, 0.0, 0.0)

        ρ = Models.model_rho(m, φ, μ, T; p_num=24, t_num=6, xi=0.0)
        # 在 μ=0 且 u/d 近似对称时数密度应接近零
        @test all(abs.(ρ) .< 0.1)
    end

    # ============================================================================
    # model_thermo
    # ============================================================================

    @testset "model_thermo 返回四元组" begin
        m = Models.create_model(:NJL)
        φ = SVector{3}(0.01, 0.02, 0.5)
        T = 0.8
        μ = SVector{3}(0.0, 0.0, 0.0)

        pressure, rho_norm, entropy, energy = Models.model_thermo(m, φ, μ, T; p_num=24, t_num=6, xi=0.0)
        @test isfinite(pressure)
        @test isfinite(rho_norm)
        @test isfinite(entropy)
        @test isfinite(energy)
    end
end
