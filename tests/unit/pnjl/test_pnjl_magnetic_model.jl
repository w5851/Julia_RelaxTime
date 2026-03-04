# PNJLMagneticModel 单元测试
#
# 测试内容：
# 1. PNJLMagneticModel 构造
# 2. 接口代理 (calculate_mass_vec, calculate_chiral, vacuum_contribution)
# 3. omega 与 solve_gap 调用

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================
# PNJLMagneticModel 构造
# ============================================================================

@testset "PNJLMagneticModel" begin

    @testset "默认构造 (eB=0)" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.0)
        @test m isa Models.AbstractPNJLModel
        @test isdefined(m, :base)
        @test m.base isa Models.PNJLModel
    end

    @testset "非零磁场构造" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.5)
        @test m isa Models.PNJLMagneticModel
    end

    # ============================================================================
    # 接口委托测试
    # ============================================================================

    @testset "calculate_mass_vec 委托 base" begin
        m_mag = Models.PNJLMagneticModel(; eB_fm2=0.0)
        m_base = m_mag.base
        φ = SVector{3}(0.01, 0.02, 0.5)

        masses_mag = Models.calculate_mass_vec(m_mag, φ)
        masses_base = Models.calculate_mass_vec(m_base, φ)
        @test masses_mag ≈ masses_base rtol=1e-14
    end

    @testset "calculate_chiral 委托 base" begin
        m_mag = Models.PNJLMagneticModel(; eB_fm2=0.0)
        φ = SVector{3}(0.01, 0.02, 0.5)

        chi_mag = Models.calculate_chiral(m_mag, φ)
        chi_base = Models.calculate_chiral(m_mag.base, φ)
        @test chi_mag ≈ chi_base rtol=1e-14
    end

    @testset "vacuum_contribution 委托 base" begin
        m_mag = Models.PNJLMagneticModel(; eB_fm2=0.0)
        φ = SVector{3}(0.01, 0.02, 0.5)
        masses = Models.calculate_mass_vec(m_mag, φ)

        vac_mag = Models.vacuum_contribution(m_mag, masses)
        vac_base = Models.vacuum_contribution(m_mag.base, masses)
        @test vac_mag ≈ vac_base rtol=1e-14
    end

    # ============================================================================
    # omega 调用
    # ============================================================================

    @testset "omega 可调用" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.0)
        x5 = SVector{5}(-1.5, -1.5, -2.0, 0.3, 0.3)
        T = 0.5
        μ = SVector{3}(0.0, 0.0, 0.0)

        Ω = Models.omega(m, x5, T, μ; p_num=24, t_num=6, xi=0.0)
        @test isfinite(Ω)
    end

    # ============================================================================
    # create_model(:PNJLMagnetic) 工厂路径
    # ============================================================================

    @testset "create_model(:PNJLMagnetic)" begin
        m = Models.create_model(:PNJLMagnetic; eB_fm2=0.0)
        @test m isa Models.PNJLMagneticModel
    end
end
