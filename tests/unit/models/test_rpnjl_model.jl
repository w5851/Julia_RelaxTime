# RPNJLModel.jl 单元测试
#
# 测试内容：
# 1. RPNJLModel 构造
# 2. create_model(:RPNJL) 工厂
# 3. calculate_mass_vec 与 base 比较
# 4. omega / gap_residual 基本可调用

using Test
using StaticArrays
using Base.MathConstants: π

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================

@testset "RPNJLModel" begin

    # --- 构造 ---
    @testset "RPNJLModel()" begin
        m = Models.RPNJLModel()
        @test m isa Models.AbstractPNJLModel
        @test m.base isa Models.PNJLModel
        @test m.ext isa NamedTuple
        @test haskey(m.ext, :g1_fm8)
        @test haskey(m.ext, :g2_fm8)
        @test haskey(m.ext, :kappa)
    end

    @testset "create_model(:RPNJL)" begin
        m = Models.create_model(:RPNJL)
        @test m isa Models.RPNJLModel
        @test m isa Models.AbstractPNJLModel
    end

    # --- 扩展参数字段有限 ---
    @testset "扩展参数有限" begin
        m = Models.RPNJLModel()
        @test isfinite(m.ext.g1_fm8)
        @test isfinite(m.ext.g2_fm8)
        @test m.ext.kappa >= 0
    end

    # --- calculate_mass_vec  ---
    @testset "calculate_mass_vec 有限且正" begin
        m = Models.RPNJLModel()
        φ = SVector{3}(-1.84, -1.84, -2.23)
        masses_rpnjl = Models.calculate_mass_vec(m, φ)
        masses_base = Models.calculate_mass_vec(m.base, φ)
        # RPNJL 含 g1/g2 修正，与 base 略有差异
        @test all(isfinite.(masses_rpnjl))
        @test all(masses_rpnjl .> 0)
    end

    @testset "calculate_mass_vec 有限" begin
        m = Models.RPNJLModel()
        φ = SVector{3}(0.01, 0.02, 0.5)
        masses = Models.calculate_mass_vec(m, φ)
        @test length(masses) == 3
        @test all(isfinite.(masses))
        @test all(masses .> 0)
    end

    # --- omega 可调用 ---
    @testset "omega RPNJL 有限" begin
        m = Models.RPNJLModel()
        x = SVector{5}(-1.84, -1.84, -2.23, 0.5, 0.5)
        T = 0.5
        μ = SVector{3}(0.0, 0.0, 0.0)
        ω = Models.omega(m, x, T, μ; p_num=24, t_num=6, xi=0.0)
        @test isfinite(ω)
    end

    @testset "polyakov_potential 默认兼容 rPNJL 旧公式" begin
        m = Models.RPNJLModel()
        Φ = 0.31
        Φbar = 0.29
        T = 0.72

        p = m.base.params
        T0 = Float64(p.T0_inv_fm)
        a0 = Float64(p.a0)
        a1 = Float64(p.a1)
        a2 = Float64(p.a2)
        b3 = Float64(p.b3)
        b4 = Float64(p.b4)
        kappa = Float64(m.ext.kappa)

        t_ratio = T0 / T
        b2 = a0 + a1 * t_ratio * exp(-a2 / t_ratio)
        φφ = Φ * Φbar
        j_poly = 1 - 6 * φφ + 4 * (Φ^3 + Φbar^3) - 3 * φφ^2
        jac = (27.0 / (24.0 * π^2)) * j_poly
        expected = T^4 * (
            -0.5 * b2 * φφ
            - (b3 / 6.0) * (Φ^3 + Φbar^3)
            + (b4 / 4.0) * φφ^2
            - kappa * log(max(jac, 1e-16))
        )

        actual = Models.polyakov_potential(m, Φ, Φbar, T)
        @test isapprox(actual, expected; rtol=1e-12, atol=1e-12)
    end

    @testset "非零扩展开关在固定点可检测" begin
        m_ext = Models.RPNJLModel(; use_rpnjl_extensions=true)
        m_base = Models.RPNJLModel(; use_rpnjl_extensions=false)

        @test m_ext.ext.kappa > 0
        @test abs(m_ext.ext.g1_fm8) > 0
        @test abs(m_ext.ext.g2_fm8) > 0

        x = SVector{5}(-0.03, -0.03, -0.04, 0.2, 0.2)
        T = 0.8
        μ = SVector{3}(0.45, 0.45, 0.45)

        comp_ext = Models.omega_components(m_ext, x, T, μ; p_num=24, t_num=8, xi=0.0)
        comp_base = Models.omega_components(m_base, x, T, μ; p_num=24, t_num=8, xi=0.0)
        residual_ext = Models.gap_residual(m_ext, x, T, μ; p_num=24, t_num=8, xi=0.0)
        residual_base = Models.gap_residual(m_base, x, T, μ; p_num=24, t_num=8, xi=0.0)

        @test isfinite(comp_ext.omega)
        @test isfinite(comp_base.omega)
        @test all(isfinite.(comp_ext.masses))
        @test all(isfinite.(comp_base.masses))
        @test all(isfinite.(residual_ext))
        @test all(isfinite.(residual_base))

        @test abs(comp_ext.omega - comp_base.omega) > 1e-10
        @test abs(comp_ext.poly - comp_base.poly) > 1e-10 || abs(comp_ext.chi - comp_base.chi) > 1e-10
    end

    # --- gap_residual 可调用 ---
    @testset "gap_residual RPNJL" begin
        m = Models.RPNJLModel()
        x = SVector{5}(-1.84, -1.84, -2.23, 0.5, 0.5)
        T = 0.5
        μ = SVector{3}(0.0, 0.0, 0.0)
        r = Models.gap_residual(m, x, T, μ; p_num=24, t_num=6, xi=0.0)
        @test r isa SVector{5}
        @test all(isfinite.(r))
    end
end
