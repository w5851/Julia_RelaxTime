# PNJLModel 单元测试
#
# 测试内容：
# 1. PNJLModel 构造
# 2. calculate_mass_vec / calculate_chiral 接口
# 3. omega_components 计算
# 4. number_densities 接口

using Test
using ForwardDiff
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================
# PNJLModel 构造
# ============================================================================

@testset "PNJLModel" begin

    @testset "默认构造" begin
        m = Models.PNJLModel()
        @test m isa Models.AbstractPNJLModel
        @test m isa Models.AbstractQCDModel
        @test isdefined(m, :params)
        @test isdefined(m, :consts)
    end

    @testset "profile 构造" begin
        m = Models.PNJLModel(; profile="default")
        @test m isa Models.PNJLModel
    end

    # ============================================================================
    # calculate_mass_vec
    # ============================================================================

    @testset "calculate_mass_vec" begin
        m = Models.PNJLModel()
        φ = SVector{3}(0.01, 0.02, 0.5)
        masses = Models.calculate_mass_vec(m, φ)
        @test length(masses) == 3
        @test all(isfinite.(masses))
        @test all(masses .> 0)  # 质量应为正

        # u/d 裸质量相近（m_ud0 相同），φ_u ≈ φ_d 时质量接近
        φ_sym = SVector{3}(0.01, 0.01, 0.5)
        masses_sym = Models.calculate_mass_vec(m, φ_sym)
        @test masses_sym[1] ≈ masses_sym[2] rtol=1e-12
    end

    # ============================================================================
    # calculate_chiral
    # ============================================================================

    @testset "calculate_chiral" begin
        m = Models.PNJLModel()
        φ = SVector{3}(0.01, 0.02, 0.5)
        chi = Models.calculate_chiral(m, φ)
        @test isfinite(chi)
    end

    # ============================================================================
    # omega_components
    # ============================================================================

    @testset "omega_components" begin
        m = Models.PNJLModel()
        x5 = SVector{5}(-1.5, -1.5, -2.0, 0.3, 0.3)
        T = 0.5
        μ = SVector{3}(0.0, 0.0, 0.0)

        comp = Models.omega_components(m, x5, T, μ; p_num=24, t_num=6, xi=0.0)

        @test all(isfinite, (comp.chi, comp.poly, comp.vac, comp.therm, comp.omega))
        @test comp.omega ≈ (comp.chi + comp.poly + comp.vac + comp.therm) rtol=1e-10
    end

    # ============================================================================
    # number_densities
    # ============================================================================

    @testset "number_densities" begin
        m = Models.PNJLModel()
        st = Models.MeanFieldState(SVector{3}(-1.5, -1.5, -2.0); Phi=0.3, PhiBar=0.3)
        T = 0.5
        μ = SVector{3}(1.0, 1.0, 1.0)

        nd = Models.number_densities(m, st, T, μ; p_num=24, t_num=6, xi=0.0)
        @test haskey(nd, :quark)
        @test haskey(nd, :antiquark)
        @test length(nd.quark) == 3
        @test length(nd.antiquark) == 3
        @test all(isfinite.(nd.quark))
        @test all(isfinite.(nd.antiquark))
    end

    @testset "自适应热积分与直接密度使用同一 RS 分布口径" begin
        m = Models.PNJLModel()
        st = Models.MeanFieldState(SVector{3}(-1.5, -1.5, -2.0); Phi=0.3, PhiBar=0.25)
        x_state = Models.state_vector(st)
        T = 0.2
        μ = SVector(0.4, 0.35, 0.2)
        xi = 0.2
        quadrature = (
            thermo_quadrature_policy=:rs_reduced_adaptive,
            thermo_quadrature_rtol=1e-9,
            thermo_quadrature_atol=1e-11,
        )

        comp0 = Models.omega_components(m, x_state, T, μ; xi=0.0, quadrature...)
        compξ = Models.omega_components(m, x_state, T, μ; xi=xi, quadrature...)
        @test compξ.chi == comp0.chi
        @test compξ.poly == comp0.poly
        @test compξ.vac == comp0.vac
        @test compξ.therm ≈ comp0.therm / sqrt(1 + xi) rtol=5e-9

        nd = Models.number_densities(m, st, T, μ; xi=xi, quadrature...)
        rho_ad = ForwardDiff.gradient(μtrial -> -Models.omega(
            m, x_state, T, μtrial; xi=xi, quadrature...,
        ), μ)
        @test rho_ad ≈ nd.quark - nd.antiquark rtol=2e-7 atol=1e-10
    end

    @testset "严格零温只支持固定态内核" begin
        m = Models.PNJLModel()
        st = Models.MeanFieldState(SVector{3}(-1.5, -1.5, -2.0); Phi=0.3, PhiBar=0.25)
        μ = SVector(2.0, 2.0, 2.0)
        nd = Models.number_densities(
            m, st, 0.0, μ;
            xi=0.1,
            thermo_quadrature_policy=:rs_reduced_adaptive,
        )
        @test all(isfinite, nd.quark)
        @test all(iszero, nd.antiquark)
        @test Models.PNJLCore.polyakov_potential(m.params, 0.3, 0.25, 0.0) == 0.0
        @test_throws ArgumentError Models.solve_gap(
            m, 0.0, μ;
            thermo_quadrature_policy=:rs_reduced_adaptive,
        )
    end

    # ============================================================================
    # Polyakov 势
    # ============================================================================

    @testset "polyakov_potential" begin
        m = Models.PNJLModel()
        Φ = 0.3
        Φbar = 0.3
        T = 0.5

        V = Models.PNJLCore.polyakov_potential(m.params, Φ, Φbar, T)
        @test isfinite(V)

        # Φ=Φbar=0 在低温退化
        V0 = Models.PNJLCore.polyakov_potential(m.params, 0.0, 0.0, 0.1)
        @test isfinite(V0)
    end
end
