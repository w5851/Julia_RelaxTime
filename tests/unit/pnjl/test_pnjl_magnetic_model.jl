# PNJLMagneticModel 单元测试
#
# 测试内容：
# 1. PNJLMagneticModel 构造
# 2. 接口代理 (calculate_mass_vec, calculate_chiral, vacuum_contribution)
# 3. omega 与 solve_gap 调用

using Test
using StaticArrays
using ForwardDiff

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

    @testset "磁场五维 residual contract" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.1)
        x5 = SVector{5}(-0.03, -0.03, -0.04, 0.2, 0.2)
        μ = SVector{3}(0.4, 0.4, 0.4)
        residual = Models.magnetic_gap_residual(m, x5, 0.7, μ;
            p_num=4, t_num=4, n_max=1, pz_max=5.0)
        @test length(residual) == 5
        @test all(isfinite, residual)
        @test_throws ArgumentError Models.magnetic_gap_residual(m, x5, 0.0, μ;
            p_num=4, t_num=4, n_max=1, pz_max=5.0)
        @test_throws ArgumentError Models.magnetic_gap_residual(m, x5, 0.7, μ;
            xi=0.1, p_num=4, t_num=4, n_max=1, pz_max=5.0)
    end

    @testset "Hessian 只提供诊断，不拒绝驻点" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.1)
        x5 = SVector{5}(-0.03, -0.03, -0.04, 0.2, 0.2)
        μ = SVector{3}(0.4, 0.4, 0.4)
        result = Models.solve_magnetic_gap(
            m,
            0.7,
            μ;
            p_num=4,
            t_num=4,
            pz_max=5.0,
            n_max=1,
            initial_guess=x5,
            include_default_seeds=false,
            iterations=8,
            residual_norm_max=1e-3,
            classify_stability=true,
        )
        @test result.stability_classified
        @test result.converged
        @test result.state !== nothing
        @test !isempty(result.candidates)
        @test all(candidate -> candidate.stability !== :not_evaluated, result.candidates)
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

    @testset "磁场质量核保留泛型实数类型" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.1)
        phi = SVector{3}(-0.03, -0.03, -0.04)
        jac = ForwardDiff.jacobian(v -> Models.calculate_mass_vec(m, SVector{3}(v)), phi)
        @test size(jac) == (3, 3)
        @test all(isfinite, jac)
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

    @testset "零场 number_densities 保持 PNJL 独立分量" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.0)
        x5 = SVector{5}(-1.5, -1.5, -2.0, 0.3, 0.3)
        μ = SVector{3}(0.1, 0.1, 0.1)
        nd = Models.number_densities(m, x5, 0.7, μ; p_num=4, t_num=4, xi=0.0)
        @test length(nd.quark) == 3
        @test length(nd.antiquark) == 3
        @test all(isfinite, nd.quark)
        @test all(isfinite, nd.antiquark)
    end

    @testset "非零磁场公共热力学入口" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.1)
        x5 = SVector{5}(-0.03, -0.03, -0.04, 0.2, 0.2)
        μ = SVector{3}(0.4, 0.4, 0.4)
        rho = Models.model_rho(m, x5, μ, 0.7; p_num=4, t_num=4, pz_max=5.0, n_max=1)
        thermo = Models.model_thermo(m, x5, μ, 0.7; p_num=4, t_num=4, pz_max=5.0, n_max=1)
        @test length(rho) == 3
        @test all(isfinite, rho)
        @test all(isfinite, thermo)
    end

    @testset "磁场配置节点与共享约束边界" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.1, p_num=4, pz_max=5.0, n_max=1)
        x5 = SVector{5}(-0.03, -0.03, -0.04, 0.2, 0.2)
        mu = SVector{3}(0.4, 0.4, 0.4)

        rho_default = Models.model_rho(m, x5, mu, 0.7; t_num=4, pz_max=5.0, n_max=1)
        rho_explicit = Models.model_rho(m, x5, mu, 0.7; p_num=4, t_num=4, pz_max=5.0, n_max=1)
        @test rho_default ≈ rho_explicit rtol=1e-12 atol=1e-12

        @test_throws ArgumentError Models.solve_constraint(
            m,
            Models.FixedMu(),
            0.7;
            μ_fm=0.4,
            seed_guess=collect(x5),
            p_num=4,
            t_num=4,
        )

        @test !Models.model_capabilities(m).supports_number_densities
        @test Models.model_capabilities(Models.PNJLMagneticModel(; eB_fm2=0.0)).supports_number_densities
        @test_throws Models.UnsupportedCapabilityError Models.require_capability(m, :number_densities)
    end

    # ============================================================================
    # create_model(:PNJLMagnetic) 工厂路径
    # ============================================================================

    @testset "create_model(:PNJLMagnetic)" begin
        m = Models.create_model(:PNJLMagnetic; eB_fm2=0.0)
        @test m isa Models.PNJLMagneticModel
    end
end
