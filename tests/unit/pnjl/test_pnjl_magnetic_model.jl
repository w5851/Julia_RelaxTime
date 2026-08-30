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

    @testset "正磁场构造合同" begin
        @test_throws ArgumentError Models.PNJLMagneticModel(; eB_fm2=0.0)
        @test_throws ArgumentError Models.PNJLMagneticModel(; eB_fm2=-0.1)
        m = Models.PNJLMagneticModel(; eB_fm2=0.1)
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

    @testset "磁场 AD residual 生产路径" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.1)
        x5 = SVector{5}(-0.03, -0.03, -0.04, 0.2, 0.2)
        μ = SVector{3}(0.4, 0.4, 0.4)
        controls = (p_num=4, t_num=4, pz_max=5.0, n_max=1, xi=0.0)

        residual = Models.magnetic_gap_residual(m, x5, 0.7, μ; controls...)
        alias = Models.magnetic_gap_residual_autodiff(m, x5, 0.7, μ; controls...)
        @test alias ≈ residual rtol=1e-12 atol=1e-12
        @test all(isfinite, residual)

        jac = ForwardDiff.jacobian(
            v -> Models.magnetic_gap_residual(m, v, 0.7, μ; controls...),
            collect(x5),
        )
        @test size(jac) == (5, 5)
        @test all(isfinite, jac)

        @test_throws ArgumentError Models.magnetic_gap_residual(
            m, x5, 0.7, μ; p_num=4, t_num=4, pz_max=5.0, xi=0.0,
        )

        result = Models.solve_magnetic_gap(
            m,
            0.7,
            μ;
            p_num=4,
            t_num=4,
            pz_max=5.0,
            n_max=1,
            initial_guess=SVector{5}(-1.84, -1.84, -2.22, 0.01, 0.01),
            include_default_seeds=false,
            fallback_method=nothing,
            iterations=8,
            residual_norm_max=1e-3,
        )
        @test result.converged
        @test result.attempt_count == 1
        @test result.state !== nothing
        @test result.candidates[1].residual_norm <= 1e-3

        resolved_nmax = Models._resolve_magnetic_attempt_nmax(
            m,
            x5,
            0.7,
            μ;
            xi=0.0,
            p_num=4,
            t_num=4,
            pz_max=5.0,
            n_max=nothing,
            cutoff_N=10,
        )
        @test resolved_nmax >= 0

        # The thermal-tail profile is a physical-point policy, not a per-seed
        # policy. Distinct seeds must therefore resolve to the same Landau
        # layer count before the solver enters NLsolve.
        alternate_seed = SVector{5}(-1.1, -1.1, -1.7, 0.4, 0.4)
        resolved_alt = Models._resolve_magnetic_attempt_nmax(
            m,
            alternate_seed,
            0.7,
            μ;
            xi=0.0,
            p_num=4,
            t_num=4,
            pz_max=5.0,
            n_max=nothing,
            cutoff_N=10,
        )
        @test resolved_alt == resolved_nmax
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
            initial_guess=SVector{5}(-1.84, -1.84, -2.22, 0.01, 0.01),
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

    @testset "磁场质量核" begin
        m_mag = Models.PNJLMagneticModel(; eB_fm2=0.1)
        φ = SVector{3}(0.01, 0.02, 0.5)

        masses_mag = Models.calculate_mass_vec(m_mag, φ)
        @test all(isfinite, masses_mag)
        @test all(>(0.0), masses_mag)
    end

    @testset "磁场质量核保留泛型实数类型" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.1)
        phi = SVector{3}(-0.03, -0.03, -0.04)
        jac = ForwardDiff.jacobian(v -> Models.calculate_mass_vec(m, SVector{3}(v)), phi)
        @test size(jac) == (3, 3)
        @test all(isfinite, jac)
    end

    @testset "磁场手征势" begin
        m_mag = Models.PNJLMagneticModel(; eB_fm2=0.1)
        φ = SVector{3}(0.01, 0.02, 0.5)

        chi_mag = Models.calculate_chiral(m_mag, φ)
        @test isfinite(chi_mag)
    end

    # ============================================================================
    # omega 调用
    # ============================================================================

    @testset "omega 可调用" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.1)
        x5 = SVector{5}(-0.03, -0.03, -0.04, 0.2, 0.2)
        T = 0.5
        μ = SVector{3}(0.0, 0.0, 0.0)

        Ω = Models.omega(m, x5, T, μ; p_num=4, t_num=4, pz_max=5.0, n_max=1, xi=0.0)
        @test isfinite(Ω)
    end

    @testset "磁场 number_densities 使用净密度合同" begin
        m = Models.PNJLMagneticModel(; eB_fm2=0.1)
        x5 = SVector{5}(-0.03, -0.03, -0.04, 0.2, 0.2)
        μ = SVector{3}(0.1, 0.1, 0.1)
        nd = Models.number_densities(m, x5, 0.7, μ; p_num=4, t_num=4, pz_max=5.0, n_max=1, xi=0.0)
        @test length(nd.quark) == 3
        @test nd.antiquark === nothing
        @test nd.net === nd.quark
        @test all(isfinite, nd.quark)
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
        @test_throws Models.UnsupportedCapabilityError Models.require_capability(m, :number_densities)
    end

    # ============================================================================
    # create_model(:PNJLMagnetic) 工厂路径
    # ============================================================================

    @testset "create_model(:PNJLMagnetic)" begin
        @test_throws ArgumentError Models.create_model(:PNJLMagnetic; eB_fm2=0.0)
        m = Models.create_model(:PNJLMagnetic; eB_fm2=0.1)
        @test m isa Models.PNJLMagneticModel
    end
end
