# transport_provider.jl 单元测试
#
# 测试内容：
# 1. energy_from_p / energy_from_p_aniso 色散关系
# 2. TransportProvider 构造
# 3. transport_provider(model) 工厂
# 4. _mass_from_quark_params / _mu_from_quark_params 分发

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================

@testset "transport_provider" begin

    # --- energy_from_p ---
    @testset "energy_from_p 色散关系" begin
        # E = sqrt(p² + m²)
        @test Models.energy_from_p(0.0, 1.0) ≈ 1.0
        @test Models.energy_from_p(3.0, 4.0) ≈ 5.0  # 3-4-5 直角三角形
        @test Models.energy_from_p(1.0, 0.0) ≈ 1.0
    end

    @testset "energy_from_p_aniso" begin
        # ξ=0 → 同 isotropic
        @test Models.energy_from_p_aniso(3.0, 4.0, 0.0, 1.0) ≈ 5.0
        # ξ>0 → E > E_iso (for cosθ≠0)
        E_aniso = Models.energy_from_p_aniso(3.0, 4.0, 1.0, 1.0)
        @test E_aniso > 5.0
        # cosθ=0 → 同 isotropic
        @test Models.energy_from_p_aniso(3.0, 4.0, 1.0, 0.0) ≈ 5.0
    end

    # --- TransportProvider 构造 ---
    @testset "TransportProvider struct" begin
        @test isdefined(Models, :TransportProvider)
    end

    # --- transport_provider(PNJLModel) ---
    @testset "transport_provider PNJLModel" begin
        m = Models.create_model(:PNJL)
        tp = Models.transport_provider(m)
        @test tp isa Models.TransportProvider
        # 分布函数字段存在
        @test tp.energy_from_p isa Function
        @test tp.quark_distribution isa Function
        @test tp.antiquark_distribution isa Function
        @test tp.quark_distribution_aniso isa Function
        @test tp.antiquark_distribution_aniso isa Function
        @test tp.mass_for_species isa Function
        @test tp.mu_for_species isa Function
    end

    # --- transport_provider 不支持的模型 ---
    @testset "transport_provider NJL 抛出 ArgumentError" begin
        m = Models.create_model(:NJL)
        @test_throws ArgumentError Models.transport_provider(m)
    end

    # --- _mass_from_quark_params ---
    @testset "质量/化学势 species 分发" begin
        qp = (m=(u=0.01, d=0.01, s=0.5), μ=(u=0.3, d=0.3, s=0.3))
        @test Models._mass_from_quark_params(:u, qp) ≈ 0.01
        @test Models._mass_from_quark_params(:s, qp) ≈ 0.5
        @test Models._mass_from_quark_params(:ubar, qp) ≈ 0.01
        @test_throws ArgumentError Models._mass_from_quark_params(:unknown, qp)

        @test Models._mu_from_quark_params(:u, qp) ≈ 0.3
        @test Models._mu_from_quark_params(:sbar, qp) ≈ 0.3
        @test_throws ArgumentError Models._mu_from_quark_params(:unknown, qp)
    end
end
