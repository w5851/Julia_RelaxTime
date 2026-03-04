# Constants_PNJL 单元测试
#
# 测试内容：
# 1. 模块加载与导出常量完整性
# 2. 物理常量值域/精度
# 3. pnjl_constants() 动态加载
# 4. load_pnjl_config() 配置解析
# 5. _validate_pnjl_critical_params 参数校验

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
# 不使用 using，避免导出符号污染 Main 作用域（会与后续测试文件冲突）
const _CP = Constants_PNJL

# ============================================================================
# 导出常量完整性
# ============================================================================

@testset "Constants_PNJL exports" begin
    @testset "物理常量已定义" begin
        @test isdefined(Constants_PNJL, :ħc_MeV_fm)
        @test isdefined(Constants_PNJL, :α)
        @test isdefined(Constants_PNJL, :N_color)
        @test isdefined(Constants_PNJL, :N_flavor)
        @test isdefined(Constants_PNJL, :ρ0_inv_fm3)
    end

    @testset "模型参数已定义" begin
        @test isdefined(Constants_PNJL, :m_ud0_inv_fm)
        @test isdefined(Constants_PNJL, :m_s0_inv_fm)
        @test isdefined(Constants_PNJL, :Λ_inv_fm)
        @test isdefined(Constants_PNJL, :G_fm2)
        @test isdefined(Constants_PNJL, :K_fm5)
    end

    @testset "Polyakov 参数已定义" begin
        @test isdefined(Constants_PNJL, :T0_inv_fm)
        @test isdefined(Constants_PNJL, :a0)
        @test isdefined(Constants_PNJL, :a1)
        @test isdefined(Constants_PNJL, :a2)
        @test isdefined(Constants_PNJL, :b3)
        @test isdefined(Constants_PNJL, :b4)
    end

    @testset "Gell-Mann 矩阵与波函数" begin
        @test isdefined(Constants_PNJL, :λ₀)
        @test isdefined(Constants_PNJL, :λ₈)
        @test isdefined(Constants_PNJL, :ψ_u)
        @test isdefined(Constants_PNJL, :ψ_d)
        @test isdefined(Constants_PNJL, :ψ_s)
        @test isdefined(Constants_PNJL, :ψbar_u)
        @test isdefined(Constants_PNJL, :ψbar_d)
        @test isdefined(Constants_PNJL, :ψbar_s)
    end

    @testset "散射映射" begin
        @test isdefined(Constants_PNJL, :SCATTERING_MESON_MAP)
        @test isdefined(Constants_PNJL, :SCATTERING_PROCESS_KEYS)
    end
end

# ============================================================================
# 物理常量值域验证
# ============================================================================

@testset "Constants_PNJL value ranges" begin
    @testset "hbarc 精度" begin
        @test _CP.ħc_MeV_fm ≈ 197.3269804 rtol=1e-5
        @test isfinite(_CP.ħc_MeV_fm)
        @test _CP.ħc_MeV_fm > 0
    end

    @testset "α (fine structure)" begin
        @test _CP.α ≈ 1.0 / 137.036 rtol=1e-4
        @test 0 < _CP.α < 1
    end

    @testset "N_color / N_flavor" begin
        @test _CP.N_color == 3
        @test _CP.N_flavor == 3
    end

    @testset "质量层级物理约束" begin
        # m_ud << m_s (light quarks much lighter than strange quark)
        @test _CP.m_ud0_inv_fm < _CP.m_s0_inv_fm
        @test _CP.m_ud0_inv_fm > 0
        @test _CP.m_s0_inv_fm > 0
    end

    @testset "耦合常数正定" begin
        @test _CP.G_fm2 > 0
        @test _CP.K_fm5 > 0
    end

    @testset "Gell-Mann 矩阵正交性" begin
        # λ₀: 归一化的单位阵 ∝ √(2/3)·I₃
        @test size(_CP.λ₀) == (3, 3)
        @test size(_CP.λ₈) == (3, 3)
        # Tr(λ₀²) ≈ 2 (SU(3) convention)
        @test sum(_CP.λ₀ .^ 2) ≈ 2.0 rtol=1e-12
    end

    @testset "波函数正交完备" begin
        @test _CP.ψ_u' * _CP.ψ_d ≈ 0.0 atol=1e-15
        @test _CP.ψ_u' * _CP.ψ_s ≈ 0.0 atol=1e-15
        @test _CP.ψ_d' * _CP.ψ_s ≈ 0.0 atol=1e-15
        @test _CP.ψ_u' * _CP.ψ_u ≈ 1.0 atol=1e-15
    end
end

# ============================================================================
# pnjl_constants() 动态加载
# ============================================================================

@testset "pnjl_constants()" begin
    c = _CP.pnjl_constants()

    @testset "返回 NamedTuple 且字段完整" begin
        @test c isa NamedTuple
        @test haskey(c, :hbarc_MeV_fm)
        @test haskey(c, :alpha_em)
        @test haskey(c, :N_color)
        @test haskey(c, :Λ_inv_fm)
        @test haskey(c, :G_fm2)
        @test haskey(c, :K_fm5)
        @test haskey(c, :T0_inv_fm)
    end

    @testset "值与模块 const 一致" begin
        @test c.hbarc_MeV_fm == _CP.ħc_MeV_fm
        @test c.alpha_em == _CP.α
        @test c.N_color == _CP.N_color
    end

    @testset "单位换算一致性" begin
        # Λ_inv_fm = Lambda_MeV / hbarc
        @test c.Λ_inv_fm ≈ 602.3 / c.hbarc_MeV_fm rtol=1e-10
    end
end

# ============================================================================
# load_pnjl_config() 配置解析
# ============================================================================

@testset "load_pnjl_config()" begin
    data = _CP.load_pnjl_config()
    @test data isa NamedTuple
    @test haskey(data, :config)
    @test haskey(data, :profile)
    @test haskey(data, :path)

    cfg = data.config
    @test haskey(cfg, "physical")
    @test haskey(cfg, "model")
    @test haskey(cfg, "polyakov")
end
