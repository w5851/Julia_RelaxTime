# RPNJLModel.jl 单元测试
#
# 测试内容：
# 1. RPNJLModel 构造
# 2. create_model(:RPNJL) 工厂
# 3. calculate_mass_vec 与 base 比较
# 4. omega / gap_residual 基本可调用

using Test
using StaticArrays

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
