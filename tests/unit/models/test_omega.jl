# omega.jl 单元测试
#
# 测试内容：
# 1. omega_components 返回完整 NamedTuple
# 2. omega ≡ grand_potential 别名
# 3. NJL 3维 / PNJL 5维 状态

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================
# omega_components
# ============================================================================

@testset "omega" begin

    @testset "omega_components NJL 返回 NamedTuple" begin
        m = Models.create_model(:NJL)
        φ = SVector{3}(0.01, 0.02, 0.5)
        T = 0.8
        μ = SVector{3}(0.0, 0.0, 0.0)
        x = SVector{5}(φ[1], φ[2], φ[3], 1.0, 1.0)

        oc = Models.omega_components(m, x, T, μ; p_num=24, t_num=6, xi=0.0)
        @test oc isa NamedTuple
        @test haskey(oc, :chi)
        @test haskey(oc, :poly)
        @test haskey(oc, :vac)
        @test haskey(oc, :therm)
        @test haskey(oc, :masses)
        @test haskey(oc, :omega)
        @test isfinite(oc.omega)
    end

    @testset "omega 标量 == omega_components.omega" begin
        m = Models.create_model(:NJL)
        φ = SVector{3}(0.01, 0.02, 0.5)
        T = 0.8
        μ = SVector{3}(0.0, 0.0, 0.0)
        x = SVector{5}(φ[1], φ[2], φ[3], 1.0, 1.0)

        ω = Models.omega(m, x, T, μ; p_num=24, t_num=6, xi=0.0)
        oc = Models.omega_components(m, x, T, μ; p_num=24, t_num=6, xi=0.0)
        @test ω ≈ oc.omega rtol=1e-14
    end

    @testset "grand_potential 别名" begin
        m = Models.create_model(:NJL)
        φ = SVector{3}(0.01, 0.02, 0.5)
        T = 0.8
        μ = SVector{3}(0.0, 0.0, 0.0)
        x = SVector{5}(φ[1], φ[2], φ[3], 1.0, 1.0)

        @test Models.grand_potential(m, x, T, μ; p_num=24, t_num=6, xi=0.0) ≈
              Models.omega(m, x, T, μ; p_num=24, t_num=6, xi=0.0) rtol=1e-14
    end

    @testset "PNJL 5维状态" begin
        Models.pnjl_module()
        m = Models.create_model(:PNJL)
        x = SVector{5}(-1.84, -1.84, -2.23, 0.5, 0.5)
        T = 0.8
        μ = SVector{3}(0.0, 0.0, 0.0)

        ω = Models.omega(m, x, T, μ; p_num=24, t_num=6, xi=0.0)
        @test isfinite(ω)
    end

    @testset "omega 有限性 (多温度)" begin
        m = Models.create_model(:NJL)
        x = SVector{5}(0.01, 0.02, 0.5, 1.0, 1.0)
        μ = SVector{3}(0.0, 0.0, 0.0)
        for T in [0.1, 0.5, 1.0, 2.0]
            ω = Models.omega(m, x, T, μ; p_num=24, t_num=6, xi=0.0)
            @test isfinite(ω)
        end
    end
end
