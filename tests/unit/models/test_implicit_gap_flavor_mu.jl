# flavor-mu implicit_gap 单元测试
#
# 测试内容：
# 1. flavor 化学势版本的隐函数求解器可构造与 forward solve
# 2. flavor 化学势版本的一阶导数接口
# 3. 对称路径退化到旧标量 μ 接口时的一致性

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

@testset "implicit_gap flavor mu" begin
    @testset "legacy flavor factory is qualified compat-only" begin
        @test !(:create_flavor_mu_implicit_gap_solver in names(Models))
    end

    @testset "compat create_flavor_mu_implicit_gap_solver" begin
        m = Models.create_model(:PNJL)
        igf = Models.create_flavor_mu_implicit_gap_solver(m; p_num=24, t_num=6)
        @test igf isa Any
    end

    @testset "forward solve flavor mu symmetric point" begin
        m = Models.create_model(:PNJL)
        igf = Models.create_flavor_mu_implicit_gap_solver(m; p_num=24, t_num=6)
        θ = [0.5, 0.0, 0.0, 0.0]
        result = igf(θ)
        x = result isa Tuple ? result[1] : result
        @test length(x) == 5
        @test all(isfinite.(x))
    end

    @testset "solve_pnjl_with_flavor_mu_derivatives symmetric consistency" begin
        T_fm = 0.5
        μ_fm = 0.2
        μ_vec = SVector(μ_fm, μ_fm, μ_fm)

        old_result = Models.solve_pnjl_with_derivatives(T_fm, μ_fm; order=1, p_num=24, t_num=6)
        new_result = Models.solve_pnjl_with_flavor_mu_derivatives(T_fm, μ_vec; order=1, p_num=24, t_num=6)
        auto_result = Models.solve_pnjl_with_flavor_mu_derivatives(T_fm, μ_vec; order=1, p_num=8, t_num=4, derivative_backend=:auto)
        td_result = Models.solve_pnjl_with_flavor_mu_derivatives(T_fm, μ_vec; order=1, p_num=8, t_num=4, derivative_backend=:taylordiff)

        @test length(new_result.x) == 5
        @test length(new_result.dx_dT) == 5
        @test size(new_result.dx_dmu_vec) == (5, 3)
        @test all(isfinite.(new_result.x))
        @test all(isfinite.(new_result.dx_dT))
        @test all(isfinite.(new_result.dx_dmu_vec))

        @test all(isapprox.(new_result.x, old_result.x; rtol=1e-7, atol=1e-9))
        @test all(isapprox.(new_result.dx_dT, old_result.dx_dT; rtol=1e-6, atol=1e-8))

        symmetric_direction = vec(sum(new_result.dx_dmu_vec; dims=2))
        @test all(isapprox.(symmetric_direction, old_result.dx_dμ; rtol=1e-6, atol=1e-8))
        @test all(isapprox.(td_result.x, auto_result.x; rtol=1e-12, atol=1e-12))
        @test all(isapprox.(td_result.dx_dT, auto_result.dx_dT; rtol=1e-12, atol=1e-12))
        @test all(isapprox.(td_result.dx_dmu_vec, auto_result.dx_dmu_vec; rtol=1e-12, atol=1e-12))
        @test_throws ArgumentError Models.solve_pnjl_with_flavor_mu_derivatives(T_fm, μ_vec; order=1, p_num=8, t_num=4, derivative_backend=:forwarddiff)
    end

    @testset "solve_pnjl_with_flavor_mu_derivatives asymmetric point" begin
        T_fm = 0.5
        μ_vec = SVector(0.18, 0.12, 0.06)

        result = Models.solve_pnjl_with_flavor_mu_derivatives(T_fm, μ_vec; order=1, p_num=24, t_num=6)

        @test length(result.x) == 5
        @test result.mu_vec == μ_vec
        @test length(result.dx_dT) == 5
        @test size(result.dx_dmu_vec) == (5, 3)
        @test all(isfinite.(result.x))
        @test all(isfinite.(result.dx_dT))
        @test all(isfinite.(result.dx_dmu_vec))
    end
end
