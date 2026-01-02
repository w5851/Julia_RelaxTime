# PNJL Core Integrals 单元测试
#
# 测试内容：
# 1. 积分节点缓存
# 2. 真空项积分
# 3. 热项积分
# 4. 能量色散关系

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "core", "Integrals.jl"))

using .Integrals

# ============================================================================
# 积分节点缓存测试
# ============================================================================

@testset "cached_nodes" begin
    @testset "默认节点数" begin
        p_mesh, cosθ_mesh, coefficients = cached_nodes(DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT)
        @test size(p_mesh) == (DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT)
        @test size(cosθ_mesh) == (DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT)
        @test size(coefficients) == (DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT)
    end
    
    @testset "自定义节点数" begin
        p_num, t_num = 24, 8
        p_mesh, cosθ_mesh, coefficients = cached_nodes(p_num, t_num)
        @test size(p_mesh) == (p_num, t_num)
        @test size(cosθ_mesh) == (p_num, t_num)
        @test size(coefficients) == (p_num, t_num)
    end
    
    @testset "缓存一致性" begin
        nodes1 = cached_nodes(48, 12)
        nodes2 = cached_nodes(48, 12)
        @test nodes1[1] === nodes2[1]  # 同一对象
        @test nodes1[2] === nodes2[2]
        @test nodes1[3] === nodes2[3]
    end
    
    @testset "节点范围" begin
        p_mesh, cosθ_mesh, _ = cached_nodes(48, 12)
        @test all(p_mesh .>= 0)  # 动量非负
        @test all(p_mesh .<= 10)  # 动量上限
        @test all(cosθ_mesh .>= 0)  # cosθ ∈ [0, 1]
        @test all(cosθ_mesh .<= 1)
    end
end

# ============================================================================
# 真空项积分测试
# ============================================================================

@testset "vacuum_integral" begin
    @testset "正质量" begin
        mass = 1.5  # fm⁻¹, ~300 MeV
        result = vacuum_integral(mass)
        @test isfinite(result)
        @test result > 0  # 真空项应为正
    end
    
    @testset "小质量" begin
        mass = 0.01  # fm⁻¹, ~2 MeV
        result = vacuum_integral(mass)
        @test isfinite(result)
        @test result > 0
    end
    
    @testset "质量单调性" begin
        # 质量越大，真空项越大
        m1 = vacuum_integral(0.5)
        m2 = vacuum_integral(1.0)
        m3 = vacuum_integral(1.5)
        @test m1 < m2 < m3
    end
    
    @testset "负质量处理" begin
        # 应该使用绝对值
        result_pos = vacuum_integral(1.0)
        result_neg = vacuum_integral(-1.0)
        @test isapprox(result_pos, result_neg; rtol=1e-10)
    end
end

@testset "calculate_energy_sum" begin
    masses = SVector{3}(1.5, 1.5, 2.5)  # u, d, s 质量
    result = calculate_energy_sum(masses)
    @test isfinite(result)
    @test result < 0  # -2Nc * sum(I) 应为负
end

# ============================================================================
# 能量色散关系测试
# ============================================================================

@testset "calculate_energy_isotropic" begin
    mass = 1.5
    p = 1.0
    E = calculate_energy_isotropic(mass, p)
    @test E ≈ sqrt(p^2 + mass^2)
    @test E > mass  # E > m
    @test E > p     # E > p
end

@testset "calculate_energy_anisotropic" begin
    mass = 1.5
    p = 1.0
    
    @testset "xi=0 退化为各向同性" begin
        E_iso = calculate_energy_isotropic(mass, p)
        E_aniso = calculate_energy_anisotropic(mass, p, 0.0, 0.5)
        @test isapprox(E_iso, E_aniso; rtol=1e-10)
    end
    
    @testset "xi>0 增加能量" begin
        cosθ = 0.5
        E_iso = calculate_energy_isotropic(mass, p)
        E_aniso = calculate_energy_anisotropic(mass, p, 0.5, cosθ)
        @test E_aniso > E_iso
    end
    
    @testset "cosθ=0 时无各向异性修正" begin
        E_iso = calculate_energy_isotropic(mass, p)
        E_aniso = calculate_energy_anisotropic(mass, p, 0.5, 0.0)
        @test isapprox(E_iso, E_aniso; rtol=1e-10)
    end
end

# ============================================================================
# 热项积分测试
# ============================================================================

@testset "calculate_log_sum" begin
    masses = SVector{3}(1.5, 1.5, 2.5)
    p_nodes, cosθ_nodes, coefficients = cached_nodes(24, 6)
    Φ = 0.3
    Φ̄ = 0.3
    mu_vec = SVector{3}(1.0, 1.0, 1.0)
    T_fm = 0.5
    xi = 0.0
    
    result = calculate_log_sum(masses, p_nodes, cosθ_nodes, coefficients, Φ, Φ̄, mu_vec, T_fm, xi)
    @test isfinite(result)
    @test result < 0  # 热项贡献通常为负
end

@testset "calculate_log_sum_derivatives" begin
    masses = SVector{3}(1.5, 1.5, 2.5)
    p_nodes, cosθ_nodes, coefficients = cached_nodes(24, 6)
    Φ = 0.3
    Φ̄ = 0.3
    mu_vec = SVector{3}(1.0, 1.0, 1.0)
    T_fm = 0.5
    xi = 0.0
    
    log_sum, d_log_sum_dmu, d_log_sum_dT = calculate_log_sum_derivatives(
        masses, p_nodes, cosθ_nodes, coefficients, Φ, Φ̄, mu_vec, T_fm, xi)
    
    @testset "返回值有限" begin
        @test isfinite(log_sum)
        @test all(isfinite.(d_log_sum_dmu))
        @test isfinite(d_log_sum_dT)
    end
    
    @testset "导数维度" begin
        @test length(d_log_sum_dmu) == 3
    end
    
    @testset "与 calculate_log_sum 一致" begin
        log_sum_direct = calculate_log_sum(masses, p_nodes, cosθ_nodes, coefficients, Φ, Φ̄, mu_vec, T_fm, xi)
        @test isapprox(log_sum, log_sum_direct; rtol=1e-10)
    end
end
