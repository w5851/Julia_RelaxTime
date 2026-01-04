# PNJL Core Thermodynamics 单元测试
#
# 测试内容：
# 1. 有效质量计算
# 2. 手征凝聚项
# 3. Polyakov loop 有效势
# 4. 热力学量计算

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "core", "Thermodynamics.jl"))

using .Thermodynamics
using .Thermodynamics.Integrals: cached_nodes

# ============================================================================
# 有效质量计算测试
# ============================================================================

@testset "calculate_mass_vec" begin
    @testset "从 φ 计算质量" begin
        φ = SVector{3}(-1.8, -1.8, -2.2)  # 典型强子相值
        masses = calculate_mass_vec(φ)
        
        @test length(masses) == 3
        @test all(isfinite.(masses))
        @test masses[1] ≈ masses[2]  # u, d 质量相等（同位旋对称）
        @test masses[3] > masses[1]  # s 质量更大
    end
    
    @testset "从 x_state 计算质量" begin
        x_state = SVector{5}(-1.8, -1.8, -2.2, 0.1, 0.1)
        masses = calculate_mass_vec(x_state)
        
        @test length(masses) == 3
        @test all(isfinite.(masses))
    end
    
    @testset "质量正定性" begin
        # 强子相：大质量
        φ_hadron = SVector{3}(-1.8, -1.8, -2.2)
        masses_hadron = calculate_mass_vec(φ_hadron)
        @test all(masses_hadron .> 0)
        
        # 夸克相：小质量
        φ_quark = SVector{3}(-0.2, -0.2, -1.8)
        masses_quark = calculate_mass_vec(φ_quark)
        @test all(masses_quark .> 0)
    end
end

# ============================================================================
# 势能项测试
# ============================================================================

@testset "calculate_chiral" begin
    φ = SVector{3}(-1.8, -1.8, -2.2)
    chi = calculate_chiral(φ)
    @test isfinite(chi)
end

@testset "calculate_U" begin
    T_fm = 0.5  # ~100 MeV
    
    @testset "禁闭相 (Φ ≈ 0)" begin
        U = calculate_U(T_fm, 0.01, 0.01)
        @test isfinite(U)
    end
    
    @testset "解禁闭相 (Φ ≈ 0.6)" begin
        U = calculate_U(T_fm, 0.6, 0.6)
        @test isfinite(U)
    end
    
    @testset "Φ = Φ̄ 对称性" begin
        U1 = calculate_U(T_fm, 0.3, 0.4)
        U2 = calculate_U(T_fm, 0.4, 0.3)
        @test isapprox(U1, U2; rtol=1e-10)
    end
    
    @testset "温度依赖性" begin
        # 高温时 Polyakov loop 势更深
        U_low = calculate_U(0.3, 0.5, 0.5)
        U_high = calculate_U(1.0, 0.5, 0.5)
        @test isfinite(U_low)
        @test isfinite(U_high)
    end
end

@testset "calculate_U_derivative_T" begin
    T_fm = 0.5
    Φ = 0.3
    Φ̄ = 0.3
    
    dU_dT = calculate_U_derivative_T(T_fm, Φ, Φ̄)
    @test isfinite(dU_dT)
    
    # 数值验证
    δT = 1e-5
    U_plus = calculate_U(T_fm + δT, Φ, Φ̄)
    U_minus = calculate_U(T_fm - δT, Φ, Φ̄)
    dU_dT_fd = (U_plus - U_minus) / (2δT)
    @test isapprox(dU_dT, dU_dT_fd; rtol=1e-4)
end

# ============================================================================
# 热力学势与压强测试
# ============================================================================

@testset "calculate_omega" begin
    x_state = SVector{5}(-1.5, -1.5, -2.1, 0.2, 0.2)
    mu_vec = SVector{3}(1.0, 1.0, 1.0)
    T_fm = 0.5
    thermal_nodes = cached_nodes(24, 6)
    xi = 0.0
    
    omega = calculate_omega(x_state, mu_vec, T_fm, thermal_nodes, xi)
    @test isfinite(omega)
end

@testset "calculate_pressure" begin
    x_state = SVector{5}(-1.5, -1.5, -2.1, 0.2, 0.2)
    mu_vec = SVector{3}(1.0, 1.0, 1.0)
    T_fm = 0.5
    thermal_nodes = cached_nodes(24, 6)
    xi = 0.0
    
    P = calculate_pressure(x_state, mu_vec, T_fm, thermal_nodes, xi)
    omega = calculate_omega(x_state, mu_vec, T_fm, thermal_nodes, xi)
    
    @test isfinite(P)
    @test isapprox(P, -omega; rtol=1e-10)
end

# ============================================================================
# 热力学量计算测试
# ============================================================================

@testset "calculate_rho" begin
    x_state = SVector{5}(-1.5, -1.5, -2.1, 0.2, 0.2)
    mu_vec = SVector{3}(1.5, 1.5, 1.5)
    T_fm = 0.5
    thermal_nodes = cached_nodes(24, 6)
    xi = 0.0
    
    rho_vec = calculate_rho(x_state, mu_vec, T_fm, thermal_nodes, xi)
    
    @test length(rho_vec) == 3
    @test all(isfinite.(rho_vec))
end

@testset "calculate_thermo" begin
    x_state = SVector{5}(-1.5, -1.5, -2.1, 0.2, 0.2)
    mu_vec = SVector{3}(1.5, 1.5, 1.5)
    T_fm = 0.5
    thermal_nodes = cached_nodes(24, 6)
    xi = 0.0
    
    P, rho_norm, s, ε = calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
    
    @test isfinite(P)
    @test isfinite(rho_norm)
    @test isfinite(s)
    @test isfinite(ε)
    @test s >= 0  # 熵非负
end
