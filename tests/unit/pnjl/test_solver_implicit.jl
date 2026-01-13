# PNJL Solver ImplicitSolver 单元测试
#
# 测试内容：
# 1. SolverResult 结构
# 2. solve 各模式
# 3. 隐函数求解器
# 4. 导数计算

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

# 避免跨测试文件导出冲突：不把符号 `using` 进 Main。
P = PNJL

const ħc = 197.327  # MeV·fm

# ============================================================================
# solve FixedMu 测试
# ============================================================================

@testset "solve FixedMu" begin
    T_MeV = 100.0
    μ_MeV = 200.0
    T_fm = T_MeV / ħc
    μ_fm = μ_MeV / ħc
    
    @testset "基本求解" begin
        result = P.solve(P.FixedMu(), T_fm, μ_fm; p_num=24, t_num=6)
        
        @test result isa P.SolverResult
        @test result.mode isa P.FixedMu
        @test result.converged
        @test length(result.solution) == 5
        @test length(result.x_state) == 5
        @test length(result.mu_vec) == 3
        @test length(result.masses) == 3
    end
    
    @testset "热力学量有限" begin
        result = P.solve(P.FixedMu(), T_fm, μ_fm; p_num=24, t_num=6)
        
        @test isfinite(result.omega)
        @test isfinite(result.pressure)
        @test isfinite(result.rho_norm)
        @test isfinite(result.entropy)
        @test isfinite(result.energy)
        @test all(isfinite.(result.masses))
    end
    
    @testset "P = -Ω" begin
        result = P.solve(P.FixedMu(), T_fm, μ_fm; p_num=24, t_num=6)
        @test isapprox(result.pressure, -result.omega; rtol=1e-10)
    end
    
    @testset "各向异性参数" begin
        result_iso = P.solve(P.FixedMu(), T_fm, μ_fm; xi=0.0, p_num=24, t_num=6)
        result_aniso = P.solve(P.FixedMu(), T_fm, μ_fm; xi=0.2, p_num=24, t_num=6)
        
        @test result_iso.converged
        @test result_aniso.converged
        @test result_iso.xi == 0.0
        @test result_aniso.xi == 0.2
    end
    
    @testset "自定义初值策略" begin
        result = P.solve(P.FixedMu(), T_fm, μ_fm; 
                      seed_strategy=P.DefaultSeed(phase_hint=:hadron),
                      p_num=24, t_num=6)
        @test result.converged
    end
end

# ============================================================================
# solve FixedRho 测试
# ============================================================================

@testset "solve FixedRho" begin
    T_MeV = 100.0
    T_fm = T_MeV / ħc
    
    @testset "基本求解" begin
        result = P.solve(P.FixedRho(1.0), T_fm; p_num=24, t_num=6)
        
        @test result isa P.SolverResult
        @test result.mode isa P.FixedRho
        @test length(result.solution) == 8
        @test length(result.x_state) == 5
        @test length(result.mu_vec) == 3
    end
    
    @testset "密度约束满足" begin
        target_rho = 1.0
        result = P.solve(P.FixedRho(target_rho), T_fm; p_num=24, t_num=6)
        
        if result.converged
            @test isapprox(result.rho_norm, target_rho; rtol=0.01)
        end
    end
    
    @testset "化学势相等" begin
        result = P.solve(P.FixedRho(1.0), T_fm; p_num=24, t_num=6)
        
        if result.converged
            @test isapprox(result.mu_vec[1], result.mu_vec[2]; rtol=1e-6)
            @test isapprox(result.mu_vec[2], result.mu_vec[3]; rtol=1e-6)
        end
    end
end

# ============================================================================
# solve FixedEntropy 测试
# ============================================================================

@testset "solve FixedEntropy" begin
    T_MeV = 150.0
    T_fm = T_MeV / ħc
    
    @testset "基本求解" begin
        result = P.solve(P.FixedEntropy(0.5), T_fm; p_num=24, t_num=6)
        
        @test result isa P.SolverResult
        @test result.mode isa P.FixedEntropy
        @test length(result.solution) == 8
    end
end

# ============================================================================
# solve FixedSigma 测试
# ============================================================================

@testset "solve FixedSigma" begin
    T_MeV = 150.0
    T_fm = T_MeV / ħc
    
    @testset "基本求解" begin
        result = P.solve(P.FixedSigma(10.0), T_fm; p_num=24, t_num=6)
        
        @test result isa P.SolverResult
        @test result.mode isa P.FixedSigma
        @test length(result.solution) == 8
    end
end

# ============================================================================
# 隐函数求解器测试
# ============================================================================

@testset "create_implicit_solver" begin
    solver = P.create_implicit_solver(xi=0.0, p_num=24, t_num=6)
    
    T_fm = 0.5
    μ_fm = 1.0
    θ = [T_fm, μ_fm]
    
    x, _ = solver(θ)
    
    @test length(x) == 5
    @test all(isfinite.(x))
end

@testset "solve_with_derivatives" begin
    T_fm = 0.5
    μ_fm = 1.0
    
    @testset "一阶导数" begin
        result = P.solve_with_derivatives(T_fm, μ_fm; order=1, p_num=24, t_num=6)
        
        @test haskey(result, :x)
        @test haskey(result, :dx_dT)
        @test haskey(result, :dx_dμ)
        @test length(result.x) == 5
        @test length(result.dx_dT) == 5
        @test length(result.dx_dμ) == 5
        @test all(isfinite.(result.x))
        @test all(isfinite.(result.dx_dT))
        @test all(isfinite.(result.dx_dμ))
    end
    
    @testset "二阶导数" begin
        result = P.solve_with_derivatives(T_fm, μ_fm; order=2, p_num=24, t_num=6)
        
        @test haskey(result, :d2x_dT2)
        @test haskey(result, :d2x_dμ2)
        @test haskey(result, :d2x_dTdμ)
        @test all(isfinite.(result.d2x_dT2))
        @test all(isfinite.(result.d2x_dμ2))
        @test all(isfinite.(result.d2x_dTdμ))
    end
end

# ============================================================================
# 物理一致性测试
# ============================================================================

@testset "physical consistency" begin
    @testset "熵非负" begin
        T_fm = 0.5
        μ_fm = 1.0
        result = P.solve(P.FixedMu(), T_fm, μ_fm; p_num=24, t_num=6)
        
        if result.converged
            @test result.entropy >= 0
        end
    end
    
    @testset "质量正定" begin
        T_fm = 0.5
        μ_fm = 1.0
        result = P.solve(P.FixedMu(), T_fm, μ_fm; p_num=24, t_num=6)
        
        if result.converged
            @test all(result.masses .> 0)
        end
    end
    
    @testset "s 质量大于 u/d" begin
        T_fm = 0.5
        μ_fm = 1.0
        result = P.solve(P.FixedMu(), T_fm, μ_fm; p_num=24, t_num=6)
        
        if result.converged
            @test result.masses[3] > result.masses[1]  # M_s > M_u
            @test result.masses[3] > result.masses[2]  # M_s > M_d
        end
    end
    
    @testset "u/d 质量相等（同位旋对称）" begin
        T_fm = 0.5
        μ_fm = 1.0
        result = P.solve(P.FixedMu(), T_fm, μ_fm; p_num=24, t_num=6)
        
        if result.converged
            @test isapprox(result.masses[1], result.masses[2]; rtol=1e-6)
        end
    end
end

# ============================================================================
# 温度/化学势扫描测试
# ============================================================================

@testset "parameter scan" begin
    @testset "温度扫描" begin
        μ_fm = 1.0
        T_values = [0.3, 0.5, 0.7]  # ~60, 100, 140 MeV
        
        for T_fm in T_values
            result = P.solve(P.FixedMu(), T_fm, μ_fm; p_num=24, t_num=6)
            @test result.converged
            @test isfinite(result.pressure)
        end
    end
    
    @testset "化学势扫描" begin
        T_fm = 0.5
        μ_values = [0.5, 1.0, 1.5]  # ~100, 200, 300 MeV
        
        for μ_fm in μ_values
            result = P.solve(P.FixedMu(), T_fm, μ_fm; p_num=24, t_num=6)
            @test result.converged
            @test isfinite(result.pressure)
        end
    end
end
