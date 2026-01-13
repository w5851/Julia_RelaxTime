# PNJL Solver SeedStrategies 单元测试
#
# 测试内容：
# 1. 默认种子策略
# 2. 多初值策略
# 3. 连续性跟踪策略
# 4. 相变感知策略
# 5. 相变感知连续性跟踪策略

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

# 避免跨测试文件的导出冲突：不把符号 `using` 进 Main，统一通过模块前缀访问。
P = PNJL

# ============================================================================
# 内置种子值测试
# ============================================================================

@testset "builtin seeds" begin
    @test length(P.HADRON_SEED_5) == 5
    @test length(P.QUARK_SEED_5) == 5
    @test length(P.HADRON_SEED_8) == 8
    @test length(P.QUARK_SEED_8) == 8
    @test length(P.MEDIUM_SEED_5) == 5
    @test length(P.HIGH_DENSITY_SEED_5) == 5
    @test length(P.HIGH_TEMP_SEED_5) == 5
    
    # 强子相：Polyakov loop 接近零
    @test abs(P.HADRON_SEED_5[4]) < 0.1
    @test abs(P.HADRON_SEED_5[5]) < 0.1
    
    # 夸克相：Polyakov loop 较大
    @test P.QUARK_SEED_5[4] > 0.3
    @test P.QUARK_SEED_5[5] > 0.3
end

# ============================================================================
# DefaultSeed 测试
# ============================================================================

@testset "DefaultSeed" begin
    @testset "构造" begin
        s = P.DefaultSeed()
        @test s isa P.SeedStrategy
        @test s.phase_hint == :auto
        
        s_hadron = P.DefaultSeed(phase_hint=:hadron)
        @test s_hadron.phase_hint == :hadron
        
        s_quark = P.DefaultSeed(phase_hint=:quark)
        @test s_quark.phase_hint == :quark
    end
    
    @testset "get_seed FixedMu" begin
        s = P.DefaultSeed(phase_hint=:hadron)
        θ = [0.5, 1.0]  # T, μ in fm⁻¹
        seed = P.get_seed(s, θ, P.FixedMu())
        @test length(seed) == 5
        @test all(isfinite.(seed))
    end
    
    @testset "get_seed FixedRho" begin
        s = P.DefaultSeed(phase_hint=:hadron)
        θ = [0.5]  # T in fm⁻¹
        seed = P.get_seed(s, θ, P.FixedRho(1.0))
        @test length(seed) == 8
        @test all(isfinite.(seed))
    end
    
    @testset "auto phase hint" begin
        s = P.DefaultSeed(phase_hint=:auto)
        
        # 低温低密度 -> 强子相
        θ_hadron = [0.3, 0.5]  # ~60 MeV, ~100 MeV
        seed_hadron = P.get_seed(s, θ_hadron, P.FixedMu())
        
        # 高温 -> 夸克相
        θ_quark = [1.0, 0.5]  # ~200 MeV, ~100 MeV
        seed_quark = P.get_seed(s, θ_quark, P.FixedMu())
        
        @test length(seed_hadron) == 5
        @test length(seed_quark) == 5
    end
end

# ============================================================================
# MultiSeed 测试
# ============================================================================

@testset "MultiSeed" begin
    @testset "构造" begin
        s = P.MultiSeed()
        @test s isa P.SeedStrategy
        # 候选初值数量可能随实现演进而调整（例如增加更多相/密度/温度候选）。
        @test length(s.candidates) >= 2
    end
    
    @testset "get_seed" begin
        s = P.MultiSeed()
        θ = [0.5, 1.0]
        seed = P.get_seed(s, θ, P.FixedMu())
        @test length(seed) == 5
    end
    
    @testset "get_all_seeds" begin
        s = P.MultiSeed()
        θ = [0.5, 1.0]
        seeds = P.get_all_seeds(s, θ, P.FixedMu())
        @test length(seeds) == length(s.candidates)
        @test all(length.(seeds) .== 5)
    end
end

# ============================================================================
# ContinuitySeed 测试
# ============================================================================

@testset "ContinuitySeed" begin
    @testset "构造" begin
        s = P.ContinuitySeed()
        @test s isa P.SeedStrategy
        @test s.previous_solution === nothing
    end
    
    @testset "初始 get_seed（无前解）" begin
        s = P.ContinuitySeed()
        θ = [0.5, 1.0]
        seed = P.get_seed(s, θ, P.FixedMu())
        @test length(seed) == 5
    end
    
    @testset "update! 和 get_seed" begin
        s = P.ContinuitySeed()
        θ = [0.5, 1.0]
        
        # 模拟一个解
        solution = [-1.5, -1.5, -2.1, 0.2, 0.2]
        P.update!(s, solution)
        
        @test s.previous_solution !== nothing
        @test length(s.previous_solution) == 5
        
        # 下一次 get_seed 应该返回前一个解
        seed = P.get_seed(s, θ, P.FixedMu())
        @test seed ≈ solution
    end
    
    @testset "reset!" begin
        s = P.ContinuitySeed()
        P.update!(s, [-1.5, -1.5, -2.1, 0.2, 0.2])
        @test s.previous_solution !== nothing
        
        P.reset!(s)
        @test s.previous_solution === nothing
    end
    
    @testset "维度扩展" begin
        s = P.ContinuitySeed()
        solution_5 = [-1.5, -1.5, -2.1, 0.2, 0.2]
        P.update!(s, solution_5)
        
        # 请求 8 维种子
        seed = P.get_seed(s, [0.5], P.FixedRho(1.0))
        @test length(seed) == 8
        @test seed[1:5] ≈ solution_5
    end
end

# ============================================================================
# PhaseAwareSeed 测试
# ============================================================================

@testset "PhaseAwareSeed" begin
    @testset "无数据构造" begin
        s = P.PhaseAwareSeed()
        @test s isa P.SeedStrategy
        @test s.boundary_data === nothing
    end
    
    @testset "从 xi 构造" begin
        # 尝试加载 xi=0.0 的数据
        s = P.PhaseAwareSeed(0.0)
        @test s isa P.SeedStrategy
        # 如果数据文件存在，boundary_data 不为空
    end
    
    @testset "get_seed 无数据" begin
        s = P.PhaseAwareSeed()
        θ = [0.5, 1.0]
        seed = P.get_seed(s, θ, P.FixedMu())
        @test length(seed) == 5
    end
end

# ============================================================================
# PhaseAwareContinuitySeed 测试
# ============================================================================

@testset "PhaseAwareContinuitySeed" begin
    @testset "构造" begin
        s = P.PhaseAwareContinuitySeed()
        @test s isa P.SeedStrategy
        @test s.previous_solution === nothing
        @test s.previous_phase == :unknown
    end
    
    @testset "从 xi 构造" begin
        s = P.PhaseAwareContinuitySeed(0.0)
        @test s isa P.SeedStrategy
    end
    
    @testset "get_seed 初始" begin
        s = P.PhaseAwareContinuitySeed()
        θ = [0.5, 1.0]
        seed = P.get_seed(s, θ, P.FixedMu())
        @test length(seed) == 5
    end
    
    @testset "update! 带相位" begin
        s = P.PhaseAwareContinuitySeed()
        solution = [-1.5, -1.5, -2.1, 0.2, 0.2]
        T_MeV = 100.0
        μ_MeV = 200.0
        
        P.update!(s, solution, T_MeV, μ_MeV)
        
        @test s.previous_solution !== nothing
        @test s.previous_solution ≈ solution
    end
    
    @testset "reset!" begin
        s = P.PhaseAwareContinuitySeed()
        P.update!(s, [-1.5, -1.5, -2.1, 0.2, 0.2], 100.0, 200.0)
        
        P.reset!(s)
        @test s.previous_solution === nothing
        @test s.previous_phase == :unknown
    end
    
    @testset "set_phase!" begin
        s = P.PhaseAwareContinuitySeed()
        P.set_phase!(s, :hadron)
        @test s.previous_phase == :hadron
        
        P.set_phase!(s, :quark)
        @test s.previous_phase == :quark
    end
end

# ============================================================================
# PhaseBoundaryData 测试
# ============================================================================

@testset "PhaseBoundaryData" begin
    @testset "load_phase_boundary" begin
        # 尝试加载数据
        data = P.load_phase_boundary(0.0)
        @test data isa P.PhaseBoundaryData
        @test data.xi == 0.0
        
        # 如果数据存在，检查基本属性
        if !isempty(data.T_values)
            @test length(data.T_values) == length(data.mu_values)
            @test issorted(data.T_values)
        end
    end
    
    @testset "interpolate_mu_c" begin
        data = P.load_phase_boundary(0.0)
        
        if !isempty(data.T_values)
            # 在数据范围内插值
            T_mid = (data.T_values[1] + data.T_values[end]) / 2
            μ_c = P.interpolate_mu_c(data, T_mid)
            @test isfinite(μ_c)
            
            # 超过 CEP 返回 NaN
            if !isnan(data.T_CEP)
                μ_c_above = P.interpolate_mu_c(data, data.T_CEP + 10)
                @test isnan(μ_c_above)
            end
        end
    end
end

# ============================================================================
# extend_seed 测试
# ============================================================================

@testset "extend_seed" begin
    base_seed = [-1.5, -1.5, -2.1, 0.2, 0.2]
    
    @testset "FixedMu" begin
        seed = P.extend_seed(base_seed, P.FixedMu())
        @test length(seed) == 5
        @test seed ≈ base_seed
    end
    
    @testset "FixedRho" begin
        seed = P.extend_seed(base_seed, P.FixedRho(1.0))
        @test length(seed) == 8
        @test seed[1:5] ≈ base_seed
        @test all(isfinite.(seed[6:8]))
    end
    
    @testset "FixedEntropy" begin
        seed = P.extend_seed(base_seed, P.FixedEntropy(0.5))
        @test length(seed) == 8
        @test seed[1:5] ≈ base_seed
    end
    
    @testset "FixedSigma" begin
        seed = P.extend_seed(base_seed, P.FixedSigma(10.0))
        @test length(seed) == 8
        @test seed[1:5] ≈ base_seed
    end
end
