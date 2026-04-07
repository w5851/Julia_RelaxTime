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

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const P = Models.pnjl_module()

# 避免跨测试文件的导出冲突：不把符号 `using` 进 Main，统一通过模块前缀访问。

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

    # 清理约束：重复候选应复用同一常量定义
    @test P.SeedStrategies.HT_GUESS_0p9_SEED_5 === P.SeedStrategies.VERY_HIGH_TEMP_SEED_5
end

@testset "PhaseAwareContinuitySeed cleanup" begin
    @test !isdefined(P, :PhaseAwareContinuitySeed)
    @test !isdefined(P, :set_phase!)
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
