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

using .PNJL.ConstraintModes
using .PNJL.SeedStrategies

# ============================================================================
# 内置种子值测试
# ============================================================================

@testset "builtin seeds" begin
    @test length(HADRON_SEED_5) == 5
    @test length(QUARK_SEED_5) == 5
    @test length(HADRON_SEED_8) == 8
    @test length(QUARK_SEED_8) == 8
    @test length(MEDIUM_SEED_5) == 5
    @test length(HIGH_DENSITY_SEED_5) == 5
    @test length(HIGH_TEMP_SEED_5) == 5
    
    # 强子相：Polyakov loop 接近零
    @test abs(HADRON_SEED_5[4]) < 0.1
    @test abs(HADRON_SEED_5[5]) < 0.1
    
    # 夸克相：Polyakov loop 较大
    @test QUARK_SEED_5[4] > 0.3
    @test QUARK_SEED_5[5] > 0.3
end

# ============================================================================
# DefaultSeed 测试
# ============================================================================

@testset "DefaultSeed" begin
    @testset "构造" begin
        s = DefaultSeed()
        @test s isa SeedStrategy
        @test s.phase_hint == :auto
        
        s_hadron = DefaultSeed(phase_hint=:hadron)
        @test s_hadron.phase_hint == :hadron
        
        s_quark = DefaultSeed(phase_hint=:quark)
        @test s_quark.phase_hint == :quark
    end
    
    @testset "get_seed FixedMu" begin
        s = DefaultSeed(phase_hint=:hadron)
        θ = [0.5, 1.0]  # T, μ in fm⁻¹
        seed = get_seed(s, θ, FixedMu())
        @test length(seed) == 5
        @test all(isfinite.(seed))
    end
    
    @testset "get_seed FixedRho" begin
        s = DefaultSeed(phase_hint=:hadron)
        θ = [0.5]  # T in fm⁻¹
        seed = get_seed(s, θ, FixedRho(1.0))
        @test length(seed) == 8
        @test all(isfinite.(seed))
    end
    
    @testset "auto phase hint" begin
        s = DefaultSeed(phase_hint=:auto)
        
        # 低温低密度 -> 强子相
        θ_hadron = [0.3, 0.5]  # ~60 MeV, ~100 MeV
        seed_hadron = get_seed(s, θ_hadron, FixedMu())
        
        # 高温 -> 夸克相
        θ_quark = [1.0, 0.5]  # ~200 MeV, ~100 MeV
        seed_quark = get_seed(s, θ_quark, FixedMu())
        
        @test length(seed_hadron) == 5
        @test length(seed_quark) == 5
    end
end

# ============================================================================
# MultiSeed 测试
# ============================================================================

@testset "MultiSeed" begin
    @testset "构造" begin
        s = MultiSeed()
        @test s isa SeedStrategy
        @test length(s.candidates) == 2
    end
    
    @testset "get_seed" begin
        s = MultiSeed()
        θ = [0.5, 1.0]
        seed = get_seed(s, θ, FixedMu())
        @test length(seed) == 5
    end
    
    @testset "get_all_seeds" begin
        s = MultiSeed()
        θ = [0.5, 1.0]
        seeds = get_all_seeds(s, θ, FixedMu())
        @test length(seeds) == 2
        @test all(length.(seeds) .== 5)
    end
end

# ============================================================================
# ContinuitySeed 测试
# ============================================================================

@testset "ContinuitySeed" begin
    @testset "构造" begin
        s = ContinuitySeed()
        @test s isa SeedStrategy
        @test s.previous_solution === nothing
    end
    
    @testset "初始 get_seed（无前解）" begin
        s = ContinuitySeed()
        θ = [0.5, 1.0]
        seed = get_seed(s, θ, FixedMu())
        @test length(seed) == 5
    end
    
    @testset "update! 和 get_seed" begin
        s = ContinuitySeed()
        θ = [0.5, 1.0]
        
        # 模拟一个解
        solution = [-1.5, -1.5, -2.1, 0.2, 0.2]
        update!(s, solution)
        
        @test s.previous_solution !== nothing
        @test length(s.previous_solution) == 5
        
        # 下一次 get_seed 应该返回前一个解
        seed = get_seed(s, θ, FixedMu())
        @test seed ≈ solution
    end
    
    @testset "reset!" begin
        s = ContinuitySeed()
        update!(s, [-1.5, -1.5, -2.1, 0.2, 0.2])
        @test s.previous_solution !== nothing
        
        reset!(s)
        @test s.previous_solution === nothing
    end
    
    @testset "维度扩展" begin
        s = ContinuitySeed()
        solution_5 = [-1.5, -1.5, -2.1, 0.2, 0.2]
        update!(s, solution_5)
        
        # 请求 8 维种子
        seed = get_seed(s, [0.5], FixedRho(1.0))
        @test length(seed) == 8
        @test seed[1:5] ≈ solution_5
    end
end

# ============================================================================
# PhaseAwareSeed 测试
# ============================================================================

@testset "PhaseAwareSeed" begin
    @testset "无数据构造" begin
        s = PhaseAwareSeed()
        @test s isa SeedStrategy
        @test s.boundary_data === nothing
    end
    
    @testset "从 xi 构造" begin
        # 尝试加载 xi=0.0 的数据
        s = PhaseAwareSeed(0.0)
        @test s isa SeedStrategy
        # 如果数据文件存在，boundary_data 不为空
    end
    
    @testset "get_seed 无数据" begin
        s = PhaseAwareSeed()
        θ = [0.5, 1.0]
        seed = get_seed(s, θ, FixedMu())
        @test length(seed) == 5
    end
end

# ============================================================================
# PhaseAwareContinuitySeed 测试
# ============================================================================

@testset "PhaseAwareContinuitySeed" begin
    @testset "构造" begin
        s = PhaseAwareContinuitySeed()
        @test s isa SeedStrategy
        @test s.previous_solution === nothing
        @test s.previous_phase == :unknown
    end
    
    @testset "从 xi 构造" begin
        s = PhaseAwareContinuitySeed(0.0)
        @test s isa SeedStrategy
    end
    
    @testset "get_seed 初始" begin
        s = PhaseAwareContinuitySeed()
        θ = [0.5, 1.0]
        seed = get_seed(s, θ, FixedMu())
        @test length(seed) == 5
    end
    
    @testset "update! 带相位" begin
        s = PhaseAwareContinuitySeed()
        solution = [-1.5, -1.5, -2.1, 0.2, 0.2]
        T_MeV = 100.0
        μ_MeV = 200.0
        
        update!(s, solution, T_MeV, μ_MeV)
        
        @test s.previous_solution !== nothing
        @test s.previous_solution ≈ solution
    end
    
    @testset "reset!" begin
        s = PhaseAwareContinuitySeed()
        update!(s, [-1.5, -1.5, -2.1, 0.2, 0.2], 100.0, 200.0)
        
        reset!(s)
        @test s.previous_solution === nothing
        @test s.previous_phase == :unknown
    end
    
    @testset "set_phase!" begin
        s = PhaseAwareContinuitySeed()
        set_phase!(s, :hadron)
        @test s.previous_phase == :hadron
        
        set_phase!(s, :quark)
        @test s.previous_phase == :quark
    end
end

# ============================================================================
# PhaseBoundaryData 测试
# ============================================================================

@testset "PhaseBoundaryData" begin
    @testset "load_phase_boundary" begin
        # 尝试加载数据
        data = load_phase_boundary(0.0)
        @test data isa PhaseBoundaryData
        @test data.xi == 0.0
        
        # 如果数据存在，检查基本属性
        if !isempty(data.T_values)
            @test length(data.T_values) == length(data.mu_values)
            @test issorted(data.T_values)
        end
    end
    
    @testset "interpolate_mu_c" begin
        data = load_phase_boundary(0.0)
        
        if !isempty(data.T_values)
            # 在数据范围内插值
            T_mid = (data.T_values[1] + data.T_values[end]) / 2
            μ_c = interpolate_mu_c(data, T_mid)
            @test isfinite(μ_c)
            
            # 超过 CEP 返回 NaN
            if !isnan(data.T_CEP)
                μ_c_above = interpolate_mu_c(data, data.T_CEP + 10)
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
        seed = extend_seed(base_seed, FixedMu())
        @test length(seed) == 5
        @test seed ≈ base_seed
    end
    
    @testset "FixedRho" begin
        seed = extend_seed(base_seed, FixedRho(1.0))
        @test length(seed) == 8
        @test seed[1:5] ≈ base_seed
        @test all(isfinite.(seed[6:8]))
    end
    
    @testset "FixedEntropy" begin
        seed = extend_seed(base_seed, FixedEntropy(0.5))
        @test length(seed) == 8
        @test seed[1:5] ≈ base_seed
    end
    
    @testset "FixedSigma" begin
        seed = extend_seed(base_seed, FixedSigma(10.0))
        @test length(seed) == 8
        @test seed[1:5] ≈ base_seed
    end
end
