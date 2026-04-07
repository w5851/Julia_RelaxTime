"""
    SeedStrategies

PNJL 求解器初值策略模块。

## 支持的策略
- `DefaultSeed`: 基于物理直觉的固定默认值
- `MultiSeed`: 多初值尝试（处理多值解）
- `HybridContinuitySeed`: 连续性优先（无历史解时回落到 continuity fallback）

## 使用示例
```julia
# 默认策略
strategy = DefaultSeed()
seed = get_seed(strategy, [T, μ], FixedMu())

# 多初值候选全评估（由上层治理规则统一选优）
result = solve(FixedMu(), T, μ; xi=xi, seed_strategy=MultiSeed())
```
"""
module SeedStrategies

# 从 Models 域导入约束模式定义
import Main.Models: ConstraintMode, FixedMu, FixedRho, FixedAsymmetricRho, FixedEntropy, FixedSigma, state_dim

export SeedStrategy, DefaultSeed, MultiSeed, HybridContinuitySeed
export get_seed, update!, extend_seed
export HADRON_SEED_5, QUARK_SEED_5, HADRON_SEED_8, QUARK_SEED_8
export MEDIUM_SEED_5, HIGH_DENSITY_SEED_5, HIGH_TEMP_SEED_5, VERY_HIGH_TEMP_SEED_5
export MEDIUM_SEED_8, HIGH_DENSITY_SEED_8

# ============================================================================
# 内置默认种子值（基于 trho_seed_table.csv 数据分析，xi=0.0）
# ============================================================================

"""
强子相典型初值（5维：φ_u, φ_d, φ_s, Φ, Φ̄）
来源：T=50 MeV, μ≈0, ρ=0 的收敛解
特征：手征凝聚完整，Polyakov loop 接近零（禁闭相）
"""
const HADRON_SEED_5 = [-1.84329, -1.84329, -2.22701, 1.0e-5, 4.0e-5]

"""
中等密度初值（5维）
来源：T=100 MeV, ρ=1.0 的收敛解
特征：部分手征恢复，Polyakov loop 中等
"""
const MEDIUM_SEED_5 = [-1.3647, -1.3647, -2.14502, 0.10594, 0.15569]

"""
高密度初值（5维）
来源：T=100 MeV, ρ=3.0 的收敛解
特征：手征对称性大部分恢复，Polyakov loop 较高
"""
const HIGH_DENSITY_SEED_5 = [-0.21695, -0.21695, -2.01372, 0.18601, 0.22333]

"""
高温初值（5维）
来源：T=200 MeV, ρ=0 的收敛解
特征：高温解禁闭相，Polyakov loop 接近 0.6
"""
const HIGH_TEMP_SEED_5 = [-0.73192, -0.73192, -1.79539, 0.60532, 0.60532]

"""
更高温初值（5维）
用途：T≈300–400 MeV 区间的数值求解更稳健的起点（更接近手征恢复 + 去禁闭）。
注：这是经验种子（并非唯一），用于提升收敛到物理解（Φ/Φ̄∈[0,1]、有效质量为正）概率。
"""
const VERY_HIGH_TEMP_SEED_5 = [-0.30, -0.30, -0.90, 0.90, 0.90]

# 额外候选：用于提高极端参数（高温/强各向异性）下的收敛鲁棒性。
# 与 scripts/pnjl/diagnose_gap_point.jl 里的人工候选保持一致。
const HT_GUESS_0p8_SEED_5 = [-0.50, -0.50, -1.20, 0.80, 0.80]
const HT_GUESS_0p9_SEED_5 = VERY_HIGH_TEMP_SEED_5
const HT_GUESS_0p95_SEED_5 = [-0.20, -0.20, -0.70, 0.95, 0.95]
const WEAK_CHIRAL_CONF_SEED_5 = [-0.50, -0.50, -1.20, 1e-3, 1e-3]

# 兼容旧接口的别名
"""夸克相典型初值（高温或高密度）"""
const QUARK_SEED_5 = HIGH_TEMP_SEED_5

"""强子相 8 维初值（含 μ_u, μ_d, μ_s，单位 fm⁻¹）"""
const HADRON_SEED_8 = [HADRON_SEED_5..., 0.22367, 0.22367, 0.22367]

"""中等密度 8 维初值"""
const MEDIUM_SEED_8 = [MEDIUM_SEED_5..., 1.70267, 1.70267, 1.70267]

"""高密度 8 维初值"""
const HIGH_DENSITY_SEED_8 = [HIGH_DENSITY_SEED_5..., 1.7516, 1.7516, 1.7516]

"""夸克相 8 维初值"""
const QUARK_SEED_8 = [QUARK_SEED_5..., 0.5, 0.5, 0.5]

# ============================================================================
# 抽象类型
# ============================================================================

"""
    SeedStrategy

初值策略的抽象基类型。所有具体策略都应继承此类型。
"""
abstract type SeedStrategy end

# ============================================================================
# 辅助函数
# ============================================================================

"""
    auto_phase_hint(T_fm, μ_fm) -> Symbol

基于温度和化学势自动判断相位。

简单启发式：高温(>150 MeV)或高化学势(>300 MeV)倾向夸克相。
"""
function auto_phase_hint(T_fm::Real, μ_fm::Real)
    ħc = 197.327
    T_mev = T_fm * ħc
    μ_mev = μ_fm * ħc
    return (T_mev > 150 || μ_mev > 300) ? :quark : :hadron
end

"""
    estimate_mu_from_rho(rho_norm) -> Float64

基于归一化密度估计化学势（自由费米气体近似）。

粗略近似：ρ/ρ₀ ≈ (μ/μ₀)³
"""
function estimate_mu_from_rho(rho_norm::Real)
    μ₀ = 1.5  # fm⁻¹，约 300 MeV
    return μ₀ * cbrt(max(rho_norm, 0.1))
end

"""
    extend_seed(base_seed, mode) -> Vector{Float64}

根据求解模式扩展基础种子。
"""
function extend_seed(base_seed::AbstractVector{<:Real}, ::FixedMu)
    return collect(Float64, base_seed[1:5])
end

function extend_seed(base_seed::AbstractVector{<:Real}, mode::FixedRho)
    seed_5 = base_seed[1:5]
    μ_guess = estimate_mu_from_rho(mode.rho_target)
    return Float64[seed_5..., μ_guess, μ_guess, μ_guess]
end

function extend_seed(base_seed::AbstractVector{<:Real}, mode::FixedAsymmetricRho)
    seed_5 = base_seed[1:5]
    μ_d = estimate_mu_from_rho(mode.rho_target)
    ratio = max(mode.ud_ratio_target, 1e-6)
    μ_u = μ_d * cbrt(ratio)
    μ_s = abs(mode.s_target) <= 1e-12 ? 0.05 : μ_d
    return Float64[seed_5..., μ_u, μ_d, μ_s]
end

function extend_seed(base_seed::AbstractVector{<:Real}, mode::FixedEntropy)
    seed_5 = base_seed[1:5]
    # 熵密度模式：μ 初值使用中等值
    μ_guess = 1.0  # fm⁻¹
    return Float64[seed_5..., μ_guess, μ_guess, μ_guess]
end

function extend_seed(base_seed::AbstractVector{<:Real}, mode::FixedSigma)
    seed_5 = base_seed[1:5]
    # 比熵模式：μ 初值使用中等值
    μ_guess = 1.0  # fm⁻¹
    return Float64[seed_5..., μ_guess, μ_guess, μ_guess]
end

# ============================================================================
# 策略1：固定默认值
# ============================================================================

"""
    DefaultSeed <: SeedStrategy

基于物理直觉的固定默认值策略。

# 字段
- `hadron_seed::Vector{Float64}`: 强子相典型值
- `quark_seed::Vector{Float64}`: 夸克相典型值
- `phase_hint::Symbol`: 相位提示 (:hadron, :quark, 或 :auto)
"""
struct DefaultSeed <: SeedStrategy
    hadron_seed::Vector{Float64}
    quark_seed::Vector{Float64}
    phase_hint::Symbol
end

"""创建默认策略，使用内置种子值"""
function DefaultSeed(; phase_hint::Symbol=:auto)
    return DefaultSeed(copy(HADRON_SEED_5), copy(QUARK_SEED_5), phase_hint)
end

function get_seed(s::DefaultSeed, θ::AbstractVector, mode::ConstraintMode)
    hint = s.phase_hint
    if hint == :auto && length(θ) >= 2
        hint = auto_phase_hint(θ[1], θ[2])
    elseif hint == :auto
        hint = :hadron  # 默认强子相
    end
    
    base = if hint == :quark
        # 高温（或高化学势）时，T=200MeV 的 quark 种子有时会收敛到非物理解；对更高温区间使用更“去禁闭/手征恢复”的种子。
        if length(θ) >= 1
            T_fm = θ[1]
            T_mev = T_fm * 197.327
            (T_mev >= 300.0) ? VERY_HIGH_TEMP_SEED_5 : s.quark_seed
        else
            s.quark_seed
        end
    else
        s.hadron_seed
    end
    return extend_seed(base, mode)
end

# ============================================================================
# 策略2：多初值尝试
# ============================================================================

"""
    MultiSeed <: SeedStrategy

多初值尝试策略，用于处理多值解。

# 字段
- `candidates::Vector{SeedStrategy}`: 候选策略列表
"""
struct MultiSeed <: SeedStrategy
    candidates::Vector{SeedStrategy}
end

"""创建多初值策略，默认尝试强子相和夸克相"""
function MultiSeed()
    candidates = [
        DefaultSeed(phase_hint=:hadron),
        DefaultSeed(phase_hint=:quark),
        # 更偏禁闭：Φ 很小，但凝聚较弱（用于避免卡到坏分支）
        DefaultSeed(copy(WEAK_CHIRAL_CONF_SEED_5), copy(WEAK_CHIRAL_CONF_SEED_5), :hadron),
        # 人工高温候选：凝聚幅度更小、Polyakov loop 更高
        DefaultSeed(copy(HT_GUESS_0p8_SEED_5), copy(HT_GUESS_0p8_SEED_5), :hadron),
        DefaultSeed(copy(HT_GUESS_0p9_SEED_5), copy(HT_GUESS_0p9_SEED_5), :hadron),
        DefaultSeed(copy(HT_GUESS_0p95_SEED_5), copy(HT_GUESS_0p95_SEED_5), :hadron),
    ]
    return MultiSeed(candidates)
end

function get_seed(s::MultiSeed, θ::AbstractVector, mode::ConstraintMode)
    # 返回第一个候选的种子
    # 实际多初值尝试在求解器层面通过 get_all_seeds 实现
    return get_seed(s.candidates[1], θ, mode)
end

"""
    get_all_seeds(s::MultiSeed, θ, mode) -> Vector{Vector{Float64}}

获取所有候选初值。
"""
function get_all_seeds(s::MultiSeed, θ::AbstractVector, mode::ConstraintMode)
    return [get_seed(c, θ, mode) for c in s.candidates]
end

export get_all_seeds

# ============================================================================
# 策略3b：连续性优先
# ============================================================================

"""
    HybridContinuitySeed <: SeedStrategy

连续性优先的混合策略：
- 优先使用上一点解（连续性）
- 若当前点无历史解，则回落到 `continuity_fallback`
"""
mutable struct HybridContinuitySeed <: SeedStrategy
    previous_solution::Union{Nothing, Vector{Float64}}
    continuity_fallback::SeedStrategy
end

"""创建混合连续性策略"""
function HybridContinuitySeed(; continuity_fallback::SeedStrategy=DefaultSeed())
    return HybridContinuitySeed(nothing, continuity_fallback)
end

function get_seed(s::HybridContinuitySeed, θ::AbstractVector, mode::ConstraintMode)
    if s.previous_solution !== nothing
        expected_dim = state_dim(mode)
        if length(s.previous_solution) == expected_dim
            return copy(s.previous_solution)
        elseif length(s.previous_solution) >= 5
            return extend_seed(s.previous_solution, mode)
        end
    end
    return get_seed(s.continuity_fallback, θ, mode)
end

function update!(s::HybridContinuitySeed, solution::AbstractVector{<:Real})
    s.previous_solution = collect(Float64, solution)
    return s
end

# ============================================================================
# 显示方法
# ============================================================================

Base.show(io::IO, s::DefaultSeed) = print(io, "DefaultSeed(phase_hint=$(s.phase_hint))")
Base.show(io::IO, s::MultiSeed) = print(io, "MultiSeed($(length(s.candidates)) candidates)")
Base.show(io::IO, s::HybridContinuitySeed) = print(io, "HybridContinuitySeed(prev=$(s.previous_solution !== nothing))")

end # module SeedStrategies

