"""
    SeedStrategies

PNJL 求解器初值策略模块。

## 支持的策略
- `DefaultSeed`: 基于物理直觉的固定默认值
- `MultiSeed`: 多初值尝试（处理多值解）
- `ContinuitySeed`: 连续性跟踪（参数扫描）
- `PhaseAwareSeed`: 基于相图的智能选择（一阶相变区域）

## 使用示例
```julia
# 默认策略
strategy = DefaultSeed()
seed = get_seed(strategy, [T, μ], FixedMu())

# 连续性跟踪
tracker = ContinuitySeed()
for μ in μ_range
    seed = get_seed(tracker, [T, μ], FixedMu())
    result = solve(...)
    update!(tracker, result.solution)
end
```
"""
module SeedStrategies

using StaticArrays

# 从父模块导入 ConstraintModes
using ..ConstraintModes: ConstraintMode, FixedMu, FixedRho, FixedEntropy, FixedSigma, state_dim

export SeedStrategy, DefaultSeed, MultiSeed, ContinuitySeed, PhaseAwareSeed
export get_seed, update!, extend_seed
export HADRON_SEED_5, QUARK_SEED_5, HADRON_SEED_8, QUARK_SEED_8
export MEDIUM_SEED_5, HIGH_DENSITY_SEED_5, HIGH_TEMP_SEED_5
export MEDIUM_SEED_8, HIGH_DENSITY_SEED_8
export default_omega_selector

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
    
    base = hint == :quark ? s.quark_seed : s.hadron_seed
    return extend_seed(base, mode)
end

# ============================================================================
# 策略2：多初值尝试
# ============================================================================

"""
    default_omega_selector(results) -> result

默认选择器：选择 Ω 最小的收敛解。
"""
function default_omega_selector(results)
    converged = filter(r -> r.converged, results)
    isempty(converged) && return first(results)  # 无收敛解时返回第一个
    return argmin(r -> r.omega, converged)
end

"""
    MultiSeed <: SeedStrategy

多初值尝试策略，用于处理多值解。

# 字段
- `candidates::Vector{SeedStrategy}`: 候选策略列表
- `selector::Function`: 结果选择器
"""
struct MultiSeed <: SeedStrategy
    candidates::Vector{SeedStrategy}
    selector::Function
end

"""创建多初值策略，默认尝试强子相和夸克相"""
function MultiSeed(; selector::Function=default_omega_selector)
    candidates = [
        DefaultSeed(phase_hint=:hadron),
        DefaultSeed(phase_hint=:quark),
    ]
    return MultiSeed(candidates, selector)
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
# 策略3：连续性跟踪
# ============================================================================

"""
    ContinuitySeed <: SeedStrategy

连续性跟踪策略，用于参数扫描。

使用前一个收敛解作为下一个点的初值。

# 字段
- `previous_solution::Union{Nothing, Vector{Float64}}`: 上一个解
- `fallback::SeedStrategy`: 回退策略
"""
mutable struct ContinuitySeed <: SeedStrategy
    previous_solution::Union{Nothing, Vector{Float64}}
    fallback::SeedStrategy
end

"""创建连续性跟踪策略"""
function ContinuitySeed(; fallback::SeedStrategy=DefaultSeed())
    return ContinuitySeed(nothing, fallback)
end

function get_seed(s::ContinuitySeed, θ::AbstractVector, mode::ConstraintMode)
    if s.previous_solution !== nothing
        # 确保维度匹配
        expected_dim = state_dim(mode)
        if length(s.previous_solution) == expected_dim
            return copy(s.previous_solution)
        elseif length(s.previous_solution) >= 5
            # 尝试扩展
            return extend_seed(s.previous_solution, mode)
        end
    end
    return get_seed(s.fallback, θ, mode)
end

"""
    update!(s::ContinuitySeed, solution)

更新连续性跟踪器的解。
"""
function update!(s::ContinuitySeed, solution::AbstractVector{<:Real})
    s.previous_solution = collect(Float64, solution)
    return s
end

"""
    reset!(s::ContinuitySeed)

重置连续性跟踪器。
"""
function reset!(s::ContinuitySeed)
    s.previous_solution = nothing
    return s
end

export reset!

# ============================================================================
# 策略4：基于相图的智能选择
# ============================================================================

"""
    PhaseAwareSeed <: SeedStrategy

基于相图的智能选择策略，用于一阶相变区域。

根据 (T, μ) 与相变线的关系选择初值。

# 字段
- `phase_boundary::Function`: 相变线函数 T -> μ_c
- `hadron_strategy::SeedStrategy`: 强子相策略
- `quark_strategy::SeedStrategy`: 夸克相策略
"""
struct PhaseAwareSeed <: SeedStrategy
    phase_boundary::Function
    hadron_strategy::SeedStrategy
    quark_strategy::SeedStrategy
end

"""
创建相图感知策略。

# 参数
- `phase_boundary`: 相变线函数，输入 T (fm⁻¹)，返回 μ_c (fm⁻¹)
"""
function PhaseAwareSeed(phase_boundary::Function)
    return PhaseAwareSeed(
        phase_boundary,
        DefaultSeed(phase_hint=:hadron),
        DefaultSeed(phase_hint=:quark),
    )
end

"""默认相变线近似（简单线性）"""
function default_phase_boundary(T_fm::Real)
    # 粗略近似：μ_c ≈ 1.5 - 0.5*T (fm⁻¹)
    # 对应约 300 MeV 在 T=0，随温度下降
    return max(0.0, 1.5 - 0.5 * T_fm)
end

"""创建使用默认相变线的策略"""
function PhaseAwareSeed()
    return PhaseAwareSeed(default_phase_boundary)
end

function get_seed(s::PhaseAwareSeed, θ::AbstractVector, mode::ConstraintMode)
    T_fm = θ[1]
    μ_fm = length(θ) >= 2 ? θ[2] : 0.0
    
    μ_c = s.phase_boundary(T_fm)
    
    if μ_fm < μ_c
        return get_seed(s.hadron_strategy, θ, mode)
    else
        return get_seed(s.quark_strategy, θ, mode)
    end
end

# ============================================================================
# 显示方法
# ============================================================================

Base.show(io::IO, s::DefaultSeed) = print(io, "DefaultSeed(phase_hint=$(s.phase_hint))")
Base.show(io::IO, s::MultiSeed) = print(io, "MultiSeed($(length(s.candidates)) candidates)")
Base.show(io::IO, s::ContinuitySeed) = print(io, "ContinuitySeed(has_previous=$(s.previous_solution !== nothing))")
Base.show(io::IO, s::PhaseAwareSeed) = print(io, "PhaseAwareSeed()")

end # module SeedStrategies
