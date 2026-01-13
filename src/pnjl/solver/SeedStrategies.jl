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
export MEDIUM_SEED_5, HIGH_DENSITY_SEED_5, HIGH_TEMP_SEED_5, VERY_HIGH_TEMP_SEED_5
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

"""
更高温初值（5维）
用途：T≈300–400 MeV 区间的数值求解更稳健的起点（更接近手征恢复 + 去禁闭）。
注：这是经验种子（并非唯一），用于提升收敛到物理解（Φ/Φ̄∈[0,1]、有效质量为正）概率。
"""
const VERY_HIGH_TEMP_SEED_5 = [-0.30, -0.30, -0.90, 0.90, 0.90]

# 额外候选：用于提高极端参数（高温/强各向异性）下的收敛鲁棒性。
# 与 scripts/pnjl/diagnose_gap_point.jl 里的人工候选保持一致。
const HT_GUESS_0p8_SEED_5 = [-0.50, -0.50, -1.20, 0.80, 0.80]
const HT_GUESS_0p9_SEED_5 = [-0.30, -0.30, -0.90, 0.90, 0.90]
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
        # 更偏禁闭：Φ 很小，但凝聚较弱（用于避免卡到坏分支）
        DefaultSeed(copy(WEAK_CHIRAL_CONF_SEED_5), copy(WEAK_CHIRAL_CONF_SEED_5), :hadron),
        # 人工高温候选：凝聚幅度更小、Polyakov loop 更高
        DefaultSeed(copy(HT_GUESS_0p8_SEED_5), copy(HT_GUESS_0p8_SEED_5), :hadron),
        DefaultSeed(copy(HT_GUESS_0p9_SEED_5), copy(HT_GUESS_0p9_SEED_5), :hadron),
        DefaultSeed(copy(HT_GUESS_0p95_SEED_5), copy(HT_GUESS_0p95_SEED_5), :hadron),
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
    PhaseBoundaryData

相变线数据结构，存储从 CSV 加载的相变线信息。

# 字段
- `T_values::Vector{Float64}`: 温度数组 (MeV)
- `mu_values::Vector{Float64}`: 相变化学势数组 (MeV)
- `T_CEP::Float64`: 临界终点温度 (MeV)
- `mu_CEP::Float64`: 临界终点化学势 (MeV)
- `xi::Float64`: 各向异性参数
"""
struct PhaseBoundaryData
    T_values::Vector{Float64}   # MeV
    mu_values::Vector{Float64}  # MeV
    T_CEP::Float64              # MeV
    mu_CEP::Float64             # MeV
    xi::Float64
end

export PhaseBoundaryData, load_phase_boundary, interpolate_mu_c

"""
    load_phase_boundary(xi; boundary_path, cep_path) -> PhaseBoundaryData

从 CSV 文件加载相变线数据。

# 参数
- `xi`: 各向异性参数
- `boundary_path`: boundary.csv 路径（默认 data/reference/pnjl/boundary.csv）
- `cep_path`: cep.csv 路径（默认 data/reference/pnjl/cep.csv）

# 返回
PhaseBoundaryData 结构
"""
function load_phase_boundary(xi::Real;
    boundary_path::String = joinpath(@__DIR__, "..", "..", "..", "data", "reference", "pnjl", "boundary.csv"),
    cep_path::String = joinpath(@__DIR__, "..", "..", "..", "data", "reference", "pnjl", "cep.csv")
)
    # 加载 CEP 数据
    T_CEP = NaN
    mu_CEP = NaN
    if isfile(cep_path)
        for line in eachline(cep_path)
            startswith(line, "xi") && continue  # 跳过表头
            parts = split(line, ',')
            length(parts) >= 3 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            abs(xi_val - xi) > 1e-6 && continue
            T_CEP = tryparse(Float64, parts[2])
            mu_CEP = tryparse(Float64, parts[3])
            break
        end
    end
    
    # 加载相变线数据
    T_values = Float64[]
    mu_values = Float64[]
    if isfile(boundary_path)
        for line in eachline(boundary_path)
            startswith(line, "xi") && continue  # 跳过表头
            parts = split(line, ',')
            length(parts) >= 3 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            abs(xi_val - xi) > 1e-6 && continue
            T_val = tryparse(Float64, parts[2])
            mu_val = tryparse(Float64, parts[3])
            (T_val === nothing || mu_val === nothing) && continue
            push!(T_values, T_val)
            push!(mu_values, mu_val)
        end
    end
    
    # 按温度排序
    if !isempty(T_values)
        order = sortperm(T_values)
        T_values = T_values[order]
        mu_values = mu_values[order]
    end
    
    return PhaseBoundaryData(T_values, mu_values, T_CEP, mu_CEP, Float64(xi))
end

"""
    interpolate_mu_c(data::PhaseBoundaryData, T_MeV) -> Float64

线性插值获取给定温度下的相变化学势。

# 返回
- μ_c (MeV)，如果 T > T_CEP 或无数据则返回 NaN
"""
function interpolate_mu_c(data::PhaseBoundaryData, T_MeV::Real)
    # 检查是否超过 CEP
    if !isnan(data.T_CEP) && T_MeV > data.T_CEP
        return NaN  # Crossover 区域，无一阶相变
    end
    
    # 检查数据是否存在
    isempty(data.T_values) && return NaN
    
    T = Float64(T_MeV)
    Ts = data.T_values
    μs = data.mu_values
    
    # 边界检查
    if T <= Ts[1]
        return μs[1]  # 外推：使用最低温度的值
    elseif T >= Ts[end]
        return μs[end]  # 外推：使用最高温度的值
    end
    
    # 线性插值
    for i in 1:(length(Ts)-1)
        if Ts[i] <= T <= Ts[i+1]
            t = (T - Ts[i]) / (Ts[i+1] - Ts[i])
            return μs[i] + t * (μs[i+1] - μs[i])
        end
    end
    
    return NaN
end

"""
    PhaseAwareSeed <: SeedStrategy

基于相图的智能选择策略，用于一阶相变区域。

根据 (T, μ) 与相变线的关系选择初值：
- μ < μ_c(T): 使用强子相初值
- μ > μ_c(T): 使用夸克相初值
- T > T_CEP: 使用自动判断（crossover 区域）

# 字段
- `boundary_data::Union{Nothing, PhaseBoundaryData}`: 相变线数据
- `hadron_strategy::SeedStrategy`: 强子相策略
- `quark_strategy::SeedStrategy`: 夸克相策略
- `crossover_strategy::SeedStrategy`: Crossover 区域策略
"""
struct PhaseAwareSeed <: SeedStrategy
    boundary_data::Union{Nothing, PhaseBoundaryData}
    hadron_strategy::SeedStrategy
    quark_strategy::SeedStrategy
    crossover_strategy::SeedStrategy
end

"""
    PhaseAwareSeed(xi; kwargs...) -> PhaseAwareSeed

创建相图感知策略。

# 参数
- `xi`: 各向异性参数
- `boundary_path`: boundary.csv 路径（可选）
- `cep_path`: cep.csv 路径（可选）

# 示例
```julia
# 使用默认路径加载 xi=0.0 的数据
strategy = PhaseAwareSeed(0.0)

# 指定自定义路径
strategy = PhaseAwareSeed(0.2; boundary_path="my_boundary.csv")
```
"""
function PhaseAwareSeed(xi::Real; kwargs...)
    boundary_data = try
        load_phase_boundary(xi; kwargs...)
    catch e
        @warn "无法加载相变线数据: $(e)，将使用默认策略"
        nothing
    end
    
    return PhaseAwareSeed(
        boundary_data,
        DefaultSeed(phase_hint=:hadron),
        DefaultSeed(phase_hint=:quark),
        DefaultSeed(phase_hint=:auto),
    )
end

"""
    PhaseAwareSeed(phase_boundary::Function) -> PhaseAwareSeed

使用自定义相变线函数创建策略。

# 参数
- `phase_boundary`: 相变线函数，输入 T (fm⁻¹)，返回 μ_c (fm⁻¹)
"""
function PhaseAwareSeed(phase_boundary::Function)
    # 创建一个包装器，将函数转换为 PhaseBoundaryData 的行为
    return _PhaseAwareSeedWithFunction(
        phase_boundary,
        DefaultSeed(phase_hint=:hadron),
        DefaultSeed(phase_hint=:quark),
    )
end

# 内部类型：使用函数的 PhaseAwareSeed
struct _PhaseAwareSeedWithFunction <: SeedStrategy
    phase_boundary::Function
    hadron_strategy::SeedStrategy
    quark_strategy::SeedStrategy
end

function get_seed(s::_PhaseAwareSeedWithFunction, θ::AbstractVector, mode::ConstraintMode)
    T_fm = θ[1]
    μ_fm = length(θ) >= 2 ? θ[2] : 0.0
    
    μ_c = s.phase_boundary(T_fm)
    
    if μ_fm < μ_c
        return get_seed(s.hadron_strategy, θ, mode)
    else
        return get_seed(s.quark_strategy, θ, mode)
    end
end

"""创建使用默认相变线的策略（无数据文件时的回退）"""
function PhaseAwareSeed()
    return PhaseAwareSeed(
        nothing,
        DefaultSeed(phase_hint=:hadron),
        DefaultSeed(phase_hint=:quark),
        DefaultSeed(phase_hint=:auto),
    )
end

const ħc_MeV_fm = 197.327  # MeV·fm

function get_seed(s::PhaseAwareSeed, θ::AbstractVector, mode::ConstraintMode)
    T_fm = θ[1]
    μ_fm = length(θ) >= 2 ? θ[2] : 0.0
    
    # 转换为 MeV
    T_MeV = T_fm * ħc_MeV_fm
    μ_MeV = μ_fm * ħc_MeV_fm
    
    # 如果没有相变线数据，使用 crossover 策略
    if s.boundary_data === nothing
        return get_seed(s.crossover_strategy, θ, mode)
    end
    
    data = s.boundary_data
    
    # 检查是否在 crossover 区域
    if !isnan(data.T_CEP) && T_MeV > data.T_CEP
        return get_seed(s.crossover_strategy, θ, mode)
    end
    
    # 获取相变化学势
    μ_c_MeV = interpolate_mu_c(data, T_MeV)
    
    if isnan(μ_c_MeV)
        # 无法确定相变点，使用 crossover 策略
        return get_seed(s.crossover_strategy, θ, mode)
    end
    
    # 根据 μ 与 μ_c 的关系选择初值
    if μ_MeV < μ_c_MeV
        return get_seed(s.hadron_strategy, θ, mode)
    else
        return get_seed(s.quark_strategy, θ, mode)
    end
end

"""
    get_phase_hint(s::PhaseAwareSeed, T_MeV, μ_MeV) -> Symbol

获取给定 (T, μ) 点的相位提示。

# 返回
- `:hadron`: 强子相
- `:quark`: 夸克相
- `:crossover`: Crossover 区域
- `:unknown`: 无法确定
"""
function get_phase_hint(s::PhaseAwareSeed, T_MeV::Real, μ_MeV::Real)
    if s.boundary_data === nothing
        return :unknown
    end
    
    data = s.boundary_data
    
    # 检查是否在 crossover 区域
    if !isnan(data.T_CEP) && T_MeV > data.T_CEP
        return :crossover
    end
    
    # 获取相变化学势
    μ_c_MeV = interpolate_mu_c(data, T_MeV)
    
    if isnan(μ_c_MeV)
        return :unknown
    end
    
    return μ_MeV < μ_c_MeV ? :hadron : :quark
end

export get_phase_hint

# ============================================================================
# 策略5：相变感知的连续性跟踪
# ============================================================================

"""
    PhaseAwareContinuitySeed <: SeedStrategy

相变感知的连续性跟踪策略。

结合连续性跟踪和相变线感知：
- 正常情况下使用前一个解作为初值（连续性跟踪）
- 当跨越相变线时，切换到对应相的默认初值
- 切换后继续连续性跟踪

适用于 T-μ 扫描，可以正确处理一阶相变区域。

# 字段
- `boundary_data::Union{Nothing, PhaseBoundaryData}`: 相变线数据
- `previous_solution::Union{Nothing, Vector{Float64}}`: 上一个解
- `previous_phase::Symbol`: 上一个点的相位 (:hadron, :quark, :crossover, :unknown)
- `hadron_seed::Vector{Float64}`: 强子相默认初值
- `quark_seed::Vector{Float64}`: 夸克相默认初值
- `fallback::SeedStrategy`: 无数据时的回退策略
"""
mutable struct PhaseAwareContinuitySeed <: SeedStrategy
    boundary_data::Union{Nothing, PhaseBoundaryData}
    previous_solution::Union{Nothing, Vector{Float64}}
    previous_phase::Symbol
    hadron_seed::Vector{Float64}
    quark_seed::Vector{Float64}
    fallback::SeedStrategy
end

export PhaseAwareContinuitySeed

"""
    PhaseAwareContinuitySeed(xi; kwargs...) -> PhaseAwareContinuitySeed

创建相变感知的连续性跟踪策略。

# 参数
- `xi`: 各向异性参数
- `boundary_path`: boundary.csv 路径（可选）
- `cep_path`: cep.csv 路径（可选）

# 示例
```julia
# 创建策略
tracker = PhaseAwareContinuitySeed(0.0)

# 扫描循环
for μ in μ_range
    seed = get_seed(tracker, [T_fm, μ_fm], FixedMu())
    result = solve(...)
    update!(tracker, result.solution, T_MeV, μ_MeV)
end

# 重置（开始新的扫描）
reset!(tracker)
```
"""
function PhaseAwareContinuitySeed(xi::Real; kwargs...)
    boundary_data = try
        load_phase_boundary(xi; kwargs...)
    catch e
        @warn "无法加载相变线数据: $(e)，将使用普通连续性跟踪"
        nothing
    end
    
    return PhaseAwareContinuitySeed(
        boundary_data,
        nothing,           # previous_solution
        :unknown,          # previous_phase
        copy(HADRON_SEED_5),
        copy(QUARK_SEED_5),
        DefaultSeed(phase_hint=:auto),
    )
end

"""无参数构造函数（无相变线数据）"""
function PhaseAwareContinuitySeed()
    return PhaseAwareContinuitySeed(
        nothing,
        nothing,
        :unknown,
        copy(HADRON_SEED_5),
        copy(QUARK_SEED_5),
        DefaultSeed(phase_hint=:auto),
    )
end

"""
    _get_current_phase(s::PhaseAwareContinuitySeed, T_MeV, μ_MeV) -> Symbol

获取当前点的相位。
"""
function _get_current_phase(s::PhaseAwareContinuitySeed, T_MeV::Real, μ_MeV::Real)
    if s.boundary_data === nothing
        return :unknown
    end
    
    data = s.boundary_data
    
    # 检查是否在 crossover 区域
    if !isnan(data.T_CEP) && T_MeV > data.T_CEP
        return :crossover
    end
    
    # 获取相变化学势
    μ_c_MeV = interpolate_mu_c(data, T_MeV)
    
    if isnan(μ_c_MeV)
        return :unknown
    end
    
    return μ_MeV < μ_c_MeV ? :hadron : :quark
end

"""
    _is_phase_transition(prev_phase, curr_phase) -> Bool

判断是否发生了相变（跨越相变线）。
"""
function _is_phase_transition(prev_phase::Symbol, curr_phase::Symbol)
    # 只有 hadron <-> quark 的转变才算相变
    return (prev_phase == :hadron && curr_phase == :quark) ||
           (prev_phase == :quark && curr_phase == :hadron)
end

function get_seed(s::PhaseAwareContinuitySeed, θ::AbstractVector, mode::ConstraintMode)
    T_fm = θ[1]
    μ_fm = length(θ) >= 2 ? θ[2] : 0.0
    
    # 转换为 MeV
    T_MeV = T_fm * ħc_MeV_fm
    μ_MeV = μ_fm * ħc_MeV_fm
    
    # 获取当前相位
    current_phase = _get_current_phase(s, T_MeV, μ_MeV)
    
    # 情况1：没有前一个解（第一个点）
    if s.previous_solution === nothing
        # 根据当前相位选择初值
        if current_phase == :hadron
            return extend_seed(s.hadron_seed, mode)
        elseif current_phase == :quark
            return extend_seed(s.quark_seed, mode)
        else
            # crossover 或 unknown，使用回退策略
            return get_seed(s.fallback, θ, mode)
        end
    end
    
    # 情况2：发生相变（跨越相变线）
    if _is_phase_transition(s.previous_phase, current_phase)
        # 切换到新相的默认初值
        if current_phase == :hadron
            return extend_seed(s.hadron_seed, mode)
        else  # :quark
            return extend_seed(s.quark_seed, mode)
        end
    end
    
    # 情况3：正常连续性跟踪
    expected_dim = state_dim(mode)
    if length(s.previous_solution) == expected_dim
        return copy(s.previous_solution)
    elseif length(s.previous_solution) >= 5
        return extend_seed(s.previous_solution, mode)
    end
    
    # 回退
    return get_seed(s.fallback, θ, mode)
end

"""
    update!(s::PhaseAwareContinuitySeed, solution, T_MeV, μ_MeV)

更新连续性跟踪器的解和相位信息。

# 参数
- `solution`: 当前点的解
- `T_MeV`: 当前温度 (MeV)
- `μ_MeV`: 当前化学势 (MeV)
"""
function update!(s::PhaseAwareContinuitySeed, solution::AbstractVector{<:Real}, 
                 T_MeV::Real, μ_MeV::Real)
    s.previous_solution = collect(Float64, solution)
    s.previous_phase = _get_current_phase(s, T_MeV, μ_MeV)
    return s
end

"""
    update!(s::PhaseAwareContinuitySeed, solution)

简化版更新（不更新相位，用于 crossover 区域）。
"""
function update!(s::PhaseAwareContinuitySeed, solution::AbstractVector{<:Real})
    s.previous_solution = collect(Float64, solution)
    # 不更新 previous_phase，保持之前的值
    return s
end

"""
    reset!(s::PhaseAwareContinuitySeed)

重置跟踪器状态。
"""
function reset!(s::PhaseAwareContinuitySeed)
    s.previous_solution = nothing
    s.previous_phase = :unknown
    return s
end

"""
    set_phase!(s::PhaseAwareContinuitySeed, phase::Symbol)

手动设置当前相位（用于特殊情况）。
"""
function set_phase!(s::PhaseAwareContinuitySeed, phase::Symbol)
    s.previous_phase = phase
    return s
end

export set_phase!

# ============================================================================
# 显示方法
# ============================================================================

Base.show(io::IO, s::DefaultSeed) = print(io, "DefaultSeed(phase_hint=$(s.phase_hint))")
Base.show(io::IO, s::MultiSeed) = print(io, "MultiSeed($(length(s.candidates)) candidates)")
Base.show(io::IO, s::ContinuitySeed) = print(io, "ContinuitySeed(has_previous=$(s.previous_solution !== nothing))")
function Base.show(io::IO, s::PhaseAwareSeed)
    if s.boundary_data !== nothing
        n = length(s.boundary_data.T_values)
        xi = s.boundary_data.xi
        print(io, "PhaseAwareSeed(xi=$xi, $n boundary points)")
    else
        print(io, "PhaseAwareSeed(no data)")
    end
end
Base.show(io::IO, s::_PhaseAwareSeedWithFunction) = print(io, "PhaseAwareSeed(custom function)")
function Base.show(io::IO, s::PhaseAwareContinuitySeed)
    has_data = s.boundary_data !== nothing
    has_prev = s.previous_solution !== nothing
    print(io, "PhaseAwareContinuitySeed(data=$has_data, prev=$has_prev, phase=$(s.previous_phase))")
end

end # module SeedStrategies
