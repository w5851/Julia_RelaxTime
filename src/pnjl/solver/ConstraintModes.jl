"""
    ConstraintModes

PNJL 求解器约束模式定义。

## 支持的模式
- `FixedMu`: 固定化学势 μ
- `FixedRho`: 固定重子数密度 ρ
- `FixedEntropy`: 固定熵密度 s
- `FixedSigma`: 固定比熵 σ = s/ρ

## 使用示例
```julia
mode = FixedMu()
mode = FixedRho(1.0)  # ρ/ρ₀ = 1.0
mode = FixedEntropy(0.5)
mode = FixedSigma(10.0)
```
"""
module ConstraintModes

export ConstraintMode, FixedMu, FixedRho, FixedEntropy, FixedSigma
export state_dim, param_dim, constraint_description

# ============================================================================
# 抽象类型
# ============================================================================

"""
    ConstraintMode

求解模式的抽象基类型。所有具体模式都应继承此类型。
"""
abstract type ConstraintMode end

# ============================================================================
# 具体模式定义
# ============================================================================

"""
    FixedMu <: ConstraintMode

固定化学势模式。

参数 θ = [T, μ]，状态变量 x = [φ_u, φ_d, φ_s, Φ, Φ̄]（5维）。
方程：5 个能隙方程。
"""
struct FixedMu <: ConstraintMode end

"""
    FixedRho <: ConstraintMode

固定重子数密度模式。

参数 θ = [T]，状态变量 x = [φ_u, φ_d, φ_s, Φ, Φ̄, μ_u, μ_d, μ_s]（8维）。
方程：5 个能隙方程 + μ_u=μ_d + μ_d=μ_s + ρ约束。

# 字段
- `rho_target::Float64`: 目标重子数密度（归一化，ρ/ρ₀）
"""
struct FixedRho <: ConstraintMode
    rho_target::Float64
end

"""
    FixedEntropy <: ConstraintMode

固定熵密度模式。

参数 θ = [T]，状态变量 x = [φ_u, φ_d, φ_s, Φ, Φ̄, μ_u, μ_d, μ_s]（8维）。
方程：5 个能隙方程 + μ_u=μ_d + μ_d=μ_s + s约束。

# 字段
- `s_target::Float64`: 目标熵密度 (fm⁻³)
"""
struct FixedEntropy <: ConstraintMode
    s_target::Float64
end

"""
    FixedSigma <: ConstraintMode

固定比熵模式（等熵线）。

参数 θ = [T]，状态变量 x = [φ_u, φ_d, φ_s, Φ, Φ̄, μ_u, μ_d, μ_s]（8维）。
方程：5 个能隙方程 + μ_u=μ_d + μ_d=μ_s + σ约束。

# 字段
- `sigma_target::Float64`: 目标比熵 σ = s/n_B
"""
struct FixedSigma <: ConstraintMode
    sigma_target::Float64
end

# ============================================================================
# 模式属性查询
# ============================================================================

"""
    state_dim(mode::ConstraintMode) -> Int

返回状态变量的维度。
"""
state_dim(::FixedMu) = 5
state_dim(::FixedRho) = 8
state_dim(::FixedEntropy) = 8
state_dim(::FixedSigma) = 8

"""
    param_dim(mode::ConstraintMode) -> Int

返回参数 θ 的维度。
"""
param_dim(::FixedMu) = 2  # [T, μ]
param_dim(::FixedRho) = 1  # [T]
param_dim(::FixedEntropy) = 1  # [T]
param_dim(::FixedSigma) = 1  # [T]

"""
    constraint_description(mode::ConstraintMode) -> String

返回模式的描述字符串。
"""
constraint_description(::FixedMu) = "Fixed chemical potential μ"
constraint_description(m::FixedRho) = "Fixed baryon density ρ/ρ₀ = $(m.rho_target)"
constraint_description(m::FixedEntropy) = "Fixed entropy density s = $(m.s_target) fm⁻³"
constraint_description(m::FixedSigma) = "Fixed specific entropy σ = s/n_B = $(m.sigma_target)"

# ============================================================================
# 显示方法
# ============================================================================

Base.show(io::IO, ::FixedMu) = print(io, "FixedMu()")
Base.show(io::IO, m::FixedRho) = print(io, "FixedRho(ρ/ρ₀=$(m.rho_target))")
Base.show(io::IO, m::FixedEntropy) = print(io, "FixedEntropy(s=$(m.s_target))")
Base.show(io::IO, m::FixedSigma) = print(io, "FixedSigma(σ=$(m.sigma_target))")

end # module ConstraintModes
