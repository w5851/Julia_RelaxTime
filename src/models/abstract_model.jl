"""abstract_model.jl

新架构的抽象类型与接口定义。

设计目标：
- 明确“能组装 Ω 所需的最小接口”
- 通过多重派发为 NJL/PNJL/rPNJL 等模型提供统一扩展点
- 接口默认抛出 MethodError，避免 silent fallback
"""

# -------------------------
# Abstract types
# -------------------------

"""所有模型的根抽象类型。"""
abstract type AbstractQCDModel end

"""NJL 家族模型的抽象类型（NJL/PNJL/rPNJL 等）。"""
abstract type AbstractNJLModel <: AbstractQCDModel end

"""PNJL 家族模型的抽象类型（包含 Polyakov 变量 Φ/Φbar，典型状态维度为 5）。

说明：
- 该抽象类型用于区分“仅 3 维 φ”的 NJL 与“包含 Φ/Φbar”的 PNJL/rPNJL。
- 便于在 gap 求解/隐函数求导等层面定义 5 维默认行为。
"""
abstract type AbstractPNJLModel <: AbstractNJLModel end

struct UnsupportedCapabilityError <: Exception
    model_kind::Symbol
    capability::Symbol
end

Base.@kwdef struct ModelCapabilities
    supports_solve_gap::Bool = true
    supports_model_thermo::Bool = true
    supports_number_densities::Bool = true
end

@inline model_kind_symbol(model::AbstractQCDModel) = Symbol(nameof(typeof(model)))

@inline function model_capabilities(::AbstractQCDModel)
    return ModelCapabilities()
end

function Base.showerror(io::IO, err::UnsupportedCapabilityError)
    print(io, "unsupported capability ", err.capability, " for model ", err.model_kind)
end

@inline function require_capability(model::AbstractQCDModel, capability::Symbol)
    caps = model_capabilities(model)
    key = Symbol("supports_", capability)
    hasproperty(caps, key) || throw(UnsupportedCapabilityError(model_kind_symbol(model), capability))
    getproperty(caps, key) || throw(UnsupportedCapabilityError(model_kind_symbol(model), capability))
    return nothing
end

# -------------------------
# Gap residual (stationarity conditions)
# -------------------------

"""gap_state_dim(model) -> Int

返回 gap 求解的未知量维度。

约定：
- NJL：3（φ_u, φ_d, φ_s）
- PNJL：5（φ_u, φ_d, φ_s, Φ, Φbar）

说明：该函数本身只是“契约入口”（用于通用 residual/solver），具体方法在求解层文件中提供。
"""
function gap_state_dim end

"""gap_residual(model, x, T, mu_vec; kwargs...) -> SVector

通用 gap 条件残差（驻点条件）。默认实现通常取 ∇ₓΩ(x)=0，
因此对 NJL(3 维) 与 PNJL(5 维) 在形式上是通用的。

说明：该函数本身只是“契约入口”（用于通用 residual/solver），具体方法在求解层文件中提供。
"""
function gap_residual end

# -------------------------
# Minimal required interface
# -------------------------

"""solve_gap(model, T, mu_vec; kwargs...) -> x_state

返回给定 (T, μ) 下的平衡平均场解 `x_state`。

阶段 0 契约：
- `mu_vec` 建议为长度 3 的向量（u,d,s），也允许 `Real`（视为 (μ,μ,μ)）。
- 返回值推荐为 `MeanFieldState`（或可被 `MeanFieldState(x_state)` 接受的形式）。

说明：这是“求解层”的入口；后续阶段会把 legacy 求解器包装成实现此接口的适配器。
"""
function solve_gap(model::AbstractQCDModel, T, mu_vec; kwargs...)
    require_capability(model, :solve_gap)
    throw(MethodError(solve_gap, (model, T, mu_vec)))
end

"""number_densities(model, x_state, T, mu_vec; thermal_nodes, xi=0.0) -> NamedTuple

返回夸克与反夸克的数密度（用于 RTA τ 等）。

约定返回：`(quark=SVector{3}(...), antiquark=SVector{3}(...))`，单位 fm⁻³。

注意：
- NJL 类模型通常不包含 Polyakov 变量，但仍可定义此接口（例如令 Φ=Φbar=1）。
- PNJL/rPNJL 通常需要 `thermal_nodes` 与各向异性参数 `xi`。
"""
function number_densities(model::AbstractQCDModel, x_state, T, mu_vec; kwargs...)
    require_capability(model, :number_densities)
    throw(MethodError(number_densities, (model, x_state, T, mu_vec)))
end

"""calculate_mass_vec(model, φ) -> SVector{3}

返回组分夸克质量向量 (Mu, Md, Ms)。
"""
function calculate_mass_vec(model::AbstractQCDModel, φ; kwargs...)
    throw(MethodError(calculate_mass_vec, (model, φ)))
end

"""calculate_chiral(model, φ) -> Real

返回平均场凝聚（手征）项 χ(φ)。
"""
function calculate_chiral(model::AbstractQCDModel, φ; kwargs...)
    throw(MethodError(calculate_chiral, (model, φ)))
end

"""vacuum_contribution(model, masses) -> Real

返回真空项贡献（通常含截断 Λ 的发散正规化）。
"""
function vacuum_contribution(model::AbstractQCDModel, masses; kwargs...)
    throw(MethodError(vacuum_contribution, (model, masses)))
end

"""thermal_contribution(model, masses, mu_vec, T) -> Real

返回热激发项贡献。

NJL: 一般依赖 (T, μ) 的费米分布；
PNJL: 需要同时依赖 Polyakov (Φ, Φbar)（此时具体模型应覆盖签名）。
"""
function thermal_contribution(model::AbstractQCDModel, masses, mu_vec, T; kwargs...)
    throw(MethodError(thermal_contribution, (model, masses, mu_vec, T)))
end

"""thermal_contribution(model, masses, Φ, Φbar, mu_vec, T) -> Real

Polyakov-aware 的热项签名：PNJL/rPNJL 等需要 (Φ, Φbar)。
NJL 可忽略 Φ/Φbar，但也建议实现该签名以便统一的 `omega` 入口组装。
"""
function thermal_contribution(model::AbstractQCDModel, masses, Φ, Φbar, mu_vec, T; kwargs...)
    throw(MethodError(thermal_contribution, (model, masses, Φ, Φbar, mu_vec, T)))
end

"""polyakov_potential(model, Φ, Φbar, T) -> Real

返回 Polyakov loop 有效势。NJL 模型应显式实现为 0（而不是依赖默认）。
"""
function polyakov_potential(model::AbstractQCDModel, Φ, Φbar, T; kwargs...)
    throw(MethodError(polyakov_potential, (model, Φ, Φbar, T)))
end
