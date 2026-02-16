"""LegacyNJLModel

阶段 1：legacy 适配器（adapter）。

说明：当前仓库中没有独立的“旧 NJL core”实现可供直接复用（不像 PNJL 有完整 legacy core+solver）。
因此这里先提供一个类型层面的适配器：它包装 `NJLModel` 并通过 Models 统一接口暴露。

好处：
- 高层工作流可以用 `create_model(:LegacyNJL)` 得到一个“可回退/可对比”的模型对象；
- 后续若补回真实 legacy NJL 公式/实现，只需替换本文件内部的委托逻辑。

当前实现：
- `calculate_mass_vec/calculate_chiral/vacuum_contribution/thermal_contribution/number_densities` 直接委托到 `NJLModel`。
- `solve_gap` 使用 `AbstractNJLModel` 的通用求解入口（见 gap_solver.jl），因此不需要在此重写。
"""

export LegacyNJLModel

"""legacy NJL 模型适配器（当前委托到 NJLModel）。"""
struct LegacyNJLModel <: AbstractNJLModel
    inner::NJLModel
end

function LegacyNJLModel(; profile::String=get(ENV, "NJL_PARAM_PROFILE", "default"))
    return LegacyNJLModel(NJLModel(; profile=profile))
end

# -----------------------------------------------------------------------------
# Unified entrypoints (delegate to NJLModel)
# -----------------------------------------------------------------------------

calculate_mass_vec(model::LegacyNJLModel, φ; kwargs...) = calculate_mass_vec(model.inner, φ; kwargs...)
calculate_chiral(model::LegacyNJLModel, φ; kwargs...) = calculate_chiral(model.inner, φ; kwargs...)
polyakov_potential(model::LegacyNJLModel, Φ, Φbar, T; kwargs...) = polyakov_potential(model.inner, Φ, Φbar, T; kwargs...)

vacuum_contribution(model::LegacyNJLModel, masses; kwargs...) = vacuum_contribution(model.inner, masses; kwargs...)
thermal_contribution(model::LegacyNJLModel, masses, Φ, Φbar, mu_vec, T; kwargs...) =
    thermal_contribution(model.inner, masses, Φ, Φbar, mu_vec, T; kwargs...)

number_densities(model::LegacyNJLModel, x_state, T, mu_vec; kwargs...) = number_densities(model.inner, x_state, T, mu_vec; kwargs...)

# Gap residual: delegate to NJLModel's implementation (finite-diff/AD stationarity of Ω).
gap_residual(model::LegacyNJLModel, x, T, mu_vec; kwargs...) = gap_residual(model.inner, x, T, mu_vec; kwargs...)
