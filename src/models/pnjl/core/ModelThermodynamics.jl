"""
    ModelThermodynamics

桥接层：基于新架构 `omega(model, ...)` 推导热力学量（P、ρ、s、ε）。

定位：放在 `src/pnjl/core/` 下，但不并入旧 `Thermodynamics` 模块，
避免把两套实现（旧 core vs 新 models）耦合在同一个模块里。

约定：
- 压强 P = -Ω
- ρ_i = ∂P/∂μ_i
- s = ∂P/∂T
- ε = -P + \\sum μ_i ρ_i + T s
"""
module ModelThermodynamics

using StaticArrays

# Include-once helper
const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

# constants
const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "Constants_PNJL.jl"))
IncludeOnce.include_once!(Main, :Constants_PNJL, _CONSTANTS_PATH)
using Main.Constants_PNJL: ρ0_inv_fm3

"""ρ0（单位 fm⁻³）。"""
const ρ0 = ρ0_inv_fm3

# new models loader (kept in Main namespace by design)
const _MODELS_ENTRY_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "models", "Models.jl"))

# Eager load Models at module load time to avoid world-age issues in nested AD.
if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :omega_components))
    Base.include(Main, _MODELS_ENTRY_PATH)
end

@inline function ensure_models_loaded()
    if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :omega_components))
        Base.include(Main, _MODELS_ENTRY_PATH)
    end
    return nothing
end

"""基于新模型计算压强：P = -Ω。"""
@inline function pressure_model(model, x_state, mu_vec, T_fm;
    p_num::Int,
    t_num::Int,
    xi)
    ensure_models_loaded()
    return Main.Models.model_pressure(model, x_state, mu_vec, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

"""基于新模型计算数密度向量：ρ_i = ∂P/∂μ_i。"""
function rho_model(model, x_state, mu_vec, T_fm;
    p_num::Int,
    t_num::Int,
    xi)
    ensure_models_loaded()
    return Main.Models.model_rho(model, x_state, mu_vec, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

"""基于新模型计算 (pressure, rho_norm, entropy, energy)。"""
function thermo_model(model, x_state, mu_vec, T_fm;
    p_num::Int,
    t_num::Int,
    xi)
    ensure_models_loaded()
    return Main.Models.model_thermo(model, x_state, mu_vec, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

"""基于新模型计算 (quark, antiquark) 数密度。

返回值约定同 legacy：`(quark=SVector{3}(...), antiquark=SVector{3}(...))`。
"""
@inline function number_densities_model(model, x_state, mu_vec, T_fm;
    thermal_nodes,
    xi)
    ensure_models_loaded()
    isdefined(Main, :Models) || error("Models not loaded; expected Main.Models")
    isdefined(Main.Models, :number_densities) || error("Models.number_densities is not defined")
    return Main.Models.number_densities(model, x_state, T_fm, mu_vec; thermal_nodes=thermal_nodes, xi=xi)
end

export ρ0
export ensure_models_loaded
export pressure_model, rho_model, thermo_model
export number_densities_model

end # module ModelThermodynamics

