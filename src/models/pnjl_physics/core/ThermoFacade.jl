if !isdefined(Main, :ThermoFacade)
    @eval module ThermoFacade

"""
    ThermoFacade

薄封装：把派生热力学量 (P, ρ/ρ0, s, ε) 的计算统一到一个入口，
内部根据 `thermo_backend=:legacy|:models` 选择实现。

设计目标：
- workflow/scan/solver 不再各写一套 if/else
- 保持 legacy/models 两后端语义一致
- 保持 lazy include，避免不必要的 load-time 耦合
"""

using StaticArrays

# PNJLCore namespace (params/formulas + thermal node helpers)
const _PNJL_CORE_PATH = normpath(joinpath(@__DIR__, "..", "PNJLCore.jl"))
if !isdefined(@__MODULE__, :PNJLCore)
    include(_PNJL_CORE_PATH)
end

# Models-based thermo bridge
const _MODEL_THERMO_PATH = normpath(joinpath(@__DIR__, "ModelThermodynamics.jl"))
if !isdefined(@__MODULE__, :ModelThermodynamics)
    include(_MODEL_THERMO_PATH)
end

export calculate_thermo_backend
export calculate_rho_backend
export calculate_pressure_backend
export calculate_omega_backend
export calculate_omega_components_backend
export calculate_thermo, calculate_rho, calculate_pressure, calculate_omega, calculate_omega_components
export calculate_mass_vec_backend
export calculate_number_densities_backend
export get_models_model
export ensure_models_loaded
export rho0
export default_momentum_count, default_theta_count, cached_nodes_legacy

const _CACHED_MODELS = Dict{Symbol, Any}()
const _CACHED_MODELS_MODULE = Ref{Any}(nothing)

@inline ensure_models_loaded() = ModelThermodynamics.ensure_models_loaded()

@inline rho0() = ModelThermodynamics.ρ0
@inline default_momentum_count() = PNJLCore.DEFAULT_MOMENTUM_COUNT
@inline default_theta_count() = PNJLCore.DEFAULT_THETA_COUNT
@inline cached_nodes_legacy(p_num::Int=default_momentum_count(), t_num::Int=default_theta_count()) =
    PNJLCore.cached_nodes(p_num, t_num)

@inline function _effective_model(model, model_kind::Symbol)
    return model === nothing ? get_models_model(model_kind) : model
end

"""get_models_model(kind::Symbol=:PNJL) -> model

返回一个缓存的 models 模型对象（默认 `:PNJL`）。
"""
@inline function get_models_model(kind::Symbol=:PNJL)
    ensure_models_loaded()
    isdefined(Main, :Models) || error("Models not loaded; expected Main.Models")
    isdefined(Main.Models, :create_model) || error("Models.create_model is not defined")

    current_models = Main.Models
    if _CACHED_MODELS_MODULE[] !== current_models
        empty!(_CACHED_MODELS)
        _CACHED_MODELS_MODULE[] = current_models
    end

    return get!(_CACHED_MODELS, kind) do
        Main.Models.create_model(kind)
    end
end

"""calculate_thermo_backend(x_state, mu_vec, T_fm; ...) -> (pressure, rho_norm, entropy, energy)

统一入口：统一走 models 派生量实现。

参数：
- `model`: 当 `thermo_backend=:models` 时必须提供；若为 `nothing` 则自动用 `get_models_model(model_kind)`。
- `thermal_nodes`: legacy 路径可显式传入节点；否则使用 `PNJLCore.cached_nodes(p_num, t_num)`。
"""
function calculate_thermo_backend(
    x_state,
    mu_vec,
    T_fm;
    model=nothing,
    model_kind::Symbol=:PNJL,
    p_num::Int=24,
    t_num::Int=8,
    thermal_nodes=nothing,
    xi::Real=0.0,
)
    m = _effective_model(model, model_kind)
    return ModelThermodynamics.thermo_model(m, x_state, mu_vec, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

"""calculate_rho_backend(x_state, mu_vec, T_fm; ...) -> rho_vec

统一入口：统一走 models 数密度向量实现。
"""
function calculate_rho_backend(
    x_state,
    mu_vec,
    T_fm;
    model=nothing,
    model_kind::Symbol=:PNJL,
    p_num::Int=24,
    t_num::Int=8,
    thermal_nodes=nothing,
    xi::Real=0.0,
)
    m = _effective_model(model, model_kind)
    return ModelThermodynamics.rho_model(m, x_state, mu_vec, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

"""calculate_pressure_backend(x_state, mu_vec, T_fm; ...) -> pressure

统一入口：统一走 models 压强实现。
"""
function calculate_pressure_backend(
    x_state,
    mu_vec,
    T_fm;
    model=nothing,
    model_kind::Symbol=:PNJL,
    p_num::Int=24,
    t_num::Int=8,
    thermal_nodes=nothing,
    xi::Real=0.0,
)
    m = _effective_model(model, model_kind)
    return ModelThermodynamics.pressure_model(m, x_state, mu_vec, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

"""calculate_omega_backend(x_state, mu_vec, T_fm; ...) -> omega

统一入口：统一走 models Ω 实现。
"""
function calculate_omega_backend(
    x_state,
    mu_vec,
    T_fm;
    model=nothing,
    model_kind::Symbol=:PNJL,
    p_num::Int=24,
    t_num::Int=8,
    thermal_nodes=nothing,
    xi::Real=0.0,
)
    m = _effective_model(model, model_kind)
    return Main.Models.omega(m, x_state, T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
end

"""calculate_omega_components_backend(x_state, mu_vec, T_fm; ...) -> NamedTuple

统一入口：统一走 models 分解实现，返回 (chi, poly, vac, therm, masses, omega)。
"""
function calculate_omega_components_backend(
    x_state,
    mu_vec,
    T_fm;
    model=nothing,
    model_kind::Symbol=:PNJL,
    p_num::Int=24,
    t_num::Int=8,
    thermal_nodes=nothing,
    xi::Real=0.0,
)
    m = _effective_model(model, model_kind)
    return Main.Models.omega_components(m, x_state, T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
end

"""calculate_mass_vec_backend(x_state; ...) -> masses

统一入口：统一走 models 质量向量实现。
"""
function calculate_mass_vec_backend(
    x_state;
    model=nothing,
    model_kind::Symbol=:PNJL,
)
    m = _effective_model(model, model_kind)
    φ = SVector{3}(x_state[1], x_state[2], x_state[3])
    return Main.Models.calculate_mass_vec(m, φ)
end

"""calculate_mass_vec_backend(φ::SVector{3}; ...) -> masses

质量向量入口（φ 版本）。
"""
function calculate_mass_vec_backend(
    φ::SVector{3};
    model=nothing,
    model_kind::Symbol=:PNJL,
)
    m = _effective_model(model, model_kind)
    return Main.Models.calculate_mass_vec(m, φ)
end

"""calculate_number_densities_backend(x_state, mu_vec, T_fm; ...) -> (quark, antiquark)

统一入口：统一走 models (quark, antiquark) 数密度实现。
"""
function calculate_number_densities_backend(
    x_state,
    mu_vec,
    T_fm;
    model=nothing,
    model_kind::Symbol=:PNJL,
    p_num::Int=24,
    t_num::Int=8,
    thermal_nodes=nothing,
    xi::Real=0.0,
)
    nodes = isnothing(thermal_nodes) ? PNJLCore.cached_nodes(p_num, t_num) : thermal_nodes

    m = _effective_model(model, model_kind)
    return ModelThermodynamics.number_densities_model(m, x_state, mu_vec, T_fm; thermal_nodes=nodes, xi=xi)
end

@inline function calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes=nothing, xi::Real=0.0; kwargs...)
    return calculate_thermo_backend(x_state, mu_vec, T_fm; thermal_nodes=thermal_nodes, xi=xi, kwargs...)
end

@inline function calculate_rho(x_state, mu_vec, T_fm, thermal_nodes=nothing, xi::Real=0.0; kwargs...)
    return calculate_rho_backend(x_state, mu_vec, T_fm; thermal_nodes=thermal_nodes, xi=xi, kwargs...)
end

@inline function calculate_pressure(x_state, mu_vec, T_fm, thermal_nodes=nothing, xi::Real=0.0; kwargs...)
    return calculate_pressure_backend(x_state, mu_vec, T_fm; thermal_nodes=thermal_nodes, xi=xi, kwargs...)
end

@inline function calculate_omega(x_state, mu_vec, T_fm, thermal_nodes=nothing, xi::Real=0.0; kwargs...)
    return calculate_omega_backend(x_state, mu_vec, T_fm; thermal_nodes=thermal_nodes, xi=xi, kwargs...)
end

@inline function calculate_omega_components(x_state, mu_vec, T_fm, thermal_nodes=nothing, xi::Real=0.0; kwargs...)
    return calculate_omega_components_backend(x_state, mu_vec, T_fm; thermal_nodes=thermal_nodes, xi=xi, kwargs...)
end

end # module ThermoFacade
end

