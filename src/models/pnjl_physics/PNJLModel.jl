"""PNJLModel

新架构下的各向异性 PNJL 模型适配器（models-native）。

当前结构：
- `PNJLCore`：参数对象 + 平均场核心公式（mass/chiral/polyakov/vacuum）
- `PNJLIntegrals`：热项积分与节点缓存
- `PNJLModel`：实现 `Models` 抽象接口并组装求解流程
"""

using StaticArrays

# Include-once helper
const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _pnjl_model_cached_nodes = PNJLCore.cached_nodes
const _pnjl_model_log_sum = PNJLCore.calculate_log_sum

export PNJLModel

"""各向异性色散 PNJL 模型（参数通过 `PNJLCore.PNJLParams` 注入）。

兼容性说明：保留 `consts::NamedTuple` 字段，供历史调用方（如 `RPNJLModel`）继续读取。
"""
struct PNJLModel <: AbstractPNJLModel
    params::PNJLCore.PNJLParams
    consts::NamedTuple
end

function PNJLModel(; profile::String=get(ENV, "PNJL_PARAM_PROFILE", "default"), physics_profile::String=get(ENV, "PHYSICS_PARAM_PROFILE", "default"))
    params = PNJLCore.pnjl_params(profile=profile, physics_profile=physics_profile)
    return PNJLModel(params)
end

function PNJLModel(params::PNJLCore.PNJLParams)
    return PNJLModel(params, PNJLCore.as_namedtuple(params))
end

function PNJLModel(consts::NamedTuple)
    params = PNJLCore.pnjl_params(consts)
    return PNJLModel(params, PNJLCore.as_namedtuple(params))
end

PNJLModel() = PNJLModel(; profile=get(ENV, "PNJL_PARAM_PROFILE", "default"), physics_profile=get(ENV, "PHYSICS_PARAM_PROFILE", "default"))

# -------------------------
# Interface implementations
# -------------------------

function calculate_mass_vec(model::PNJLModel, φ::SVector{3, T}; kwargs...) where {T}
    return PNJLCore.calculate_mass_vec(model.params, φ)
end

function calculate_chiral(model::PNJLModel, φ::SVector{3, T}; kwargs...) where {T}
    return PNJLCore.chiral_potential(model.params, φ)
end

function polyakov_potential(model::PNJLModel, Φ, Φbar, T_fm; kwargs...)
    return PNJLCore.polyakov_potential(model.params, Φ, Φbar, T_fm)
end

function vacuum_contribution(model::PNJLModel, masses::SVector{3, T}; kwargs...) where {T}
    Λ = convert(T, model.params.Λ_inv_fm)
    Nc = convert(T, model.params.N_color)

    total = zero(T)
    @inbounds for i in 1:3
        total += PNJLCore.vacuum_integral_with_cutoff(masses[i], Λ)
    end
    return -2 * Nc * total
end

function thermal_contribution(
    ::PNJLModel,
    masses::SVector{3, T},
    Φ,
    Φbar,
    mu_vec,
    T_fm;
    p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT,
    t_num::Int=PNJLCore.DEFAULT_THETA_COUNT,
    xi=0.0,
    kwargs...
) where {T}
    p_mesh, cosθ_mesh, coefficients = _pnjl_model_cached_nodes(p_num, t_num)
    return _pnjl_model_log_sum(masses, p_mesh, cosθ_mesh, coefficients, Φ, Φbar, mu_vec, T_fm, xi)
end

"""返回 (quark, antiquark) 数密度（单位 fm⁻³）。

注意：此处不再委托 legacy `Thermodynamics.calculate_number_densities`，
而是在 models 侧直接按分布函数进行数值积分（便于后续彻底解耦 legacy core）。
"""
function number_densities(
    model::PNJLModel,
    x_state,
    T_fm,
    mu_vec;
    thermal_nodes=nothing,
    p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT,
    t_num::Int=PNJLCore.DEFAULT_THETA_COUNT,
    xi=0.0,
    kwargs...
)
    st = x_state isa MeanFieldState ? x_state : MeanFieldState(x_state)
    Φ = Float64(st.Phi)
    Φbar = Float64(st.PhiBar)
    T_val = Float64(T_fm)
    xi_val = Float64(xi)

    φ = SVector{3, Float64}(Float64(st.phi[1]), Float64(st.phi[2]), Float64(st.phi[3]))
    masses = calculate_mass_vec(model, φ)

    nodes = isnothing(thermal_nodes) ? _pnjl_model_cached_nodes(p_num, t_num) : thermal_nodes
    thermal_p_mesh, cosθ_mesh, thermal_coefficients = nodes
    pref = 2 * model.params.N_color

    acc_q = MVector{3, Float64}(0.0, 0.0, 0.0)
    acc_aq = MVector{3, Float64}(0.0, 0.0, 0.0)

    mu_vec = normalize_mu_vec(mu_vec)

    @inbounds for i in 1:3
        mass_i = Float64(masses[i])
        mu_i = Float64(mu_vec[i])
        total_q = 0.0
        total_aq = 0.0
        for idx in eachindex(thermal_p_mesh)
            p = Float64(thermal_p_mesh[idx])
            cosθ = Float64(cosθ_mesh[idx])
            w = Float64(thermal_coefficients[idx])
            total_q += w * pref * pnjl_quark_distribution_aniso(p, mass_i, mu_i, T_val, Φ, Φbar, xi_val, cosθ)
            total_aq += w * pref * pnjl_antiquark_distribution_aniso(p, mass_i, mu_i, T_val, Φ, Φbar, xi_val, cosθ)
        end
        acc_q[i] = total_q
        acc_aq[i] = total_aq
    end

    return (quark=SVector{3}(acc_q), antiquark=SVector{3}(acc_aq))
end

"""solve_gap(::PNJLModel, T, mu_vec; kwargs...) -> MeanFieldState

统一通过 Models 求解链路完成平衡态求解。
"""
function solve_gap(
    model::PNJLModel,
    T_fm,
    mu_vec;
    solver_backend::Symbol=:legacy,
    fallback_legacy_on_failure::Bool=true,
    xi::Real=0.0,
    p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT,
    t_num::Int=PNJLCore.DEFAULT_THETA_COUNT,
    kwargs...
)
    legacy_kwargs = Dict{Symbol,Any}(pairs(kwargs))
    for key in (:solver, :initial_guess, :residual_norm_max)
        pop!(legacy_kwargs, key, nothing)
    end

    effective_backend = solver_backend === :auto ? :models : solver_backend
    if effective_backend !== :models && effective_backend !== :legacy
        throw(ArgumentError("unknown solver_backend=$solver_backend (expected :legacy, :models or :auto)"))
    end

    try
        return invoke(
            solve_gap,
            Tuple{AbstractPNJLModel, Any, Any},
            model,
            T_fm,
            mu_vec;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            kwargs...,
        )
    catch err
        if !fallback_legacy_on_failure
            rethrow(err)
        end
        return invoke(
            solve_gap,
            Tuple{AbstractPNJLModel, Any, Any},
            model,
            T_fm,
            mu_vec;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            legacy_kwargs...,
        )
    end
end
