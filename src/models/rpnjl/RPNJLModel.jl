"""RPNJLModel

阶段 6 MVP：rPNJL 的最小可运行模型壳层。

设计目标：
- 先打通 `create_model(:RPNJL)` 与 `omega`/`solve_gap` 主链；
- 复用 PNJL 的已稳定数值路径，避免一次性引入八夸克/Vandermonde 细节；
- 预留 `g1/g2/kappa` 参数位，为后续物理特化提供挂载点。
"""

using StaticArrays
using Base.MathConstants: π

const _CONSTANTS_PNJL_PATH = normpath(joinpath(@__DIR__, "..", "..", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PNJL_PATH)
end
using Main.Constants_PNJL: pnjl_constants, load_config

export RPNJLModel

const _RPNJL_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "models", "rpnjl"))

const _RPNJL_DEFAULT_EXTENSION = (
    g1_fm8=0.0,
    g2_fm8=0.0,
    kappa=0.0,
    profile="default",
    path=nothing,
)

@inline function _rpnjl_profile_default()
    return get(ENV, "RPNJL_PARAM_PROFILE", get(ENV, "PNJL_PARAM_PROFILE", "default"))
end

@inline function _rpnjl_env_flag(name::String, default::Bool=false)
    s = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    return s in ("1", "true", "yes", "y", "on")
end

const _RPNJL_DEFAULT_CONFIG = Dict{String, Any}(
    "rpnjl" => Dict{String, Any}(
        "g1_MeV_inv8" => 0.0,
        "g2_MeV_inv8" => 0.0,
        "kappa" => 0.0,
    ),
)

function _load_rpnjl_config(; profile::String, log_config::Bool=_rpnjl_env_flag("RPNJL_CONFIG_LOG", false))
    cfg = load_config(_RPNJL_CONFIG_DIR, _RPNJL_DEFAULT_CONFIG; profile=profile)
    if log_config
        @info "RPNJL config resolved" profile=profile path=cfg.path
    end
    return cfg
end

function _validate_rpnjl_extension(sec::Dict{String, Any}, path)
    g1 = Float64(get(sec, "g1_MeV_inv8", 0.0))
    g2 = Float64(get(sec, "g2_MeV_inv8", 0.0))
    kappa = Float64(get(sec, "kappa", 0.0))

    isfinite(g1) || error("invalid rpnjl config$(path === nothing ? "" : " at $(path)"): rpnjl.g1_MeV_inv8 must be finite")
    isfinite(g2) || error("invalid rpnjl config$(path === nothing ? "" : " at $(path)"): rpnjl.g2_MeV_inv8 must be finite")
    isfinite(kappa) || error("invalid rpnjl config$(path === nothing ? "" : " at $(path)"): rpnjl.kappa must be finite")
    kappa >= 0 || error("invalid rpnjl config$(path === nothing ? "" : " at $(path)"): rpnjl.kappa must be >= 0")
    return nothing
end

function _load_rpnjl_extension(; profile::String, hbarc_MeV_fm::Float64, log_config::Bool=_rpnjl_env_flag("RPNJL_CONFIG_LOG", false))
    cfg = _load_rpnjl_config(profile=profile, log_config=log_config)
    path = cfg.path
    path === nothing && return _RPNJL_DEFAULT_EXTENSION

    sec = get(cfg.config, "rpnjl", Dict{String, Any}())
    _validate_rpnjl_extension(sec, path)

    g1_MeV_inv8 = Float64(get(sec, "g1_MeV_inv8", 0.0))
    g2_MeV_inv8 = Float64(get(sec, "g2_MeV_inv8", 0.0))
    kappa = Float64(get(sec, "kappa", 0.0))

    g1_fm8 = g1_MeV_inv8 * hbarc_MeV_fm^8
    g2_fm8 = g2_MeV_inv8 * hbarc_MeV_fm^8

    return (
        g1_fm8=g1_fm8,
        g2_fm8=g2_fm8,
        kappa=kappa,
        profile=profile,
        path=path,
    )
end

function _apply_rpnjl_polyakov_overrides(base_params::PNJLCore.PNJLParams, extension_dict::Dict{String, Any})
    has_poly_override = any(k -> haskey(extension_dict, k), (
        "scheme",
        "T0_MeV",
        "a0",
        "a1",
        "a2",
        "a3",
        "b3",
        "b4",
        "fukushima_a_MeV",
        "fukushima_b_MeV3",
    ))
    has_poly_override || return base_params

    hbarc = Float64(base_params.hbarc_MeV_fm)
    scheme = Symbol(lowercase(String(get(extension_dict, "scheme", base_params.polyakov_scheme))))
    fukushima_a_inv_fm = Float64(get(extension_dict, "fukushima_a_MeV", Float64(base_params.fukushima_a_inv_fm) * hbarc)) / hbarc
    fukushima_b_inv_fm3 = Float64(get(extension_dict, "fukushima_b_MeV3", Float64(base_params.fukushima_b_inv_fm3) * hbarc^3)) / hbarc^3

    return PNJLCore.PNJLParams(
        ;
        PNJLCore.as_namedtuple(base_params)...,
        T0_inv_fm=Float64(get(extension_dict, "T0_MeV", Float64(base_params.T0_inv_fm) * hbarc)) / hbarc,
        polyakov_scheme=scheme,
        a0=Float64(get(extension_dict, "a0", base_params.a0)),
        a1=Float64(get(extension_dict, "a1", base_params.a1)),
        a2=Float64(get(extension_dict, "a2", base_params.a2)),
        a3=Float64(get(extension_dict, "a3", base_params.a3)),
        b3=Float64(get(extension_dict, "b3", base_params.b3)),
        b4=Float64(get(extension_dict, "b4", base_params.b4)),
        fukushima_a_inv_fm=fukushima_a_inv_fm,
        fukushima_b_inv_fm3=fukushima_b_inv_fm3,
    )
end

function _apply_rpnjl_model_overrides(base_params::PNJLCore.PNJLParams, extension_dict::Dict{String, Any})
    has_model_override = any(k -> haskey(extension_dict, k), (
        "Lambda_MeV",
        "m_ud0_MeV",
        "m_s0_MeV",
        "G_over_Lambda2",
        "K_over_Lambda5",
    ))
    has_model_override || return base_params

    hbarc = Float64(base_params.hbarc_MeV_fm)

    Lambda_MeV = Float64(get(extension_dict, "Lambda_MeV", Float64(base_params.Λ_inv_fm) * hbarc))
    m_ud0_MeV = Float64(get(extension_dict, "m_ud0_MeV", Float64(base_params.m_ud0_inv_fm) * hbarc))
    m_s0_MeV = Float64(get(extension_dict, "m_s0_MeV", Float64(base_params.m_s0_inv_fm) * hbarc))

    Λ_inv_fm = Lambda_MeV / hbarc
    m_ud0_inv_fm = m_ud0_MeV / hbarc
    m_s0_inv_fm = m_s0_MeV / hbarc

    G_over_Lambda2 = Float64(get(extension_dict, "G_over_Lambda2", Float64(base_params.G_fm2) * Float64(base_params.Λ_inv_fm)^2))
    K_over_Lambda5 = Float64(get(extension_dict, "K_over_Lambda5", Float64(base_params.K_fm5) * Float64(base_params.Λ_inv_fm)^5))

    G_fm2 = G_over_Lambda2 / (Λ_inv_fm^2)
    K_fm5 = K_over_Lambda5 / (Λ_inv_fm^5)

    return PNJLCore.PNJLParams(
        ;
        PNJLCore.as_namedtuple(base_params)...,
        Λ_inv_fm=Λ_inv_fm,
        m_ud0_inv_fm=m_ud0_inv_fm,
        m_s0_inv_fm=m_s0_inv_fm,
        G_fm2=G_fm2,
        K_fm5=K_fm5,
    )
end

struct RPNJLModel <: AbstractPNJLModel
    base::PNJLModel
    ext::NamedTuple
    use_extensions::Bool
end

function RPNJLModel(
    ;
    profile::String=_rpnjl_profile_default(),
    physics_profile::String=get(ENV, "PHYSICS_PARAM_PROFILE", "default"),
    use_rpnjl_extensions::Bool=true,
    log_config::Bool=_rpnjl_env_flag("RPNJL_CONFIG_LOG", false),
)
    base_params = PNJLCore.pnjl_params(; profile=profile, physics_profile=physics_profile, log_config=log_config)

    cfg_data = use_rpnjl_extensions ? _load_rpnjl_config(profile=profile, log_config=log_config) : nothing
    path = (cfg_data === nothing) ? nothing : cfg_data.path
    ext_cfg = (cfg_data === nothing || path === nothing) ? Dict{String, Any}() : get(cfg_data.config, "rpnjl", Dict{String, Any}())
    use_rpnjl_extensions && _validate_rpnjl_extension(ext_cfg, path)
    merged_params = if use_rpnjl_extensions
        model_overridden = _apply_rpnjl_model_overrides(base_params, ext_cfg)
        _apply_rpnjl_polyakov_overrides(model_overridden, ext_cfg)
    else
        base_params
    end

    extension = use_rpnjl_extensions ?
        _load_rpnjl_extension(profile=profile, hbarc_MeV_fm=Float64(base_params.hbarc_MeV_fm), log_config=log_config) :
        merge(_RPNJL_DEFAULT_EXTENSION, (profile=profile, path=nothing))
    return RPNJLModel(PNJLModel(merged_params), extension, use_rpnjl_extensions)
end

RPNJLModel() = RPNJLModel(; profile=_rpnjl_profile_default(), physics_profile=get(ENV, "PHYSICS_PARAM_PROFILE", "default"))

@inline _base_params(model::RPNJLModel) = model.base.params

function calculate_mass_vec(model::RPNJLModel, φ::SVector{3, T}; kwargs...) where {T}
    params = _base_params(model)
    m_ud0 = convert(T, Float64(params.m_ud0_inv_fm))
    m_s0 = convert(T, Float64(params.m_s0_inv_fm))
    G = convert(T, Float64(params.G_fm2))
    K = convert(T, Float64(params.K_fm5))

    g1 = convert(T, Float64(model.ext.g1_fm8))
    g2 = convert(T, Float64(model.ext.g2_fm8))

    σ2 = sum(φ .^ 2)

    base_u = m_ud0 - 2 * G * φ[1] + (K / 4) * φ[2] * φ[3]
    base_d = m_ud0 - 2 * G * φ[2] + (K / 4) * φ[1] * φ[3]
    base_s = m_s0 - 2 * G * φ[3] + (K / 4) * φ[1] * φ[2]

    δM_u = -2 * g1 * φ[1] * σ2 - 4 * g2 * φ[1]^3
    δM_d = -2 * g1 * φ[2] * σ2 - 4 * g2 * φ[2]^3
    δM_s = -2 * g1 * φ[3] * σ2 - 4 * g2 * φ[3]^3

    return SVector{3, T}(
        base_u + δM_u,
        base_d + δM_d,
        base_s + δM_s,
    )
end

function calculate_chiral(model::RPNJLModel, φ::SVector{3, T}; kwargs...) where {T}
    params = _base_params(model)
    G = convert(T, Float64(params.G_fm2))
    K = convert(T, Float64(params.K_fm5))

    g1 = convert(T, Float64(model.ext.g1_fm8))
    g2 = convert(T, Float64(model.ext.g2_fm8))

    σ2 = sum(φ .^ 2)
    base_chiral = G * σ2 - (K / 2) * prod(φ)
    octet_term = (3 * g1 / 2) * σ2^2 + 3 * g2 * sum(φ .^ 4)
    return base_chiral + octet_term
end

@inline function _rpnjl_safe_log(x; min_val=1e-16)
    min_x = one(x) * min_val
    x <= zero(x) && return log(min_x)
    return x < min_x ? log(min_x) : log(x)
end

function polyakov_potential(model::RPNJLModel, Φ, Φbar, T_fm; kwargs...)
    model.use_extensions || return polyakov_potential(model.base, Φ, Φbar, T_fm; kwargs...)

    base_value = PNJLCore.polyakov_potential(_base_params(model), Φ, Φbar, T_fm)
    kappa = model.ext.kappa
    kappa == 0 && return base_value

    TT = promote_type(typeof(Φ), typeof(Φbar), typeof(T_fm), typeof(kappa))
    ΦT = convert(TT, Φ)
    ΦbarT = convert(TT, Φbar)
    φφ = ΦT * ΦbarT
    j_poly = 1 - 6 * φφ + 4 * (ΦT^3 + ΦbarT^3) - 3 * φφ^2
    j_pref = convert(TT, 27.0 / (24.0 * π^2))
    jac = j_pref * j_poly
    return convert(TT, base_value) - convert(TT, kappa) * _rpnjl_safe_log(jac)
end

function vacuum_contribution(model::RPNJLModel, masses::SVector{3, T}; kwargs...) where {T}
    params = _base_params(model)
    Λ = convert(T, Float64(params.Λ_inv_fm))
    Nc = convert(T, Int(params.N_color))

    total = zero(T)
    @inbounds for i in 1:3
        total += PNJLCore.vacuum_integral_with_cutoff(masses[i], Λ)
    end
    return -2 * Nc * total
end

function thermal_contribution(
    ::RPNJLModel,
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

@inline number_densities(model::RPNJLModel, x_state, T, mu_vec; kwargs...) = number_densities(model.base, x_state, T, mu_vec; kwargs...)
