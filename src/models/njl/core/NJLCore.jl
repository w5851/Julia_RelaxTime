"""
    NJLCore

三味 NJL 模型的可复用 core（参数加载 + 单位换算 + 平均场核心公式）。

来源：docs/reference/formula/models/njl/NJL_core.md

单位约定：
- 配置文件使用 MeV（如 `m_ud0_MeV`, `Lambda_MeV`），`hbarc` 默认为 MeV·fm。
- 代码内部默认使用 fm⁻¹ / fmⁿ：
  - 质量、温度、化学势：fm⁻¹
  - 截断：fm⁻¹
  - G：fm²
  - K：fm⁵

核心公式：
- 组分质量（式 2-13）：
  M_i = m_i - 4G φ_i + 2K φ_j φ_k,  i≠j≠k
- 平均场凝聚项（对应式 2-12 的常见写法）：
    χ(φ) = 2G \\sum_i φ_i^2 - 4K φ_u φ_d φ_s

注意：这里的 φ 采用与项目现有 PNJL 路径一致的“凝聚变量”记号。
"""
module NJLCore

using StaticArrays

# Shared config loader (reusable across models)
const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    include(_CONFIG_LOADER_PATH)
end
using .ConfigLoader: deep_merge, load_config

export NJLParams
export load_njl_config, njl_params
export mevinv_to_fminv, mevinv2_to_fm2, mevinv5_to_fm5
export calculate_mass_vec, chiral_potential

"""MeV → fm⁻¹：x_fm⁻¹ = x_MeV / ħc"""
@inline mevinv_to_fminv(x_MeV::Real, hbarc_MeV_fm::Real) = float(x_MeV) / float(hbarc_MeV_fm)

"""MeV⁻² → fm²：x_fm² = x_MeV⁻² * ħc²"""
@inline mevinv2_to_fm2(x_MeV_inv2::Real, hbarc_MeV_fm::Real) = float(x_MeV_inv2) * float(hbarc_MeV_fm)^2

"""MeV⁻⁵ → fm⁵：x_fm⁵ = x_MeV⁻⁵ * ħc⁵"""
@inline mevinv5_to_fm5(x_MeV_inv5::Real, hbarc_MeV_fm::Real) = float(x_MeV_inv5) * float(hbarc_MeV_fm)^5

"""NJL 参数（内部单位以 fm 体系存储）。"""
Base.@kwdef struct NJLParams
    # meta
    label::String = "njl"

    # dof
    N_color::Int = 3
    N_flavor::Int = 3

    # scales
    hbarc_MeV_fm::Float64 = 197.3269804
    rho0_fm3::Float64 = 0.16

    # model params (internal fm units)
    Lambda_inv_fm::Float64
    G_fm2::Float64
    K_fm5::Float64
    m_ud0_inv_fm::Float64
    m_s0_inv_fm::Float64
end

const DEFAULT_PROFILE = "default"

const PHYSICS_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "physics"))
const NJL_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "models", "njl"))

const DEFAULT_PHYSICS_CONFIG = Dict{String, Any}(
    "physical" => Dict(
        "hbarc" => 197.3269804,
        "alpha_em" => 0.0072973525664,
    ),
)

const DEFAULT_NJL_MODEL_CONFIG = Dict{String, Any}(
    "model" => Dict(
        "label" => "njl-standard",
        "version" => "1.0",
        "N_color" => 3,
        "N_flavor" => 3,
        "rho0_fm3" => 0.16,
        "Lambda_MeV" => 602.3,
        "G_over_Lambda2" => 1.835,
        "K_over_Lambda5" => 12.36,
        "m_ud0_MeV" => 5.5,
        "m_s0_MeV" => 140.7,
    ),
)

"""读取配置并返回解析后的 Dict。

使用方案B：
- physics: `config/physics/<profile>.toml`
- model: `config/models/njl/<profile>.toml`
"""
function load_njl_config(; profile::String=get(ENV, "NJL_PARAM_PROFILE", DEFAULT_PROFILE))
    physics_profile = get(ENV, "PHYSICS_PARAM_PROFILE", DEFAULT_PROFILE)

    physics_data = load_config(PHYSICS_CONFIG_DIR, DEFAULT_PHYSICS_CONFIG; profile=physics_profile)
    shared_model = get(physics_data.config, "model_shared", Dict{String, Any}())
    inherited_model_cfg = isempty(shared_model) ? Dict{String, Any}[] : [Dict("model" => shared_model)]

    model_data = load_config(
        NJL_CONFIG_DIR,
        DEFAULT_NJL_MODEL_CONFIG;
        profile=profile,
        inherited_configs=inherited_model_cfg,
    )

    return deep_merge(physics_data.config, model_data.config)
end

"""从配置文件构造 `NJLParams`（默认 profile=default）。"""
function njl_params(; profile::String=get(ENV, "NJL_PARAM_PROFILE", "default"))
    cfg = load_njl_config(profile=profile)

    physical = get(cfg, "physical", Dict{String, Any}())
    model = get(cfg, "model", Dict{String, Any}())

    hbarc = Float64(get(physical, "hbarc", 197.3269804))
    label = String(get(model, "label", "njl"))

    N_color = Int(get(model, "N_color", 3))
    N_flavor = Int(get(model, "N_flavor", 3))
    rho0 = Float64(get(model, "rho0_fm3", 0.16))

    Lambda_MeV = Float64(get(model, "Lambda_MeV", 602.3))
    G_over_Lambda2 = Float64(get(model, "G_over_Lambda2", 1.835))
    K_over_Lambda5 = Float64(get(model, "K_over_Lambda5", 12.36))

    m_ud0_MeV = Float64(get(model, "m_ud0_MeV", 5.5))
    m_s0_MeV = Float64(get(model, "m_s0_MeV", 140.7))

    # unit conversions
    Lambda_inv_fm = mevinv_to_fminv(Lambda_MeV, hbarc)

    # GΛ², KΛ⁵ are dimensionless; first recover MeV^{-2/-5}, then convert to fm^{2/5}
    G_MeV_inv2 = G_over_Lambda2 / (Lambda_MeV^2)
    K_MeV_inv5 = K_over_Lambda5 / (Lambda_MeV^5)

    G_fm2 = mevinv2_to_fm2(G_MeV_inv2, hbarc)
    K_fm5 = mevinv5_to_fm5(K_MeV_inv5, hbarc)

    return NJLParams(
        label=label,
        N_color=N_color,
        N_flavor=N_flavor,
        hbarc_MeV_fm=hbarc,
        rho0_fm3=rho0,
        Lambda_inv_fm=Float64(Lambda_inv_fm),
        G_fm2=Float64(G_fm2),
        K_fm5=Float64(K_fm5),
        m_ud0_inv_fm=mevinv_to_fminv(m_ud0_MeV, hbarc),
        m_s0_inv_fm=mevinv_to_fminv(m_s0_MeV, hbarc),
    )
end

"""组分夸克质量（式 2-13），返回 (Mu, Md, Ms)。

输入 φ = (φu, φd, φs)。
"""
function calculate_mass_vec(p::NJLParams, φ::SVector{3, T}) where {T}
    if p.N_flavor != 3
        throw(ArgumentError("calculate_mass_vec currently supports N_flavor=3, got $(p.N_flavor)"))
    end

    φu, φd, φs = φ
    G = convert(T, p.G_fm2)
    K = convert(T, p.K_fm5)

    m_ud0 = convert(T, p.m_ud0_inv_fm)
    m_s0 = convert(T, p.m_s0_inv_fm)

    Mu = m_ud0 - 4 * G * φu + 2 * K * φd * φs
    Md = m_ud0 - 4 * G * φd + 2 * K * φu * φs
    Ms = m_s0 - 4 * G * φs + 2 * K * φu * φd

    return SVector{3, T}(Mu, Md, Ms)
end

"""平均场凝聚项 χ(φ) = 2G Σ φ_i² - 4K φu φd φs。"""
function chiral_potential(p::NJLParams, φ::SVector{3, T}) where {T}
    G = convert(T, p.G_fm2)
    K = convert(T, p.K_fm5)
    return 2 * G * sum(φ .^ 2) - 4 * K * prod(φ)
end

end # module NJLCore
