"""
    NJL2Core

两味 NJL 模型 core（参数加载 + 单位换算 + 平均场核心公式）。

公式依据：docs/reference/formula/models/njl/NJL_core.md（两味补充）
- 质量方程：M_i = m_i - 4G φ_i, i∈{u,d}
- 手征项：χ(φ)=2G(φ_u^2+φ_d^2)

单位约定与 NJLCore 一致：内部使用 fm 体系。
"""
module NJL2Core

using StaticArrays

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    include(_CONFIG_LOADER_PATH)
end
using .ConfigLoader: deep_merge, load_config

export NJL2Params
export load_njl2_config, njl2_params
export mevinv_to_fminv, mevinv2_to_fm2
export calculate_mass_vec, chiral_potential

@inline mevinv_to_fminv(x_MeV::Real, hbarc_MeV_fm::Real) = float(x_MeV) / float(hbarc_MeV_fm)
@inline mevinv2_to_fm2(x_MeV_inv2::Real, hbarc_MeV_fm::Real) = float(x_MeV_inv2) * float(hbarc_MeV_fm)^2

Base.@kwdef struct NJL2Params
    label::String = "njl-2f"
    N_color::Int = 3
    N_flavor::Int = 2
    hbarc_MeV_fm::Float64 = 197.327
    rho0_fm3::Float64 = 0.16

    Lambda_inv_fm::Float64
    G_fm2::Float64
    m_ud0_inv_fm::Float64
end

const DEFAULT_PROFILE = "default"

const PHYSICS_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "physics"))
const NJL2_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "config", "models", "njl2"))

const DEFAULT_PHYSICS_CONFIG = Dict{String, Any}(
    "physical" => Dict(
        "hbarc" => 197.327,
        "alpha_em" => 1.0 / 137.035999084,
    ),
)

const DEFAULT_NJL2_MODEL_CONFIG = Dict{String, Any}(
    "model" => Dict(
        "label" => "njl-2f",
        "version" => "1.0",
        "N_color" => 3,
        "N_flavor" => 2,
        "rho0_fm3" => 0.16,
        "Lambda_MeV" => 631.0,
        "G_over_Lambda2" => 2.14,
        "m_ud0_MeV" => 5.5,
    ),
)

function load_njl2_config(; profile::String=get(ENV, "NJL2_PARAM_PROFILE", DEFAULT_PROFILE))
    physics_profile = get(ENV, "PHYSICS_PARAM_PROFILE", DEFAULT_PROFILE)

    physics_data = load_config(PHYSICS_CONFIG_DIR, DEFAULT_PHYSICS_CONFIG; profile=physics_profile)
    model_data = load_config(NJL2_CONFIG_DIR, DEFAULT_NJL2_MODEL_CONFIG; profile=profile)

    return deep_merge(physics_data.config, model_data.config)
end

function njl2_params(; profile::String=get(ENV, "NJL2_PARAM_PROFILE", DEFAULT_PROFILE))
    cfg = load_njl2_config(profile=profile)

    physical = get(cfg, "physical", Dict{String, Any}())
    model = get(cfg, "model", Dict{String, Any}())

    hbarc = Float64(get(physical, "hbarc", 197.327))
    label = String(get(model, "label", "njl-2f"))

    N_color = Int(get(model, "N_color", 3))
    N_flavor = Int(get(model, "N_flavor", 2))
    rho0 = Float64(get(model, "rho0_fm3", 0.16))

    Lambda_MeV = Float64(get(model, "Lambda_MeV", 631.0))
    G_over_Lambda2 = Float64(get(model, "G_over_Lambda2", 2.14))
    m_ud0_MeV = Float64(get(model, "m_ud0_MeV", 5.5))

    Lambda_inv_fm = mevinv_to_fminv(Lambda_MeV, hbarc)
    G_MeV_inv2 = G_over_Lambda2 / (Lambda_MeV^2)
    G_fm2 = mevinv2_to_fm2(G_MeV_inv2, hbarc)

    return NJL2Params(
        label=label,
        N_color=N_color,
        N_flavor=N_flavor,
        hbarc_MeV_fm=hbarc,
        rho0_fm3=rho0,
        Lambda_inv_fm=Float64(Lambda_inv_fm),
        G_fm2=Float64(G_fm2),
        m_ud0_inv_fm=mevinv_to_fminv(m_ud0_MeV, hbarc),
    )
end

"""两味组分质量，返回 (Mu, Md, Ms_stub=0)。"""
function calculate_mass_vec(p::NJL2Params, φ::SVector{3, T}) where {T}
    if p.N_flavor != 2
        throw(ArgumentError("calculate_mass_vec expects N_flavor=2, got $(p.N_flavor)"))
    end

    φu, φd = φ[1], φ[2]
    G = convert(T, p.G_fm2)
    m0 = convert(T, p.m_ud0_inv_fm)

    Mu = m0 - 4 * G * φu
    Md = m0 - 4 * G * φd

    return SVector{3, T}(Mu, Md, zero(T))
end

"""两味平均场凝聚项 χ(φ)=2G(φ_u²+φ_d²)。"""
function chiral_potential(p::NJL2Params, φ::SVector{3, T}) where {T}
    G = convert(T, p.G_fm2)
    return 2 * G * (φ[1]^2 + φ[2]^2)
end

end # module NJL2Core
