"""
# Constants_PNJL.jl

集中维护 RelaxTime 项目中使用的共用常量。
"""
module Constants_PNJL

using TOML

# Shared config loader (reusable across models)
const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    include(_CONFIG_LOADER_PATH)
end
using .ConfigLoader: deep_merge, load_config

# Transport/scattering process constants (decoupled from PNJL model parameters)
const _TRANSPORT_CONSTANTS_PATH = normpath(joinpath(@__DIR__, "constants", "TransportConstants.jl"))
if !isdefined(@__MODULE__, :TransportConstants)
    include(_TRANSPORT_CONSTANTS_PATH)
end
using .TransportConstants: SCATTERING_MESON_MAP, SCATTERING_PROCESS_KEYS

export ħc_MeV_fm, α
export N_color, N_flavor, ρ0_inv_fm3, m_ud0_inv_fm, m_s0_inv_fm, Λ_inv_fm, G_fm2, K_fm5, thermal_p_max_inv_fm
export T0_inv_fm, polyakov_scheme, a0, a1, a2, a3, b3, b4, fukushima_a_inv_fm, fukushima_b_inv_fm3
export λ₀, λ₈, ψ_u, ψ_d, ψ_s, ψbar_u, ψbar_d, ψbar_s
export PNJL_PROFILE, PNJL_CONFIG_PATH, load_pnjl_config
export pnjl_constants
export SCATTERING_MESON_MAP, SCATTERING_PROCESS_KEYS

const PHYSICS_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "config", "physics"))
const PNJL_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "config", "models", "pnjl"))
const DEFAULT_PROFILE = "default"

const DEFAULT_PHYSICS_CONFIG = Dict{String, Any}(
    "physical" => Dict(
        "hbarc" => 197.3269804,
        "alpha_em" => 0.0072973525664,
    ),
)

const DEFAULT_PNJL_MODEL_CONFIG = Dict{String, Any}(
    "model" => Dict(
        "N_color" => 3,
        "N_flavor" => 3,
        "rho0_fm3" => 0.16,
        "Lambda_MeV" => 602.3,
        "G_over_Lambda2" => 1.835,
        "K_over_Lambda5" => 12.36,
        "m_ud0_MeV" => 5.5,
        "m_s0_MeV" => 140.7,
        "thermal_p_max_MeV" => 1973.269804,
    ),
    "polyakov" => Dict(
        "scheme" => "log",
        "T0_MeV" => 210.0,
        "a0" => 3.51,
        "a1" => -2.47,
        "a2" => 15.2,
        "a3" => 0.0,
        "b3" => -1.75,
        "b4" => 7.555,
        "fukushima_a_MeV" => 664.0,
        "fukushima_b_MeV3" => 7529536.0,
    ),
)

@inline function _env_flag(name::String, default::Bool=false)
    s = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    return s in ("1", "true", "yes", "y", "on")
end

@inline function _require_finite_positive(value::Float64, label::String)
    isfinite(value) || error("invalid config: $(label) must be finite, got $(value)")
    value > 0 || error("invalid config: $(label) must be > 0, got $(value)")
    return value
end

@inline function _require_finite_nonnegative(value::Float64, label::String)
    isfinite(value) || error("invalid config: $(label) must be finite, got $(value)")
    value >= 0 || error("invalid config: $(label) must be >= 0, got $(value)")
    return value
end

@inline function _require_finite_nonzero(value::Float64, label::String)
    isfinite(value) || error("invalid config: $(label) must be finite, got $(value)")
    value != 0 || error("invalid config: $(label) must be != 0, got $(value)")
    return value
end

function _validate_pnjl_critical_params(
    ;
    hbarc_MeV_fm::Float64,
    alpha_em::Float64,
    N_color::Int,
    N_flavor::Int,
    rho0_fm3::Float64,
    Lambda_MeV::Float64,
    m_ud0_MeV::Float64,
    m_s0_MeV::Float64,
    G_over_Lambda2::Float64,
    K_over_Lambda5::Float64,
    thermal_p_max_MeV::Float64,
    T0_MeV::Float64,
)
    _require_finite_positive(hbarc_MeV_fm, "physical.hbarc")
    _require_finite_positive(alpha_em, "physical.alpha_em")
    N_color >= 1 || error("invalid config: model.N_color must be >= 1, got $(N_color)")
    N_flavor >= 1 || error("invalid config: model.N_flavor must be >= 1, got $(N_flavor)")
    _require_finite_positive(rho0_fm3, "model.rho0_fm3")
    _require_finite_positive(Lambda_MeV, "model.Lambda_MeV")
    _require_finite_nonnegative(m_ud0_MeV, "model.m_ud0_MeV")
    _require_finite_nonnegative(m_s0_MeV, "model.m_s0_MeV")
    _require_finite_positive(G_over_Lambda2, "model.G_over_Lambda2")
    _require_finite_nonzero(K_over_Lambda5, "model.K_over_Lambda5")
    _require_finite_positive(thermal_p_max_MeV, "model.thermal_p_max_MeV")
    _require_finite_positive(T0_MeV, "polyakov.T0_MeV")
    return nothing
end

function load_pnjl_config(
    ;
    profile::String=get(ENV, "PNJL_PARAM_PROFILE", DEFAULT_PROFILE),
    physics_profile::String=get(ENV, "PHYSICS_PARAM_PROFILE", DEFAULT_PROFILE),
    log_config::Bool=_env_flag("PNJL_CONFIG_LOG", false),
)

    physics_data = load_config(PHYSICS_CONFIG_DIR, DEFAULT_PHYSICS_CONFIG; profile=physics_profile)
    shared_model = get(physics_data.config, "model_shared", Dict{String, Any}())
    inherited_model_cfg = isempty(shared_model) ? Dict{String, Any}[] : [Dict("model" => shared_model)]

    model_data = load_config(
        PNJL_CONFIG_DIR,
        DEFAULT_PNJL_MODEL_CONFIG;
        profile=profile,
        inherited_configs=inherited_model_cfg,
    )

    merged = deep_merge(physics_data.config, model_data.config)

    if log_config
        @info "PNJL config resolved" profile=profile physics_profile=physics_profile model_path=model_data.path physics_path=physics_data.path
    end

    return (config=merged, profile=profile, path=model_data.path)
end

"""pnjl_constants(; profile, physics_profile) -> NamedTuple

返回一组 PNJL 相关“常量/参数”（单位已换算到项目内部 fm 体系），但**不依赖模块加载时的 `const` 缓存**。

用途：
- 单测中用 `withenv` 或直接传参验证 profile 切换行为
- 新代码逐步从全局常量迁移到 params 注入时的过渡入口
"""
function pnjl_constants(
    ;
    profile::String=get(ENV, "PNJL_PARAM_PROFILE", DEFAULT_PROFILE),
    physics_profile::String=get(ENV, "PHYSICS_PARAM_PROFILE", DEFAULT_PROFILE),
    validate::Bool=true,
    log_config::Bool=_env_flag("PNJL_CONFIG_LOG", false),
)
    data = load_pnjl_config(profile=profile, physics_profile=physics_profile, log_config=log_config)
    cfg = data.config

    physical_cfg = get(cfg, "physical", Dict{String, Any}())
    model_cfg = get(cfg, "model", Dict{String, Any}())
    polyakov_cfg = get(cfg, "polyakov", Dict{String, Any}())

    hbarc_MeV_fm = Float64(get(physical_cfg, "hbarc", 197.3269804))
    alpha_em = Float64(get(physical_cfg, "alpha_em", 0.0072973525664))

    N_color = Int(get(model_cfg, "N_color", 3))
    N_flavor = Int(get(model_cfg, "N_flavor", 3))
    rho0_fm3 = Float64(get(model_cfg, "rho0_fm3", 0.16))

    Lambda_MeV = Float64(get(model_cfg, "Lambda_MeV", 602.3))
    Λ_inv_fm = Lambda_MeV / hbarc_MeV_fm

    m_ud0_MeV = Float64(get(model_cfg, "m_ud0_MeV", 5.5))
    m_s0_MeV = Float64(get(model_cfg, "m_s0_MeV", 140.7))
    m_ud0_inv_fm = m_ud0_MeV / hbarc_MeV_fm
    m_s0_inv_fm = m_s0_MeV / hbarc_MeV_fm

    G_over_Lambda2 = Float64(get(model_cfg, "G_over_Lambda2", 1.835))
    K_over_Lambda5 = Float64(get(model_cfg, "K_over_Lambda5", 12.36))
    G_fm2 = G_over_Lambda2 / (Λ_inv_fm^2)
    K_fm5 = K_over_Lambda5 / (Λ_inv_fm^5)
    thermal_p_max_MeV = Float64(get(model_cfg, "thermal_p_max_MeV", 1973.269804))
    thermal_p_max_inv_fm = thermal_p_max_MeV / hbarc_MeV_fm

    T0_MeV = Float64(get(polyakov_cfg, "T0_MeV", 210.0))
    T0_inv_fm = T0_MeV / hbarc_MeV_fm
    polyakov_scheme = Symbol(lowercase(String(get(polyakov_cfg, "scheme", "log"))))

    a0 = Float64(get(polyakov_cfg, "a0", 3.51))
    a1 = Float64(get(polyakov_cfg, "a1", -2.47))
    a2 = Float64(get(polyakov_cfg, "a2", 15.2))
    a3 = Float64(get(polyakov_cfg, "a3", 0.0))
    b3 = Float64(get(polyakov_cfg, "b3", -1.75))
    b4 = Float64(get(polyakov_cfg, "b4", 7.555))
    fukushima_a_MeV = Float64(get(polyakov_cfg, "fukushima_a_MeV", 664.0))
    fukushima_b_MeV3 = Float64(get(polyakov_cfg, "fukushima_b_MeV3", 7529536.0))
    fukushima_a_inv_fm = fukushima_a_MeV / hbarc_MeV_fm
    fukushima_b_inv_fm3 = fukushima_b_MeV3 / (hbarc_MeV_fm^3)

    if validate
        _validate_pnjl_critical_params(
            hbarc_MeV_fm=hbarc_MeV_fm,
            alpha_em=alpha_em,
            N_color=N_color,
            N_flavor=N_flavor,
            rho0_fm3=rho0_fm3,
            Lambda_MeV=Lambda_MeV,
            m_ud0_MeV=m_ud0_MeV,
            m_s0_MeV=m_s0_MeV,
            G_over_Lambda2=G_over_Lambda2,
            K_over_Lambda5=K_over_Lambda5,
            thermal_p_max_MeV=thermal_p_max_MeV,
            T0_MeV=T0_MeV,
        )
    end

    return (
        profile=profile,
        physics_profile=physics_profile,
        path=data.path,
        hbarc_MeV_fm=hbarc_MeV_fm,
        alpha_em=alpha_em,
        N_color=N_color,
        N_flavor=N_flavor,
        rho0_fm3=rho0_fm3,
        Λ_inv_fm=Λ_inv_fm,
        m_ud0_inv_fm=m_ud0_inv_fm,
        m_s0_inv_fm=m_s0_inv_fm,
        G_fm2=G_fm2,
        K_fm5=K_fm5,
        thermal_p_max_inv_fm=thermal_p_max_inv_fm,
        T0_inv_fm=T0_inv_fm,
        polyakov_scheme=polyakov_scheme,
        a0=a0,
        a1=a1,
        a2=a2,
        a3=a3,
        b3=b3,
        b4=b4,
        fukushima_a_inv_fm=fukushima_a_inv_fm,
        fukushima_b_inv_fm3=fukushima_b_inv_fm3,
    )
end

const PNJL_CONFIG_DATA = load_pnjl_config()
const PNJL_PROFILE = PNJL_CONFIG_DATA.profile
const PNJL_CONFIG_PATH = PNJL_CONFIG_DATA.path
const PNJL_CONFIG = PNJL_CONFIG_DATA.config

# 基本物理常量-------------------------------------
const physical_cfg = PNJL_CONFIG["physical"]
const ħc_MeV_fm = Float64(physical_cfg["hbarc"])  # MeV·fm
const α::Float64 = Float64(physical_cfg["alpha_em"])  # 精细结构常数

# PNJL模型参数-------------------------------------
const model_cfg = PNJL_CONFIG["model"]
const N_color = Int(model_cfg["N_color"])  # 夸克颜色数
const N_flavor = Int(model_cfg["N_flavor"])  # 夸克味道数
const ρ0_inv_fm3 = Float64(model_cfg["rho0_fm3"])  # 核子数密度, 单位fm⁻³
const Λ_inv_fm = Float64(model_cfg["Lambda_MeV"]) / ħc_MeV_fm  # 截断参数
const m_ud0_inv_fm = Float64(model_cfg["m_ud0_MeV"]) / ħc_MeV_fm  # u,d夸克裸质量
const m_s0_inv_fm = Float64(model_cfg["m_s0_MeV"]) / ħc_MeV_fm  # s夸克裸质量
const G_fm2 = Float64(model_cfg["G_over_Lambda2"]) / Λ_inv_fm^2  # NJL四夸克耦合
const K_fm5 = Float64(model_cfg["K_over_Lambda5"]) / Λ_inv_fm^5  # NJL六夸克耦合
const thermal_p_max_inv_fm = Float64(get(model_cfg, "thermal_p_max_MeV", 1973.269804)) / ħc_MeV_fm  # 热项动量积分上限

# Polyakov环有效势参数-------------------------------------
const polyakov_cfg = PNJL_CONFIG["polyakov"]
const T0_inv_fm = Float64(polyakov_cfg["T0_MeV"]) / ħc_MeV_fm  # Polyakov有效势参数
const polyakov_scheme = Symbol(lowercase(String(get(polyakov_cfg, "scheme", "log"))))
const a0 = Float64(polyakov_cfg["a0"])
const a1 = Float64(polyakov_cfg["a1"])
const a2 = Float64(polyakov_cfg["a2"])
const a3 = Float64(get(polyakov_cfg, "a3", 0.0))
const b3 = Float64(polyakov_cfg["b3"])
const b4 = Float64(polyakov_cfg["b4"])
const fukushima_a_inv_fm = Float64(get(polyakov_cfg, "fukushima_a_MeV", 664.0)) / ħc_MeV_fm
const fukushima_b_inv_fm3 = Float64(get(polyakov_cfg, "fukushima_b_MeV3", 7529536.0)) / (ħc_MeV_fm^3)

# Gell-Mann矩阵(SU(3)味对称性)-------------------------------------
# λ₀: 味单位矩阵(归一化)
const λ₀ = [
    1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
] * sqrt(2/3)

# λ₈: Gell-Mann第8矩阵(u,d对称,s不同)
const λ₈ = [
    1.0  0.0  0.0
    0.0  1.0  0.0
    0.0  0.0 -2.0
] / sqrt(3)

# 夸克味波函数(列向量)-------------------------------------
# ψ_u: u夸克波函数
const ψ_u = [1.0, 0.0, 0.0]

# ψ_d: d夸克波函数
const ψ_d = [0.0, 1.0, 0.0]

# ψ_s: s夸克波函数
const ψ_s = [0.0, 0.0, 1.0]

# 夸克味波函数(行向量/1×3矩阵)-------------------------------------
# ψbar_u: u夸克共轭波函数
const ψbar_u = [1.0 0.0 0.0]

# ψbar_d: d夸克共轭波函数
const ψbar_d = [0.0 1.0 0.0]

# ψbar_s: s夸克共轭波函数
const ψbar_s = [0.0 0.0 1.0]

"""散射过程与介子映射常量已迁移至 `src/constants/TransportConstants.jl`。

`Constants_PNJL` 保留同名导出作为兼容层。
"""

end # module Constants_PNJL
