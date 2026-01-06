"""
# Constants_PNJL.jl

集中维护 RelaxTime 项目中使用的共用常量。
"""
module Constants_PNJL

using TOML

export ħc_MeV_fm, α
export N_color, N_flavor, ρ0_inv_fm3, m_ud0_inv_fm, m_s0_inv_fm, Λ_inv_fm, G_fm2, K_fm5
export T0_inv_fm, a0, a1, a2, b3, b4
export λ₀, λ₈, ψ_u, ψ_d, ψ_s, ψbar_u, ψbar_d, ψbar_s
export PNJL_PROFILE, PNJL_CONFIG_PATH, load_pnjl_config
export SCATTERING_MESON_MAP, SCATTERING_PROCESS_KEYS

const CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "config", "pnjl"))
const DEFAULT_PROFILE = "default"

const DEFAULT_CONFIG = Dict{String, Any}(
    "physical" => Dict(
        "hbarc" => 197.327,
        "alpha_em" => 1.0 / 137.035999084,
    ),
    "model" => Dict(
        "N_color" => 3,
        "N_flavor" => 3,
        "rho0_fm3" => 0.16,
        "Lambda_MeV" => 602.3,
        "G_over_Lambda2" => 1.835,
        "K_over_Lambda5" => 12.36,
        "m_ud0_MeV" => 5.5,
        "m_s0_MeV" => 140.7,
    ),
    "polyakov" => Dict(
        "T0_MeV" => 210.0,
        "a0" => 3.51,
        "a1" => -2.47,
        "a2" => 15.2,
        "b3" => -1.75,
        "b4" => 7.555,
    ),
)

function deep_merge(base::Dict{String, Any}, override::Dict{String, Any})
    result = deepcopy(base)
    for (k, v) in override
        if v isa Dict{String, Any} && haskey(result, k) && result[k] isa Dict{String, Any}
            result[k] = deep_merge(result[k]::Dict{String, Any}, v)
        else
            result[k] = v
        end
    end
    return result
end

function load_pnjl_config(; profile::String=get(ENV, "PNJL_PARAM_PROFILE", DEFAULT_PROFILE))
    path = joinpath(CONFIG_DIR, string(profile, ".toml"))
    cfg = deepcopy(DEFAULT_CONFIG)
    if isfile(path)
        try
            parsed = TOML.parsefile(path)
            cfg = deep_merge(cfg, parsed)
        catch err
            @warn "解析 PNJL 配置失败，使用内置默认值" profile path exception = (err, catch_backtrace())
        end
    else
        if profile != DEFAULT_PROFILE
            @warn "未找到指定 PNJL 配置文件，退回默认参数" profile path
        end
    end
    return (config=cfg, profile=profile, path=isfile(path) ? path : nothing)
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

# Polyakov环有效势参数-------------------------------------
const polyakov_cfg = PNJL_CONFIG["polyakov"]
const T0_inv_fm = Float64(polyakov_cfg["T0_MeV"]) / ħc_MeV_fm  # Polyakov有效势参数
const a0 = Float64(polyakov_cfg["a0"])
const a1 = Float64(polyakov_cfg["a1"])
const a2 = Float64(polyakov_cfg["a2"])
const b3 = Float64(polyakov_cfg["b3"])
const b4 = Float64(polyakov_cfg["b4"])

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

# 散射过程到介子种类的映射表-------------------------------------
"""
散射过程介子映射表

映射关系：散射过程 → 散射类型 + 各散射道的介子列表

# 数据结构
```julia
Dict(
    :process_name => Dict(
        :type => :qq 或 :qqbar,  # 散射类型（qq有t/u道，qqbar有t/s道）
        :channels => Dict(
            :t => Dict(:simple => [...], :mixed_P => true/false, :mixed_S => true/false),
            :u => Dict(...),  # 仅qq散射有u道
            :s => Dict(...)   # 仅qqbar散射有s道
        )
    )
)
```

# 介子类型说明
- `:simple` 列表：一般介子（:pi, :K, :sigma_pi, :sigma_K）
- `:mixed_P => true`：存在赝标量混合介子（η/η'整体）
- `:mixed_S => true`：存在标量混合介子（σ/σ'整体）

# 参考文献
doc/formula/散射过程所有可能.md
"""
const SCATTERING_MESON_MAP = Dict{Symbol, Dict}(
    # ========== 表1：夸克-夸克散射过程（4个，有t道和u道）==========

    # u d → u d
    :ud_to_ud => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => [:pi, :sigma_pi], :mixed_P => false, :mixed_S => false)
        )
    ),

    # u u → u u
    :uu_to_uu => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true)
        )
    ),
    
    # u s → u s
    :us_to_us => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => [:K, :sigma_K], :mixed_P => false, :mixed_S => false)
        )
    ),

    # s s → s s
    :ss_to_ss => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true)
        )
    ),

    # ========== 夸克-夸克的电荷共轭：反夸克-反夸克散射（用于弛豫时间）==========

    # ū đ → ū đ （与 u d → u d 等价）
    :ubardbar_to_ubardbar => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => [:pi, :sigma_pi], :mixed_P => false, :mixed_S => false)
        )
    ),

    # ū ū → ū ū （与 u u → u u 等价）
    :ubarubar_to_ubarubar => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true)
        )
    ),

    # ū s̄ → ū s̄ （与 u s → u s 等价）
    :ubarsbar_to_ubarsbar => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => [:K, :sigma_K], :mixed_P => false, :mixed_S => false)
        )
    ),

    # s̄ s̄ → s̄ s̄ （与 s s → s s 等价）
    :sbarsbar_to_sbarsbar => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true)
        )
    ),
    
    # ========== 表2：夸克-反夸克散射过程（7个 + 电荷共轭2个，有t道和s道）==========
    
    # u đ → u đ
    :udbar_to_udbar => Dict(
        :type => :qqbar,
        :channels => Dict(
            # t 道包含赝标量混合介子(η/η')与标量混合介子(σ/σ')：与 C++ debug 输出的 heavy(P/S) (k00/k88/k08, detK) 对齐
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :s => Dict(:simple => [:pi, :sigma_pi], :mixed_P => false, :mixed_S => false)
        )
    ),
    
    # d ū → d ū (电荷共轭过程，与u đ等价)
    :dubar_to_dubar => Dict(
        :type => :qqbar,
        :channels => Dict(
            # 电荷共轭过程：与 udbar_to_udbar 采用相同的介子组合
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :s => Dict(:simple => [:pi, :sigma_pi], :mixed_P => false, :mixed_S => false)
        )
    ),
    
    # u ū → u ū
    :uubar_to_uubar => Dict(
        :type => :qqbar,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :s => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true)
        )
    ),
    
    # u ū → d đ
    :uubar_to_ddbar => Dict(
        :type => :qqbar,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => false, :mixed_S => false),
            :s => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true)
        )
    ),

    # u s̄ → u s̄  (修正：原表格误写为s s̄)
    :usbar_to_usbar => Dict(
        :type => :qqbar,
        :channels => Dict(
            :t => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true),
            :s => Dict(:simple => [:K, :sigma_K], :mixed_P => false, :mixed_S => false)
        )
    ),

    # s ū → s ū (电荷共轭过程，与u s̄等价)
    :subar_to_subar => Dict(
        :type => :qqbar,
        :channels => Dict(
            :t => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true),
            :s => Dict(:simple => [:K, :sigma_K], :mixed_P => false, :mixed_S => false)
        )
    ),
    
    # u ū → s s̄
    :uubar_to_ssbar => Dict(
        :type => :qqbar,
        :channels => Dict(
            :t => Dict(:simple => [:K, :sigma_K], :mixed_P => false, :mixed_S => false),
            :s => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true)
        )
    ),
    
    # s s̄ → u ū
    :ssbar_to_uubar => Dict(
        :type => :qqbar,
        :channels => Dict(
            :t => Dict(:simple => [:K, :sigma_K], :mixed_P => false, :mixed_S => false),
            :s => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true)
        )
    ),
    
    # s s̄ → s s̄
    :ssbar_to_ssbar => Dict(
        :type => :qqbar,
        :channels => Dict(
            :t => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true),
            :s => Dict(:simple => Symbol[], :mixed_P => true, :mixed_S => true)
        )
    )
)

"""散射过程 key 的固定列表。

用于避免在其它模块里重复维护 process 列表（例如弛豫时间的 REQUIRED_PROCESSES）。
顺序与 `SCATTERING_MESON_MAP` 的插入顺序一致。
"""
const SCATTERING_PROCESS_KEYS = Tuple(keys(SCATTERING_MESON_MAP))

end # module Constants_PNJL
