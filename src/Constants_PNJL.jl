"""
# Constants_PNJL.jl

集中维护 RelaxTime 项目中使用的共用常量。
"""
module Constants_PNJL
export ħc_MeV_fm, α
export N_color, N_flavor, ρ0_inv_fm3, m_ud0_inv_fm, m_s0_inv_fm, Λ_inv_fm, G_fm2, K_fm5
export T0_inv_fm, a0, a1, a2, b3, b4
export λ₀, λ₈, ψ_u, ψ_d, ψ_s, ψbar_u, ψbar_d, ψbar_s
# 基本物理常量-------------------------------------
const ħc_MeV_fm = 197.327  # 197.327 MeV·fm
const α::Float64 = 1.0 / 137.035999084  # 精细结构常数

# PNJL模型参数-------------------------------------
const N_color = 3  # 夸克颜色数(红、绿、蓝)
const N_flavor = 3  # 夸克味道数(u, d, s)
const ρ0_inv_fm3 = 0.16  # 核子数密度, 单位fm⁻³ (ρ0 ≈ 0.16 fm⁻³)
const m_ud0_inv_fm = 5.5 / ħc_MeV_fm  # u, d 夸克裸质量，单位 fm⁻¹ (m_ud0 ≈ 5.5 MeV)
const m_s0_inv_fm = 140.7 / ħc_MeV_fm  # s 夸克裸质量，单位 fm⁻¹ (m_s0 ≈ 140.7 MeV)
const Λ_inv_fm = 602.3 / ħc_MeV_fm  # 截断参数，单位 fm⁻¹ (Λ ≈ 602.3 MeV)
const G_fm2 = 1.835 / Λ_inv_fm^2  # NJL 四夸克相互作用耦合常数，单位 fm² (G ≈ 1.835Λ²)
const K_fm5 = 12.36 / Λ_inv_fm^5  # NJL 六夸克相互作用耦合常数(Kobayashi-Maskawa-'t Hooft)，单位 fm⁵ (K ≈ 12.36Λ⁵)

# Polyakov环有效势参数-------------------------------------
const T0_inv_fm = 210.0 / ħc_MeV_fm  # Polyakov有效势参数 ，单位 fm⁻¹ (T0 ≈ 210 MeV)
const a0 = 3.51
const a1 = -2.47
const a2 = 15.2
const b3 = -1.75
const b4 = 7.555

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
doc/formula/散射过程所有可能.md 表1和表2
"""
const SCATTERING_MESON_MAP = Dict{Symbol, Dict}(
    # ========== 表1：夸克-夸克散射过程（4个，有t道和u道）==========
    
    # u u → u u
    :uu_to_uu => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true)
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
    
    # u d → u d
    :ud_to_ud => Dict(
        :type => :qq,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :u => Dict(:simple => [:pi, :sigma_pi], :mixed_P => false, :mixed_S => false)
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
    
    # ========== 表2：夸克-反夸克散射过程（7个，有t道和s道）==========
    
    # u đ → u đ
    :udbar_to_udbar => Dict(
        :type => :qqbar,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :s => Dict(:simple => [:pi, :sigma_pi], :mixed_P => false, :mixed_S => false)
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
    
    # d ū → d ū (电荷共轭过程，与u đ等价)
    :dubar_to_dubar => Dict(
        :type => :qqbar,
        :channels => Dict(
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :s => Dict(:simple => [:pi, :sigma_pi], :mixed_P => false, :mixed_S => false)
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

end # module Constants_PNJL
