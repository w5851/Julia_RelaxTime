"""TransportConstants

与输运/散射过程相关的“物理过程枚举与映射”常量。

动机：
- 这些常量属于输运/散射振幅层，而不是 PNJL 模型本身
- 后续不同模型（NJL/PNJL/rPNJL）若共享同一套散射过程结构，可直接复用

注意：
- 当前数据来源与原 `Constants_PNJL.SCATTERING_MESON_MAP` 保持一致
- 该模块不依赖 PNJL 参数文件，仅提供过程结构信息
"""
module TransportConstants

export SCATTERING_MESON_MAP, SCATTERING_PROCESS_KEYS

# 散射过程到介子种类的映射表-------------------------------------
"""散射过程介子映射表

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
            # t 道包含赝标量混合介子(η/η')与标量混合介子(σ/σ')
            :t => Dict(:simple => [:pi, :sigma_pi], :mixed_P => true, :mixed_S => true),
            :s => Dict(:simple => [:pi, :sigma_pi], :mixed_P => false, :mixed_S => false)
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

    # u s̄ → u s̄
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

end # module TransportConstants
