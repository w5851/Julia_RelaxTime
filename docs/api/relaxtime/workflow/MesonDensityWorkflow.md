# MesonDensityWorkflow

将 meson workflow 的返回值后处理为介子数密度。

本页是 meson density workflow 的领域细节页，重点说明它如何与现有 meson workflow 对接，以及为什么它被约束为后处理层。

如果你只是想判断“应该从哪个 `Models` 入口开始”，优先阅读 [docs/api/models/workflows/MesonDensityWorkflow.md](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/docs/api/models/workflows/MesonDensityWorkflow.md)。

## 入口

### `solve_meson_density_from_meson_point`

```julia
solve_meson_density_from_meson_point(
    meson_point;
    pi_channel=:pi,
    k_channel=:K,
    μ_pi=0.0,
    μ_K=0.0,
    d_pi=3,
    d_K=4,
    qmax_pi=nothing,
    qmax_K=nothing,
    num_q_nodes=256,
)
```

它要求输入对象至少包含：

- `thermo_params.T`
- `thermo_params.ξ`
- `meson_results[:pi].mass`
- `meson_results[:K].mass`

workflow 会直接从这些字段读取温度与介子质量，然后调用 `Main.MesonDensity.stable_meson_number_density` 完成数密度计算。

### `solve_gap_and_meson_density_point`

```julia
solve_gap_and_meson_density_point(T_fm, mu_fm; density_kwargs=(;), kwargs...)
```

它内部先调用 `solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)`，再执行：

```julia
solve_meson_density_from_meson_point(meson_point; density_kwargs...)
```

因此它的物理主链固定为：

```text
gap / equilibrium
  -> meson mass / threshold / gap
  -> meson density
```

而不是：

```text
script 手工取温度与质量
  -> 自行拼数密度
```

## 当前设计原则

1. workflow 是唯一计算入口；
2. 脚本只能消费 workflow 返回值；
3. 数密度层不重写 meson 求解链；
4. 后续 BW / BU 扩展仍应在此后处理边界内继续演进。
