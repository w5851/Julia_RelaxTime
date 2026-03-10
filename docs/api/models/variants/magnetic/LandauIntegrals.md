# Landau 能级与积分接口

本页说明 magnetic 主题中与 Landau 能级离散化直接相关的低层公开接口。实现位于 [src/models/pnjl_physics/core/MagneticIntegrals.jl](../../../../src/models/pnjl_physics/core/MagneticIntegrals.jl#L1)。

## 主要导出

- `QUARK_CHARGE_ABS`
- `alpha_n`
- `energy_landau`
- `smooth_cutoff`
- `resolve_nmax_from_cutoff`
- `omega0_flavor_landau`
- `omegat_flavor_landau`
- `density_flavor_landau`

## `QUARK_CHARGE_ABS`

三味电荷绝对值常量，按 flavor 顺序组织。它属于公开常量，但更偏实现与推导辅助，而不是首页首选入口。

## `alpha_n`

Landau 简并度：

- `n = 0 -> 1`
- `n > 0 -> 2`

## `energy_landau`

定义单个 flavor、单个 Landau 层上的能谱：

$$E_{f,n} = \sqrt{2 n |q_f| eB + p_z^2 + M_f^2}$$

## `smooth_cutoff`

平滑截断因子，用于 Landau 真空项积分中的 UV 抑制。若只想做常规业务调用，不建议直接把它当作入口。

## `resolve_nmax_from_cutoff`

根据质量、化学势、电荷与磁场强度估计 `n_max`。它的职责是给出收敛起点，而不是替代收敛验证；真正的主题级收敛检查仍应通过 `magnetic_nmax_convergence_report` 完成。

## `omega0_flavor_landau` / `omegat_flavor_landau`

分别给出单味 Landau 真空项与热项。它们是 thermodynamics 组装层的低级输入，不建议多数用户直接从这里开始。

## `density_flavor_landau`

提供单味数密度近似路径。该接口通常与 `calculate_magnetic_number_densities` 配合理解，而不是单独当作主题主入口。