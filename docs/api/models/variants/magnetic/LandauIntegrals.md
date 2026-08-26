# Landau 能级与积分接口

本页说明 magnetic 主题中与 Landau 能级离散化直接相关的低层公开接口。实现位于 [src/models/pnjl_physics/core/MagneticIntegrals.jl](../../../../src/models/pnjl_physics/core/MagneticIntegrals.jl#L1)。

## 主要导出

- `QUARK_CHARGE_ABS`
- `alpha_n`
- `energy_landau`
- `smooth_cutoff`
- `resolve_nmax_from_cutoff`
- `resolve_magnetic_nmax`
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

## `resolve_magnetic_nmax` 与 `resolve_nmax_from_cutoff`

`resolve_magnetic_nmax(T_fm, mu_vec, magnetic; masses=nothing)` 是生产 profile 的
统一入口。显式 `MagneticConfig.n_max` 优先；默认 `:thermal_tail` 使用
`E_tail=max(abs.(mu_vec)) + thermal_tail_factor*T_fm` 估计热尾，并在
`n_max_floor`/`n_max_cap` 预算内返回一次共享层数。`solve_magnetic_gap` 对同一个
`(T,mu,eB)` 点的全部 seed、primary 和 fallback attempt 复用该结果。

`resolve_nmax_from_cutoff` 仍根据质量、化学势、电荷与磁场强度给出 legacy
`:vacuum_cutoff` 策略的起点；它不是默认 MFIR 真空项的求和上限，也不替代收敛验证。
真正的主题级收敛检查仍应通过 `magnetic_nmax_convergence_report` 完成。

## `omega0_flavor_landau` / `omegat_flavor_landau`

分别给出单味 Landau 真空项与热项。它们是 thermodynamics 组装层的低级输入，不建议多数用户直接从这里开始。

## `density_flavor_landau`

提供单味数密度近似路径。该接口通常与 `calculate_magnetic_number_densities` 配合理解，而不是单独当作主题主入口。
