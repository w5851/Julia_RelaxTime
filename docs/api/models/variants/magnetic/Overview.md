# Magnetic 主题总览

本页优先回答“如何在 `Models` 入口下使用 magnetic 变体”，而不是先展开 Landau 求和与热力学细节。

## 何时使用该主题

当你需要以下任一能力时，应从本主题开始：

- 构造外磁场 PNJL 模型对象
- 使用 magnetic 配置对象计算磁场下的 `omega`、压力与数密度
- 检查 `n_max` 截断是否足够收敛
- 对固定点或固定 `eB` 扫描脚本建立文档化入口

## 首选公开入口

主题主入口来自 [src/models/Models.jl](../../../../src/models/Models.jl#L25) 以及 [src/models/pnjl_physics/PNJLMagneticModel.jl](../../../../src/models/pnjl_physics/PNJLMagneticModel.jl#L1)。

优先关注：

- `PNJLMagneticModel`
- `MagneticConfig`
- `default_magnetic_config`
- `calculate_magnetic_omega`
- `calculate_magnetic_pressure`
- `calculate_magnetic_number_densities`
- `magnetic_nmax_convergence_report`

完整导出基线见 [generated/Exports.md](generated/Exports.md)。

## 最短使用路径

```julia
using .Models

model = Models.PNJLMagneticModel(; eB_fm2=0.10)
config = Models.default_magnetic_config(; eB_fm2=0.10)

x_state = @SVector [0.0, 0.0, 0.0, 0.2, 0.2]
mu_vec = @SVector [0.0, 0.0, 0.0]
T_fm = 150.0 / Models.Constants_PNJL.ħc_MeV_fm

omega_val = Models.calculate_magnetic_omega(x_state, mu_vec, T_fm, config)
pressure_val = Models.calculate_magnetic_pressure(x_state, mu_vec, T_fm, config)
densities = Models.calculate_magnetic_number_densities(x_state, mu_vec, T_fm, config)
report = Models.magnetic_nmax_convergence_report(x_state, mu_vec, T_fm, config)
```

## 输入与单位口径

- `T_fm`、`mu_vec`、质量、动量使用 `fm^-1`
- `eB_fm2` 使用 `fm^-2`
- `x_state` 的 5 维顺序为 `(phi_u, phi_d, phi_s, Phi, PhiBar)`
- `mu_vec` 的 3 维顺序为 `(mu_u, mu_d, mu_s)`

参数与配置细节见 [ModelAndConfig.md](ModelAndConfig.md)。

## 非首页首选入口

以下接口属于进阶或实现相关入口：

- `alpha_n`
- `energy_landau`
- `smooth_cutoff`
- `resolve_nmax_from_cutoff`
- `omega0_flavor_landau`
- `omegat_flavor_landau`
- `density_flavor_landau`

这些能力的职责边界见 [CoreConcepts.md](CoreConcepts.md) 与 [LandauIntegrals.md](LandauIntegrals.md)。