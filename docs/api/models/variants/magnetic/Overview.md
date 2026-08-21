# Magnetic 主题总览

本页优先回答“如何在 `Models` 入口下使用 magnetic 变体”，而不是先展开 Landau 求和与热力学细节。

## 何时使用该主题

当你需要以下任一能力时，应从本主题开始：

- 构造外磁场 PNJL 模型对象
- 使用 magnetic 配置对象计算磁场下的 `omega`、压力与数密度
- 检查 `n_max` 截断是否足够收敛
- 对固定点或固定 `eB` 扫描脚本建立文档化入口
- 运行完整五维 `FixedMu` 的 `(T,mu,eB)` 生产扫描

## 首选公开入口

主题主入口来自 [src/models/Models.jl](../../../../src/models/Models.jl#L25) 以及 [src/models/pnjl_physics/PNJLMagneticModel.jl](../../../../src/models/pnjl_physics/PNJLMagneticModel.jl#L1)。

优先关注：

- `PNJLMagneticModel`
- `solve_magnetic_gap`
- `magnetic_gap_residual`
- `MagneticConfig`
- `default_magnetic_config`
- `calculate_magnetic_omega`
- `calculate_magnetic_pressure`
- `calculate_magnetic_number_densities`
- `magnetic_nmax_convergence_report`
- `Models.run_magnetic_scan`

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

# 非零磁场的完整五维 branch-aware equilibrium
model = Models.PNJLMagneticModel(; eB_fm2=0.10)
result = Models.solve_magnetic_gap(model, T_fm, mu_vec; p_num=24, t_num=8)
state = result.state
all_branches = result.candidates
```

## 输入与单位口径

- `T_fm`、`mu_vec`、质量、动量使用 `fm^-1`
- `eB_fm2` 使用 `fm^-2`
- `x_state` 的 5 维顺序为 `(phi_u, phi_d, phi_s, Phi, PhiBar)`
- `mu_vec` 的 3 维顺序为 `(mu_u, mu_d, mu_s)`
- `solve_magnetic_gap` 要求 `T_fm > 0`，当前只允许 `xi=0`；`result.candidates` 才是分支审计的完整输出。未启用稳定性分类时，`result.state` 按已找到候选中的最低 `Omega` 选择，但不证明该根局部稳定。启用 `classify_stability=true` 只增加有限差分 Hessian 标签和显式诊断信息；该标签不是 PNJL 默认生产过滤条件，也不会替代候选集合或分支选择策略。
- 非零 `eB` 时，`calculate_magnetic_number_densities` 返回 `net` 和历史兼容字段 `quark`（二者均为 `q-qbar`），并将 `antiquark` 明确设为 `nothing`；磁场路线没有独立的夸克/反夸克数密度输出，不应把该结果当作普通 PNJL 输运输入。`PNJLMagneticModel.number_densities` 在 `eB≈0` 时退化到普通 PNJL 的独立 `quark/antiquark` 结果，以保持零场兼容。
- 非零 `eB` 的模型 capability 会将通用 `number_densities` 标为不支持，因为该通用合同要求独立 `quark/antiquark`；磁场密度请使用上面的专用 API。`eB≈0` 时 capability 和适配器都恢复普通 PNJL 的独立密度合同。
- 共享 `solve_constraint`/ProblemSpec 约束入口暂不承载磁场专用 Omega，会显式拒绝 `PNJLMagneticModel`；FixedMu 磁场平衡态请使用上面的 `solve_magnetic_gap`。

参数与配置细节见 [ModelAndConfig.md](ModelAndConfig.md)。

生产扫描的外部单位、selected/candidates 双 CSV 合同和 CLI 用法见
[MagneticScan.md](../../scans/MagneticScan.md)。

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
