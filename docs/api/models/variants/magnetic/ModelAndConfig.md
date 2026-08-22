# 模型与配置合同

本页集中说明 magnetic 主题的模型对象、配置对象与参数模板。

## `PNJLMagneticModel`

定义位于 [src/models/pnjl_physics/PNJLMagneticModel.jl](../../../../src/models/pnjl_physics/PNJLMagneticModel.jl#L1)。

结构：

- `base::PNJLModel`
- `magnetic::MagneticConfig`

它的职责是：

- 验证正磁场配置并把磁场 `Omega` 接入五维、多 seed 的 branch-aware gap solver
- 把 magnetic 配置与磁场热力学能力接入 `Models` 聚合表面

## `MagneticIMCParams`

定义位于 [src/models/pnjl_physics/core/MagneticThermodynamics.jl](../../../../src/models/pnjl_physics/core/MagneticThermodynamics.jl#L45)。

字段：

- `a`
- `b`
- `c`
- `d`
- `Lambda_QCD_MeV`

默认构造器：

- `default_imc_params()`

## `MagneticConfig`

定义位于 [src/models/pnjl_physics/core/MagneticThermodynamics.jl](../../../../src/models/pnjl_physics/core/MagneticThermodynamics.jl#L53)。

关键字段：

- `eB_fm2`
- `n_max`
- `p_num`
- `pz_max`
- `cutoff_N`
- `imc`

默认构造器：

- `default_magnetic_config(; eB_fm2=...)`；必须显式提供不低于门槛的 `eB_fm2`

## 单位与向量合同

- `x_state`：5 维，`(phi_u, phi_d, phi_s, Phi, PhiBar)`
- `mu_vec`：3 维，`(mu_u, mu_d, mu_s)`
- `T_fm`、`mu_vec`、质量、动量：`fm^-1`
- `eB_fm2`：`fm^-2`
- 外部生产扫描的 `eB` 使用 `MeV^2`；当前最小值为 `MAGNETIC_EB_MIN_MEV2 = 100`
- 内部最小值为 `MAGNETIC_EB_MIN_FM2 = 100 / hbarc^2`

`PNJLMagneticModel` 与 `MagneticConfig` 对 `eB_fm2 < MAGNETIC_EB_MIN_FM2` 抛出
`ArgumentError`；通过验证后始终使用 Landau magnetic route。`solve_magnetic_gap`
要求有限且 `T_fm > 0`，当前只允许 `xi=0`。`MagneticGapResult.candidates` 保存去重后的
多分支候选；未启用稳定性分类时标准 `solve_gap` 按已找到候选中的最低 `Omega`
返回 `MeanFieldState`，但不提供局部稳定性证明。启用 `classify_stability=true` 只增加
有限差分 Hessian 诊断标签；该标签不是 PNJL 默认生产过滤条件，所有可行候选仍应保留供
分支策略使用。

## AD residual 诊断路径

`magnetic_gap_residual_autodiff(model, x_state, T_fm, mu_vec; ...)` 使用
`ForwardDiff.gradient` 对固定 `n_max` 的 Landau `Omega` 构造五维 stationarity
residual。`solve_magnetic_gap(...; residual_method=:forward)` 会把该 residual
交给 `NLsolve` 的 `autodiff=:forward`，从而对 residual 再求 Jacobian；这条路径的
目的只是为 finite-difference baseline 提供可重复的 A/B 诊断，不改变默认的
`residual_method=:finite`。

AD 路径必须显式提供 `n_max`（或在 `MagneticConfig` 中固定它），因为自动估计包含
离散的 `floor` 截断规则。它只对固定 cutoff 的同一 Landau kernel 做导数；不宣称
自动微分一定更快，也不改变五维状态、branch candidate 或 Hessian 诊断合同。
建议使用
[`scripts/perf/pnjl/magnetic_autodiff_probe.jl`](../../../../scripts/perf/pnjl/magnetic_autodiff_probe.jl)
在低节点、单 seed 条件下比较 wall time、`f_calls`、`g_calls`、iterations 和残差，
再决定是否扩大实验或改变生产默认值。

磁场模型的 `model_capabilities(model).supports_number_densities` 恒为 `false`：通用
模型接口要求独立的 `quark`/`antiquark` 密度，而磁场适配器只提供
`calculate_magnetic_number_densities(...).net` 净密度；这不影响磁场专用密度 API 的可调用性。

## 配置模板

推荐配置模板： [config/models/pnjl/magnetic_default.toml](../../../../config/models/pnjl/magnetic_default.toml)

关键段落：

- `[magnetic]`
- `[magnetic.imc]`

典型使用顺序：

1. 用 `default_magnetic_config` 构造最小配置
2. 必要时从 TOML 模板同步 `eB_fm2`、`n_max`、`p_num`、`pz_max`
3. 运行后用 `magnetic_nmax_convergence_report` 检查截断质量
