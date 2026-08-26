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
- `route`（默认 `:mfir`；可显式指定 `:landau_legacy` 做历史诊断）
- `zeta_num`（MFIR Hurwitz-zeta 特殊函数的固定积分节点数）
- `imc`
- `n_max_policy`（默认 `:thermal_tail`；也可显式使用 `:vacuum_cutoff` 做 legacy 诊断）
- `thermal_tail_factor`、`n_max_floor`、`n_max_cap`（自动热 Landau 截断预算）

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
`ArgumentError`；通过验证后默认使用“零场三动量截断 + MFIR/Hurwitz-zeta 磁场真空
修正 + Landau 热项”路线。只有显式设置 `route=:landau_legacy` 才使用历史完整
Landau 真空项与平滑截断；该路线不属于生产或 acceptance。`solve_magnetic_gap`
要求有限且 `T_fm > 0`，当前只允许 `xi=0`。`MagneticGapResult.candidates` 保存去重后的
多分支候选；未启用稳定性分类时标准 `solve_gap` 按已找到候选中的最低 `Omega`
返回 `MeanFieldState`，但不提供局部稳定性证明。启用 `classify_stability=true` 只增加
有限差分 Hessian 诊断标签；该标签不是 PNJL 默认生产过滤条件，所有可行候选仍应保留供
分支策略使用。`n_max` 显式给定时优先使用该值；未给定时，默认
`resolve_magnetic_nmax` 按 `E_tail=|mu|_max + thermal_tail_factor*T` 和最小夸克电荷
估计热项所需的 Landau 层数，并应用 `n_max_floor`/`n_max_cap`。同一个物理点的全部
seed、primary 和 fallback attempt 共用一次解析结果，避免候选之间比较不同截断。

## AD residual 生产路径

`magnetic_gap_residual(model, x_state, T_fm, mu_vec; ...)` 使用
`ForwardDiff.gradient` 对固定 `n_max` 的 Landau `Omega` 构造五维 stationarity
residual；`magnetic_gap_residual_autodiff` 仅保留为同一实现的兼容别名。
`solve_magnetic_gap` 对 primary/fallback 都固定使用 `NLsolve(autodiff=:forward)`。

离散 Landau 截断不能在 ForwardDiff 内动态改变：如果调用方没有显式给出
`n_max`，solver 会在进入 seed loop 前为该物理点解析一次，并在所有 seed 的
primary/fallback residual、Jacobian、收敛后复核和分支候选中复用它。直接调用 residual
时则必须显式给出 `n_max`，以避免把离散截断混入导数。`resolve_magnetic_nmax` 是
可单独复核该 profile 解析结果的公开入口。

`finite_difference_step` 只属于可选 Hessian 稳定性诊断的外层差分，不再参与
stationarity residual 或 NLsolve Jacobian。五维状态、branch candidate 和
“Hessian 不是默认生产筛选条件”的合同保持不变。稳态诊断入口见
[`scripts/perf/pnjl/magnetic_autodiff_probe.jl`](../../../../scripts/perf/pnjl/magnetic_autodiff_probe.jl)。

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
2. 生产默认使用 `magnetic_default.toml`；跨求解器 source-parity 诊断可显式选择
   `profile="magnetic_source_parity"`（固定 `n_max=79`、`p_num=128`、`pz_max=40`、
   `zeta_num=256`）
3. 必要时显式覆盖 `n_max` 或热尾项参数；显式值优先于自动 profile
4. 运行后用 `magnetic_nmax_convergence_report` 检查代表点的截断质量
