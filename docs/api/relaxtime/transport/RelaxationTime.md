# 弛豫时间 API

`relaxation_times(quark_params, thermo_params, K_coeffs; densities, ...)`
返回六种粒子的 `tau`、`tau_inv` 和可复用的通道 `rates`。
支持 `QuarkParams/ThermoParams` 及等价 NamedTuple，内部使用自然单位制。

## 聚合关系

`tau_i⁻¹ = Σ_j n_j * w_ij`，其中数密度单位为 fm⁻³，平均速率为 fm²，τ 为 fm。
同位旋对称下 d/dbar 与 u/ubar 共用结果。各道的重数和粒子/反粒子密度组合保持原有公式，
见 [公式说明](../../../reference/formula/relaxtime/transport/RelaxationTime_FromAverageScatteringRate.md)。

- `relaxation_rates(densities, rates)`：只做已给定密度、通道速率的代数组合。
- `compute_average_rates(...; existing_rates, cs_caches, ...)`：补全缺失的通道。
- `relaxation_times(...; densities, existing_rates, cs_caches, ...)`：完整入口。
- `relaxation_rate_contribution_rows(densities,rates)`：逐通道贡献诊断。

`REQUIRED_PROCESSES` 给出兼容的基本过程清单；完整 τ 聚合和生产建表使用
`REQUIRED_RATE_KEYS_FOR_TAU` 的 17 道顺序，其中包含反粒子过程。

`existing_rates` 足够时直接复用。过程别名与 `rate_lookup -> rate_value` 的既有兼容关系保留；
它们与被移除的截面采样算法是不同的接口。

## 采样默认值

缺失通道现在仅使用 `ThermalWeightedCache`。按 `REQUIRED_RATE_KEYS_FOR_TAU` 的顺序处理；
新建每道缓存前回收内存并重置极化缓存，避免跨道的有效位量化复用改变已验收路径。
已经存在的缓存强制校验物理上下文、节点、权重和误差预算。

默认外层为半无穷 p48/cos36/phi36，t 积分 24 点，截面的物理支持使用 Λ。
`sampling_config=ThermalSamplingConfig(rtol=1e-3)` 控制固定外层求积上的插值预算。
`sigma_cache_policy` 只接受 `:thermal_weighted`，`interpolation_mode` 接受 `:auto|:thermal_weighted`。
旧 w0cdf 网格、扣除/加回参数、`CrossSectionCache` 和无身份的两列 CSV 缓存加载已移除。

`densities` 由调用方提供，须与背景及尾项约定相符。低层默认值不会自动配置背景、A 或末端输运；
生产使用 [TransportSamplingProfiles](TransportSamplingProfiles.md) 统一配置这些积分层。
独立控制末端 η、ζ、σ 的动量积分，不要把它们的上限与截面物理截断混用。

## 示例

```julia
result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
result.tau
result.tau_inv
result.rates
```

可用 `cs_caches=Dict(process=>cache, ...)` 复用已认证缓存。
`p_grid/p_w`、`cos_grid/cos_w` 和 `phi_grid/phi_w` 支持成对的自定义求积，
但改变节点后必须重建缓存并另做外层收敛检查。

默认传播子与热分布使用同一 ξ；诊断选项 `propagator_xi_policy=:isotropic` 仅改变传播子上下文，
可配合显式 `propagator_quark_params`。细节与失败语义见 [AverageScatteringRate](AverageScatteringRate.md)。

旧 Fortran/C++ 的有限积分与旧截面表不能当作新算法的收敛参考。
保留同点分层比较；外部读图的 5% 标准与新配置最终量的 1% 收敛目标分别记录。
