# `Models.run_tmu_scan`

本页从 `Models` 统一入口说明 T-μ 扫描，而不是从 `PNJL.TmuScan` 子模块出发。

## 定位

`Models.run_tmu_scan` 是 T-μ 网格扫描的首选业务入口。它适合：

- 批量生成 `(T, μ, ξ)` 网格结果
- 以 phase-aware 或 continuity 方式沿扫描路径推进
- 为后续分析、作图或回归脚本提供稳定 CSV 产物

实现转发位于 `src/models/entrypoints.jl`，底层执行位于 `src/models/scans/TmuScan.jl`。

## 入口形态

从当前实现看，核心关键字参数包括：

- `T_values`
- `mu_values`
- `xi_values`
- `output_path`
- `overwrite`
- `resume`
- `use_phase_aware`
- `bootstrap_multiseed`
- `continuation`（`:auto`、`:legacy`、`:single_pass`、`:palc`）与 `palc_options`
- `solver_backend`
- `p_num`, `t_num`
- `thermo_quadrature_policy`（`:tensor_gauss` 或 PNJL 标量热核的 `:rs_reduced_adaptive`）
- `thermo_quadrature_rtol`, `thermo_quadrature_atol`, `thermo_quadrature_maxevals`
- `progress_cb`

`:rs_reduced_adaptive` 当前明确限定为 `model_kind=:PNJL`；其他模型会提前报错，不会把配置写入 provenance 后仍静默使用固定网格。

返回值为命名元组，包含：

- `total`
- `success`
- `failure`
- `skipped`
- `output`

## 何时优先使用它

优先于直接调用底层模块的情形：

- 你只关心完成扫描任务，而不想自己管理模块路径
- 你希望调用口径与 `Models` 其他主题保持一致
- 你需要稳定的 CSV 结果与断点续扫行为

## 与 phase-aware 策略的关系

`use_phase_aware` 是兼容参数。`:legacy` 沿每条等温线使用连续性种子和全候选选择。
`:single_pass` 展开候选种子并去重后逐点求解；`:palc` 从多起点双向延拓，校正到请求的
每个 μq 后选择通过硬约束的最大压力根。`:auto` 当前沿用 `:legacy`，待生产对照资格确认。
显式新策略支持 PNJL/RPNJL、tensor 积分与默认 ground_state selector；不支持的组合抛出
`ArgumentError`。`palc_options` 传入 `PALCContinuation`，例如 `(ds=0.02, dsmax=0.08)`。
自定义 `solver` 或路径未支持的求解器选项由 `:auto` 保留在旧链路；显式新策略在扫描前报错。

## 路径依赖与扫描顺序

T-μ 扫描的连续性跟踪是路径依赖的：策略对象只知道“上一个点”的收敛解，因此你按什么顺序遍历网格，就沿那条路径延续。

推荐做法：

- 固定一条 `μ` 扫描线，沿 `T` 推进，并在换线时 `reset!`
- 或固定一条 `T` 扫描线，沿 `μ` 推进，并在换线时 `reset!`

不推荐把二维网格当作任意蛇形路径处理，否则分支跟踪结果会随遍历顺序改变。

## 输出合同

T-μ 扫描的稳定跨模块产物是 CSV，而不是 Julia 内部结构本身。当前 header 为：

```text
T_MeV,mu_MeV,xi,pressure_fm4,rho,entropy_fm3,energy_fm4,phi_u,phi_d,phi_s,Phi1,Phi2,M_u_MeV,M_d_MeV,M_s_MeV,iterations,residual_norm,converged,message
```

读取约束：

- 失败点也会写入 CSV，并通过 `converged=false` 与 `message` 标记
- 下游应按列名而不是列位置解析
- 若将来扩展 schema，应优先追加列而不是重命名现有列
- `T_values` 和 `mu_values` 使用 MeV，μ 是共同夸克化学势 μq；CSV 的 `rho` 与
  callback 中历史命名的 `result.rho_norm` 均为重子密度 fm⁻³。通用 `solve_path` 的
  `SolverResult.rho_norm` 则为 ρ/ρ₀，扫描适配层显式换算。
- `path_diagnostics_cb((T, xi), path_result)` 在 PALC 每条等温线结束后调用。
  结果仍按输入网格顺序输出；resume 跳过已有点，若一条线部分完成，PALC 重建该线后只写缺点。

## 相关主题

- [Overview.md](Overview.md)
- [Algorithms.md](Algorithms.md)
- [TrhoScan.md](TrhoScan.md)
