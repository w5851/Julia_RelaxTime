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
- `solver_backend`
- `p_num`, `t_num`
- `progress_cb`

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

当 `use_phase_aware=true` 时，入口会尝试建立 `PhaseAwareContinuitySeed`，并在首点可选地配合 `bootstrap_multiseed=true` 做稳定分支自举。

这意味着：

- 该入口并不只是把参数循环交给求解器
- 它同时承担“路径组织 + 分支跟踪 + 输出合同”三层职责

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

## 相关主题

- [Overview.md](Overview.md)
- [Algorithms.md](Algorithms.md)
- [TrhoScan.md](TrhoScan.md)