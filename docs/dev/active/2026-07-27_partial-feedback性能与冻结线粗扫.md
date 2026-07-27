# Partial-feedback 性能基准与冻结线粗扫

更新日期：2026-07-27

当前状态：in progress。本文档只记录当前 `x_min_cut` partial-feedback 的成本与定性趋势诊断，不是 strict-support、节点收敛或 production 任务单。

## 目标与固定边界

- 在排除 Julia 冷启动/JIT 后估计一个固定 `(T,mu_B)` partial-feedback 单点的平均耗时。
- 沿仓库现有 `default` / `baseline_freezeout` 参数化稀疏扫描 `K^+/pi^+` 与 `K^-/pi^-`。
- 保留 `rho_Q/rho_B=0.4`、`rho_S=0` 的 quark-only 内层约束，并在外层加入当前四个 charged BU 密度。
- 使用 `p_num=8`、`t_num=4`、`qmax=4`、`q_nodes=4`、`omega=[0.05,3]`、`omega_nodes=8`、`density_policy=:x_min_cut`、`bose_x_min=0.05` 的低成本诊断设置。
- outer 从 quark-only seed 起步时允许最多 12 次迭代、50 次 evaluator calls；残差门限仍为 `2e-3`，不通过放宽 residual tolerance 获取收敛。
- 不修改稳定 `Models.run_freezeout_meson_density_scan`、正式 BU kernel、transport 或既有 regression baseline。

## 运行入口

1. 热启动 benchmark：

   `julia --project=. scripts/perf/relaxtime/bench_meson_conserved_charge_partial_feedback.jl`

   首次完整点作为 warm-up；默认 3 个新 cache 样本写入 `data/outputs/perf/`。

2. 冻结线粗扫：

   `julia --project=. scripts/analysis/relaxtime/scan_meson_conserved_charge_feedback_freezeout.jl`

   默认读取 benchmark summary 的中位耗时，并在约 600 s 预算内从 7/5/3 点网格中选择最大网格；计算从高 `sqrt(s_NN)` 向低能 continuation，CSV 最终按能量升序写出。

## 结果记录

运行完成后在本节追加：

- warm-up 与各热启动样本耗时、median/mean/min/max；
- 选用的 `sqrt(s_NN)` 网格及每点状态；
- `K^+/pi^+`、`K^-/pi^-`、四个组成密度、守恒荷残差和 `m_M-mu_M`/Bose status；
- 只讨论趋势与数量级，明确 `x_min_cut` 与低节点限制。

原始 CSV 和 perf JSON 保持为本地 diagnostic 产物，不作为 production artifact 提交。

## 2026-07-27 运行结果

### 热启动 benchmark

- warm-up：`138.787 s`，包含本进程首次 JIT/初始化，不计入统计；
- 3 个 fresh-cache 样本：`min=0.196 s`、`median=0.248 s`、`mean=0.254 s`、`max=0.318 s`；
- 每个样本均收敛，`outer residual_norm=1.9200e-3`，`outer_iterations=9`，`outer_evaluations=28`；
- baseline 约 `0.057--0.098 s`；候选 evaluator 合计约 `0.135--0.209 s`，其中 gap 约 `0.026--0.048 s`、四通道 BU 约 `0.107--0.158 s`；
- benchmark JSON：`data/outputs/perf/relaxtime/meson_conserved_charge_partial_feedback/hot_start_summary.json`。

### 冻结线 sparse scan

benchmark median 与 600 s 预算选择 7 点：
`sqrt(s_NN)=[3,4.5,7.7,11.5,20,62.4,200] GeV`。

| `sqrt(s_NN)` [GeV] | `T` [MeV] | `K+/pi+` | `K-/pi-` | residual norm |
|---:|---:|---:|---:|---:|
| 3.0 | 79.96 | 0.958 | 0.952 | 1.83e-3 |
| 4.5 | 111.82 | 0.954 | 0.949 | 1.15e-4 |
| 7.7 | 139.61 | 0.961 | 0.941 | 3.42e-4 |
| 11.5 | 151.59 | 0.970 | 0.938 | 2.05e-4 |
| 20.0 | 160.21 | 0.977 | 0.941 | 4.01e-4 |
| 62.4 | 165.27 | 0.971 | 0.951 | 3.16e-4 |
| 200.0 | 165.92 | 0.965 | 0.958 | 3.97e-4 |

7 个点均报告 `partial_feedback_converged`，四个 charged density status 均为 `ok`。最终重跑中 `200 GeV` 首点耗时约 `49.8 s`，是该扫描进程的冷启动/JIT 影响；其余点耗时约 `0.064--0.785 s`。

当前口径下两个 ratio 都近似平坦在 `O(1)`，没有显著的 `0.3` 数量级或 horn-like 结构。该结论受 `x_min_cut`、低节点和未纳入 `Omega_M` stationarity 的 partial feedback 限制，只用于决定后续 kernel/反馈研究优先级，不是实验拟合或 production 结论。

扫描 CSV/README：`data/outputs/results/relaxtime/analysis/meson_conserved_charge_feedback_freezeout_sparse/`。
