# PNJL CEP 窄窗口 pilot PR

## 范围与不变项

- 基于 `main@246a72cc0cffe24abbe8bca17ec83c726e943e3d` 的隔离分支，仅比较 C2 dense baseline、rho-support cascade 与 high-resolution oracle。
- 覆盖 `xi=-0.5,0,0.5` 的旧 canonical CEP 邻域；完整数值计算只允许由独立 GitHub Actions workflow 执行。
- 不启动全相图 production、C3/O1、transport，不修改或晋升 `data/reference/pnjl`、既有 C0/C1/C2 artifact 和 production/figure 目录。
- `SolverResult v1` 保持不变；solver telemetry 是请求作用域可选对象，不使用全局计数器。

## 固定口径

- `p_num=24`、`t_num=8`、`rs_reduced_adaptive`、`rtol=1e-8`、`atol=1e-10`、`maxevals=10^7`、`iterations=80`、ground-state pressure governance。
- 初始温区为 canonical CEP `±4 MeV`，必要时扩至 `±8 MeV`，温度 bracket 最细 `0.125 MeV`。
- baseline 为 C2 level-3 `Δrho=0.00625` dense secant/Maxwell；cascade 使用 rho-support、previous-T prior 与每切片最多 12 个 targeted 点；oracle 在 `Δrho=0.00625` 后以 `Δrho=0.003125` 复核最终 bracket 两侧。

## 交付物与验收

- [x] `SolverWorkTelemetry` 与 FixedRho joint solver instrumentation。
- [x] 独立 runner、exact `(T,xi,rho)` memoization、9-job Actions matrix 与 failed-job 重跑契约。
- [x] collector、schema/gate smoke、CSV/manifest/claim ledger/plot 脚本。
- [x] Actions 9 jobs 完成并下载到 `D:\Desktop\Julia_RelaxTime_issue130_artifacts\cep_narrow_pilot_20260726`（最终 run `30199406478`）。
- [x] 审计 NaN/Inf、重复键、失败点、SHA/provenance、oracle 双分辨率和 baseline/cascade 成本。
- [x] 根据自动 diagnostic verdict 更新 evidence package、Issue #130，并停在作者物理判断与后续算法决策处。

## 最终诊断状态（2026-07-26）

- calculation SHA：`7e8b91adc0eaa2a01e1e50428870eb562febab70`；9/9 matrix job 与 aggregate success。
- 自动 gate：`diagnostic_only`，`oracle_inconclusive`；`xi=0.5` 未形成 CEP bracket，`xi=0` oracle 的 `ΔT=+0.6875 MeV` 超过预声明 0.5 MeV 容差。
- cascade 的 solver 工作量、residual/Jacobian 调用和 runner 秒数均低于 dense baseline；oracle 因双分辨率复核更贵。
- 证据已导入 `docs/analysis/pnjl/cep_maxwell/narrow_pilot/cep_narrow_pilot_v1/`；不晋升 reference，不启动 transport，等待作者物理判断。

自动 gate 只产生 `pass` 或 `diagnostic_only` 诊断，不替代作者对 CEP 物理正确性的判断，也不触发 reference promotion。
