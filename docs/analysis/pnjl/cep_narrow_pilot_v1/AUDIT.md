# PNJL CEP 窄窗口 pilot v1 审计

## Provenance

- calculation SHA：`7e8b91adc0eaa2a01e1e50428870eb562febab70`
- Actions run：[#30199406478](https://github.com/w5851/Julia_RelaxTime/actions/runs/30199406478)
- source baseline：`main@246a72cc0cffe24abbe8bca17ec83c726e943e3d`
- 矩阵：9 个 `(xi, method)` job，`xi ∈ {-0.5, 0, 0.5}`；全部 job 与 merged artifact 使用同一 calculation SHA。
- Actions critical path：625 s；原始总 job 时间：3440 s；向上取整 runner-minutes：64（含 aggregate job）。

完整数值只在上述 GitHub Actions run 执行；本地仅重新收集、绘图和审计 artifact。

## 数值完整性

- 9/9 matrix key 唯一且齐全；所有 job `finite_and_converged=true`。
- equilibrium/fixed-rho 解没有失败点、Inf 或非有限解；所有 targeted slice 的最大补点数为 9，满足每切片 ≤12。
- CSV 中出现的 NaN 只属于派生诊断量：单调 slice 的 `mu_transition`、Maxwell 未解析的 `near_critical` 量，或未找到 CEP bracket 的 CEP 字段；不代表采用了非有限 equilibrium 解。
- oracle 的 `0.00625 → 0.003125` 复核在 `xi=-0.5` 和 `xi=0` 的 slice classification 稳定；`xi=0.5` 在 ±8 MeV 窗口内未形成 bracket，因此 oracle 证据不完整。

## CEP 结果摘要

| xi | 方法 | T_CEP [MeV] | μq_CEP [MeV] | ΔT vs canonical | Δμq vs canonical | 备注 |
|---:|---|---:|---:|---:|---:|---|
| -0.5 | baseline/oracle | 147.1572266 | 202.24668 | +0.0625 | -0.1475 | oracle 双分辨率稳定 |
| -0.5 | cascade | 147.1572266 | 202.24669 | +0.0625 | -0.1475 | 位于同一 bracket |
| 0 | baseline/oracle | 131.0244141 | 291.26170 | +0.6875 | -1.1191 | `|ΔT|` 超过预声明 0.5 MeV |
| 0 | cascade | 131.1494141 | — | +0.8125 | — | 低侧为 near-critical，Maxwell transition 未解析 |
| 0.5 | 三方法 | — | — | — | — | ±8 MeV 内未找到 bracket，oracle inconclusive |

因此自动 gate verdict 为 `diagnostic_only`，不是 reference 晋升结论。当前证据不能支持把新算法或新 CEP 结果标记为 production-grade。

## Solver cost

累计 `method_costs.csv` 显示：

| 方法 | equilibrium requests | unique solves | point requests | cache hits | residual calls | Jacobian calls | runner seconds |
|---|---:|---:|---:|---:|---:|---:|---:|
| C2 dense baseline | 26,281 | 26,281 | 26,281 | 0 | 288,466 | 231,681 | 843.6 |
| rho-support cascade | 3,536 | 3,536 | 7,072 | 3,536 | 38,696 | 31,119 | 157.7 |
| high-resolution oracle | 28,841 | 28,841 | 31,405 | 2,564 | 326,287 | 263,700 | 1,066.6 |

在本窗口内 cascade 的 solver 工作量和 runner 时间均显著低于 dense baseline；oracle 因双分辨率复核而更贵。分类层耗时不是 gate blocker。fallback/rescue/retry 计数已保留在 `method_costs.csv`，不应从本 pilot 推断生产路径的全温区性能。

## Gate 与后续边界

自动 gate 为：

- `oracle_status=oracle_inconclusive`；`cascade_status=outside_oracle`；`status=diagnostic_only`。
- 失败原因：`xi=0.5` 缺少 oracle bracket/复核；`xi=0` oracle 的 `ΔT=+0.6875 MeV` 超出预声明 canonical 容差；cascade 不能在全部 xi 上落入 oracle bracket。

本 PR 不修改 `data/reference/pnjl`，不覆盖 C0/C1/C2 artifact，不导入 production/reference，也不启动 C3/O1/transport。后续应由作者先做物理判断，再决定是否需要修正 S-shape/Maxwell/CEP bracket 算法或重新设计 pilot；不应依据本自动 gate 直接晋升。

## 失败尝试记录

- run `30198397023`（SHA `49dfae49…`）：cascade job 在 artifact 收尾时因 JSON3 拒绝 NaN 失败；已修复为 JSON `null`，未重算其余仍运行 job。
- run `30198755883`（SHA `53bb44b3…`）：序列化修复后完成，但暴露 `minimum(...; init=NaN)` 使 bracket 永远丢失以及 collector 把 job 总补点误当作 slice 上限；仅作为诊断证据保留。
- 本审计使用修复后的 run `30199406478`，未改动旧 reference。
