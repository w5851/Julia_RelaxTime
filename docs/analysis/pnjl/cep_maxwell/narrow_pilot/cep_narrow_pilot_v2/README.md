# PNJL CEP 窄窗口 pilot v2

这是 Issue #130 的三态 CEP diagnostic candidate。它只验证局部 CEP 区间、两层
`rho-support cascade` 证据和 solver 总成本；不覆盖 `data/reference/pnjl`，不晋升
reference，不接入 production、全温区 phase-reference 或 transport。旧 canonical
只作为历史叠图和差异说明，不参与自动正确性 gate。

## Provenance

- 三态合同基线：`main@0603e4df656ac610e12d51f9c810ec93be6e8f14`。
- 数值 calculation SHA：`9d325b62ef366f1de70cd3b2a3555e71bd12664a`。
- 数值 source run：[#30339094251](https://github.com/w5851/Julia_RelaxTime/actions/runs/30339094251)。
- collector/replay 修复 SHA：`d2f8b1d0c2fb7a6cec7d21a1076a1d24639c9ca5`。
- aggregate replay：[#30342409727](https://github.com/w5851/Julia_RelaxTime/actions/runs/30342409727)。
- replay 只重新收集 source run 的数值 artifacts；没有重算任何 PNJL equilibrium 点。

`manifest.json` 保留 source run 的 calculation SHA、文件哈希、gate 和 Actions
成本；`validation_windows.json` 是冻结的 bounds/common-anchor 输入。完整的
`curve_points.csv`（174001 行）保留在 Actions artifact/本地审计目录，仓库只导入
`curve_index.csv`、代表性曲线图和可追溯的 source hash，以避免把原始曲线副本提交进仓库。

## 自动 gate（不替代作者物理判断）

- 总 verdict：`pass`。
- oracle：`stable`；三层 xi 均有两层 rho 的双端证据，端点分类和低端 Maxwell
  `muq` 满足预声明收敛阈值。
- cascade：`within_oracle`；公共非 ambiguous anchors 分类一致，区间相交且没有
  排除 oracle 端点超过一个 `0.125 MeV` resolution step。
- fallback/retry 风险：未触发；cascade `0.001256`，dense `0.001617`。

`pass` 仅表示本 pilot 的数值/成本/provenance contract 通过，不表示 CEP 已被解析、
不替代作者物理审核，也不构成 production/reference 晋升许可。`cep_accuracy.csv`
中的 9 个方法/xi 行仍是 `ambiguous`，而 `slice_metrics.csv` 的切片级证据包含
`confirmed_first_order`、`confirmed_monotone` 和 `ambiguous_near_critical` 三态。

## 区间摘要

下表首先给出 CSV 的组合 CEP 区间（`T_last_first_order_MeV` 到
`T_first_monotone_MeV`）；`fine_*` 两层 rho endpoint 是 oracle 收敛 gate 使用的
独立证据。两者不能混写：例如 `xi=0` 的组合区间为 `130.5869→131.0869 MeV`，
而两层 fine endpoint 为 `130.9619→131.0869 MeV`。

| xi | 方法 | 组合区间 [MeV] | fine 两层端点 [MeV] | 状态 |
|---:|---|---:|---:|---|
| -0.5 | cascade | 139.9697→147.2197 | 147.0947→147.2197 | ambiguous |
| -0.5 | dense/oracle | 147.0947→147.2197 | 147.0947→147.2197 | ambiguous |
| 0 | cascade | 130.5869→131.0869 | 130.9619→131.0869 | ambiguous |
| 0 | dense/oracle | 130.5869→131.0869 | 130.9619→131.0869 | ambiguous |
| 0.5 | cascade | 101.4600→107.0850 | 106.9600→107.0850 | ambiguous |
| 0.5 | dense/oracle | 106.9600→107.0850 | 106.9600→107.0850 | ambiguous |

实际 ambiguity interval 可以大于 `temperature_resolution_target_MeV=0.125 MeV`；
这正是三态合同不把 near-critical 强制改写成 monotone 的边界。

## Solver 成本

| 方法 | unique solves | residual | Jacobian | Newton | runner s |
|---|---:|---:|---:|---:|---:|
| rho-support cascade | 10,353 | 112,439 | 90,632 | 36,156 | 343.526 |
| C2 dense baseline | 35,255 | 377,225 | 302,222 | 120,429 | 1,008.904 |
| high-resolution oracle | 70,455 | 754,168 | 604,150 | 240,728 | 1,681.359 |

Actions 聚合关键路径为 `1247 s`，原始 job 时间为 `5173 s`，向上取整
runner-minutes 为 `91`。cascade 的 solver 工作量、runner 时间和 fallback/retry
率均不高于 memoized dense baseline；分类层微秒耗时不是本 pilot blocker。

## 图和后续边界

- `figures/cep_intervals.png`：三态区间及两层 endpoint。
- `figures/slice_states.png`：切片三态分类。
- `figures/cost_frontier.png`：solver cost frontier。
- `figures/rho_mu_xi_*.png`：每个 xi 的代表性 rho–mu 曲线叠图。
- `curve_index.csv`：每个 xi/方法的最后 Maxwell、代表性 ambiguous、首个
  monotone 曲线及 source artifact/hash 定位。

当前停在作者对 CEP 区间、原始曲线和性能证据的物理审核。通过后才另开
production integration PR；随后才允许 phase-reference replay、reference 审核和
transport。
