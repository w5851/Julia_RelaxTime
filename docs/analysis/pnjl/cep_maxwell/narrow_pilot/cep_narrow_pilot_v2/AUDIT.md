# PNJL CEP 窄窗口 pilot v2 审计

## 1. 范围与 provenance

本包审计 `cep_narrow_pilot_v2` 的三阶段 Actions 数值证据：3 个 cascade
discovery、冻结 validation windows、6 个独立 dense/oracle validation job，以及
aggregate replay。完整数值只在 Actions 执行；本地只下载、收集、绘图和审计。

- calculation SHA：`9d325b62ef366f1de70cd3b2a3555e71bd12664a`。
- source numerical run：[#30339094251](https://github.com/w5851/Julia_RelaxTime/actions/runs/30339094251)，
  14 个 job，`completed/success`。
- replay collector SHA：`d2f8b1d0c2fb7a6cec7d21a1076a1d24639c9ca5`。
- aggregate replay：[#30342409727](https://github.com/w5851/Julia_RelaxTime/actions/runs/30342409727)，
  `completed/success`；只下载 source numerical artifacts，没有重新求解。
- source artifact curve hash：
  `curve_points.csv` = `3294e66227e4a1b36ffad8abaf254035dcae613805a3e0939ea6ab0c0191852f`。
- validation windows 的 calculation SHA 与上述 source SHA 一致，label policy 为
  `bounds_and_common_anchors_only`，不传播 cascade 标签。

## 2. 完整性检查

- job matrix：9 个 key 唯一且齐全（3 xi × 1 discovery cascade + 2 validation
  methods）；所有 job schema 为 `cep_narrow_pilot_v2`，calculation SHA 一致。
- `curve_points.csv`：174001 行，所有 equilibrium 点 `finite=true`、`converged=true`；
  没有重复 `(xi, method, stage, T, rho_level, rho)` key。
- `slice_metrics.csv`：172 个唯一切片 key；最大 cascade targeted additions 为 9，
  满足每切片两层累计 ≤12。
- `method_costs.csv`、`cep_accuracy.csv`、`actions_costs.csv`：无重复 key，所有计数
  和耗时字段有限。
- CSV 中检查到的非有限值只出现在 118 个派生 `mu_last_first_order_MeV=NaN`
  诊断字段（monotone/geometry 不适用）；原始 equilibrium 曲线和最终收敛字段没有
  NaN/Inf。
- collector manifest 的 `workflow_contract_errors`、`diagnostic_errors`、
  `oracle_errors`、`cascade_errors` 和 `performance_errors` 均为空。

## 3. 三态证据

切片级计数为：54 `confirmed_first_order`、63 `confirmed_monotone`、55
`ambiguous_near_critical`。9 个 CEP 聚合行保留 `result_status=ambiguous`，这是合同
要求：在两条证据边界之间没有同一点的 resolved `(T, μ)`，不能把区间中点静默写成
CEP。

Oracle 的 fine 两层 endpoint（`0.00625→0.003125`）为：

| xi | fine last first-order | fine first monotone | fine ΔT |
|---:|---:|---:|---:|
| -0.5 | 147.0947266 | 147.2197266 | 0.125 |
| 0 | 130.9619141 | 131.0869141 | 0.125 |
| 0.5 | 106.9599609 | 107.0849609 | 0.125 |

`cep_accuracy.csv` 的组合字段同时保留所有搜索证据。例如 `xi=0` 的
`T_last_first_order_MeV→T_first_monotone_MeV` 为 `130.5869141→131.0869141 MeV`，
而 fine 两层 gate endpoint 为 `130.9619141→131.0869141 MeV`。这两个字段分别表示
组合搜索边界和 rho-resolution 证书，不能互相替换。

## 4. 自动 gate

- `status=pass`、`oracle_status=stable`：三 xi 都有双端证据；coarse/fine 分类一致，
  endpoint 温差和低端 Maxwell `muq` 差满足阈值。
- `cascade_status=within_oracle`：公共非 ambiguous anchors 分类一致，ambiguity
  interval 相交，oracle 两端没有被排除超过一个 `0.125 MeV` step。
- `fallback/retry` rate：cascade `0.0012556747`，dense `0.0016167919`；没有触发
  25% 相对恶化风险。

这些 gate 只验证本 pilot 的数值证据、成本和 provenance；不证明旧 canonical 或新
CEP 物理解的优劣，也不授权 production/reference promotion。

## 5. Solver/Actions 成本

| 方法 | solves | residual | Jacobian | Newton | fallback | runner s |
|---|---:|---:|---:|---:|---:|---:|
| rho-support cascade | 10353 | 112439 | 90632 | 36156 | 13 | 343.526 |
| C2 dense baseline | 35255 | 377225 | 302222 | 120429 | 57 | 1008.904 |
| high-resolution oracle | 70455 | 754168 | 604150 | 240728 | 110 | 1681.359 |

Actions 关键路径 `1247 s`，raw total `5173 s`，runner-minutes `91`。三种方法都
使用独立 exact `(T, xi, rho)` memoization；cascade 使用 `0.05→0.025`，dense 使用
`0.0125→0.00625`，oracle 独立使用 `0.00625→0.003125`。因此成本差异没有被跨方法
共享解掩盖。

## 6. 证据边界与下一步

本包不导入 52 MB 原始 `curve_points.csv`，避免仓库膨胀；`curve_index.csv`、图像和
source hash 指向本地审计目录及 Actions artifact。作者审核时应同时检查原始 rho–mu
曲线、`ambiguous_near_critical` 的几何失败原因和 canonical 历史叠图。

当前 verdict 为 diagnostic candidate。只有作者确认区间和曲线物理解释后，才另开
production integration PR；随后才按既定顺序进行全温区 phase-reference replay、
reference 审核和 transport。任何自动 gate `pass` 都不能跳过这些审核。
