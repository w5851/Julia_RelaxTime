# PNJL Maxwell endpoint-limit contract v1

本分析包固定读取 GitHub Actions run `30980094983` 的单点 endpoint refinement
artifact，并在不调用 equilibrium/Maxwell solver 的前提下形成端点极限证书。

## 结论

- source artifact 原 verdict 保持 `candidate_only_endpoint_inconclusive`，没有改写历史产物。
- 13 个 refinement levels 均保持唯一三交点 Maxwell candidate。
- 左侧 bracket 从 `[0, 0.003125]` 确定性收缩至
  `[0, 7.62939453125e-07]`。
- 最终 `mu_M=331.5739309844038 MeV`，
  `rho_q=2.8043699741472983`，面积残差
  `-2.9237516382062305e-07`；geometry 与 solver-work gate 均通过。
- 作者已确认该有界端点区间足以作为一阶相变证据。因此 derived verdict 为
  `endpoint_limited_first_order_candidate`，三态分类为 `confirmed_first_order`。

## 合同语义

`endpoint_limited_first_order` 不声称存在严格正的 `rho_h` 下界。兼容 scalar 字段在后续
production 实现中应以 `rho_hadron=0.0` 表示端点极限，同时把真实区间
`[0, rho_hadron_upper_bound]` 写入 diagnostics/CSV。该证书不增加第四种 CEP 物理状态；
最终仍映射到 `confirmed_first_order`。

## 证据边界

本包只覆盖 `(xi,T)=(-0.5,5 MeV)`，不构成 CEP resolved、24-anchor shadow 通过、
phase-reference 晋升或 transport 输入。公共 Maxwell 核心和 hybrid production 尚未修改；
完整原始 `curve_points.csv` 仅保留在 Actions/local artifact，通过 provenance/hash 追溯。
