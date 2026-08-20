# Shadow audit

## 输入与 provenance

- source workflow：`30527097410`
- calculation/head SHA：`d1e386deb950169b6d5db9cc3dd0617845ee60e8`
- postprocess replay：由 `d6ffa6b172d022157bea1ebb5118d5c500528dfc` 的
  collector/plotter 生成
- source workflow 最终状态：`completed/success`
- reference write：`false`

`manifest.json` 中的 `aggregation_corrections` 保留历史首-anchor cache 计数与
全 anchor 归一化计数的 raw/normalized 对照；完整原始曲线不提交仓库，仅由
`curve_index.csv` 和 source artifact hash 追溯。

原始九个 job artifact 未被改写；聚合 `method_costs.csv` 对历史 runner 的首-anchor
cache 计数进行了可追溯重建，并在 `manifest.json` 的
`aggregation_corrections` 中保留 raw/normalized 对照。

## 数值与物理边界

- 所有 24 个 oracle anchor 的最终 rho–mu 曲线 finite/converged，且无重复曲线键。
- cascade 相比 memoized dense 明显减少 fixed-rho、residual/Jacobian 和 runner work。
- `xi=0.5,T=5` 的 cascade confirmed-first-order 没有 oracle 同类确定证书，属于
  unsupported confirmation。
- `xi=-0.5` 的 cascade 缺少 first-order/monotone 双端证据，不能形成可靠 CEP bracket。

因此 verdict 为 `hybrid_required`。这不是 workflow failure，也不是 reference gate；
下一阶段需要三级 hybrid shadow，并继续由作者审核代表性 rho–mu 曲线。
