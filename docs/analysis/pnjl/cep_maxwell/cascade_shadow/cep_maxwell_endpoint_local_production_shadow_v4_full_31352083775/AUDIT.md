# Full endpoint-local shadow audit

输入是 Actions aggregate replay 的 final evidence。校验结果：

- slice rows: `72`；raw curve rows（仅用于导入前审计）: `54659`；
  raw key 唯一且所有点 finite/converged；
- 派生表中的 `NaN` 仅出现在 schema 明确允许的“不适用字段”（例如单调切片没有
  Maxwell 端点）；所有 `Inf` 及未声明的 NaN 均拒绝；
- 三个 endpoint certificates 的 support envelope 均覆盖初始 left bracket 与 anchor；
- 代表图的右侧面板使用独立的 `rho–mu` 局部纵轴，并标出 Maxwell/spinodal/coexistence
  位置；一阶图的 rho 窗口遵循 `max(4% of phase/support envelope, 0.01 rho)`，
  mu 轴留白遵循 `max(8% of local mu span, 0.0002 MeV)`；这些图层规则只改变
  后处理，不改变原始曲线或数值 gate；
- 没有有限 Maxwell/spinodal/support 的 monotone 切片使用原始 production 曲线最长低斜率
  区间作为显示窗口；该窗口不参与物理分类或 gate；
- gate verdict: `full_hybrid_candidate`；所有 oracle/classification/
  endpoint/coverage/performance errors 为空；
- 完整 raw curve 不提交仓库，外部 SHA: `bec96fd7c9a1887a63729355db47b55632126c1a85a3c8ec8a8297a8e34d8c63`。

边界：这是 diagnostic candidate，不自动晋升 production/reference，不启动 C0/C1/C2、
phase-reference 或 transport。作者需审阅 figures/ 和 tables/claim_ledger.csv 后再决定。
