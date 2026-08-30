# C2 convergence audit v1

这是 C0/C1/C2 的 solver-free 后处理审计包。输入固定为 Actions 已完成 artifact；脚本不调用 Julia equilibrium solver，也不修改数值产物。

## 结论

- C1→C2 的主 verdict 为 `classification_regression`。
- 9 个 public first-order anchor 在 C2 退为 ambiguous；它们是分类回归，不被解释为 monotone。
- 16 个 CEP bracket 超过 `0.1 MeV`，并单独记录 endpoint search resolution 与实际 ambiguity width。
- crossover 只在双方共同的 monotone/crossover 物理温区、物理 `mu` 并集上插值比较；剩余风险集中在 `xi=0.2875` 的高 `mu` 局部点。
- C2 保持 `diagnostic-only`，不创建 production PR，不启动新的 C0/C1/C2、reference 或 transport。

## 输入与追溯

完整曲线和全量逐点比较留在外部 Actions/local artifact；本目录只保留 anchor/geometry/CEP 表、comparison summary、例外行和 `xi=0.2875` crossover 绘图子集。全量比较仍参与 replay verdict，但不作为仓库内原始曲线替代物。所有输入文件 SHA、calculation SHA、postprocess SHA、source run 和 generator SHA 见 `manifest.json`。

## 作者检查项

后续 feasibility 必须在固定 cap=12 下验证 density certificate、CEP bracket 和 crossover refinement；不得通过增加 cap 或放宽现有容差掩盖上述回归。
