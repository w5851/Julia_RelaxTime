# PNJL CEP 窄窗口 pilot v1

这是一个候选诊断产物，不覆盖 `data/reference/pnjl`，也不代表 reference 晋升。
每个矩阵 job 固定同一 calculation SHA、`p_num=24`、`t_num=8`、
`rs_reduced_adaptive` 和 ground-state pressure governance；完整数值计算只在
GitHub Actions 执行。

## 自动 gate（作者物理判断仍然必需）

- 总 verdict：`diagnostic_only`
- oracle verdict：`oracle_inconclusive`
- cascade verdict：`outside_oracle`
- 诊断信息：`oracle refinement evidence is missing；oracle 0.00625 -> 0.003125 classification is not stable；oracle does not satisfy the predeclared canonical CEP tolerance；cascade CEP is outside the oracle bracket/uncertainty`
- Actions critical path：625.0 s
- Actions runner-minutes（向上取整）：64

`method_costs.csv` 区分 equilibrium request、unique solve、memo cache hit、
residual/Jacobian、Newton/trust-region、fallback/rescue/retry；`slice_metrics.csv`
保留每个 `(T, xi, method)` 切片的局部计数。`claim_ledger.csv` 明确观察、gate
与不可推出的物理结论边界。
