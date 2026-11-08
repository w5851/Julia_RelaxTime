# Algorithmic Feasibility Evidence

本目录收纳 PNJL 算法可行性分析。`criticality_feasibility_v1` 使用解析 cusp 曲线族审计临界判据；
`palc_matched_20261108_analysis` 使用真实 PNJL 小窗口比较延拓与种子扫描，两者各自保留适用范围。

## Evidence Package

| 逻辑角色 | 当前路径 | verdict / 用途 |
| --- | --- | --- |
| criticality feasibility v1 | [`criticality_feasibility_v1/`](criticality_feasibility_v1/) | `diagnostic_candidate`；受控解析反例上的算法合同审计 |
| PALC matched pilot | [`palc_matched_20261108_analysis/`](palc_matched_20261108_analysis/) | 两窗口同密度点、cold/warm 与工作量对照；隔离诊断及候选根选择反例 |

## Boundary

- 解析 cusp 真值只用于 feasibility 和 deterministic tests，不作为 PNJL CEP 或 phase production 结论。
- `performance_rejected` 保持为性能审计 verdict；不得外推为真实 PNJL solver 成本或 production replacement。
- 本次 namespace 迁移只改变目录路径；包内 README、manifest 和表格保持生成时 provenance，不重算或重绘证据。
