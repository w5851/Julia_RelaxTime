# Algorithmic Feasibility Evidence

本目录收纳独立的 PNJL 临界判据算法可行性审计。它使用解析 cusp 曲线族验证 secant、rho-support cascade、拟合稳定性、CEP 外推和性能记录的合同；不代表 PNJL 数值结果，也不属于 phase-reference promotion 证据。

## Evidence Package

| 逻辑角色 | 当前路径 | verdict / 用途 |
| --- | --- | --- |
| criticality feasibility v1 | [`criticality_feasibility_v1/`](criticality_feasibility_v1/) | `diagnostic_candidate`；受控解析反例上的算法合同审计 |

## Boundary

- 解析 cusp 真值只用于 feasibility 和 deterministic tests，不作为 PNJL CEP 或 phase production 结论。
- `performance_rejected` 保持为性能审计 verdict；不得外推为真实 PNJL solver 成本或 production replacement。
- 本次 namespace 迁移只改变目录路径；包内 README、manifest 和表格保持生成时 provenance，不重算或重绘证据。
