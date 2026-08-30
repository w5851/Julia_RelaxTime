# PNJL Hybrid Stage-C certificate feasibility v2

verdict: `oracle_inconclusive`。本目录使用 source run `30601857594` 的 final aggregate
artifact 做 solver-free replay；没有调用 equilibrium solver，不修改 production、reference 或历史 v1 evidence。

v2 直接复用 Julia `PhaseCore`、Maxwell construction 和 phase geometry comparison。Stage-C
路由使用 Stage-B 曲线的 drawdown、density support、turning-point 和 Maxwell 特征；oracle
状态和 geometry 只用于比较与 gate，不能选择补点。固定 slope margin 不再作为弱 S 的必要条件。

- tested caps: `12, 24, 48, 96, 160`
- dense unique-solve reference: `15384.0`
- selected cap: `nothing`
- solver called: `false`

`author_adjudication` 只适用于两个作者确认的一阶点；高温 monotone 点标记为
`three_method_consensus`。如果 verdict 不是 `feasible_candidate`，本结果不授权 production
集成或 reference 重放。
