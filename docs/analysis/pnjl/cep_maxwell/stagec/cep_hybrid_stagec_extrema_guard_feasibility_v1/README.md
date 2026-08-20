# PNJL Stage-C discrete-extrema guard feasibility v1

verdict: `feasible_candidate`。

本目录是对 source run `30737739707`（calculation SHA
`467be1fce847a9c991ec362c3335be07fccbe604`）的 solver-free replay。Stage-B 使用完整
`0.00625` 全域曲线；Stage-C 只从 Stage-B 派生的特征排序中选择 guard 内的
`0.003125` 已有曲线点。没有调用 equilibrium solver，也没有修改 production、
reference 或历史 v1/v2 evidence。

Guard 定义为两个 μ 极值外侧的首个严格采样点：相等点跳过，不插值、不二分、
不使用固定 padding。缺少任一外侧点、多 S 形或 topology 不稳定时保持
`ambiguous_near_critical`。

- tested caps: `12, 24, 48, 96, 160`
- dense unique-solve reference: `15384.0`
- selected cap: `12`
- solver called: `false`
- five-point revalidation run: `30730990835` (`recorded`)

`author_adjudication` 只作为显式 provenance 使用，不把其它视觉判断升级为自动
证书。即使 feasibility gate 通过，它也只授权下一步 production focused CI，
不等于 shadow、reference 或正式 production 已通过。
