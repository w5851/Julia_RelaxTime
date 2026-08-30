# Issue #130 phase-reference layers v1 审计

- numerical Maxwell source run: `32354095831`；replay run: `32451053476`。
- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`。
- strict Maxwell rows: `7162`（C2/v6 原生行加 276 个真实补点）。
- derived rows: Maxwell `12537`、crossover `3135`、spinodal `11989`、CEP `161`。
- replay 是 `aggregate_replay`，已验证 276/276 materialized、无错误且 `solver_called=false`；replay manifest 的 source run 由 workflow 输入显式传递并保存在 `replay_provenance.json`。
- strict 层保留旧 unresolved geometry；derived 层只在相邻 xi 的共同 axis support 内局部插值，非原生行标记为 `interpolated_noncertified`，不跨缺失 support 外推。
- render 层使用 ordered xi-slice lines，不三角化跨越 masked cell；图像连续性不改变数据层状态。

## Verdict

这些层是 `diagnostic/author-review candidate`，不是 phase-reference promotion。在作者审核和独立 promotion gate 完成前，不写入 `data/reference`，不启动 RS transport。
