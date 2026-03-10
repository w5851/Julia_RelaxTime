# relaxtime API 主题总览

本目录承接 `relaxtime` 领域链路的 API 文档，重点覆盖：

- transport：输运系数、弛豫时间、平均散射率与 workflow 细节
- scattering：散射幅、微分截面、总截面
- propagator：介子传播子、Mott 阈值与有效耦合
- polarization：极化函数与缓存

## 推荐阅读入口

如果你的目标是“直接算输运系数”或“理解 transport provider 契约”，优先从以下页面开始：

- `transport/README.md`
- `transport/Overview.md`
- `transport/CoreConcepts.md`
- `transport/generated/Exports.md`
- `workflow/TransportWorkflow.md`

如果你是从 `Models.solve_gap_and_transport` 之类的统一入口进入，建议按以下顺序阅读：

1. `docs/api/models/workflows/TransportWorkflow.md`
2. `workflow/TransportWorkflow.md`
3. `transport/Overview.md`
4. `transport/CoreConcepts.md`

## 当前已落地的 relaxtime 子主题

- `ParticleSymbols.md`
- `transport/TransportCoefficients.md`
- `transport/RelaxationTime.md`
- `transport/AverageScatteringRate.md`
- `workflow/TransportWorkflow.md`
- `scattering/*`
- `propagator/*`
- `polarization/*`

## 目录迁移原则

本目录现在不再只是预留命名空间。后续若继续把旧 `docs/api/*.md` 中的 relaxtime 相关页面迁入本目录，仍遵循：

- 新路径先建立稳定主题入口；
- 旧路径保留兼容跳转说明；
- 不为维持旧结构而阻止新主题成形。
