# relaxtime API 主题总览

本目录承接 `relaxtime` 领域链路的 API 文档，重点覆盖：

- transport：输运系数、弛豫时间、平均散射率与 workflow 细节
- scattering：散射幅、微分截面、总截面
- propagator：介子传播子、Mott 阈值与有效耦合
- polarization：极化函数与缓存
- meson_density：稳定粒子极限数密度与 `K/\pi` 比值主线
- meson_thermo：介子 pressure 与最小 EOS 组合量

## 推荐阅读入口

如果你的目标是“直接算输运系数”或“理解 transport provider 契约”，优先从以下页面开始：

- `transport/README.md`
- `transport/Overview.md`
- `transport/CoreConcepts.md`
- `transport/generated/Exports.md`
- `workflow/TransportWorkflow.md`
- `../../reference/formula/relaxtime/transport/Transport_EndToEnd_Pipeline.md`
- `workflow/MesonDensityWorkflow.md`
- `meson_density/MesonDensity.md`
- `meson_density/BUPhaseGates.md`（strict BU 测度、高能相位锚点与 Levinson/Mott 门禁）
- `meson_density/PhaseNormalization.md`（散射相位、S-matrix 与 BU 测度的纯代数闭合）
- `meson_density/ChargedPhaseBackend.md`（strict charged phase/BU 诊断后端与合成路径合同）
- `meson_thermo/MesonThermodynamics.md`
- `propagator/MesonInteractionKernel.md`（Phase 1.5/2 完整 KMT 相互作用核与中性 RPA 代数后端）
- `propagator/MesonRPA.md`（中性 `(0,3,8)` RPA 矩阵代数）
- `propagator/MesonRPAAdapter.md`（Phase 3：`PolarizationAniso` 到中性 RPA 的诊断桥接）
- `propagator/ChargedRPAKernel.md`（Phase B：charged RPA 归一化/味道顺序契约，diagnostic）
- `propagator/ChargedRPAProvider.md`（Phase C：charged `Pi_us/Pi_su` provider 适配器，diagnostic）

如果你的目标是“从统一 workflow 直接生成 `n_pi(T)` / `n_K(T)` / `K/\pi(T)` 的扫描输出”，优先结合以下页面与脚本：

- `../models/workflows/MesonDensityWorkflow.md`
- `workflow/MesonDensityWorkflow.md`
- `meson_density/MesonDensity.md`
- `scripts/relaxtime/run_meson_density_scan.jl`

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
- `workflow/MesonDensityWorkflow.md`
- `meson_thermo/MesonThermodynamics.md`
- `propagator/MesonInteractionKernel.md`
- `propagator/MesonRPA.md`
- `propagator/MesonRPAAdapter.md`
- `propagator/ChargedRPAKernel.md`
- `propagator/ChargedRPAProvider.md`
- `meson_density/BUPhaseGates.md`
- `meson_density/PhaseNormalization.md`
- `meson_density/ChargedPhaseBackend.md`
- `scattering/*`
- `propagator/*`
- `polarization/*`

## 目录迁移原则

本目录现在不再只是预留命名空间。后续若继续把旧 `docs/api/*.md` 中的 relaxtime 相关页面迁入本目录，仍遵循：

- 新路径先建立稳定主题入口；
- 旧路径保留兼容跳转说明；
- 不为维持旧结构而阻止新主题成形。
