# Models 工作流主题 API

本目录是 `Models` 统一入口下工作流能力的主题主目录，面向 transport workflow、meson workflow、模块访问器与参数适配层的统一治理。

## 阅读顺序

- 用户入口：见 [Overview.md](docs/api/models/workflows/Overview.md)
- 职责核心：见 [CoreConcepts.md](docs/api/models/workflows/CoreConcepts.md)
- 导出 API 全集：见 [generated/Exports.md](docs/api/models/workflows/generated/Exports.md)

## 覆盖范围

本主题覆盖以下公开入口：

- `Models.solve_gap_and_transport`
- `Models.solve_transport_from_equilibrium`
- `Models.solve_gap_and_meson_point`
- `Models.solve_meson_density_from_meson_point`
- `Models.solve_gap_and_meson_density_point`
- `Models.transport_workflow_module`
- `Models.meson_workflow_module`
- `Models.meson_density_workflow_module`
- `Models.workflow_param_adapters_module`
- `Models.pnjl_module`

## 目录职责

- [Overview.md](docs/api/models/workflows/Overview.md)：给首次使用者的统一入口页
- [CoreConcepts.md](docs/api/models/workflows/CoreConcepts.md)：解释 workflow 输入分层、职责边界与模块访问器定位
- [generated/Exports.md](docs/api/models/workflows/generated/Exports.md)：自动生成的公开导出索引

## 旧领域文档的关系

- [docs/api/relaxtime/workflow/TransportWorkflow.md](docs/api/relaxtime/workflow/TransportWorkflow.md) 保留为 transport 领域细节页
- [docs/api/pnjl/MesonMassWorkflow.md](docs/api/pnjl/MesonMassWorkflow.md) 保留为 meson 领域细节页
- [docs/api/relaxtime/workflow/MesonDensityWorkflow.md](docs/api/relaxtime/workflow/MesonDensityWorkflow.md) 保留为介子数密度 workflow 领域细节页

本目录不试图替代这些领域文档，而是从 `Models` 统一入口视角重新组织它们。
