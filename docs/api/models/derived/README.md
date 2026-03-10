# Models 衍生量主题

本目录承接 `Models` 公开表面中不适合直接平铺到 `docs/api/models/` 一级目录的“衍生量/响应量”主题。

当前已落地主题：

- [susceptibility/README.md](susceptibility/README.md)：守恒荷 susceptibility、cumulant 与 baryon 组合量主题。
- [derivatives/README.md](derivatives/README.md)：热力学导数、质量导数与体粘滞组合导数主题。

设计原则：

- 衍生量主题优先从 `Models` 聚合导出视角组织，而不是继续沿用旧 `PNJL.*` 页面做主导航。
- 每个主题仍遵循三层视图：面向用户入口、职责核心、导出 API 全集。
- 旧 `pnjl` derivatives 兼容页已退出主导航；susceptibility 与 derivatives 均以当前主题目录为准。

当前主题分工：

- `susceptibility/`：守恒荷 susceptibility、cumulant、`baryon_*` 组合量
- `derivatives/`：`thermo_derivatives`、`mass_derivatives`、`bulk_derivative_coeffs`、`dP_dT`、`dP_dmu`、`bulk_viscosity_coefficients`、`compute_B_bracket`

说明：

- `legacy_transport_c_p` 属于 derivatives 主题中的兼容层导出，不应被误读为首选用户入口。