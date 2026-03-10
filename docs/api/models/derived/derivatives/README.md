# Derivatives 主题 API

本目录是 `Models` 衍生量层中 derivatives 主题的主入口，面向热力学导数、质量导数与体粘滞组合导数相关公开能力的统一治理。

推荐阅读顺序：

1. [Overview.md](Overview.md)：先看 `thermo_derivatives`、`mass_derivatives` 与 `bulk_derivative_coeffs` 的最短使用路径
2. [CoreConcepts.md](CoreConcepts.md)：理解隐函数求导、链式法则、AD 口径与兼容层边界
3. [MassDerivatives.md](MassDerivatives.md)：质量导数结果结构与使用场景
4. [ThermoDerivatives.md](ThermoDerivatives.md)：热力学导数主接口与高阶 `dP_dT` / `dP_dmu`
5. [BulkViscosityDerivatives.md](BulkViscosityDerivatives.md)：体粘滞组合导数、B bracket 与兼容入口
6. [generated/Exports.md](generated/Exports.md)：自动生成的公开导出全集

本主题覆盖的 `Models` 公开表面包括：

- `mass_derivatives`
- `thermo_derivatives`
- `bulk_derivative_coeffs`
- `dP_dT`
- `dP_dmu`
- `bulk_viscosity_coefficients`
- `compute_B_bracket`
- `legacy_transport_c_p`

说明：

- `legacy_transport_c_p` 仅作为 compatibility / legacy 导出保留，不应作为 derivatives 主题首页首选入口；
- 旧 `pnjl` derivatives 兼容页不再承担本主题主说明职责。