# GasLiquid 变体主题 API

本目录是 `Models` 变体层中 gas_liquid 主题的主入口，面向气液（RMF/Walecka）变体模型与单点 workflow 的公开能力。

推荐阅读顺序：

1. [Overview.md](Overview.md)：先看 `GasLiquidModel`、`solve_gas_liquid_point`、RMF 单点和扫描入口
2. [CoreConcepts.md](CoreConcepts.md)：理解模型适配层与气液核心方程层职责边界
3. [generated/Exports.md](generated/Exports.md)：自动生成的公开导出全集

本主题覆盖的 `Models` 公开表面包括：

- `GasLiquidModel`
- `solve_gas_liquid_point`
- `solve_gas_liquid_rmf_point`
- `run_gas_liquid_tmu_scan` / `run_gas_liquid_trho_scan`
- `build_gas_liquid_result_row` / `build_gas_liquid_manifest`

说明：

- 气液核心实现位于 `src/models/variants/gas_liquid/core/`，该层面向实现维护者。
- 面向调用方的首入口应优先通过 `Models` 聚合导出访问。
