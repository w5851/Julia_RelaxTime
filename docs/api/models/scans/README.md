# Models 扫描主题 API

本主题从 `Models` 统一入口视角组织扫描相关 API，而不是继续把 `PNJL.*` 兼容层页面作为主导航。

适用对象：

- 需要批量生成 T-μ 或 T-ρ 网格结果的调用方
- 需要决定采样模板、扫描顺序与输出契约的维护者
- 需要把扫描结果继续接到相图、分析或回归链路的开发者

推荐阅读顺序：

1. [Overview.md](Overview.md)：先判断应该使用哪一个 `Models` 入口
2. [Algorithms.md](Algorithms.md)：理解采样策略、扫描顺序与契约边界
3. [TmuScan.md](TmuScan.md)：T-μ 扫描入口与参数口径
4. [TrhoScan.md](TrhoScan.md)：T-ρ 扫描入口与参数口径
5. [FreezeoutScan.md](FreezeoutScan.md)：freeze-out 路径扫描入口与路径契约
6. [FreezeoutMesonDensityScan.md](FreezeoutMesonDensityScan.md)：freeze-out 路径上的介子数密度 workflow
7. [SamplingGrid.md](SamplingGrid.md)：默认密度网格与加密策略
8. [generated/Exports.md](generated/Exports.md)：公开导出全集与覆盖检查

本主题覆盖的 `Models` 公开导出为：

- `Models.run_tmu_scan`
- `Models.run_trho_scan`
- `Models.run_freezeout_fixedmu_scan`
- `Models.run_freezeout_meson_density_scan`
- `Models.build_default_rho_grid`

本主题已直接吸收以下内容，不再依赖旧 `pnjl` 页面承载主体说明：

- T-μ / T-ρ 扫描的关键字参数口径
- CSV 输出契约与字段边界
- smoke / standard / phase-diagram 三档采样模板
- `reverse_rho=true`、`use_phase_aware=true` 等推荐与禁用组合

旧 `docs/api/pnjl/` 扫描页面后续仅保留迁移说明与历史兼容定位。
