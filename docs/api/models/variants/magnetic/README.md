# Magnetic 变体主题 API

本目录是 `Models` 变体层中 magnetic 主题的主入口，面向外磁场 PNJL 与 Landau 能级离散化相关公开能力的统一治理。

推荐阅读顺序：

1. [Overview.md](Overview.md)：先看 `PNJLMagneticModel`、配置对象与最短使用路径
2. [CoreConcepts.md](CoreConcepts.md)：理解模型变体定位、Landau 离散化与 `n_max` 收敛边界
3. [ModelAndConfig.md](ModelAndConfig.md)：模型与配置合同
4. [LandauIntegrals.md](LandauIntegrals.md)：Landau 能级与积分相关接口
5. [Thermodynamics.md](Thermodynamics.md)：磁场热力学主接口与收敛报告
6. [磁场扫描 API 细账](../../scans/MagneticScan.md)：完整五维 `FixedMu` 的 `(T,mu,eB)` 生产扫描
7. [generated/Exports.md](generated/Exports.md)：自动生成的公开导出全集

本主题覆盖的 `Models` 公开表面包括：

- `PNJLMagneticModel`
- `MagneticIMCParams`
- `MagneticConfig`
- `default_imc_params`
- `default_magnetic_config`
- `alpha_n`
- `energy_landau`
- `smooth_cutoff`
- `resolve_nmax_from_cutoff`
- `resolve_magnetic_nmax`
- `omega0_flavor_landau`
- `omegat_flavor_landau`
- `density_flavor_landau`
- `coupling_GB`
- `calculate_magnetic_omega_components`
- `calculate_magnetic_omega`
- `calculate_magnetic_pressure`
- `calculate_magnetic_rho`
- `calculate_magnetic_number_densities`
- `magnetic_nmax_convergence_report`
- `Models.run_magnetic_scan`

旧 `pnjl` magnetic 兼容页已退出主说明路径；当前主题页是唯一现行入口。
