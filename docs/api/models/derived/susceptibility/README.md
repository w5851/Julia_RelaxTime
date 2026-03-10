# Susceptibility 主题 API

本目录是 `Models` 衍生量层中 susceptibility 主题的主入口，面向守恒荷静态涨落、累积量与 baryon 组合量相关公开能力的统一治理。

推荐阅读顺序：

1. [Overview.md](Overview.md)：先看何时使用 `chi_*`、`chi_BQS` 与 `cumulant_*`
2. [CoreConcepts.md](CoreConcepts.md)：理解 flavor 化学势、BQS 变换与 `T` 缩放口径
3. [ConservedChargeSusceptibilities.md](ConservedChargeSusceptibilities.md)：守恒荷 susceptibility 主接口
4. [Cumulants.md](Cumulants.md)：累积量与 `baryon_*` 组合量
5. [generated/Exports.md](generated/Exports.md)：自动生成的公开导出全集

本主题覆盖的 `Models` 公开表面包括：

- `conserved_charge_susceptibility`
- `chi_BQS`
- `chi_B`, `chi1_B`, `chi2_B`, `chi3_B`, `chi4_B`
- `chi_Q`, `chi1_Q`, `chi2_Q`, `chi3_Q`, `chi4_Q`
- `chi_S`, `chi1_S`, `chi2_S`, `chi3_S`, `chi4_S`
- `chi11_BQ`, `chi11_BS`, `chi11_QS`
- `cumulant_B`
- `cumulant_BQS`
- `baryon_Ssigma`
- `baryon_kappa_sigma2`

以下导出暂不纳入本主题主体说明，而留给后续 derivatives 子主题：

- `flavor_pressure_derivatives`
- `mass_derivatives`
- `bulk_derivative_coeffs`
- `dP_dT`
- `dP_dmu`
- `thermo_derivatives`
- `compute_B_bracket`