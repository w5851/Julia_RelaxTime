# Cumulant 与 Baryon 组合量接口

本页说明 susceptibility 主题中与 cumulant 和组合量直接相关的公开接口。

## 主要导出

- `cumulant_BQS`
- `cumulant_B`
- `baryon_Ssigma`
- `baryon_kappa_sigma2`

## `cumulant_BQS`

统一的守恒荷 cumulant 入口：

- `cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(i, j, k))`

当前采用：

- `C_{ijk}^{BQS} = V * T^3 * chi_{ijk}^{BQS}`

因此它与 `chi_BQS` 的关系是主题级核心合同，而不是单独一条实现细节。

## `cumulant_B`

这是 baryon 方向的常用快捷入口：

- `cumulant_B(T_fm, muB_fm, V; order=n)`

当你只关心 `C_n^B` 时，优先使用它，而不是手动拼 `orders=(n, 0, 0)`。

## `baryon_Ssigma`

返回：

- `Ssigma = chi_3^B / chi_2^B`

它直接表达了 baryon 方向的常见实验/表型组合量，因此在主题导航中应与 `chi_n^B`、`C_n^B` 放在一起理解。

## `baryon_kappa_sigma2`

返回：

- `kappa_sigma2 = chi_4^B / chi_2^B`

与 `baryon_Ssigma` 一样，这个接口更偏结果读取，而不是底层导数构造。

## 使用建议

- 需要体积依赖 cumulant 时，优先用 `cumulant_BQS` 或 `cumulant_B`
- 需要无量纲比值量时，优先用 `baryon_Ssigma` 与 `baryon_kappa_sigma2`
- 若结果对 `chi_2^B` 接近零的区域敏感，应在脚本侧额外检查数值稳定性