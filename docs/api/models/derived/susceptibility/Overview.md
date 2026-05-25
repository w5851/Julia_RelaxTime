# Susceptibility 主题总览

本页优先回答“如何在 `Models` 入口下使用守恒荷 susceptibility 与 cumulant 接口”，而不是先展开底层 AD 或 flavor 变换细节。

## 何时使用该主题

当你需要以下任一能力时，应从本主题开始：

- 在给定 `(T, mu_B, mu_Q, mu_S)` 下计算 `chi_n^B`、`chi_n^Q`、`chi_n^S`
- 在 baryon / charge / strangeness 方向直接读取四阶接口：`chi4_B`、`chi4_Q`、`chi4_S`
- 计算 mixed susceptibilities，例如 `chi11_BQ`、`chi11_BS`、`chi11_QS`，以及 `chi_BQS(...; orders=(2,1,0))` 这类高阶 mixed 组合
- 使用统一的 `chi_BQS` / `cumulant_BQS` 接口表达任意已支持的 BQS 组合
- 计算 `C_n^B`、`Ssigma`、`kappa_sigma2` 等 baryon 方向组合量

## 首选公开入口

主题主入口来自 `src/models/Models.jl` 对 `src/models/derivatives/ConservedChargeSusceptibilities.jl` 的聚合导出。

优先关注：

- `chi_B`
- `chi_BQS`
- `chi11_BQ`
- `chi4_B`
- `cumulant_B`
- `cumulant_BQS`
- `baryon_Ssigma`
- `baryon_kappa_sigma2`

完整导出基线见 [generated/Exports.md](generated/Exports.md)。

## 最短使用路径

```julia
using .Models

T_fm = 150.0 / Models.Constants_PNJL.ħc_MeV_fm
muB_fm = 0.0
muQ_fm = 0.0
muS_fm = 0.0
V_fm3 = 1000.0

chi2B = Models.chi2_B(T_fm, muB_fm)
chi4B = Models.chi_B(T_fm, muB_fm; order=4, derivative_backend=:taylordiff)
chi11BQ = Models.chi11_BQ(T_fm, muB_fm, muQ_fm, muS_fm)
chi210BQ = Models.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(2, 1, 0))
c2B = Models.cumulant_B(T_fm, muB_fm, V_fm3; order=2)
c4B = Models.cumulant_B(T_fm, muB_fm, V_fm3; order=4)
ks2 = Models.baryon_kappa_sigma2(T_fm, muB_fm)
```

## 输入与单位口径

- `T_fm`、`muB_fm`、`muQ_fm`、`muS_fm` 使用 `fm^-1`
- `V` 使用 `fm^3`
- `order` 表示对 `(mu_X / T)` 的导数阶数，其中 `X` 可以是 `B`、`Q`、`S`
- susceptibility 归一化采用主题级缩放 `chi_n ~ T^(n - 4) * d^n P / d(mu_X)^n`
- 默认 `derivative_backend=:auto`：纯 B/Q/S 单方向走单变量 TaylorDiff fast path；mixed BQS 组合走内部 multivariate Taylor jet；`ForwardDiff` 只保留为显式 reference/fallback
- cumulant 口径保持 `C_n = V * T^3 * chi_n`，因此四阶接口 `C_4^B` 可直接由 `cumulant_B(...; order=4)` 读取
- 当前统一入口以静态涨落为主，不包含传输系数或体粘滞导数主线

## 不属于本主题首轮范围的接口

以下接口虽然同属 `Models` 衍生量簇，但不应作为本主题主入口：

- `flavor_pressure_derivatives`
- `mass_derivatives`
- `thermo_derivatives`
- `bulk_derivative_coeffs`

这些能力与 derivatives 子主题的导数口径关系更紧，应在 `derived/derivatives/` 子主题收口。
