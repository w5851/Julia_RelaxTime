# Derivatives 主题职责核心

本页说明 derivatives 为什么应作为 `Models` 的“衍生量主题”，以及隐函数求导、链式法则、体粘滞组合导数和兼容层之间的职责边界。

## 1. 为什么它不是普通一级主题

derivatives 不是模型主流程，而是平衡态求解之后的派生分析层：

- 上游仍是 `Models` 的 gap 求解与热力学计算
- 本主题主对象是导数、组合导数与下游 transport 公式需要的系数
- 用户关心的是“如何读取派生量”，而不是“如何组织求解流程”

因此，它应归入 `derived` 总层，而不是与 `phase`、`solver` 同层平铺。

## 2. 实现主线是隐函数求导

当前 derivatives 主实现位于 `src/models/derivatives/ThermoDerivatives.jl`，核心路线是：

- 用隐函数求解器获得平衡态 `x_state(θ)`
- 令 `θ = (T, mu)`
- 对 `x_state(θ)`、热力学量和质量应用 AD 与链式法则

这意味着 derivatives 主题的职责核心不是单个公式，而是“平衡态求解 + 导数传播”的统一逻辑。

## 3. 热力学导数与质量导数是同一条链上的两个读出层

在当前实现中：

- `mass_derivatives` 提供质量及其对 `T` / `mu` 的导数
- `thermo_derivatives` 提供压力、能量、熵、数密度及其导数
- `bulk_derivative_coeffs` 则从 `thermo_derivatives` 结果中提取组合导数

所以这三类接口应被文档化为一条连续链，而不是三个互不相关的零散页面。

## 4. `dP_dT` / `dP_dmu` 是高阶导数补充入口

`dP_dT` 与 `dP_dmu` 提供总压强对温度或化学势的高阶导数。它们和 `thermo_derivatives` 的关系是：

- `thermo_derivatives` 更适合完整读出一阶热力学量
- `dP_dT` / `dP_dmu` 更适合单独关注某个标量导数并控制阶数

文档上不应把它们写成完全独立主题，而应作为 `ThermoDerivatives.md` 中的补充层。

## 5. bulk viscosity 主线是下游公式接口，而不是前台入门页

`bulk_viscosity_coefficients` 与 `compute_B_bracket` 服务的是体粘滞公式装配：

- 前者汇总 `v_n_sq`、`dμB_dT_sigma`、`dM_dT`、`dM_dμB` 等量
- 后者进一步计算体粘滞公式中的 `B` 项

这条链比 `thermo_derivatives` 更偏下游使用，因此应在主题中后置，而不是放在首页最前。

## 6. `legacy_transport_c_p` 必须被显式降级为兼容层

`legacy_transport_c_p` 当前只是从 `bulk_viscosity_coefficients(...).c_p` 取值的兼容导出：

- 它应该保留在 generated 导出全集中
- 但不应和 `thermo_derivatives`、`bulk_viscosity_coefficients` 并列成首选入口

这一点必须在人工维护页中明确写出，避免把历史兼容接口误导为长期主 API。

## 7. 与 susceptibility 的边界已经固定

当前边界应保持为：

- susceptibility：守恒荷 static fluctuations、cumulants、`baryon_*`
- derivatives：热力学导数、质量导数、bulk viscosity 组合导数

即使两个主题都依赖导数，也不应再回到旧 `ThermoDerivatives` 大页的混合叙事方式。