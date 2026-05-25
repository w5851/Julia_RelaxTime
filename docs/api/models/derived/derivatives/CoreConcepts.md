# Derivatives 主题职责核心

本页说明 derivatives 为什么应作为 `Models` 的“衍生量主题”，以及 Taylor-series gap 路线、链式读出、体粘滞组合导数和兼容层之间的职责边界。

## 1. 为什么它不是普通一级主题

derivatives 不是模型主流程，而是平衡态求解之后的派生分析层：

- 上游仍是 `Models` 的 gap 求解与热力学计算
- 本主题主对象是导数、组合导数与下游 transport 公式需要的系数
- 用户关心的是“如何读取派生量”，而不是“如何组织求解流程”

因此，它应归入 `derived` 总层，而不是与 `phase`、`solver` 同层平铺。

## 2. 实现主线是显式 Taylor-series 隐函数求解

当前 derivatives 主实现位于 `src/models/derivatives/ThermoDerivatives.jl`，PNJL 默认核心路线是：

- 先求 primal 平衡态 `x0`
- 在单变量 Taylor 代数中构造 `T(δ)`、`mu(δ)` 与 `x(δ)`
- 显式迭代求解 `F(x(δ), T(δ), mu(δ)) = 0`
- 从压强、质量或状态量 series 中提取所需导数

这意味着 derivatives 主题的职责核心不是单个公式，而是“平衡态求解 + 导数传播”的统一逻辑。ThermoDerivatives 生产路径不再提供 `ForwardDiff + ImplicitDifferentiation` fallback；迁移时使用默认 `:auto` 或显式 `:taylordiff`。

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

## 6. `c_p` 的唯一读取路径

`c_p` 统一通过 `bulk_viscosity_coefficients(...).c_p` 读取。

- 不再保留 `legacy_transport_c_p` 兼容导出
- 避免在同一语义上长期维持双 API 入口

## 7. 与 susceptibility 的边界已经固定

当前边界应保持为：

- susceptibility：守恒荷 static fluctuations、cumulants、`baryon_*`
- derivatives：热力学导数、质量导数、bulk viscosity 组合导数

即使两个主题都依赖导数，也不应再回到旧 `ThermoDerivatives` 大页的混合叙事方式。
