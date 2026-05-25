# Susceptibility 主题职责核心

本页说明 susceptibility 为什么应作为 `Models` 的“衍生量主题”，以及 BQS 组合、flavor 化学势导数和 cumulant 之间的职责边界。

## 1. 为什么它不是普通一级主题

susceptibility 不是模型主流程，而是模型热力学状态继续派生出的响应量族：

- 上游仍是 `Models` 的平衡态与压力计算
- 主题主对象是静态涨落、广义 susceptibility 与 cumulant
- 用户真正关心的是“如何读出守恒荷响应”，而不是“如何解 gap”

因此，它更接近 `Models` 的衍生量主题，而不是与 `phase`、`solver` 平级的流程主题。

## 2. 实现主线是 flavor 压强导数与 baryon Taylor-mode

当前实现位于 `src/models/derivatives/ConservedChargeSusceptibilities.jl`，主路线采用 AD，并优先围绕：

- `P(T, mu_u, mu_d, mu_s)`

来构造 `B/Q/S` 方向的 susceptibility。也就是说，BQS 不是独立重新求解的一套对象，而是 flavor 化学势导数经过线性变换后的主题级表达。

迁移说明：高阶单方向 `chi_B` 不再把高阶导数默认压在嵌套 ForwardDiff Dual 上。`derivative_backend=:auto` 现在对纯 B/Q/S 单方向都走 TaylorDiff 单变量 fast path：

1. 先在目标点求 primal gap 解 `x0`；
2. 用 ForwardDiff 在 primal 点构造 `J0 = dF/dx`，这里仍保留 ForwardDiff 作为低阶 Jacobian 工具；
3. 令单一方向变量与 `x(δ)` 携带 TaylorDiff univariate series；
4. 逐阶求解 Taylor 系数，使 `F(x(δ), μ(δ))=0`；
5. 从 pressure series 提取 `d^nP/dδ^n` 并乘以 `T^(n - 4)`。

mixed BQS 高阶导数走单独的内部 multivariate Taylor jet backend。它复用相同的 primal gap solve、`J0 = dF/dx`、逐阶线性求解、残差检查和 pressure extraction，但多变量系数布局只在非零 B/Q/S 轴数量超过 1 时启用。单方向 `B/Q/S` 不会退回通用 `D=1` jet，而是继续委派给单变量 TaylorDiff fast path。

ForwardDiff 路径没有移除：它仍作为显式 `:forwarddiff` reference/fallback，以及 TaylorDiff/jet 路径中的 primal Jacobian/gradient。

## 3. BQS 与 flavor 的边界

本主题需要显式区分两层：

- flavor 层：`mu_u`, `mu_d`, `mu_s`
- 物理读出层：`mu_B`, `mu_Q`, `mu_S`

统一接口 `conserved_charge_susceptibility` 与 `chi_BQS` 隐含了 BQS 到 flavor 的映射，因此多数用户不需要直接处理 flavor Jacobian；但这也意味着当前支持范围由底层映射和导数阶数共同限定。

## 4. 当前支持范围

按照当前实现，稳定支持范围是：

- 纯单轴 `B/Q/S` 方向在 TaylorDiff backend 下可扩展到更高单方向阶数；`ForwardDiff` fallback 仍保留 `1..4`
- mixed BQS susceptibilities 在 `:auto` / `:taylordiff` / `:mixedjet` 下由内部 multivariate jet 支持，阶数由 `sum(orders)` 决定
- `ForwardDiff` fallback 对 mixed 组合只保留总二阶 reference：`(1,1,0)`、`(1,0,1)`、`(0,1,1)`

因此，统一入口可以表达高阶 mixed BQS，但成本会随 jet 变量数和总阶数增长；性能敏感的纯单方向仍应保留默认 `:auto`。

## 5. `T` 缩放是主题级合同

实现不是直接对 `P / T^4` 做高阶 AD，而是在压力导数完成后施加：

`T^(n - 4)`

这一缩放口径。对于主题说明而言，这比函数签名更重要，因为它决定了：

- `chi_n` 的量纲与归一化语义
- `cumulant_BQS = V * T^3 * chi_BQS` 的构造关系
- `baryon_Ssigma`、`baryon_kappa_sigma2` 等组合量的理解方式

## 6. derivatives 仍是相邻但独立的主题

旧 `ThermoDerivatives` 页面里的高价值内容主要是：

- `thermo_derivatives`
- `dP_dT` / `dP_dmu`
- `bulk_derivative_coeffs`

这些接口与 susceptibility 主题共享“导数”背景，但不应在第一轮被混入同一主题。当前更合理的边界是：

- 守恒荷涨落、cumulant、baryon 组合量进入 susceptibility
- 热力学导数、质量导数、体粘滞组合导数留给 derivatives
