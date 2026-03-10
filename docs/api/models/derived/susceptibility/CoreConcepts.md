# Susceptibility 主题职责核心

本页说明 susceptibility 为什么应作为 `Models` 的“衍生量主题”，以及 BQS 组合、flavor 化学势导数和 cumulant 之间的职责边界。

## 1. 为什么它不是普通一级主题

susceptibility 不是模型主流程，而是模型热力学状态继续派生出的响应量族：

- 上游仍是 `Models` 的平衡态与压力计算
- 主题主对象是静态涨落、广义 susceptibility 与 cumulant
- 用户真正关心的是“如何读出守恒荷响应”，而不是“如何解 gap”

因此，它更接近 `Models` 的衍生量主题，而不是与 `phase`、`solver` 平级的流程主题。

## 2. 实现主线是 flavor 压强导数

当前实现位于 `src/models/derivatives/ConservedChargeSusceptibilities.jl`，主路线采用 AD，并优先围绕：

- `P(T, mu_u, mu_d, mu_s)`

来构造 `B/Q/S` 方向的 susceptibility。也就是说，BQS 不是独立重新求解的一套对象，而是 flavor 化学势导数经过线性变换后的主题级表达。

## 3. BQS 与 flavor 的边界

本主题需要显式区分两层：

- flavor 层：`mu_u`, `mu_d`, `mu_s`
- 物理读出层：`mu_B`, `mu_Q`, `mu_S`

统一接口 `conserved_charge_susceptibility` 与 `chi_BQS` 隐含了 BQS 到 flavor 的映射，因此多数用户不需要直接处理 flavor Jacobian；但这也意味着当前支持范围由底层映射和导数阶数共同限定。

## 4. 当前支持范围不是“任意阶任意组合”

按照当前实现，首轮稳定支持范围是：

- 纯单轴 `B/Q/S` 方向 `1..4` 阶
- 总二阶 mixed susceptibilities：`(1,1,0)`、`(1,0,1)`、`(0,1,1)` 以及纯二阶 `(2,0,0)`、`(0,2,0)`、`(0,0,2)`

因此，文档必须把“当前支持范围”写清楚，避免把统一入口误读为任意高阶泛化接口。

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