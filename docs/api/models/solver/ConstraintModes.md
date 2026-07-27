# 约束模式与约束求解入口

本页说明 `Models` 中“问题是怎么被定义的”。在 solver 主题里，约束模式不是附属概念，而是决定状态维度、参数维度和求解目标的核心合同。

## `ConstraintModes`

当前核心模式包括：

- `FixedMu()`
- `FixedRho(rho_target)`
- `FixedAsymmetricRho(rho_target, ud_ratio_target, s_target)`
- `FixedMuBConservedCharges(muB_fm, charge_to_baryon_ratio=0.4, strangeness_density_target=0.0)`
- `FixedEntropy(s_target)`
- `FixedSigma(sigma_target)`

它们回答的是：给定什么约束，求解器应该把哪些量当作未知量。

## 状态维度与参数维度

当前约定：

- `FixedMu()`：状态维度 5，参数维度 2
- 其它固定密度/熵/比熵/非对称/守恒荷模式：状态维度 8，参数维度 1

这个维度合同会直接影响：

- residual 的构造
- 初值如何扩展
- `solve` / `solve_constraint` 应该返回什么结构

补充（Plan-B / Gate A 迁移期）：

- `ModelStateSchema` 已作为状态字段布局的显式合同引入；
- 推荐通过 `schema_for_model(model_kind)` + `flatten_state` / `unflatten_state` 处理状态向量映射；
- 约束组装层开始支持 schema 驱动入口（例如 `build_conditions(mode, params, schema; mu_dim=...)`），用于逐步替代固定切片假设。

## `solve` 与 `solve_constraint`

对多数调用方：

- 固定化学势可直接使用 `solve(FixedMu(), ...)`
- 其它约束模式同样可通过 `solve(mode, ...)` 使用标准求解入口

而 `solve_constraint` 是进阶统一入口：

- `solve_constraint`

> Plan-B / Gate A 当前主路径为 `solve_constraint`；旧 fixed-* compat 入口已从公开 API 中移除。

这些接口更适合“我明确要控制某个约束求解子路径”的场景，而不是首页首屏主推荐入口。

## `GapParams`、`build_conditions` 与 `build_residual!`

约束模式不会单独工作，它们需要配合条件构建层：

- `GapParams`：打包温度、积分节点、各向异性和显式热积分策略/误差控制等残差构造上下文；复制到邻近温度时保留同一数值口径
- `build_conditions`：构造条件/方程对象
- `build_residual!`：生成可直接交给求解器的残差函数
- `gap_state_dim` / `gap_residual`：提供模型级的状态维度与残差入口

迁移收敛补充（PR-D / Residual Spine）：

- `gap_core_residual!`：作为 solver 主链中的统一残差内核（2/3/5 维）。
- `ConstraintSolver` 在约束求解中计算 gap-norm 时，统一复用 `gap_core_residual!` 路径，不再在各 mode 重复直调 `gap_residual` 计算。
- `gap_residual` 的定位是“模型级公开残差入口（API/适配层）”；`gap_core_residual!` 的定位是“solver 内核主链（装配与一致性）”。

这类接口多数更偏维护者或算法开发者，但必须在主题中说明清楚，因为它们定义了 `solve` 系列入口背后的问题形式。

## `ProblemSpec` 与组件化约束主链（R1）

在 R1 迭代中，`FixedRho` 已支持通过 `ProblemSpec` 进入主链：

- `build_problem_spec(mode)` 提供 mode 对应的默认求解契约
- `solve_constraint(...; problem_spec=spec)` 可显式指定契约
- `FixedRho` 的 `ProblemSpec.forward_solve` 已使用统一候选治理规则

当前默认行为：

- `solve_constraint` 对 `FixedRho` / `FixedEntropy` / `FixedSigma` / `FixedAsymmetricRho` / `FixedMuBConservedCharges` 默认走 `ProblemSpec` 主链
- `use_problem_spec` / `allow_legacy_path` / `warn_on_legacy_path` 兼容参数已移除；调用时传入会抛 `ArgumentError`
- `problem_spec` 仍作为可选覆盖项（`ProblemSpec` 或 `nothing`），用于显式契约注入

配套组件层通过 `build_constraint_components(mode)` 暴露 mode 到约束职责的映射，用于固定“语义同构”口径。当前核心组件语义包括：

- `StationarityComponent`
- `EqualMuComponent`
- `FixedBaryonDensityComponent`
- `AsymmetricDensityComponent`
- `FixedMuBComponent`
- `ConservedChargeDensityComponent`
- `FixedEntropyComponent`
- `FixedSigmaComponent`

## 固定 `mu_B` 的 B/Q/S 守恒荷模式

`FixedMuBConservedCharges` 用于在固定 `(T,mu_B)` 下联立求解平均场和三味 chemical potentials。输入采用自然单位：

- `muB_fm`：fm^-1；
- `charge_to_baryon_ratio`：无量纲，重离子碰撞常用场景为 `0.4`；
- `strangeness_density_target`：fm^-3，首个物理场景为 `0`。

三味 chemical potentials 与守恒荷 chemical potentials 的约定为：

```text
mu_u = mu_B/3 + 2 mu_Q/3
mu_d = mu_B/3 - mu_Q/3
mu_s = mu_B/3 - mu_Q/3 - mu_S
```

反演关系为 `mu_B=mu_u+2mu_d`、`mu_Q=mu_u-mu_d`、`mu_S=mu_d-mu_s`。稳定 helper 为：

- `flavor_mu_from_bqs`
- `conserved_mu_from_flavor`
- `conserved_densities_from_flavor`

密度 helper 的输入必须是 flavor **净密度** `(rho_u,rho_d,rho_s)`，不是 `number_densities(...).quark` 的粒子密度。若从粒子/反粒子入口读取，先使用 `quark-antiquark`。

该 mode 的三个非 gap residual 为：

```text
mu_u + 2mu_d - mu_B
(rho_Q - charge_to_baryon_ratio*rho_B)/rho0
(rho_S - strangeness_density_target)/rho0
```

其中 `rho_B=(rho_u+rho_d+rho_s)/3`、`rho_Q=(2rho_u-rho_d-rho_s)/3`、`rho_S=-rho_s`。charge residual 使用 affine form，不在 `rho_B≈0` 时除以密度。

```julia
model = Models.create_model(:PNJL)
mode = Models.FixedMuBConservedCharges(240 / 197.3269804, 0.4, 0.0)
result = Models.solve(model, mode, 170 / 197.3269804; p_num=8, t_num=4)

mu_bqs = Models.conserved_mu_from_flavor(result.mu_vec...)
rho_flavor = Models.model_rho(model, result.x_state, result.mu_vec, 170 / 197.3269804)
rho_bqs = Models.conserved_densities_from_flavor(rho_flavor)
```

当前实现是 quark/mean-field-only 守恒荷平衡，不包含介子数密度对 `rho_Q`、`rho_S` 的反馈。介子守恒荷应在显式外层修正或完整热力学反馈阶段加入，不能把本 mode 的 `rho_S=0` 误写成“全体系奇异数严格为零”。

## 为什么这部分不应继续只放在旧 `pnjl` 页面

约束模式本质上已经是 `Models` 的公共合同，而不只是 legacy PNJL 的内部文档主题。继续把它们主要留在旧页，会让新用户误以为 `Models` 只有流程 façade，没有自己的求解合同。
