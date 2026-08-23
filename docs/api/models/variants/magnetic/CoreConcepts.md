# Magnetic 主题职责核心

本页说明 magnetic 为什么应作为 `Models` 的“模型变体主题”，以及 Landau 离散化、IMC 耦合和 `n_max` 收敛报告之间的职责边界。

## 1. 为什么它不是普通一级主题

magnetic 不是一个单纯流程主题，而是对现有 PNJL 模型族的一条物理变体支路：

- 零磁场基线仍由 `PNJLModel` 提供
- 非零磁场通过 `solve_magnetic_gap` 对磁场 `Omega` 做完整五维驻点求解；普通 `solve_gap` 可返回已找到候选中按 `Omega` 选出的 convenience state，但不替代 branch-aware 候选输出。Hessian 分类是可选诊断，不是 PNJL 默认生产过滤条件；分支 API 始终保留全部候选
- 在热力学层切换到 Landau 能级离散化与磁场相关公式

因此，它更接近“模型变体”而不是又一个与 `phase`、`scans` 平级的流程目录。

## 2. 两层实现边界

当前磁场能力大致分成两层：

- 模型适配层：`PNJLMagneticModel`
- 变体核心层：`MagneticIntegrals.jl` 与 `MagneticThermodynamics.jl`

前者负责把 magnetic 变体纳入 `Models` 聚合语义；后者负责 Landau 积分、IMC 耦合和磁场热力学计算。

磁场 solver 的正式边界是 `T_fm > 0`、`xi=0`。每个候选包含 seed、五维 residual、`Omega`、物理性、求解方法、迭代数和 `n_max`；`classify_stability=true` 时额外用局部 Hessian 标记 `local_minimum` 或 `saddle_or_maximum`。这些标签用于诊断和显式研究策略，不构成默认生产拒绝条件；即使候选被标为 `saddle_or_maximum`，也必须保留在 branch-aware 结果中。Hessian 诊断不是默认求解成本的一部分。

## 3. `eB -> 0` 退化是主题级合同

`calculate_magnetic_omega_components` 在 `|eB_fm2| <= 1e-14` 时自动退化到零磁场 PNJL 路径。这不是实现小技巧，而是主题级兼容合同：

- 保证 magnetic 变体与零磁场主链可比较
- 避免用户为零磁场点强行切换另一套接口
- 使 fixed-point 或回归脚本在小磁场极限下保持稳定行为

## 4. `n_max` 收敛治理不是附属功能

Landau 求和的核心风险不是语法使用，而是截断是否足够：

- `resolve_nmax_from_cutoff` 给出自动估计
- `magnetic_nmax_convergence_report` 比较 `n_base` 与 `n_base + delta_n` 的相对差

这意味着 magnetic 主题必须把“如何检查收敛”放进主说明，而不是仅列函数签名。

## 5. IMC 与热力学组件的边界

`coupling_GB` 负责 IMC 耦合 `G(B)`；热力学主接口负责把：

- 手征项
- Polyakov 势
- Landau 真空项
- Landau 热项

组装成完整 `omega`。因此 `coupling_GB` 不应被误解为独立业务入口，而应放在 thermodynamics 的核心概念链上理解。

## 6. 旧文档吸收原则

旧 `pnjl` magnetic 兼容页中最有价值的，是单位约定、参数口径、脚本与 baseline 指引。这些内容已应被新主题直接吸收，而不是继续让新主题依赖旧页做主体补充。
