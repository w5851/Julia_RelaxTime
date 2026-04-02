# 约束模式与约束求解入口

本页说明 `Models` 中“问题是怎么被定义的”。在 solver 主题里，约束模式不是附属概念，而是决定状态维度、参数维度和求解目标的核心合同。

## `ConstraintModes`

当前核心模式包括：

- `FixedMu()`
- `FixedRho(rho_target)`
- `FixedAsymmetricRho(rho_target, ud_ratio_target, s_target)`
- `FixedEntropy(s_target)`
- `FixedSigma(sigma_target)`

它们回答的是：给定什么约束，求解器应该把哪些量当作未知量。

## 状态维度与参数维度

当前约定：

- `FixedMu()`：状态维度 5，参数维度 2
- 其它固定密度/熵/比熵/非对称模式：状态维度 8，参数维度 1

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

- `GapParams`：打包温度、积分节点、各向异性等残差构造上下文
- `build_conditions`：构造条件/方程对象
- `build_residual!`：生成可直接交给求解器的残差函数
- `gap_state_dim` / `gap_residual`：提供模型级的状态维度与残差入口

这类接口多数更偏维护者或算法开发者，但必须在主题中说明清楚，因为它们定义了 `solve` 系列入口背后的问题形式。

## `ProblemSpec` 与组件化约束主链（R1）

在 R1 迭代中，`FixedRho` 已支持通过 `ProblemSpec` 进入主链：

- `build_problem_spec(mode)` 提供 mode 对应的默认求解契约
- `solve_constraint(...; problem_spec=spec)` 可显式指定契约
- `FixedRho` 的 `ProblemSpec.forward_solve` 已使用统一候选治理规则

配套组件层通过 `build_constraint_components(mode)` 暴露 mode 到约束职责的映射，用于固定“语义同构”口径。当前核心组件语义包括：

- `StationarityComponent`
- `EqualMuComponent`
- `FixedBaryonDensityComponent`
- `AsymmetricDensityComponent`
- `FixedEntropyComponent`
- `FixedSigmaComponent`

## 为什么这部分不应继续只放在旧 `pnjl` 页面

约束模式本质上已经是 `Models` 的公共合同，而不只是 legacy PNJL 的内部文档主题。继续把它们主要留在旧页，会让新用户误以为 `Models` 只有流程 façade，没有自己的求解合同。
