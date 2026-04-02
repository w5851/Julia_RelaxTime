# 核心概念与职责边界

本页是 solver 主题的第二层。它不按源码顺序复述实现，而是解释状态合同、约束模式、seed 策略和求解器工厂之间的边界。

## 1. 三类入口的职责分工

这一主题里的公开接口大致分成三类：

- 用户入口：`create_model`、`solve_gap`、`solve`、`solve_multi`
- 合同与策略：`MeanFieldState`、`ConstraintModes`、`SeedStrategy` 家族
- 求解器工厂：`create_implicit_gap_solver`、`create_pnjl_implicit_solver`、`create_flavor_mu_implicit_gap_solver`

理解这三类边界，比记忆某个文件名更重要。

## 2. 状态合同先于算法细节

`Models` 体系首先统一的是“状态与输入如何表示”：

- `x_state` 统一收敛到 `MeanFieldState`
- `mu_vec` 统一收敛到三味 `SVector{3}`
- `state_vector` 负责在对象表示和 5 维向量表示之间切换

这意味着后续 `omega`、`solve_gap`、`solve`、implicit AD、workflow 都共享同一状态合同，而不是各自维护一套向量口径。

## 3. `solve_gap` 与 `solve` 不是重复接口

它们服务的层次不同：

- `solve_gap`：给出平衡态状态对象，偏模型层
- `solve`：给出带热力学载荷和统计信息的 `SolverResult`，偏约束求解层

如果你只是要一个平衡解状态，优先使用 `solve_gap`；如果你要在 `FixedMu` / `FixedRho` / `FixedEntropy` / `FixedSigma` / `FixedAsymmetricRho` 下拿统一结果结构，使用 `solve`。

## 4. 约束模式定义的是问题形态

`ConstraintModes` 的价值不在于名字本身，而在于它把“参数维度”和“未知量维度”固定下来。

例如：

- `FixedMu()`：5 维状态，参数是 `[T, μ]`
- `FixedRho(...)`：8 维状态，参数是 `[T]`
- `FixedEntropy(...)` / `FixedSigma(...)`：同样属于 8 维约束问题

因此，约束模式不是一个装饰标签，而是直接决定 residual 构造和求解形态的合同对象。

## 4.1 ProblemSpec 与 ConstraintComponents 的职责边界

R1 收敛后，`FixedRho` 主链不再只依赖“单个大 residual 函数”，而是通过“问题规格 + 组件映射”组织：

- `ProblemSpec`：描述某个 mode 的求解契约（`conditions` / `forward_solve` / selector / hard_rules）
- `ConstraintComponents`：描述约束由哪些职责组件拼接（stationarity、equal-mu、macro targets 等）

这种分层的作用是把“语义定义”与“数值策略”分开：

- 语义定义由 `ConstraintModes + ConstraintComponents` 冻结
- 数值策略由 `forward_solve` 与候选治理规则演进

在当前实现中，`FixedRho` 的 `ProblemSpec.forward_solve` 已进入显式候选池 + 统一 hard constraints + selector 路径，避免语义漂移被隐藏在特化分支内。

## 5. Seed strategy 负责“怎么开始”，不负责“怎么收尾”

`SeedStrategy` 家族的职责是给求解器提供初值，而不是代替求解器决定所有分支逻辑。

这几个层次要分清：

- `DefaultSeed`：固定经验初值
- `MultiSeed`：显式尝试多初值候选
- `ContinuitySeed` / `HybridContinuitySeed`：扫描路径上的状态延续
- `PhaseAwareSeed` / `PhaseAwareContinuitySeed`：在相变语境下做更稳妥的初值选择

它们与 `solve` / `solve_multi` 的配合关系，是本主题最关键的维护知识之一。

## 6. 求解器工厂面向导数链路

`create_implicit_gap_solver`、`create_pnjl_implicit_solver`、`create_flavor_mu_implicit_gap_solver`、`create_implicit_solver` 不服务普通单点求解，而是服务隐函数求导与涨落/导数链路。

因此它们的职责边界是：

- 普通求解：`solve_gap` / `solve`
- 需要 `dx/dθ`：implicit solver factory，包括更通用的 `create_implicit_solver`
- 需要导数组合：`solve_with_derivatives`、`solve_pnjl_with_derivatives` 等高阶接口

## 7. 为什么这个主题必须吸收旧页内容

旧 `pnjl` 文档里已经沉淀了很多高价值信息：

- `solve_multi` 与 `MultiSeed` 的关系
- `PhaseAwareContinuitySeed` 的首点自举语义
- 约束模式对应的状态维度与参数维度
- `ImplicitFunction` 路径和 AD 边界

如果这些信息继续只留在旧页，新 `models/solver` 主题就只是一个空壳入口。这个任务的目标正是把这些说明吸收进来，让新主题独立成立。
