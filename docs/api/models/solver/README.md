# Models 核心求解与约束主题 API

本主题组织 `Models` 统一入口下最核心的求解公共表面：模型创建、状态合同、约束模式、初值策略与隐式求解器工厂。

它回答的不是“某个流程怎么跑完”，而是更基础的问题：

- `Models` 的最小稳定求解入口是什么
- 平衡态状态 `x_state` 与 `mu_vec` 的统一合同是什么
- `solve_gap`、`solve`、`solve_multi` 分别适合什么场景
- 约束模式、seed strategy、implicit solver factory 各自承担什么职责

推荐阅读顺序：

1. [Overview.md](Overview.md)：先理解 `create_model`、`solve_gap`、`solve`、`solve_multi` 怎么选
2. [CoreConcepts.md](CoreConcepts.md)：理解状态合同、约束模式、seed 与 solver factory 的职责边界
3. [StateContract.md](StateContract.md)：平均场状态与化学势输入口径
4. [ConstraintModes.md](ConstraintModes.md)：固定化学势/密度/熵/比熵等约束入口
5. [SeedStrategies.md](SeedStrategies.md)：默认、多初值、连续性与相变感知策略
6. [ImplicitSolvers.md](ImplicitSolvers.md)：隐式求解器工厂与导数接口
7. [generated/Exports.md](generated/Exports.md)：公开导出全集与覆盖检查

本主题优先覆盖的 `Models` 公开表面包括：

- `Models.create_model`
- `Models.solve_gap`
- `Models.solve`
- `Models.solve_multi`
- `Models.solve_constraint`
- `Models.ProblemSpec` / `Models.build_problem_spec`
- `Models.AbstractConstraintComponent` / `Models.build_constraint_components`
- `Models.MeanFieldState` / `Models.meanfield_state` / `Models.state_vector`
- `Models.ConstraintModes`
- `Models.SeedStrategy` 家族
- `Models.create_implicit_gap_solver` / `Models.create_pnjl_implicit_solver` / `Models.create_flavor_mu_implicit_gap_solver`
- `Models.solve_with_derivatives`

本主题已经吸收旧 `docs/api/pnjl/` 求解器相关页面中的主要价值。旧页后续应只承担迁移说明或兼容层定位，不再作为新主题主说明页。

## 维护者阅读顺序（ProblemSpec 主链）

当你需要维护 `ProblemSpec` 路径时，建议按以下顺序读源码：

1. `src/models/solver/ProblemSpec.jl`：只看契约核心（`ExtraConstraints`、`ProblemSpec`、`build_problem_spec` 绑定入口）
2. `src/models/solver/ProblemSpecOrchestrator.jl`：看 attempt plan 生成、governed attempt 执行编排、mode forward_solve 组织
3. `src/models/solver/SolverDiagnostics.jl`：看 `diagnostic_level=:none/:summary/:full` 的摘要/全量候选诊断拼装逻辑
4. `src/models/solver/Solver.jl`：看入口层如何将 runtime options 传入 `ProblemSpec.forward_solve`
