# Models 核心求解与约束主题 API

本主题组织 `Models` 统一入口下最核心的求解公共表面：模型创建、状态合同、约束模式、初值策略与导数入口。

它回答的不是“某个流程怎么跑完”，而是更基础的问题：

- `Models` 的最小稳定求解入口是什么
- 平衡态状态 `x_state` 与 `mu_vec` 的统一合同是什么
- `solve_gap`、`solve`、`solve_multi` 分别适合什么场景
- 约束模式、seed strategy、导数入口和 residual problem builders 各自承担什么职责

推荐阅读顺序：

1. [Overview.md](Overview.md)：先理解 `create_model`、`solve_gap`、`solve`、`solve_multi` 怎么选
2. [CoreConcepts.md](CoreConcepts.md)：理解状态合同、约束模式、seed 与 solver factory 的职责边界
3. [StateContract.md](StateContract.md)：平均场状态与化学势输入口径
4. [ResultDiagnosticErrorContracts.md](ResultDiagnosticErrorContracts.md)：稳定 Result/Diagnostic/Error 契约与版本策略
5. [ConstraintModes.md](ConstraintModes.md)：固定化学势/密度/熵/比熵等约束入口
6. [SeedStrategies.md](SeedStrategies.md)：默认、多初值、连续性与相变感知策略
7. [ImplicitSolvers.md](ImplicitSolvers.md)：TD-first 导数接口与 residual problem builders
8. [generated/Exports.md](generated/Exports.md)：公开导出全集与覆盖检查

本主题优先覆盖的 `Models` 公开表面包括：

- `Models.create_model`
- `Models.solve_gap`
- `Models.solve`
- `Models.solve_multi`
- `Models.solve_constraint`
- `Models.coerce_solver_result` / `Models.solver_result_view`
- `Models.coerce_solver_diagnostic_public_view` / `Models.to_public_namedtuple`
- `Models.ProblemSpec` / `Models.build_problem_spec`
- `Models.AbstractConstraintComponent` / `Models.build_constraint_components`
- `Models.ImplicitProblem`
- `Models.build_pnjl_fixedmu_problem` / `Models.build_pnjl_flavor_mu_problem` / `Models.build_njl_problem`
- `Models.MeanFieldState` / `Models.meanfield_state` / `Models.state_vector`
- `Models.ConstraintModes`
- `Models.SeedStrategy` 家族
- `Models.solve_with_derivatives`
- `Models.solve_pnjl_with_derivatives` / `Models.solve_pnjl_with_flavor_mu_derivatives`
- `Models.ThermoDiffContext` / `Models.ParamSpec` / `Models.DiffTarget`
- `Models.build_thermo_diff_context` / `Models.diff_target` / `Models.jacobian`
- `Models.build_pilot_diff_context` / `Models.eval_pilot_derivatives`（Issue #81 试点统一导数服务）

本主题已经吸收旧 `docs/api/pnjl/` 求解器相关页面中的主要价值。旧页后续应只承担迁移说明或兼容层定位，不再作为新主题主说明页。

## 维护者阅读顺序（ProblemSpec 主链）

当你需要维护 `ProblemSpec` 路径时，建议按以下顺序读源码：

1. `src/models/solver/ProblemSpec.jl`：只看契约核心（`ExtraConstraints`、`ProblemSpec`、`build_problem_spec` 绑定入口）
2. `src/models/solver/ProblemSpecOrchestrator.jl`：看 attempt plan 生成、governed attempt 执行编排、mode forward_solve 组织
3. `src/models/solver/SolverDiagnostics.jl`：看 `diagnostic_level=:none/:summary/:full` 的摘要/全量候选诊断拼装逻辑
4. `src/models/solver/Solver.jl`：看入口层如何将 runtime options 传入 `ProblemSpec.forward_solve`

## 维护者阅读顺序（治理主引擎）

PR-4 后，solver 主链的 attempt/fallback/selection 采用单治理引擎语义；维护时建议按以下顺序阅读：

1. `src/models/solver/CandidateGovernance.jl`：`execute_attempt_pool`、`governance_quality_tag`、`execute_governance_selector`（主治理入口）
2. `src/models/solver/ConstraintSolverCommon.jl`：`select_pressure_max_candidate` / `select_residual_min_candidate`（统一排序规则）
3. `src/models/solver/ProblemSpecOrchestrator.jl`：模式求解如何接入主治理入口与诊断拼装
4. `src/models/solver/GenericRootEngine.jl`：保留为通用根求解兼容层，不再承担 solver 主链治理语义
