# Models Solver 文件归位与隐式最小适配实施单

## 背景与目标

- 基线：承接 PR #56 与 PR #57 的主线收敛结果（`ProblemSpec`/`GenericRootEngine`/`Solver` 主链已建立，`ImplicitSolver` 已退役）。
- 问题：`src/models/` 根目录仍暴露 `gap_solver.jl`、`constraint_solver.jl`、`implicit_gap.jl` 三个 solver 相关实现，分层边界不够清晰。
- 目标：
  - 将三者归位到 `src/models/solver/`；
  - 建立“隐式问题 + 求解配置”的最小通用构建层；
  - 对导数层保持稳定 API，不引入物理公式变更。

## 范围与非目标

### 本轮范围

- 文件归位（目录层级收口）。
- 新增最小隐式求解器抽象（接口级，不重写数值内核）。
- 通过薄兼容层维持现有导出 API 与调用路径。

### 非目标

- 不修改物理模型公式与热力学定义。
- 不引入新 mode 语义或业务工作流分支。
- 不在本轮做大规模 API 清理（保留历史导出名，后续再退役）。

## 实施任务清单

### P1 文件归位（不改语义）

- [x] 移动 `src/models/gap_solver.jl` -> `src/models/solver/GapSolver.jl`
- [x] 移动 `src/models/constraint_solver.jl` -> `src/models/solver/ConstraintSolver.jl`
- [x] 移动 `src/models/implicit_gap.jl` -> `src/models/solver/ImplicitGapLegacy.jl`
- [x] 更新 `src/models/Models.jl` include 路径，保持与现有加载顺序等价

### P2 新建隐式最小通用层

- [x] 新增 `src/models/solver/ImplicitProblem.jl`
- [x] 新增 `src/models/solver/ImplicitBuilder.jl`
- [x] 在 `Models.jl` 中纳入 include，并按需导出最小构建接口

### P3 以现有内核组装第一版 problem

- [x] 新增 `src/models/solver/ImplicitAdapters.jl`
- [x] 实现 `build_pnjl_fixedmu_problem(model; ...)`
- [x] 实现 `build_pnjl_flavor_mu_problem(model; ...)`
- [x] 实现 `build_njl_problem(model; ...)`
- [x] 内部仅复用 `solve_gap`/`gap_residual`，禁止重复实现求解逻辑

### P4 兼容层与导数层接线

- [x] 将 `ImplicitGapLegacy.jl` 改为薄转调层（`create_*` / `solve_*_with_derivatives` / `derive_*` 对外保持兼容）
- [x] 确认 `src/models/derivatives/ThermoDerivatives.jl` 仅依赖稳定入口，不感知旧文件位置

### P5 最小回归验证

- [x] `julia --project=. -e 'include("tests/unit/models/test_implicit_gap.jl")'`
- [x] `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl")'`
- [x] 若失败，记录 root cause 与最小修复，不扩展需求边界

## 最小接口草图（拟落地）

```julia
struct ImplicitProblem{F,C}
    forward_solve::F      # θ -> (x, meta)
    conditions::C         # (θ, x, meta) -> residual vector
    x_dim::Int
    theta_dim::Int
end

Base.@kwdef struct ImplicitSolverConfig
    linear_solver = DirectLinearSolver()
    representation = MatrixRepresentation()
end

function build_implicit_solver(
    problem::ImplicitProblem,
    cfg::ImplicitSolverConfig=ImplicitSolverConfig(),
)
    return ImplicitFunction(
        problem.forward_solve,
        problem.conditions;
        linear_solver=cfg.linear_solver,
        representation=cfg.representation,
    )
end
```

## 适配器接口草图（第一版）

```julia
build_pnjl_fixedmu_problem(model; xi=0.0, p_num=64, t_num=8, kwargs...)
build_pnjl_flavor_mu_problem(model; xi=0.0, p_num=64, t_num=8, kwargs...)
build_njl_problem(model; xi=0.0, p_num=64, t_num=8, kwargs...)
```

约束：

- `forward_solve` 允许 primal 路径使用 `Float64` 转换。
- `conditions` 必须 Dual-friendly，不允许将 AD 路径强制压回 `Float64`。

## 验收标准（DoD）

- 三个 solver 相关文件已从 `src/models/` 根目录迁出。
- 新增最小隐式通用层并可由 legacy 接口转调。
- 单元与集成最小回归通过。
- 对外导出 API 不破坏（本轮保兼容）。

## 风险与缓解

- 风险：一次性重写导致导数回归。
  - 缓解：采用“新建最小层 + 旧层薄转调”的渐进策略。
- 风险：include 顺序变化引入加载时依赖错误。
  - 缓解：先做“路径替换不变序”迁移，再做接口收敛。
- 风险：AD 路径被不当类型转换破坏。
  - 缓解：代码审查中强制检查 `conditions` 路径的类型保真。

## 执行记录

- 2026-04-07：已完成 P1-P5 本轮目标并通过最小回归。
- 验证结果：
  - `julia --project=. -e 'include("tests/unit/models/test_implicit_gap.jl")'` -> PASS (17/17)
  - `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl")'` -> PASS (4/4 + 8/8)

## 参考

- `docs/decisions/0004-solver-three-layer-contract-and-vector-kernel.md`
- `src/models/Models.jl`
- `src/models/derivatives/ThermoDerivatives.jl`
- `src/models/solver/ProblemSpec.jl`
- `src/models/solver/GenericRootEngine.jl`
