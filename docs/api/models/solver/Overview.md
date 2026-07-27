# Models 核心求解入口总览

本页优先回答“我应该调用哪个入口”，而不是先展开底层实现细节。

## 最小稳定入口

对大多数新调用方，核心入口顺序如下：

1. `Models.create_model`：先选择模型对象
2. `Models.solve_gap`：直接求平衡态平均场状态
3. `Models.solve`：在 PNJL 约束模式下求完整 `SolverResult`
4. `Models.solve_multi`：当你明确需要多初值筛选物理解时使用
5. `Models.coerce_solver_result` / `Models.solver_result_view`：当你需要稳定序列化视图或版本化契约检查时使用

Phase 3 收口说明：`fixedmu_problem_spec_active` 兼容字段与 `fixedmu_use_problem_spec` 开关均已移除，`solve_constraint` 统一走 ProblemSpec 主链。

## 什么时候用 `create_model`

当你明确在 `Models` 多重派发体系下工作时，先构造模型对象：

```julia
model = Models.create_model(:PNJL)
st = Models.solve_gap(model, T_fm, mu_vec)
```

适合场景：

- 你需要直接调用 `omega`、`omega_components`、`model_thermo` 等模型方法
- 你要在 NJL / PNJL / RPNJL 之间切换实现
- 你希望后续工作流沿模型对象接口继续展开

## 什么时候用 `solve_gap`

`Models.solve_gap` 适合“我只想拿到平衡态平均场状态”的场景。

典型输出是 `MeanFieldState` 风格的状态对象，可继续用于：

- `omega` / `omega_components`
- `state_vector`
- 自定义后处理或 TD-first 导数链路

它是最轻量的统一入口，不承担 CSV、扫描或完整 workflow 打包职责。

## 什么时候用 `solve`

`Models.solve` 适合“我已经确定约束模式，希望拿到完整的求解结果载荷”的场景。

典型形式：

```julia
res = Models.solve(Models.FixedMu(), T_fm, mu_fm)
```

适合场景：

- 你要同时拿到 `x_state`、`mu_vec`、`omega`、`pressure`、`rho_norm`、`entropy`、`energy`、`masses`
- 你在固定密度、固定熵等约束模式下工作
- 你需要一个标准化 `SolverResult` 而不是仅状态对象

## 什么时候用 `solve_multi`

`Models.solve_multi` 适合一阶相变附近或多值解区域，需要显式尝试多组初值并按物理判据选解的场景。

模式支持说明：`solve_multi` 当前支持 `FixedMu | FixedRho | FixedAsymmetricRho | FixedMuBConservedCharges`，不覆盖 `FixedEntropy` / `FixedSigma`。

它不是常规入口的默认替代品，而是“我明确知道这里有分支竞争”的增强入口。

## Models 边界签名与参数合同（Plan-A）

在 `Models` 边界，`solve*` 系列当前稳定公共表面为：

- `solve(model, mode, ...)` 与 `solve(mode, ...)`
- `solve_multi(model, mode, ...)` 与 `solve_multi(mode, ...)`
- `solve_constraint(model, mode, T_fm; ...)`
- `solve_vec(model, mode, theta_vec; ...)`
- `solve_named(model, mode, theta_named; ...)`

关键参数合同：

- `solve_constraint` 统一走 `ProblemSpec` 主链；`use_problem_spec` / `allow_legacy_path` / `warn_on_legacy_path` / `fixedmu_use_problem_spec` / `legacy_fallback_plugin` 已移除。
- `solve_vec`：`FixedMu` 约定 `theta_vec=[T_fm, μ_fm]`（长度 2）；其它模式约定 `theta_vec=[T_fm]`（长度 1）。`FixedMuBConservedCharges` 的 `muB_fm` 固定在 mode 中，不重复放入 `theta_vec`。
- `solve_named`：`FixedMu` 需要 `(:T_fm, :μ_fm)`；其它模式需要 `:T_fm`。

`solve_vec` / `solve_named` 在本页按签名与参数契约层覆盖；逐符号细节以自动生成导出索引为准（见 `docs/api/models/solver/generated/Exports.md`）。

## 用户入口分层建议

### 新使用者

- `create_model`
- `solve_gap`
- `solve`

### 进阶使用者

- `solve_multi`
- `solve_constraint`
- `MeanFieldState` / `state_vector` / `normalize_mu_vec`

### 维护者或算法开发者

- `ConstraintModes`
- `SeedStrategy` 家族
- `build_conditions` / `build_residual!`
- `solve_pnjl_with_derivatives` / `solve_pnjl_with_flavor_mu_derivatives`
- residual problem builders：`build_njl_problem` / `build_pnjl_fixedmu_problem` / `build_pnjl_flavor_mu_problem`

## 与其它主题的关系

- phase 主题、workflow 主题、scans 主题都建立在本主题定义的求解合同之上
- 如果你先理解本主题，再去看 [../workflows/README.md](../workflows/README.md) 或 [../scans/README.md](../scans/README.md)，会更容易看清哪些是“业务流程入口”，哪些是“核心求解公共表面”

## 迁移期契约说明（Plan-B / Gate A）

- `SolverResult` 已从固定 `SVector{5}/SVector{3}` 扩展为向量泛型字段，结果消费侧不应再假设状态维度恒为 5、化学势维度恒为 3。
- scan 结果适配层（`TmuScan` / `TrhoScan`）已按动态向量透传 `x_state` 与 `mu_vec`，为后续 schema-driven 主链路铺路。
- 当需要显式状态布局时，优先使用 `ModelStateSchema` 系列 API（`schema_for_model`、`flatten_state`、`unflatten_state`）。
- `SolverResult` 稳定契约版本字段为 `contract_version=:v1`；上层消费建议统一通过 `solver_result_view` 取公共视图。
- 诊断契约稳定公共视图与内部调试视图已分层；上层默认应消费 `to_public_namedtuple` / `coerce_solver_diagnostic_public_view`。

## 导数契约骨架（Phase-1 / Issue #79）

solver/diff 当前提供最小稳定骨架，先冻结“上下文 + 目标 + 参数 + Jacobian 形状”约定，再在后续阶段接入真实导数数值实现：

- `build_thermo_diff_context(result; mode, model, theta, spec_override=nothing)`
  - 构建 `ThermoDiffContext`；`theta` 采用 NamedTuple 口径（首批支持 `T_fm` / `mu_fm(兼容 μ_fm)` / `xi`，并归一化到 `mu_fm`）。
- `ParamSpec([:T_fm, :mu_fm, :xi])`
  - 约束外参导数输入维度与顺序；非法参数名会抛 `ArgumentError`。
- `diff_target(:pressure | :entropy | :rho_norm | :energy)`
  - 目标注册入口；未知目标名抛 `ArgumentError`。
- `jacobian(ctx, target, params)` / `jacobian(ctx, targets, params)`
  - 始终返回矩阵；单目标单参数固定 `1x1`，多目标单参数返回 `N x 1`。
  - 当 `ParamSpec` 包含多个参数时，要求 backend 返回 `N x P` 矩阵（`N=目标数`, `P=参数数`）。

Issue #80 已接入默认数值 Jacobian（当前主路径聚焦 `FixedMu`），并保持 `N x P` 形状契约与输入校验语义。

## 试点统一导数服务（Issue #81）

在不改动 `ProblemSpec` 主契约前提下，solver/diff 提供两条试点服务入口，供脚本链与分析链以统一方式接入导数：

- `build_pilot_diff_context(result; mode, model, theta, spec_override=nothing, jacobian_backend=nothing)`
  - 对 `build_thermo_diff_context` 的薄封装；统一 `mu_fm`/`μ_fm` 别名语义。
- `eval_pilot_derivatives(ctx; target_names, param_names)`
  - 统一执行 `DiffTarget + ParamSpec + jacobian`，并返回 `jacobian` 与 `by_name` 命名映射。

该层定位为“试点评估服务”，用于验证可维护性、稳定性与性能边界；是否回收进 `ProblemSpec` 由后续 Phase-2 决策确定。
