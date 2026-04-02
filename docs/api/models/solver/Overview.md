# Models 核心求解入口总览

本页优先回答“我应该调用哪个入口”，而不是先展开底层实现细节。

## 最小稳定入口

对大多数新调用方，核心入口顺序如下：

1. `Models.create_model`：先选择模型对象
2. `Models.solve_gap`：直接求平衡态平均场状态
3. `Models.solve`：在 PNJL 约束模式下求完整 `SolverResult`
4. `Models.solve_multi`：当你明确需要多初值筛选物理解时使用

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
- 自定义后处理或隐式导数链路

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

它不是常规入口的默认替代品，而是“我明确知道这里有分支竞争”的增强入口。

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
- `create_implicit_gap_solver` 与相关导数接口

## 与其它主题的关系

- phase 主题、workflow 主题、scans 主题都建立在本主题定义的求解合同之上
- 如果你先理解本主题，再去看 [../workflows/README.md](../workflows/README.md) 或 [../scans/README.md](../scans/README.md)，会更容易看清哪些是“业务流程入口”，哪些是“核心求解公共表面”

## 迁移期契约说明（Plan-B / Gate A）

- `SolverResult` 已从固定 `SVector{5}/SVector{3}` 扩展为向量泛型字段，结果消费侧不应再假设状态维度恒为 5、化学势维度恒为 3。
- scan 结果适配层（`TmuScan` / `TrhoScan`）已按动态向量透传 `x_state` 与 `mu_vec`，为后续 schema-driven 主链路铺路。
- 当需要显式状态布局时，优先使用 `ModelStateSchema` 系列 API（`schema_for_model`、`flatten_state`、`unflatten_state`）。
