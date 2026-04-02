# Models 入口工作流总览

本文档是 `Models` 统一工作流主题的首选入口，优先面向“想通过 `Models` 快速完成一个业务工作流调用”的使用者。

## 何时使用本主题

当你希望通过 `Models` 聚合入口完成以下任一任务时，应从这里开始：

- 先求平衡态，再计算输运系数
- 基于已有平衡态结果，直接计算输运系数
- 先求平衡态，再计算介子质量、宽度与 Mott 阈值

## 首选公开入口

相关导出位于 [src/models/entrypoints.jl](src/models/entrypoints.jl#L13) 与 [src/models/Models.jl](src/models/Models.jl#L80)。

优先使用的 API：

- `Models.solve_gap_and_transport`
- `Models.solve_transport_from_equilibrium`
- `Models.solve_gap_and_meson_point`

进阶接口：

- `Models.transport_workflow_module`
- `Models.meson_workflow_module`
- `Models.workflow_module_for`
- `Models.workflow_param_adapters_module`
- `Models.pnjl_module`

完整导出基线见 [generated/Exports.md](docs/api/models/workflows/generated/Exports.md)。

## 三条典型路径

### 1. 平衡态 + 输运一体化

`Models.solve_gap_and_transport` 适合“给定 `(T_fm, mu_fm, xi)`，直接得到平衡态与 `(eta, zeta, sigma)`”的场景。

更详细的输入合同与性能提示见 [docs/api/relaxtime/workflow/TransportWorkflow.md](docs/api/relaxtime/workflow/TransportWorkflow.md)。

### 2. 基于已知平衡态直接算输运

`Models.solve_transport_from_equilibrium` 适合你已经有 `equilibrium` 或已经整理好 `quark_params` 与 `thermo_params`，不想重复求解 gap 的场景。

### 3. 平衡态 + 介子质量/Mott 阈值一体化

`Models.solve_gap_and_meson_point` 适合“给定一个点，直接求出多个介子通道的质量、宽度、阈值与 gap”的场景。

更详细的通道与返回结构见 [docs/api/pnjl/MesonMassWorkflow.md](docs/api/pnjl/MesonMassWorkflow.md)。

## 最短示例

### `solve_gap_and_transport`

```julia
res = Models.solve_gap_and_transport(
    0.15,
    0.0;
    xi=0.0,
    compute_tau=false,
    compute_bulk=false,
)
```

### `solve_gap_and_meson_point`

```julia
res = Models.solve_gap_and_meson_point(
    0.15,
    0.0;
    xi=0.0,
    mesons=(:pi, :K, :eta, :eta_prime),
)
```

## 不是首选入口的导出

以下导出虽然属于公开表面，但不应作为多数用户的首页首选：

- `transport_workflow_module`
- `meson_workflow_module`
- `workflow_param_adapters_module`
- `pnjl_module`

它们主要用于：

- 获取底层 workflow 模块本体
- 明确参数归一化职责（`as_relaxtime_inputs`）
- 做聚合入口边界检查或调试

这些能力的职责边界见 [CoreConcepts.md](docs/api/models/workflows/CoreConcepts.md)。
