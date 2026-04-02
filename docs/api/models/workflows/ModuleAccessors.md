# Module Accessors 与参数适配层

本文档说明 `Models` 统一入口中那些“公开但非首选业务入口”的导出：模块访问器与参数适配层。

## 涉及导出

- `transport_workflow_module`
- `meson_workflow_module`
- `workflow_param_adapters_module`
- `pnjl_module`

这些导出位于 [src/models/entrypoints.jl](src/models/entrypoints.jl#L61)。

## `transport_workflow_module`

该访问器返回 transport workflow 模块本体，并在返回前校验关键 API 是否存在。

作用：

- 让调试脚本或高级调用方显式拿到底层模块
- 在聚合入口层集中做 API 存在性检查
- 避免调用方自己散落地依赖 include 逻辑

## `meson_workflow_module`

该访问器返回 meson workflow 模块本体，并确保 `solve_gap_and_meson_point` 已定义。

作用与 transport 访问器相同，但目标是 meson workflow 模块。

## `pnjl_module`

该访问器直接返回当前 `Models` 模块本体。它更偏向聚合入口/调试用途，而不是常规业务接口。

## `workflow_param_adapters_module`

该访问器返回参数适配模块，对应 [src/models/workflows/WorkflowParamAdapters.jl](src/models/workflows/WorkflowParamAdapters.jl#L1)。

其中最核心的能力是：

- `normalize_quark_params`
- `normalize_thermo_params`
- `as_relaxtime_inputs`

它的职责是把 `QuarkParams` 与 `ThermoParams` 校验并归一化为 workflow 内部稳定使用的 NamedTuple 结构。

当前约定：

- 参数适配层统一以 `as_relaxtime_inputs` 为单一入口，不再保留 `as_legacy_inputs`。

## 为什么这些接口仍然应该写进主题文档

原因不是鼓励所有用户去直接调用它们，而是为了避免两类误解：

- 误把它们当成普通业务入口
- 完全忽略它们的存在，导致高级调用方绕过统一入口自己维护一套模块查找与参数适配逻辑

因此，本页的目标是把它们从“未被解释的公开导出”变成“有边界说明的进阶接口”。
