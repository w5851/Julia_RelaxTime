# Models 工作流职责核心

本文档面向维护者或需要理解 `Models` 统一入口边界的读者，解释 transport workflow、meson workflow、模块访问器与参数适配层的职责分工。

## 为什么需要单独的 workflow 主题

`src/models/entrypoints.jl` 做的不只是简单 re-export。它把“调用方应该使用什么入口”“底层 workflow 模块在哪里”“哪些参数需要先归一化”这些问题集中到 `Models` 统一表面上。

因此，这个主题的核心不是某个单一算法，而是职责分层。

## 三层职责分工

### 1. 业务工作流入口

这一层是大多数用户真正该调用的入口：

- `solve_gap_and_transport`
- `solve_transport_from_equilibrium`
- `solve_gap_and_meson_point`

它们的共同特征：

- 面向完整任务闭环，而不是单个局部算法
- 从 `Models` 聚合入口直接可用
- 返回可继续消费的结构化结果，而不是只暴露底层中间量

### 2. 模块访问器

这一层负责把聚合入口与底层 workflow 模块连接起来：

- `transport_workflow_module`
- `meson_workflow_module`
- `pnjl_module`

这些导出本质上是边界与调试接口，不是首页业务入口。它们的价值主要在于：

- 让调用方或调试脚本显式拿到底层模块
- 在入口层做必需 API 是否存在的检查
- 减少散落的动态查找或 include 逻辑

### 3. 参数适配层

`workflow_param_adapters_module` 对应 [src/models/workflows/WorkflowParamAdapters.jl](src/models/workflows/WorkflowParamAdapters.jl#L1)。

这层的职责是：

- 校验 `QuarkParams` 与 `ThermoParams`
- 转换成 workflow 内部稳定使用的 NamedTuple 结构
- 提供 `as_legacy_inputs` 这类已弃用的兼容桥接能力

兼容说明：

- `as_legacy_inputs` 仅为历史调用方保留，已进入弃用窗口。
- workflow 主路径与常规测试已切换到 `as_relaxtime_inputs`。
- 该旧入口不再属于推荐使用的公开表面，只应在显式兼容场景下访问。

它的意义不在于“再提供一个业务入口”，而在于把参数归一化从 workflow 主逻辑里抽离出来。

## transport workflow 的边界

transport workflow 的实现位于 [src/models/workflows/TransportWorkflow.jl](src/models/workflows/TransportWorkflow.jl#L1)。

应重点理解三层输入：

- 平衡层：`T_fm`、`mu_fm`、`xi`、`seed_state`、`solver_kwargs`
- τ 计算层：`compute_tau`、`K_coeffs`、`tau_kwargs`
- 输运层：`tau`、`transport_kwargs`、`prefer_energy_aniso`

这解释了为什么 `solve_gap_and_transport` 和 `solve_transport_from_equilibrium` 需要并存：

- 前者负责完整闭环
- 后者负责复用已有平衡态，避免重复求解

### provider 桥接边界

`Models` 还导出了：

- `transport_provider`
- `TransportProvider`
- `prepare_transport_provider`

但它们的主题角色并不是“再提供一层业务入口”，而是把 `Models` 平衡态结果桥接到 transport 积分层可消费的 provider 契约。

因此：

- 本主题只说明它们为什么存在于 `Models` 导出表面；
- 具体字段、回退规则、`prefer_energy_aniso` 路径与 `ctx` 语义，主说明放在 `../../relaxtime/transport/CoreConcepts.md`；
- workflow 页只保留“何时需要 provider、provider 如何在 workflow 中被预备”的桥接说明。

## meson workflow 的边界

meson workflow 的实现位于 [src/models/workflows/MesonMassWorkflow.jl](src/models/workflows/MesonMassWorkflow.jl#L1)。

其职责是把：

- 平衡态求解
- 介子质量/宽度求解
- Mott 阈值与 gap 计算

串联成一个单点工作流。

这个入口的价值主要在于：

- 统一返回 `equilibrium`、`quark_params`、`thermo_params`
- 在多个 meson 通道之间复用状态
- 将“物理通道扫描”组织成可直接消费的结果字典

## 为什么模块访问器必须写文档

在当前入口导出中，真正未被人工文档吸收的恰恰是：

- `transport_workflow_module`
- `meson_workflow_module`
- `workflow_param_adapters_module`

如果不单独说明，读者会面临两个误区：

- 误以为这些导出是普通业务入口
- 完全不知道这些导出存在，进而在调试和扩展时绕过统一入口自己做一套动态加载逻辑

因此，本主题必须显式把它们定位为“公开但偏进阶的边界接口”。