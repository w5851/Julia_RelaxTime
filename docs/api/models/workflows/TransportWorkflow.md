# Models Transport Workflow

本文档从 `Models` 统一入口视角描述 transport workflow。底层实现位于 [src/models/workflows/TransportWorkflow.jl](src/models/workflows/TransportWorkflow.jl#L1)，领域细节页保留在 [docs/api/relaxtime/workflow/TransportWorkflow.md](docs/api/relaxtime/workflow/TransportWorkflow.md)。

## 首选入口

- `Models.solve_gap_and_transport`
- `Models.solve_transport_from_equilibrium`

这两个入口都通过 [src/models/entrypoints.jl](src/models/entrypoints.jl#L49) 转发到 workflow 模块。

## `solve_gap_and_transport`

这是完整闭环入口：先求平衡态，再进入 τ、bulk 与输运系数计算。

实现签名见 [src/models/workflows/TransportWorkflow.jl](src/models/workflows/TransportWorkflow.jl#L453)。

它适合以下场景：

- 你只有 `(T_fm, mu_fm, xi)`，想直接得到平衡态与 `(eta, zeta, sigma)`
- 你希望 workflow 内部自行决定是否计算 `tau` 与 `bulk_coeffs`
- 你希望把平衡求解与输运后处理保持在同一入口下

返回结果通常包括：

- `equilibrium`
- `quark_params`, `thermo_params`
- `masses`
- `densities`
- `tau`, `tau_inv`, `rates`
- `bulk_coeffs`
- `transport`
- `reproducibility`（`physics_profile`, `physics_config_path`）

## `solve_transport_from_equilibrium`

这是后处理入口：在你已经拥有平衡态结果的前提下，只执行输运侧计算。

实现签名见 [src/models/workflows/TransportWorkflow.jl](src/models/workflows/TransportWorkflow.jl#L522)。

它适合：

- 避免重复 gap 求解
- 在同一平衡态上做多组 `tau`、`transport_config` 或 provider 对比
- 将平衡态求解与输运后处理解耦到不同批次或不同脚本中

## 输入分层

从 `Models` 入口角度，transport workflow 的关键输入可分三层：

- 平衡层：`T_fm`、`mu_fm`、`xi`、`seed_state`、`solver_kwargs`
- τ 计算层：`compute_tau`、`K_coeffs`、`tau_kwargs`
- 输运层：`tau`、`transport_config`、`transport_kwargs`、`prefer_energy_aniso`

这也是为什么本主题把 transport workflow 放进“职责核心”而不是只给一个函数列表。

本页只保留到“知道有哪些输入层级”为止；参数优先级、workflow 内部调用顺序与返回结构细节交给 `docs/api/relaxtime/workflow/TransportWorkflow.md`。

## 何时看原领域页

如果你需要更细的内容，例如：

- `prefer_energy_aniso` 的优先级细节
- `TransportIntegrationConfig` 的积分配置
- provider 字段与 `prepare_transport_provider` 的契约
- bulk 计算的性能提示

应继续阅读：

- [docs/api/relaxtime/workflow/TransportWorkflow.md](docs/api/relaxtime/workflow/TransportWorkflow.md)
- [docs/api/relaxtime/transport/CoreConcepts.md](docs/api/relaxtime/transport/CoreConcepts.md)

本页只负责说明它在 `Models` 公开入口中的定位。
