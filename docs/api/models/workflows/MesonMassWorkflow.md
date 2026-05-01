# Models Meson Workflow

本文档从 `Models` 统一入口视角描述 meson workflow。底层实现位于 [src/models/workflows/MesonMassWorkflow.jl](src/models/workflows/MesonMassWorkflow.jl#L1)，领域细节页保留在 [docs/api/pnjl/MesonMassWorkflow.md](docs/api/pnjl/MesonMassWorkflow.md)。

## 首选入口

- `Models.solve_gap_and_meson_point`

该入口通过 [src/models/entrypoints.jl](src/models/entrypoints.jl#L57) 转发到 workflow 模块。

## `solve_gap_and_meson_point`

这是一个单点工作流入口：先求平衡态，再对若干介子通道求质量、宽度、Mott 阈值与 gap。

实现签名见 [src/models/workflows/MesonMassWorkflow.jl](src/models/workflows/MesonMassWorkflow.jl#L169)。

它适合：

- 对单个 `(T_fm, mu_fm, xi)` 点做介子谱分析
- 在多个 meson 通道间复用同一个平衡态
- 需要统一拿到 `equilibrium`、`quark_params`、`thermo_params` 与 `meson_results`

## 返回结构

该 workflow 返回的核心结果包括：

- `equilibrium`
- `quark_params`, `thermo_params`
- `meson_results::Dict{Symbol, NamedTuple}`

其中 `meson_results` 会按通道给出：

- `mass`
- `gamma`
- `converged`
- `residual`
- 非混合通道的 `threshold`, `gap`
- 混合通道的 `threshold=(uu, ss, min)`, `gaps=(uu, ss, min)`

此外，workflow 返回值还包含：

- `continuation_state`

它打包了后续温度/参数续算所需的最小状态：

- `equilibrium_seed_state`
- `meson_seed_state`
- `mixed_seed_tracking_state`

扫描脚本应优先传递这一对象，而不是在脚本层分别维护多份 seed。

## 默认通道与复用接口

workflow 模块同时导出：

- `DEFAULT_MESONS`
- `build_equilibrium_params`

这意味着调用方既可以直接走 `solve_gap_and_meson_point`，也可以先求平衡态，再把参数转换为 meson 计算所需的稳定输入结构。

## 何时看原领域页

如果你需要更细的内容，例如：

- 默认通道列表
- 与 Fortran 对照的历史说明
- 扫描脚本位置与数值健壮性细节

应继续阅读 [docs/api/pnjl/MesonMassWorkflow.md](docs/api/pnjl/MesonMassWorkflow.md)。本页只负责说明它在 `Models` 公开入口中的角色。
