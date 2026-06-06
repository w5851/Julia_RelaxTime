# Models Meson Workflow

本文档从 `Models` 统一入口视角描述 meson workflow。底层实现位于 [src/models/workflows/MesonMassWorkflow.jl](src/models/workflows/MesonMassWorkflow.jl#L1)，领域细节页保留在 [docs/api/pnjl/MesonMassWorkflow.md](docs/api/pnjl/MesonMassWorkflow.md)。

## 首选入口

- `Models.solve_gap_and_meson_point`
- `Models.solve_meson_point_from_equilibrium`

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
- `root_sign_flipped`：若求解器落到负号等价 pole，workflow 已按物理约定输出 `mass >= 0`，并将该字段置为 `true`
- 非混合通道的 `threshold`, `gap`
- 混合通道的 `threshold=(uu, ss, min)`, `gaps=(uu, ss, min)`

此外，workflow 返回值还包含：

- `continuation_state`

它打包了后续温度/参数续算所需的最小状态：

- `equilibrium_seed_state`
- `meson_seed_state`
- `mixed_seed_tracking_state`

扫描脚本应优先传递这一对象，而不是在脚本层分别维护多份 seed。

对于需要热力学稳定相而不是沿上一温度点追踪亚稳分支的扫描，可传入
`equilibrium_seed_strategy=Models.MultiSeed()`。此时平衡态 FixedMu 求解会在每个点执行多初值压力选优；
`continuation_state` 仍可用于介子根的连续性，但不会强制 equilibrium 沿上一点分支延拓。

## `solve_meson_point_from_equilibrium`

这是 equilibrium-to-meson adapter：调用方先通过 `Models.solve(...)` 或其他
稳定 equilibrium source 得到 `SolverResult`，再把它交给该入口继续求介子质量、
宽度、Mott 阈值与 gap。adapter 不重新运行 gap solver。

典型用途：

- `FixedMu` 路径的 parity / diagnostic 复用
- `FixedAsymmetricRho` 这类 density-constrained equilibrium source 的介子后处理
- 扫描脚本中把“上游路径策略”和“介子数密度口径”分开组合

最小调用形态：

```julia
model = Models.create_model(:PNJL)
eq = Models.solve(model, Models.FixedAsymmetricRho(0.05, 0.876, 0.0), T_fm)
meson_point = Models.solve_meson_point_from_equilibrium(
    eq,
    T_fm;
    mesons=(:pi_plus, :K_plus),
)
```

## 默认通道与复用接口

workflow 模块同时导出：

- `DEFAULT_MESONS`
- `build_equilibrium_params`
- `solve_meson_point_from_equilibrium`

这意味着调用方既可以直接走 `solve_gap_and_meson_point`，也可以先求平衡态，再把参数转换为 meson 计算所需的稳定输入结构。

## 何时看原领域页

如果你需要更细的内容，例如：

- 默认通道列表
- 与 Fortran 对照的历史说明
- 扫描脚本位置与数值健壮性细节

应继续阅读 [docs/api/pnjl/MesonMassWorkflow.md](docs/api/pnjl/MesonMassWorkflow.md)。本页只负责说明它在 `Models` 公开入口中的角色。
