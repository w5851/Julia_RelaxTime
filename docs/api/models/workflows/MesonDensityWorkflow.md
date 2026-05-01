# Models Meson Density Workflow

本文档从 `Models` 统一入口视角描述 meson density workflow。底层实现位于 [src/models/workflow_apps/MesonDensityWorkflow.jl](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/src/models/workflow_apps/MesonDensityWorkflow.jl)，领域细节页保留在 [docs/api/relaxtime/workflow/MesonDensityWorkflow.md](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/docs/api/relaxtime/workflow/MesonDensityWorkflow.md)。

## 首选入口

- `Models.solve_meson_density_from_meson_point`
- `Models.solve_gap_and_meson_density_point`

## `solve_meson_density_from_meson_point`

这是后处理入口。它只消费 `Models.solve_gap_and_meson_point` 的返回值，不重复求平衡态或介子极点。

它适合：

- 已经拿到 `meson_results`，想在同一个点上继续算 `n_pi`、`n_K`、`K/π`
- 保持“介子数密度挂在介子质量 workflow 之后”的主链结构
- 避免脚本层自行拼接 `thermo_params`、介子质量与数密度公式

## `solve_gap_and_meson_density_point`

这是完整闭环入口：先求平衡态与介子谱，再在其返回值上追加 `meson_density` 字段。

它适合：

- 调用方只有 `(T_fm, mu_fm, xi)`，想直接得到介子质量与稳定粒子极限数密度
- 希望脚本只走 `Models` 统一入口，而不是直接 include workflow 实现文件

## 返回结构

`solve_gap_and_meson_density_point` 返回值保留 `solve_gap_and_meson_point` 的核心字段：

- `equilibrium`
- `quark_params`, `thermo_params`
- `meson_results`

并新增：

- `meson_density`
  - `m_pi`, `m_K`
  - `n_pi`, `n_K`
  - `kpi_ratio`
  - 简并因子、化学势与积分节点配置

## 当前边界

当前 workflow 只覆盖稳定粒子极限：

- `π/K`
- `μ_B = 0` 主线默认口径
- `K/π` 比值

BW 与 BU 相移扩展仍应沿同一 workflow 链继续后接，而不是回到脚本层重组流程。
