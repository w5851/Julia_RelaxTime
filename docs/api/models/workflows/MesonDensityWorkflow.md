# Models Meson Density Workflow

本文档从 `Models` 统一入口视角描述 meson density workflow。底层实现位于 [src/models/workflow_apps/MesonDensityWorkflow.jl](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/src/models/workflow_apps/MesonDensityWorkflow.jl)，领域细节页保留在 [docs/api/relaxtime/workflow/MesonDensityWorkflow.md](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/docs/api/relaxtime/workflow/MesonDensityWorkflow.md)。

## 首选入口

- `Models.solve_meson_density_from_meson_point`
- `Models.solve_gap_and_meson_density_point`
- `Models.solve_strict_bw_meson_density_from_meson_point`
- `Models.solve_gap_and_strict_bw_meson_density_point`
- `Models.solve_phase_shift_meson_density_from_meson_point`
- `Models.solve_gap_and_phase_shift_meson_density_point`

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

## `solve_strict_bw_meson_density_from_meson_point`

这是 strict BW workflow 的后处理入口。它继续只消费
`Models.solve_gap_and_meson_point` 的返回值，但把稳定粒子极限替换为 reduced strict BW 双积分后处理。

它适合：

- 已经拿到 `meson_results`，想在同一点上继续算有限宽度极点近似下的 `n_pi`、`n_K`、`K/π`
- 保持 workflow 唯一入口，不让 strict BW 又退回到脚本层私有实现

当前口径：

- `E(q)=sqrt(q^2+m^2)`
- `Gamma(q)=Gamma(q=0)`
- `omega = E(q) + Delta omega`

当前也支持：

- `stage = :stage1_reduced`
- `stage = :stage2_qpole`

其中 `stage2_qpole` 会在 `q` 网格上逐点重解复极点 `z_p(q)`，再回填 strict BW 积分。

## `solve_gap_and_strict_bw_meson_density_point`

这是 reduced strict BW 版本的完整闭环入口：先求平衡态与介子谱，再追加
`strict_bw_meson_density` 字段。

## `solve_phase_shift_meson_density_from_meson_point`

这是当前 Phase E3 最小 BU 相移双积分的后处理入口。它仍然只消费
`Models.solve_gap_and_meson_point` 的返回值，但把稳定粒子极限替换为相移双积分后处理。

它适合：

- 已经拿到 `meson_results`，想在同一点上继续算当前最小口径下的 `n_pi`、`n_K`、`K/π`
- 保持 workflow 唯一入口，不让分析脚本直接拼接 `quark_params` / `thermo_params` / 相移积分

当前默认口径：

- `qmax = 12`
- `q_nodes = 48`
- `omega_max = 10`
- `omega_nodes = 48`

## `solve_gap_and_phase_shift_meson_density_point`

这是相移双积分版本的完整闭环入口：先求平衡态与介子谱，再追加
`phase_shift_meson_density` 字段。

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
- `strict_bw_meson_density`
  - `m_pi`, `m_K`
  - `gamma_pi`, `gamma_K`
  - `n_pi`, `n_K`
  - `kpi_ratio`
  - `qmax`, `q_nodes`, `omega_max`, `omega_nodes`
- `phase_shift_meson_density`
  - `m_pi`, `m_K`
  - `n_pi`, `n_K`
  - `kpi_ratio`
  - `pi_density`, `k_density`
  - `qmax`, `q_nodes`, `omega_max`, `omega_nodes`, `eta`

## 当前边界

当前 workflow 已覆盖：

1. 稳定粒子极限
2. strict BW Stage1 reduced
3. strict BW Stage2 q-pole
3. Phase E3 最小相移双积分

其中 strict BW 与 phase-shift 入口当前仍处于严格受限状态：

- reduced strict BW：
  - 只消费 workflow 当前点给出的 `mass/gamma`
  - 对应 `stage = :stage1_reduced`
- q-pole strict BW：
  - 在 `q` 网格上逐点调用介子极点方程
  - 当前依赖 continuation seed 串行续算
  - 对应 `stage = :stage2_qpole`
- Phase E3 phase-shift：
  - 仅支持 `xi = 0`
  - 仅支持 `π/K` 聚合通道
  - 积分方案固定为 GL + 硬截断

后续 full strict BW 与更完整的 BU 扩展仍应沿同一 workflow 链继续后接，而不是回到脚本层重组流程。
