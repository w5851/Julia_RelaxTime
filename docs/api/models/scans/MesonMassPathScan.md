# `Models.run_meson_mass_path_scan`

本页说明沿特征线路径运行 meson-mass workflow 的统一入口。

## 定位

`Models.run_meson_mass_path_scan` 负责：

- 接收已经排好序的路径点列
- 维护路径推进时的 `MesonMassWorkflow.continuation_state`
- 在 fixed-`μ_B` 与 fixed-`σ` 两类点之间选择最小必要求解链
- 统一写出带路径元数据的介子质量 CSV

它是 path/workflow 层正式入口，不引入新的 meson 物理公式。

## 当前支持的路径家族

第一版正式支持两类路径：

- `freezeout`
  - 输入点列直接提供 `T_MeV` 与 `muB_MeV`
  - 典型调用来自 `Models.run_freezeout_meson_mass_scan`
- `isentropic`
  - 输入点列提供 `T_MeV` 与 `sigma_target`
  - 入口先用 `FixedSigma` 求出平衡点，再复用 meson workflow
  - 典型调用来自 `Models.run_isentropic_meson_mass_scan`

当前没有把所有 path family 统一抽成更大的共享框架；复用范围只收敛到 meson-mass path scan 本身。

## 入口形态

核心关键字包括：

- `points`
- `xi`
- `mesons`
- `output_path`
- `overwrite`
- `p_num`, `t_num`
- `max_iter`
- `solver_backend`
- `solver_kwargs`
- `mass_kwargs`
- `continuation=:auto|:legacy|:palc`
- `path_diagnostics_cb(segment_indices, path_result)`（可选）

等熵路径的 `:auto` 对至少 21 点、固定 `sigma_target`、标准 Models 后端且未提供
`solver_kwargs` 的连续分段使用 PALC。温度顺序保持输入顺序，不跨 path group、segment
或 sigma 目标拼接。PALC 求得的平衡态直接传给介子后处理，避免再次求 gap。
短路径及自定义求解参数沿用逐点链路；`:legacy` 可显式选择原流程。
返回摘要的 `palc_paths` 记录实际使用 PALC 的分段数。单个未收敛点继续输出失败行。
这里的 `max_iter` 仍用于介子步骤，PALC 的平衡态校正保持旧约束求解器的最低 300 次上限。

`points` 的最小 contract 为：

- 公共字段：
  - `T_MeV`
  - `path_family`
  - `path_profile`
  - `path_segment`
  - `path_point_index`
  - `path_order_key`
  - `path_label`
- fixed-`μ_B` 路径额外需要：
  - `muB_MeV`
- fixed-`σ` 路径额外需要：
  - `sigma_target`

## 输出合同

CSV 固定先写路径元数据：

- `path_family`
- `path_profile`
- `path_segment`
- `path_point_index`
- `path_order_key`
- `path_label`
- `freezeout_profile`
- `sigma_target`
- `sqrt_s_NN_GeV`

再写热力学点与平衡态诊断：

- `T_MeV`
- `muB_MeV`
- `muq_MeV`
- `xi`
- `equilibrium_converged`
- `equilibrium_iterations`
- `equilibrium_residual_norm`
- `Phi`
- `Phibar`
- `m_u`, `m_d`, `m_s`

随后复用现有 meson-mass 网格扫描口径：

- `M_*`
- `Gamma_*`
- `converged_*`
- `residual_*`
- `threshold_* / gap_*`
- 混合道的 `thr_uu / thr_ss / thr_min` 与 `gap_uu / gap_ss / gap_min`

## 便捷入口

当前正式暴露两个 path-family 包装入口：

- `Models.run_freezeout_meson_mass_scan`
- `Models.run_isentropic_meson_mass_scan`

它们分别负责：

- freeze-out profile / path profile 到点列的展开
- fixed-`σ` profile 与温度点列的展开

再统一落到 `Models.run_meson_mass_path_scan`。
