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
