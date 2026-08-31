# `Models.run_freezeout_meson_density_scan`

本页说明沿 freeze-out baseline 路径运行介子数密度 workflow 的统一入口。

## 定位

`Models.run_freezeout_meson_density_scan` 负责：

- 接收 `\sqrt{s_{NN}}` 路径参数
- 通过 freeze-out profile 映射到 `(T,\mu_B)`
- 通过 path profile 决定是否在 baseline 之上施加 proxy / stitched 路径策略
- 通过 flavor chemical profile 固定 `\mu_u`、`\mu_d`、`\mu_s`
- 通过 meson chemical profile 固定 `\mu_\pi`、`\mu_K`、`d_\pi`、`d_K`
- 复用 `MesonDensityWorkflow` 的既有物理核，输出路径扫描结果

它仍是 **path/workflow 层薄封装**，不是新的独立 solver。

## 当前支持的物理口径

当前统一入口已支持以下 regime：

- `:stable`
- `:strict_bw_stage1`
- `:strict_bw_stage2`
- `:phase_shift_current`
- `:phase_shift_gbu`

它们分别复用：

- `Models.solve_gap_and_meson_density_point`
- `Models.solve_gap_and_strict_bw_meson_density_point`
- `Models.solve_gap_and_phase_shift_meson_density_point`

## profile 分层

这条入口当前显式拆成四层 profile：

1. freeze-out baseline profile  
   负责 `\sqrt{s_{NN}} -> (T,\mu_B)`
2. path profile  
   负责：
   - baseline freeze-out 是否直接使用
   - 或是否进入更高层 proxy / stitched path
3. flavor chemical profile  
   负责：
   - `\mu_u`
   - `\mu_d`
   - `\mu_s`
4. meson chemical profile  
   负责：
   - `\mu_\pi`
   - `\mu_K`
   - `d_\pi`
   - `d_K`
   - `pi_label` / `k_label`

这意味着当前入口已经能表达：

- 总 `K/\pi`
- charge-resolved 但 `\mu_\pi = 0`
- charge-resolved 且固定 `\mu_\pi = 100, 134.5 MeV`
- strange flavor phenomenology 的最小显式 profile

## continuity 契约

当前连续性机制不再复用 `TmuScan` 的 seed cache，而是复用：

- `MesonMassWorkflow` 返回的 `continuation_state`

这保证 freeze-out meson-density 路径仍然通过统一 workflow 入口推进，而不是在脚本层平行拼装。

## 入口形态

核心关键字包括：

- `sqrt_s_NN_values`
- `xi_values`
- `freezeout_profile_name`
- `path_profile_name`
- `flavor_chemical_profile_name`
- `meson_chemical_profile_name`
- `regime`
- `output_path`
- `traversal`
- `overwrite`
- `resume`
- `p_num`, `t_num`
- `solver_kwargs`
- `mass_kwargs`

并按 regime 选择性接受：

- `stable_num_q_nodes`
- strict BW 的 `qmax / q_nodes / omega_max / omega_nodes`
- phase-shift 的 `qmax / q_nodes / omega_min / omega_max / omega_nodes / eta`
- phase-shift 的 `phase_shift_phase_anchor` 与 `phase_shift_omega_measure`
  - 分别控制高能相位锚定和 `domega/pi` / `domega/(2pi)` 显式测度
- phase-shift 的 `phase_shift_density_policy`
  - 默认 `:strict_normal_domain`，用于新策略下显式标记 `omega <= mu_M` 的 unsafe Bose domain
  - 历史 freeze-out 数值基线或兼容性检查可显式使用 `:excitation_only_E_gt_mu`

## 输出合同

当前 CSV 至少固定以下信息：

- 路径键：
  - `sqrt_s_NN_GeV`
  - `muB_MeV`
  - `xi`
  - `freezeout_profile`
  - `path_profile`
  - `path_segment`
  - `flavor_chemical_profile`
  - `meson_chemical_profile`
  - `regime`
- 热力学点：
  - `T_MeV`
  - `muq_MeV`
  - `mu_u_MeV`
  - `mu_d_MeV`
  - `mu_s_MeV`
- 带电/聚合口径：
  - `pi_label`
  - `k_label`
  - `charge_resolved`
  - `mu_pi_MeV`
  - `mu_K_MeV`
  - `d_pi`
  - `d_K`
- 数值结果：
  - `m_pi_MeV`
  - `m_K_MeV`
  - `gamma_pi_MeV`
  - `gamma_K_MeV`
  - `n_pi`
  - `n_K`
  - `kpi_ratio`

## 当前边界

当前入口已经足够作为 charged / freeze-out 子链的最小统一 workflow 入口，但还没有覆盖：

- stitched critical + constant-`T` path
- 文献 Figure 4 级别的完整 path reconstruction

但它已经具备显式 path-profile 层，后续可以在不改 workflow 主核的前提下继续扩展这些策略。

## 相关主题

- [FreezeoutScan.md](FreezeoutScan.md)
- [Overview.md](Overview.md)
- [../../../reference/formula/relaxtime/meson_density/MesonDensity_Freezeout路径参数化.md](../../../reference/formula/relaxtime/meson_density/MesonDensity_Freezeout%E8%B7%AF%E5%BE%84%E5%8F%82%E6%95%B0%E5%8C%96.md)
- [../../../reference/formula/relaxtime/meson_density/MesonDensity_Flavor化学势Profile.md](../../../reference/formula/relaxtime/meson_density/MesonDensity_Flavor%E5%8C%96%E5%AD%A6%E5%8A%BFProfile.md)
