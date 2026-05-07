# `Models.run_crossover_meson_density_scan`

本页说明沿手征 crossover line 运行介子数密度 workflow 的统一入口。

## 定位

`Models.run_crossover_meson_density_scan` 负责：

- 先调用 `Models.build_crossover_line` 生成 `mu_q -> T_crossover` 路径
- 在这条路径上复用统一 `MesonDensityWorkflow`
- 通过 flavor chemical profile 固定 `\mu_u`、`\mu_d`、`\mu_s`
- 通过 meson chemical profile 固定 `\mu_\pi`、`\mu_K`、`d_\pi`、`d_K`
- 输出与 freeze-out / external-path 入口对齐的 meson-density CSV

它的角色不是新的 solver，而是“**内部 crossover 路径驱动**”的正式 workflow 入口。

## 何时使用它

适用于：

- 目标是 `Friesen 2019` 一类 charged `K/\pi` crossover-line reproduction
- 你希望由仓库内部统一的 crossover locator 先生成路径
- 需要让 stable / strict BW / current BU / generalized BU 共用同一条内部 path shell

若路径已经由外部 CSV 明确给出，应改用：

- `Models.run_external_path_meson_density_scan`

若路径来自 freeze-out baseline 或其 proxy path，应改用：

- `Models.run_freezeout_meson_density_scan`

## 当前支持的物理口径

当前统一入口已支持：

- `:stable`
- `:strict_bw_stage1`
- `:strict_bw_stage2`
- `:phase_shift_current`
- `:phase_shift_gbu`

它们共用同一条 crossover 路径，仅切换下游 meson-density kernel。

## 入口形态

核心关键字包括：

- `mu_min_MeV`
- `mu_max_MeV`
- `T_min_MeV`
- `T_max_MeV`
- `n_mu`
- `method`
- `variable`
- `flavor_chemical_profile_name`
- `meson_chemical_profile_name`
- `regime`
- `output_path`
- `overwrite`
- `p_num`, `t_num`, `max_iter`

并按 regime 选择性接受：

- `stable_num_q_nodes`
- strict BW 的 `qmax / q_nodes / omega_max / omega_nodes`
- phase-shift 的 `qmax / q_nodes / omega_min / omega_max / omega_nodes / eta`

## continuity 契约

这条入口当前同样复用：

- `MesonMassWorkflow` 返回的 `continuation_state`

因此 crossover 路径上的点推进与 freeze-out / external-path 入口保持同一治理方式，而不是在脚本层单独维护种子状态。

## 输出合同

当前 CSV 至少固定以下信息：

- 路径与定位字段：
  - `muq_MeV`
  - `muB_MeV`
  - `T_MeV`
  - `T_over_muB`
  - `crossover_method`
  - `crossover_variable`
  - `crossover_converged`
  - `crossover_derivative`
- chemical profile 字段：
  - `flavor_chemical_profile`
  - `meson_chemical_profile`
  - `mu_u_MeV`
  - `mu_d_MeV`
  - `mu_s_MeV`
  - `mu_pi_MeV`
  - `mu_K_MeV`
- 结果字段：
  - `m_pi_MeV`
  - `m_K_MeV`
  - `gamma_pi_MeV`
  - `gamma_K_MeV`
  - `n_pi`
  - `n_K`
  - `kpi_ratio`

## 当前边界

这条入口当前只解决：

- 内部 crossover 路径如何进入统一 meson-density workflow
- charged/flavor chemical profile 如何在该路径上稳定复用
- 不同 meson-density regime 如何共享同一条 path shell

它不负责：

- 外部文献路径提取
- freeze-out baseline 参数化
- target comparison / plotting

这些步骤应分别交给 external-path、freeze-out 或 analysis 层处理。

## 相关主题

- [Overview.md](Overview.md)
- [FreezeoutMesonDensityScan.md](FreezeoutMesonDensityScan.md)
- [ExternalPathMesonDensityScan.md](ExternalPathMesonDensityScan.md)
