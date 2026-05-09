# MesonThermoWorkflow

## 定位

`MesonThermoWorkflow` 是 `Models` 统一入口下的介子热力学 workflow。

当前第一版目标很窄：

- 复用 `Models.solve_gap_and_meson_point` 的平衡态与介子质量主链
- 在双通道 meson set 上补齐 mesonic pressure；当前已验证 `π/K` 与 `π/sigma_pi`
- 给出用户可消费的最小 EOS 合同：
  - `P_meson`
  - `P_meson_qp`
  - `P_meson_ld`
  - `P_quark_meanfield`
  - `P_total`
  - `entropy`
  - `epsilon`
  - `trace_anomaly`

当前派生热力学量有两种实现状态：

- `stable` / `strict BW`：已切到 `Omega_total -> Models.model_thermo -> ForwardDiff`
- `phase-shift`：已接入同一总线入口，当前默认直接走 `omega_total_ad`

## 公开入口

- `Models.solve_meson_thermo_from_meson_point`
- `Models.solve_gap_and_meson_thermo_point`
- `Models.solve_strict_bw_meson_thermo_from_meson_point`
- `Models.solve_gap_and_strict_bw_meson_thermo_point`
- `Models.solve_phase_shift_meson_thermo_from_meson_point`
- `Models.solve_gap_and_phase_shift_meson_thermo_point`
- `Models.build_meson_thermo_contract_row`
- `Models.meson_thermo_workflow_module`

## 口径

当前实现支持三类 pressure 口径：

1. `stable_meson_pressure`
2. `strict_bw_stage1_reduced_pressure`
3. `phase_shift_current` / `phase_shift_gbu_reference`

其中：

- `stable` 与 `strict BW` 直接消费 `meson_results[:pi/:K].mass/gamma`
- `phase-shift` 复用现有 propagator / polarization / weighted-phase 链
- workflow 已验证 `π/K` 与 `π/sigma_pi` 两套双通道组合；当前仍不主动扩到 full nonet

## 返回结构

完整 `solve_gap_and_*_meson_thermo_point` 返回值会在原 `meson_point` 上追加：

- `meson_thermo`
- `strict_bw_meson_thermo`
- `phase_shift_meson_thermo`

核心字段：

- `T_fm`
- `workflow`
- `channel_set`
- `primary_channel`
- `secondary_channel`
- `P_meson`
- `P_meson_qp`
- `P_meson_ld`
- `P_quark_meanfield`
- `P_total`
- `delta_P_vs_no_meson`
- `P_meson_over_P_total`
- `P_pi_qp`
- `P_pi_ld`
- `P_K_qp`
- `P_K_ld`
- `P_primary`
- `P_secondary`
- `P_primary_qp`
- `P_primary_ld`
- `P_secondary_qp`
- `P_secondary_ld`
- `entropy`
- `epsilon`
- `trace_anomaly`
- `equilibrium_converged`
- `phase_structure`
- `phase_shift_variant`
- `ld_cutoff`
- `ld_cutoff_mode`
- `ld_threshold_mode`
- `thermo_derivation_mode`

## CSV 合同

`Models.build_meson_thermo_contract_row` 会输出第一版最小平铺字段：

- `T_MeV`
- `muB_MeV`
- `workflow`
- `channel_set`
- `primary_channel`
- `secondary_channel`
- `P_meson`
- `P_meson_qp`
- `P_meson_ld`
- `P_total`
- `P_quark_meanfield`
- `epsilon`
- `entropy`
- `trace_anomaly`
- `P_meson_over_P_total`
- `delta_P_vs_no_meson`
- `P_pi_qp`
- `P_pi_ld`
- `P_K_qp`
- `P_K_ld`
- `P_primary`
- `P_secondary`
- `P_primary_qp`
- `P_primary_ld`
- `P_secondary_qp`
- `P_secondary_ld`
- `equilibrium_converged`
- `phase_structure`
- `phase_shift_variant`
- `ld_cutoff`
- `ld_cutoff_mode`
- `ld_threshold_mode`
- `thermo_derivation_mode`

## Canonical Scan

仓库内现已提供首个 canonical `mu_B = 0` 扫描脚本：

- `scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl`

当前脚本约定：

- 默认 case：`phase_shift_current + pi/sigma_pi + mu_B=0 + xi=0`
- 输出目录内至少包含：
  - `scan.csv`
  - `README.md`
  - `effective_config.json`
  - `run_manifest.json`
- `scan.csv` 的核心平铺字段来自 `Models.build_meson_thermo_contract_row`
- 为避免 CSV 歧义，脚本落盘时会把 `channel_set` 编码为无逗号形式；解释双通道结果时仍优先使用 `primary_channel / secondary_channel`

## 当前边界

- `phase_structure` 当前固定写为 `unknown`，尚未与 phase pipeline 做正式联动
- `QP / LD` 目前只在 phase-shift pressure 口径下显式输出；stable / strict BW 仍为空值
- 兼容字段 `P_pi/P_K` 仍保留，但当第二通道切到 `sigma_pi` 时，更应使用 `primary/secondary` 字段解读合同
- 当前三类 meson thermo workflow 的总热力学派生量都应标记为 `omega_total_ad`
- canonical `mu_B = 0` 脚本与最小结果目录已落地，但图资产仍未在脚本内自动生成
- 尚未沉淀 regression baseline
- 尚未做 channel 扩张决策
