# MesonThermoWorkflow

## 定位

`MesonThermoWorkflow` 是 `Models` 统一入口下的介子热力学 workflow。

当前第一版目标很窄：

- 复用 `Models.solve_gap_and_meson_point` 的平衡态与介子质量主链
- 在 `π/K` 通道上补齐 mesonic pressure
- 给出用户可消费的最小 EOS 合同：
  - `P_meson`
  - `P_quark_meanfield`
  - `P_total`
  - `entropy`
  - `epsilon`
  - `trace_anomaly`

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
- 第一版只承诺 `π/K` 通道，不主动扩到 full nonet

## 返回结构

完整 `solve_gap_and_*_meson_thermo_point` 返回值会在原 `meson_point` 上追加：

- `meson_thermo`
- `strict_bw_meson_thermo`
- `phase_shift_meson_thermo`

核心字段：

- `T_fm`
- `workflow`
- `channel_set`
- `P_meson`
- `P_quark_meanfield`
- `P_total`
- `delta_P_vs_no_meson`
- `P_meson_over_P_total`
- `entropy`
- `epsilon`
- `trace_anomaly`
- `equilibrium_converged`
- `phase_structure`
- `phase_shift_variant`

## CSV 合同

`Models.build_meson_thermo_contract_row` 会输出第一版最小平铺字段：

- `T_MeV`
- `muB_MeV`
- `workflow`
- `channel_set`
- `P_meson`
- `P_total`
- `P_quark_meanfield`
- `epsilon`
- `entropy`
- `trace_anomaly`
- `P_meson_over_P_total`
- `delta_P_vs_no_meson`
- `equilibrium_converged`
- `phase_structure`
- `phase_shift_variant`

## 当前边界

- `phase_structure` 当前固定写为 `unknown`，尚未与 phase pipeline 做正式联动
- 尚未提供 canonical temperature scan 脚本与结果目录落盘
- 尚未沉淀 regression baseline
- 尚未做 channel 扩张决策
