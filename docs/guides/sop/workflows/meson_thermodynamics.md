# 介子热力学计算 SOP

## 1. 目的与适用范围

本 SOP 用于通过 `Models` 统一 workflow 计算介子 pressure 及其与 mean-field 热力学量组合后的 EOS，并运行 canonical `mu_B=0` phase-shift temperature scan。

覆盖三类内部口径：

- stable meson pressure；
- strict BW Stage1 reduced pressure；
- phase-shift current / GBU reference pressure，含 QP/LD 分解。

当前 canonical 执行链固定为 `phase_shift_current + pi/sigma_pi + mu_B=0 + xi=0`。

## 2. 非适用范围

- 不覆盖 full pseudoscalar/scalar nonet。
- 不把 `pi/sigma_pi` 结果解释为所有介子通道总和。
- 不提供已经成立的 external literature validation；当前只有内部 regression 与 plot-review 证据。
- 不负责介子数密度或 `K/pi` 数密度比。
- 不把 `phase_structure=unknown` 当作相结构判定结果。

## 3. 权威入口

统一 API：

- `Models.solve_meson_thermo_from_meson_point`
- `Models.solve_gap_and_meson_thermo_point`
- `Models.solve_strict_bw_meson_thermo_from_meson_point`
- `Models.solve_gap_and_strict_bw_meson_thermo_point`
- `Models.solve_phase_shift_meson_thermo_from_meson_point`
- `Models.solve_gap_and_phase_shift_meson_thermo_point`
- `Models.build_meson_thermo_contract_row`

专题 canonical CLI：

- `scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl`

该脚本是具有 integration/regression 覆盖的专题入口，不属于 Phase A 核心稳定 CLI 白名单。脚本层只负责编排、落盘和 provenance，不应重新拼接总热力学公式。

## 4. 物理口径、单位与参数约束

- CLI `tmin/tmax/tstep`：MeV。
- canonical chemical potential：`mu_B=0 MeV`，内部 `mu_q=0 fm^-1`。
- `xi`：无量纲。
- `qmax/omega_min/omega_max`：`fm^-1`。
- `P_meson`、`P_quark_meanfield`、`P_total`、`epsilon`：自然单位 `fm^-4`。
- `entropy`：自然单位 `fm^-3`。
- `P/T^4`、`P_meson/P_total` 与 trace-anomaly 归一化量按输出合同解释。

双通道解释优先使用：

- `primary_channel/secondary_channel`；
- `P_primary/P_secondary` 及其 QP/LD 分量。

当第二通道为 `sigma_pi` 时，兼容字段 `P_K_*` 只是第二通道 legacy 槽位，不能按 kaon 物理解读。

Pressure reference：

- `raw_absolute`：保留原始 pressure 零点；
- `vacuum_subtracted_mu0`：只平移 mean-field/total pressure 零点，不改变 `P_meson`。

两种 reference 是口径选择，不是数值收敛层级。

## 5. 输入配置及优先级

当前 canonical CLI 以 `ScanOptions` 默认值与显式 CLI 为合同，没有独立专题 TOML：

1. 脚本默认值定义 canonical case；
2. 显式 CLI 覆盖默认值；
3. 底层模型和物理常数继续来自 `config/models/pnjl/` 与 `config/physics/`；
4. 最终有效参数写入 `effective_config.json` 和 `run_manifest.json`。

默认治理重点：

- `scheme=current`；
- `mesons=pi,sigma_pi`；
- `ld_cutoff_mode=match_model_lambda`；
- `ld_threshold_mode=omega_lt_q`；
- `allow_legacy_fd_fallback=false`。

`match_qmax` 只作为 legacy/cutoff sensitivity 对照；启用 legacy FD fallback 必须显式记录。

## 6. 环境与版本冻结

在仓库根目录记录：

```powershell
julia --version
git rev-parse HEAD
git status --short
```

运行要求：

- 使用根 `--project=.`；
- 稳定环境优先经 `run_with_sysimage.ps1/.sh` wrapper；
- manifest 必须记录 Julia、git commit、argv、config hash 和 artifacts；
- 正式 case 不得依赖未记录的本地 pressure reference 或临时代码修改。

## 7. Smoke 预检

结构和 canonical 数值 smoke：

```powershell
julia --project=. tests/integration/relaxtime/test_phase_shift_meson_thermo_scan_smoke.jl
```

等价的独立 diagnostic 命令：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl --outdir data/outputs/results/sop_smoke/meson_thermo --overwrite --tmin 210 --tmax 210 --tstep 5 --q-nodes 8 --omega-nodes 8 --qmax 4 --omega-max 3 --p-num 8 --t-num 4 --max-iter 20 --pressure-reference vacuum-subtracted-mu0 --reference-t-mev 5
```

该命令验证 equilibrium、meson point、phase-shift pressure、AD 派生量、CSV contract 和 sidecars；低积分分辨率不能用于正式 EOS 结论。

## 8. 收敛性验证

正式 temperature scan 前至少比较：

- equilibrium `p_num/t_num/max_iter`；
- phase-shift `qmax/q_nodes`；
- `omega_min/omega_max/omega_nodes`；
- `eta`；
- LD cutoff 及 QP/LD 分区稳定性；
- T 网格步长；
- vacuum reference temperature 与 fallback 是否触发。

核心比较量：

- `P_meson/P_meson_qp/P_meson_ld`；
- `P_total`；
- `entropy/epsilon/trace_anomaly`；
- primary/secondary channel 分量；
- equilibrium convergence 和路径连续性。

`current` 与 `gbu_reference` 是不同物理权重口径，不能作为彼此的网格收敛层。`raw_absolute` 与 `vacuum_subtracted_mu0` 也不能用作误差估计。

## 9. 正式计算命令

正式 case 使用独立命名目录和收敛证据支持的参数：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl --outdir data/outputs/results/relaxtime/meson_thermo/<case_id> --scheme current --mesons pi,sigma_pi --xi 0 --tmin <T_min_MeV> --tmax <T_max_MeV> --tstep <T_step_MeV> --p-num <converged_p_num> --t-num <converged_t_num> --max-iter <converged_max_iter> --qmax <converged_qmax> --q-nodes <converged_q_nodes> --omega-min <omega_min> --omega-max <converged_omega_max> --omega-nodes <converged_omega_nodes> --pressure-reference vacuum-subtracted-mu0 --reference-t-mev <T_ref_MeV>
```

本 SOP 不给出通用 production 节点数。实际口径以 `effective_config.json` 为准，并必须链接对应 convergence evidence。

## 10. 输出目录与产物合同

canonical case 至少包含：

- `scan.csv`；
- `README.md`；
- `effective_config.json`；
- `run_manifest.json`。

`scan.csv` 必须保留 temperature/chemical-potential、workflow/channel、pressure/EOS、QP/LD、equilibrium、reference 和积分配置字段。

Sidecar 语义：

- `effective_config.json`：本次实际 solver、reference 和 phase-shift 参数；
- `run_manifest.json`：git、Julia、argv、config hash、summary 与 artifacts；
- README：case 口径、成功/失败计数和兼容字段解释。

## 11. Regression / Validation 验收

最小层级：

- Unit kernel：`tests/unit/relaxtime/test_meson_thermodynamics.jl`；
- Unit workflow：`tests/unit/models/test_meson_thermo_workflow.jl`；
- Integration：`tests/integration/relaxtime/test_phase_shift_meson_thermo_scan_smoke.jl`；
- Fixedpoint regression：`tests/regression/relaxtime/test_meson_thermo_fixedpoint_regression.jl`；
- Canonical path regression：`tests/regression/relaxtime/test_meson_thermo_canonical_muB0_path_regression.jl`；
- Plot-review regression：`tests/regression/relaxtime/test_meson_thermo_plot_review_case_regression.jl`。

当前没有 meson thermo external validation gate。内部 regression 通过只能证明当前项目口径未漂移，不能证明与某篇文献定量一致。

## 12. 失败点、断点续算与重跑

- 默认 resume 按 `(T_MeV, muB_MeV, xi)` 跳过已有点；
- 参数、channel、scheme 或 reference 改变时必须使用新 case，不能向旧 CSV 追加；
- `--overwrite` 会重建 case 内核心 sidecars，应在明确重跑时使用；
- reference state 不收敛时不得自动退回 raw pressure；应保留异常并检查候选 reference temperature；
- AD 失败默认应终止，不能静默启用 legacy FD fallback；
- 单点失败不应从路径中手工删除后继续宣称完整 scan。

## 13. Diagnostic 与 Formal Production 的边界

Diagnostic-only：

- 低节点 smoke；
- `allow_legacy_fd_fallback=true`；
- legacy `ld_cutoff_mode=match_qmax` 对照；
- 未完成 q/omega/LD/T-grid 收敛；
- channel 或 pressure reference 仍在比较。

Internal production/regression-grade 至少需要：

1. 完整 sidecars；
2. equilibrium 与 phase-shift convergence evidence；
3. fixedpoint/path regression；
4. QP/LD 与 AD 派生量合同审计；
5. 明确 channel set 和 pressure reference。

即使满足内部 production 条件，在 external reference 的通道、归一化、模型差异和 target 未审计前，也不得写成 literature-validated。

## 14. 后处理与作图

图像后处理应消费冻结的 `scan.csv`。现有 plot-review 工具：

- `scripts/relaxtime/build_meson_thermo_plot_review_case.py`
- `scripts/relaxtime/build_meson_thermo_paperlike_pi_compare.py`

图像写入 `data/outputs/figures/relaxtime/meson_thermo/<case_id>/`，不得在绘图脚本中重新计算或修改 pressure 数据。

## 15. 关联公式、API 和测试

- [Models MesonThermoWorkflow](../../../api/models/workflows/MesonThermoWorkflow.md)
- [MesonThermodynamics API](../../../api/relaxtime/meson_thermo/MesonThermodynamics.md)
- [公共科学计算生命周期](../common_scientific_run.md)
- `scripts/relaxtime/run_phase_shift_meson_thermo_scan.jl`
- `tests/regression/relaxtime/test_meson_thermo_fixedpoint_regression.jl`

## 16. 最后验证记录

- 验证日期：2026-07-10
- 验证范围：canonical CLI、CSV/sidecar contract 和现有 regression 映射
- 执行命令：见第 7 节及 authority map
- 状态：通过；canonical 数值与 sidecar testset 30/30，非法 scheme 拒绝 testset 1/1
- 备注：external validation 仍未建立
