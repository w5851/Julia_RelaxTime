# scripts/analysis/

探索性分析与诊断脚本目录（非 CI 路径，仅手动运行）。

## 当前内容

- `convergence/`：数值收敛性分析相关脚本与说明
- `relaxtime/export_phase_guided_publication_clean_v15_pdf.py`：以 `--authorization <record.json>` 将已接受的 v15 PNG 包补充为 76 张 PDF；校验 1867 点来源、同字节 PNG 重绘和字体/矢量合同。使用新 `--case`，拒绝覆盖；`--font-sources <map.json>` 仅允许同字节历史字体。保留原 PNG 路径，不复制审阅图及源码树。
- `scan_mott_meson_vs_xi_mu0.jl`：在 `mu=0` 条件下扫描 Mott 点与介子质量随 `xi` 的变化
- `mott_reference_mapping.jl`：生成/检查 Mott 参考映射数据
- `relaxtime/bu_kernel_gate_utils.jl`：纯数值 finite-window BU 恒等式 helper，显式分解 bulk、上下边界项与 closure；由 unit 测试直接覆盖，不是正式数密度入口
- `relaxtime/compare_bu_derivative_vs_byparts_e5.jl`：在同一 smooth/unwrapped 相移分支上比较 derivative、by-parts bulk、有限窗口 boundary 与 total；旧 `density_*_byparts` 列只是 bulk-only 历史映射，`*_closure_*` 才是门禁量
- `relaxtime/audit_bu_meson_density_literature_alignment.jl`：面向 BU 介子数密度文献对齐审查的 fixed-point / phase normalization / charged `μ_K` 规则诊断；charged 输出在同一点覆盖 `bu2020_mu_s_0p2`、`friesen2019_mu_s_0p55` 以及 `K^+/π^+`、`K^-/π^-`
- `relaxtime/meson_conserved_charge_feedback_utils.jl`：外层介子守恒荷修正的纯数值 helper，包含 charged B/Q/S 映射、affine residual 与带 evaluation budget 的阻尼有限差分 Newton；不直接求 PNJL 或 BU
- `relaxtime/meson_conserved_charge_feedback_runtime.jl`：partial-feedback analysis runtime；集中候选 evaluator、缓存与 gap/BU/outer 分项计时，不是稳定公共入口
- `relaxtime/meson_conserved_charge_outer_feedback_spike.jl`：固定 `(T,mu_B)` 的 quark-only→partial-feedback 单点诊断；每个 `(mu_Q,mu_S)` 候选重求五平均场，并同时计算 `pi+/-`、`K+/-` current-BU 密度
- `relaxtime/scan_meson_conserved_charge_feedback_freezeout.jl`：读取 hot-start benchmark 中位耗时，沿 `default/baseline_freezeout` 选择 7/5/3 点稀疏扫描 `K+/pi+`、`K-/pi-`；输出仅为 diagnostic

输运 v1–v13 旧图层、旧 xi001 数值及其专用生成器已退役，历史恢复见
[输运专题入口](../../docs/analysis/relaxtime/phase_guided_transport/README.md)。
当前物理分析使用新采样结果，入口为 [Mott 分析](../../docs/analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/mott_local_bumps_20261028.md)。

## 使用说明

- 这些脚本用于研究与诊断，不作为测试准入门禁。
- 建议从仓库根目录执行并显式指定项目环境。
- finite-window kernel gate 不认证 folded/no-anomalous 等非光滑 display 口径，也不修改正式 `phase_shift_meson_number_density` 语义。
- charged 同点消融默认使用显式 `density_policy=:x_min_cut` 诊断延拓，并同时输出 `strict_requested_window_status`、`m_M-μ_M` 与 effective lower bound；不得把该诊断直接视为守恒荷平衡解或 production 结果。
- 外层守恒荷 spike 同样固定为显式 `density_policy=:x_min_cut`，且介子尚未进入 mean-field stationarity；即使 verdict 为 `full-feedback-candidate`，也只表示值得拉起完整热力学反馈研究，不是自洽平衡解认证。
- 可用 `BU_AUDIT_RUN_NEUTRAL=false` 只运行 charged 消融；`BU_AUDIT_Q_NODES`、`BU_AUDIT_OMEGA_NODES`、`BU_AUDIT_QMAX`、`BU_AUDIT_OMEGA_MIN/MAX` 和 `BU_AUDIT_BOSE_X_MIN` 可控制低成本诊断网格。

```bash
julia --project=. scripts/analysis/scan_mott_meson_vs_xi_mu0.jl
julia --project=. scripts/analysis/mott_reference_mapping.jl
julia --project=. scripts/analysis/relaxtime/compare_bu_derivative_vs_byparts_e5.jl
julia --project=. scripts/analysis/relaxtime/audit_bu_meson_density_literature_alignment.jl
julia --project=. scripts/analysis/relaxtime/meson_conserved_charge_outer_feedback_spike.jl
julia --project=. scripts/analysis/relaxtime/scan_meson_conserved_charge_feedback_freezeout.jl
```
