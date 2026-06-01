# run 脚本功能目录（全量）

本文件记录仓库内当前 `run_*.jl` 脚本的功能、定位与迁移建议，用作“脚本清理 -> 前端固化 -> CI 手动触发固化”的统一索引。

说明：

- 白名单稳定入口仍以 `docs/guides/scripts/README.md` 为准。
- 本目录是全量盘点，不等同于“都推荐直接对外使用”。
- `Status` 字段用于治理，不代表脚本立即可删除。

## 状态口径

- `core-candidate`：核心入口候选，优先进入前端/CI 固化。
- `domain-candidate`：专题入口候选，稳定但偏专题研究或后处理。
- `gate`：回归/门禁/分析脚本，优先用于 CI 或维护流程。
- `deprecate-candidate`：疑似过时候选，建议先迁移替代后再删除。
- `internal-module-candidate`：命名为 `run_`，但实质偏内部模块。
- `analysis-experimental`：分析/一次性实验脚本，默认不作为稳定入口。
- `internal-authoritative`：内部主链路依赖脚本（可被 include），不直接作为外部入口。

## 目录

| Script | Function | Status | Suggested Replacement / Notes |
| --- | --- | --- | --- |
| `scripts/models/run_unified_scan.jl` | `Models` 主链统一扫描/工作流入口（`scan tmu|trho`, `workflow phase`） | `core-candidate` | 当前 authoritative 的 PNJL 扫描 CLI 入口 |
| `scripts/relaxtime/run_gap_transport_scan.jl` | 平衡求解 + tau + RTA 输运主扫描 | `core-candidate` | 参数最全，适合作为工作流后端能力 |
| `scripts/relaxtime/run_gap_meson_mass_scan.jl` | 平衡求解 + 介子质量/宽度 + Mott 阈值扫描 | `core-candidate` | 支撑 Mott 相关产线 |
| `scripts/relaxtime/run_relaxtime_orchestrator.jl` | relaxtime 工作流编排入口 | `core-candidate` | 更适合 CI 手动触发入口 |
| `scripts/relaxtime/run_mott_phase_scan.jl` | Mott 相扫描主入口（含 manifest） | `core-candidate` | 可与派生 CSV / plot 脚本联动 |
| `scripts/relaxtime/run_combined_meson_density_scan.jl` | 组合式介子数密度扫描入口，按 scan path × density regime 输出同一批点的 stable / strict BW / phase-shift 多口径 CSV、README 与 SVG | `domain-candidate` | 当前实现 `--path tmu`，支持多个固定 `mu_q` 路径与 `phase_display=fold_0_pi` FIG3-like 诊断；用作四口径生产/审计桥接入口，后续可扩展 freezeout/crossover/external path |
| `scripts/analysis/relaxtime/render_combined_meson_density_fig3_like.py` | 从组合式介子数密度 CSV 渲染 FIG3-like heatmap PNG | `analysis-helper` | 只负责绘图，不重新计算物理量；正式计算仍以 `run_combined_meson_density_scan.jl` 为入口 |
| `scripts/relaxtime/run_bu2020_meson_density_audit_scan.jl` | BU2020/temp7 介子数密度主线审计扫描，输出 stable / strict BW / phase-shift / no-anomalous 多口径 CSV 与 README | `domain-candidate` | 研究审计入口；不进入稳定白名单，默认通过 `Models` workflow 而非复制 temp7 |
| `scripts/pnjl/run_conserved_charge_susceptibilities.jl` | 守恒荷广义磁化率与累积量入口 | `domain-candidate` | 已在脚本指南稳定入口中出现 |
| `scripts/pnjl/run_aniso_phase_template.jl` | 各向异性相图实验模板（扫描 + 相结构 + 绘图） | `domain-candidate` | 偏实验模板 |
| `scripts/pnjl/run_magnetic_point.jl` | 磁场单点计算 | `domain-candidate` | 磁场专题 |
| `scripts/pnjl/run_magnetic_eb_scan.jl` | 磁场 eB 轴扫描 | `domain-candidate` | 磁场专题 |
| `scripts/pnjl/run_magnetic_stability_scan.jl` | 磁场稳定性批量检查 | `domain-candidate` | 磁场专题 |
| `scripts/relaxtime/run_mott_phase_derived_csv.jl` | Mott 扫描派生字段 CSV 生成 | `domain-candidate` | 后处理脚本 |
| `scripts/relaxtime/run_mott_phase_plot_modes.jl` | Mott 图像模式产出（mode A/B） | `domain-candidate` | 绘图包装脚本 |
| `scripts/relaxtime/run_offline_transport_patch.jl` | 对质量异常点离线补算并输出 patch | `domain-candidate` | 维护/后处理工具 |
| `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl` | 手动组合产物工作流（cross_section/temperature_scan_muB0_xi0/fixed_temperature_xi_scan_muB0；兼容旧别名 `plan_a`/`plan_b`） | `domain-candidate` | 可作为迁移期桥接入口 |
| `scripts/relaxtime/run_scan.jl` | relaxtime 单入口分发器（`gap-transport` / `tau-vs-t` / `manual-workflow`） | `domain-candidate` | 作为 relaxtime 子命令包装层使用 |
| `scripts/relaxtime/run_xi_smoothness_batch.jl` | xi smoothness 批处理扫描入口 | `analysis-experimental` | 研究分析资产，不作为稳定用户白名单 |
| `scripts/analysis/convergence/run_convergence.jl` | 收敛性分析扫描 | `gate` | 分析脚本 |
| `scripts/analysis/relaxtime/meson_mixed/run_globalmin_window_experiment.jl` | mixed-meson global-min 窗口实验脚本 | `analysis-experimental` | 研究实验资产，不作为稳定用户白名单 |
| `scripts/dev/run_prune_wave_gate.jl` | 裁剪波次门禁（快照 + smoke + 回滚） | `gate` | 维护治理脚本 |
| `scripts/relaxtime/run_outputs_2026_02_05.jl` | 历史日期任务的一次性产物脚本 | `analysis-experimental` | deprecate-candidate，替代方向：`run_manual_relaxation_scan_workflow.jl` / orchestrator |

## 维护规则

- 当脚本状态变更（如 `deprecate-candidate -> removed`）时，先更新本目录，再执行代码删除。
- 当新增稳定入口时，同时更新 `docs/guides/scripts/README.md` 白名单与本目录。
