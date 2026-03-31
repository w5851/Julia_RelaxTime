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
| `scripts/pnjl/run_tmu_scan.jl` | PNJL T-mu 参数扫描（断点续扫、CSV） | `core-candidate` | 可作为前端/CI 首批固化对象 |
| `scripts/pnjl/run_dense_trho_scan.jl` | 高密度 T-rho 加密扫描 | `core-candidate` | 与 `run_adaptive_trho_scan.jl` 形成基础 + 自适应补点 |
| `scripts/pnjl/run_adaptive_trho_scan.jl` | 基于既有结果的自适应补点扫描 | `core-candidate` | 依赖已有 Trho 扫描结果 |
| `scripts/relaxtime/run_gap_transport_scan.jl` | 平衡求解 + tau + RTA 输运主扫描 | `core-candidate` | 参数最全，适合作为工作流后端能力 |
| `scripts/relaxtime/run_gap_meson_mass_scan.jl` | 平衡求解 + 介子质量/宽度 + Mott 阈值扫描 | `core-candidate` | 支撑 Mott 相关产线 |
| `scripts/relaxtime/run_relaxtime_orchestrator.jl` | relaxtime 工作流编排入口 | `core-candidate` | 更适合 CI 手动触发入口 |
| `scripts/relaxtime/run_mott_phase_scan.jl` | Mott 相扫描主入口（含 manifest） | `core-candidate` | 可与派生 CSV / plot 脚本联动 |
| `scripts/pnjl/run_conserved_charge_susceptibilities.jl` | 守恒荷广义磁化率与累积量入口 | `domain-candidate` | 已在脚本指南稳定入口中出现 |
| `scripts/pnjl/run_aniso_phase_template.jl` | 各向异性相图实验模板（扫描 + 相结构 + 绘图） | `domain-candidate` | 偏实验模板 |
| `scripts/pnjl/run_magnetic_point.jl` | 磁场单点计算 | `domain-candidate` | 磁场专题 |
| `scripts/pnjl/run_magnetic_eb_scan.jl` | 磁场 eB 轴扫描 | `domain-candidate` | 磁场专题 |
| `scripts/pnjl/run_magnetic_stability_scan.jl` | 磁场稳定性批量检查 | `domain-candidate` | 磁场专题 |
| `scripts/relaxtime/run_mott_phase_derived_csv.jl` | Mott 扫描派生字段 CSV 生成 | `domain-candidate` | 后处理脚本 |
| `scripts/relaxtime/run_mott_phase_plot_modes.jl` | Mott 图像模式产出（mode A/B） | `domain-candidate` | 绘图包装脚本 |
| `scripts/relaxtime/run_offline_transport_patch.jl` | 对质量异常点离线补算并输出 patch | `domain-candidate` | 维护/后处理工具 |
| `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl` | 手动组合产物工作流（cross_section/plan_a/plan_b） | `domain-candidate` | 可作为迁移期桥接入口 |
| `scripts/pnjl/run_scan_fixedpoint_regression.jl` | PNJL 扫描固定点回归 | `gate` | 回归门禁，不建议前端按钮化 |
| `scripts/pnjl/run_magnetic_fixedpoint_regression.jl` | 磁场固定点回归 | `gate` | 回归门禁 |
| `scripts/relaxtime/run_transport_fixedpoint_regression.jl` | 输运固定点回归（legacy/models 对照） | `gate` | 回归门禁 |
| `scripts/relaxtime/run_total_cross_section_fixedpoint_regression.jl` | 总截面固定点回归 | `gate` | 回归门禁 |
| `scripts/analysis/convergence/run_convergence.jl` | 收敛性分析扫描 | `gate` | 分析脚本 |
| `scripts/dev/run_prune_wave_gate.jl` | 裁剪波次门禁（快照 + smoke + 回滚） | `gate` | 维护治理脚本 |
| `scripts/relaxtime/run_outputs_2026_02_05.jl` | 历史日期任务的一次性产物脚本 | `analysis-experimental` | deprecate-candidate，替代方向：`run_manual_relaxation_scan_workflow.jl` / orchestrator |
| `scripts/pnjl/run_solver_smoke.jl` | PNJL 单点/小样本 smoke 自检 | `deprecate-candidate` | 替代方向：单元 smoke + `run_tmu_scan.jl` |
| `scripts/pnjl/run_solver_tmu.jl` | PNJL T-mu 小网格稳定性自检 | `deprecate-candidate` | 替代方向：`run_tmu_scan.jl` |
| `scripts/relaxtime/workflow/cross_section_orchestrated.jl` | orchestrator 依赖的截面扫描内部模块 | `internal-authoritative` | 已迁移到 workflow 子目录，作为内部模块保留 |

## 维护规则

- 当脚本状态变更（如 `deprecate-candidate -> removed`）时，先更新本目录，再执行代码删除。
- 当新增稳定入口时，同时更新 `docs/guides/scripts/README.md` 白名单与本目录。
