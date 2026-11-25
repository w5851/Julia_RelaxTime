# 1c 残留：`issue130_phase_reference_v1` 与 v2/v3/v9 渲染链退役

> 状态：`deferred_task`，独立治理重构任务。本批次**不做**，理由是它的引用面远超一次收尾改动的范围，且会触碰现行 C2 正式数据源的消费者。

## 为什么单独拆出来

2026-11-25 图形退役批次（Batch 1c 候选）在盘点时发现该残留的引用面为 **60+ 处**，包含测试、审计证据树与构建脚本三层，混入收尾批次会把风险引到已定版的相位参考数据上。按 Gemini 裁定：先保 C2 正式数据源安全，再择机单开专题。

## 现状盘点（2026-11-25，`main = b22ceecd`）

### 待退役对象

| 对象 | 规模 | 说明 |
| --- | --- | --- |
| `data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v1/` | 3 文件 / 14.07 MB | 仅图件，无对应 `results/` 目录 |
| `scripts/analysis/pnjl/render_issue130_phase_surface_v2.py` | 17.3 KB | 渲染链 |
| `scripts/analysis/pnjl/render_issue130_phase_surface_v3.py` | 18.9 KB | 渲染链 |
| `scripts/analysis/pnjl/render_issue130_phase_surface_v9.py` | 34.4 KB | 渲染链（v8 与之并存，先一并判定） |

### 引用方分层

1. **测试（5 个，硬阻断）**：`test_pnjl_phase_surface_{v2,v3,v9}.py`、`test_issue130_phase_reference_v2.py`、`test_pnjl_phase_reference_{adapter,import}.py`
2. **证据树（整棵）**：`docs/analysis/pnjl/phase_reference/` 下 4 轮 legacy audit、promotion gate、runtime parity v1/v2、physical deletion 记录、`docs/analysis/pnjl/phase_surface_series/figure_layer/`
3. **构建与审计脚本**：`audit_issue130_phase_reference_{legacy_retirement,promotion,runtime_consumers}.py`、`build_issue130_phase_reference_v2.py`、`import_issue130_phase_reference.py`、`audit_issue130_rs_runtime_parity.jl`
4. **治理与配置**：`config/governance/task_tracks.toml`（`issue130-phase-reference-import` 条目的证据）、`.gitattributes`

## 步骤（按此顺序，不并行）

1. **消费者矩阵落盘**：用现有 `issue130_phase_reference_legacy_audit_v4/tables/consumer_occurrences.csv` 重新扫描一次，产出现状矩阵，标注「可删 / 需迁移 / 保留作历史证据」。
2. **判定 v1 是否仍被运行时读取**：核对 `import_issue130_phase_reference.py` 与 adapter 测试的目标指向；若 v1 只被测试与证据引用，即可进入退役。
3. **迁移而非删除**：把证据树中的 v1 引用改写为 `note:` 墓碑（记录删除提交），沿用本仓 2026-11-25 台账墓碑写法。
4. **退役测试**：删除 3 个 v2/v3/v9 渲染测试前，先确认对应脚本是否还有其他消费方；`test_pnjl_phase_reference_{adapter,import}.py` 若服务于现行 v2 数据源则保留。
5. **物理删除**：按 `issue130_phase_reference_physical_deletion_v1` 既有流程执行，补 `.gitattributes -text` 清理与归档策略条目。
6. **回归**：`test_active_vector_packages.py`、`check_data_output_path_guard.jl`、`check_task_ledger.jl`、`check_active_docs_governance.jl` 全绿。

## 完成判据

- 全仓 `git grep` 对 `issue130_phase_reference_v1` 与 `render_issue130_phase_surface_v[239]` 只剩墓碑与历史记录；
- 现行 C2 相位参考（`data/outputs/results/pnjl/phase_reference/...accepted` 与 `native_issue130_c2_*` 图包）的哈希链不受影响；
- 治理三闸门与 Python/Julia 相关测试全绿；
- 台账证据无 `file does not exist`。

## 不要做的事

- 不要顺手改动 `phase_reference` v2/v3 现行数据与 `data/inputs/pnjl/...` 输入；
- 不要在未完成步骤 1–2 的矩阵之前删除任何 v1 图件；
- 不要把本任务与输运基线刷新（#381 讨论）混在同一分支。
