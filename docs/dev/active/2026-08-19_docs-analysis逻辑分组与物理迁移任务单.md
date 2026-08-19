# `docs/analysis` 逻辑分组与物理迁移任务单

创建日期：2026-08-19

状态：`active_batch_migration`。本任务只整理诊断分析树的 namespace 和入口，不改变数值语义、正式产物或 production promotion 状态。

## 1. Scope Lock

允许修改：

- `docs/analysis/` 下的分组 README、索引和显式目录迁移；
- 与迁移直接相关的分析脚本默认输出路径、active task evidence 和治理 ledger 路径；
- 每批迁移的文件集合、路径和 hash 审阅记录。

明确不修改：

- solver、Maxwell、C1/C2/reference、transport、CSV/JSON 数值内容和正式 figures；
- 诊断 candidate 的 verdict、phase-reference promotion gate 或 RS transport production；
- 历史生成包的 manifest、checksum、execution log 和生成时 provenance。

## 2. Completed Batches

- [x] 建立 `docs/analysis` 和 PNJL 分组索引。
- [x] 将 C1/C2 surface views、C2 audits/follow-ups 和 CEP/Maxwell evidence 按逻辑线分组。
- [x] 将 phase-reference 决策证据迁入 `docs/analysis/pnjl/phase_reference/`，提交 `8c41919a`。
- [x] 将独立算法可行性审计迁入 `docs/analysis/pnjl/algorithmic_feasibility/`。
- [x] 将 phase-guided transport v1/v2 连续证据线迁入 `docs/analysis/relaxtime/phase_guided_transport/`，新增分组入口并同步 live 脚本路径。
- [x] 将独立 PNJL/Mott 证据线从 `docs/analysis/pnjl_mott/` 迁入 `docs/analysis/pnjl/mott/`，新增分组入口并同步当前索引。
- [x] 将 figure asset registry 治理包从 `docs/analysis/figure_asset_registry_v1/` 迁入 `docs/analysis/governance/figure_asset_registry_v1/`，同步 live plotting 入口和当前文档引用。

每批迁移均使用显式路径、单独审阅、`R100`/hash 核对和独立提交；`raw_curve_archive_v1/` 继续作为独立外部归档指针保留在 PNJL 根目录。

### Batch review: phase-guided transport

- source roots：`phase_guided_transport_p128_xi001_analysis/`（22 files，1,447,033 bytes）和 `phase_guided_transport_v2_pole_sensitive_rendering/`（53 files，8,352,342 bytes）；
- pre-migration inventory SHA-256：`caac52b4cd9a7c503dce218298bf1389b90d756bd0aadcb0670e9829a7231fc4`、`f34b9b095af6ee35bb142913626f72ca9f731749e4d7a3b70eb7f9fc9b4c6c67`；
- migration boundary：仅改变物理 namespace 和 live 入口，包内生成时 manifest、图 manifest、旧路径快照、CSV/JSON/PNG 和 provenance 不重写；
- pre-existing metadata note：v1 root `manifest.json` 已有 `outputs` 2/21 hash mismatch；本批不修复，已登记到 metadata repair follow-up。

### Batch review: PNJL/Mott evidence line

- source root：`docs/analysis/pnjl_mott/`；destination root：`docs/analysis/pnjl/mott/`；既有 payload 为 13 files、1,157,737 bytes；迁移前 source-root inventory SHA-256：`cefabf25274504eaaa0608f2a2b94eb2d91c1c8e29d913ce75980c52302bbe0e`；
- logical boundary：`xi` 依赖、Mott 温度、介子谱、复极点机制和文献定位；与 Issue #130 phase-reference 决策线分离；
- migration boundary：既有 13 个分析/图像/CSV 文件整体迁移；新增 `pnjl/mott/README.md` 作为分组入口；仅修正 supporting note 中两个指向旧 namespace 的路径文字；不重算、不重绘、不改 CSV/PNG 科学内容；
- post-migration verification：destination payload 仍为 13 files、1,157,738 bytes；12 个既有文件为 `R100`，supporting note 为 `R097` 且仅有 namespace 文字差异；destination-root inventory SHA-256：`1787eac6aa611765d77fb05f0e7138c2b1dab12b868e379d652b5c858d84c2dd`；
- historical boundary：`docs/dev/archived/**` 审批记录和 `docs/analysis/governance/figure_asset_registry_v1/` 生成/清理快照中的旧路径不改写；
- metadata boundary：本批不处理任何 manifest/checksum mismatch；metadata repair 继续作为独立 follow-up。

### Batch review: figure asset registry governance

- source root：`docs/analysis/figure_asset_registry_v1/`；destination root：`docs/analysis/governance/figure_asset_registry_v1/`；既有 payload 为 7 files、1,432,945 bytes；迁移前 source-root inventory SHA-256：`83ec8a59933e3e4cf8877d1080f366e8650f1306a954fa9349ecb73525cfe51f`；
- logical boundary：历史 figure inventory、作者审核、cleanup preflight、retirement 和 relocation provenance；这是治理材料，不与科学分析包合并；
- migration boundary：7 个既有治理文件整体迁移；新增 `docs/analysis/governance/README.md`；README 只更新当前重建命令的 canonical path；三个 plotting live 入口、SOP、authority map、当前历史分析 README 和活动清理任务单同步路径；
- historical boundary：`asset_registry.json`、`cleanup_preflight_v1.json`、`retirement_execution_v1.json`、`relocation_execution_v1.json` 和作者审核记录的生成时路径、hash、状态均不重写；
- post-migration verification：destination payload 仍为 7 files、1,432,967 bytes；6 个既有文件为 `R100`，README 为 `R94` 且仅有 canonical path 文字变化；destination-root inventory SHA-256：`4ff83974b7c709b404f624f62bb8925fbc5e21849cea338a2d0db3287b527df3`；不执行删除、重绘或数值重算；
- metadata boundary：本批不处理任何 manifest/checksum mismatch；metadata repair 继续作为独立 follow-up。

## 3. Required Follow-up: Analysis Metadata Repair

- [ ] 单独建立 metadata 修复批次，处理 phase-reference 包及 phase-guided transport v1 包中既有的 `manifest.json` / `checksums.sha256` 与当前文件字节不一致问题。

当前已知范围（2026-08-18 迁移前审阅）：

- `phase_reference_current_state_freeze_v1`：`output_files` 8/8 不匹配，`checksums.sha256` 9/9 不匹配；
- `phase_reference_limited_evidence_audit_v1`：`output_files` 11/15 不匹配；
- `phase_reference_manual_overlay_promotion_audit_v1`：`output_files` 10/10 不匹配。

修复批次的验收条件：先冻结当前文件字节和旧 hash，明确是否重算 metadata；只允许修改 manifest/checksum 等元数据，不重生成 CSV/JSON/PNG，不改变历史 verdict；修复后重新验证文件集合、输入/输出 hash、路径引用和 provenance，并单独提交。该 follow-up 不因目录迁移完成而自动关闭。

## 4. Validation Minimum

- 迁移前后文件集合和 SHA-256 对照；
- JSON 解析、脚本 AST/入口 smoke 和 `git diff --check`；
- `check_task_ledger.jl`、`check_docs_consistency.jl`、`check_active_docs_governance.jl` 和 `check_script_entrypoints.jl`；
- 每批只暂存本批路径、索引、脚本和 task evidence，并使用符合近期历史的 `docs:` commit。
