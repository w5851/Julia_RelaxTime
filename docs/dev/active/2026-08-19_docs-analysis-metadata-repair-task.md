# `docs/analysis` metadata 修复任务单

创建日期：2026-08-19

状态：`triaged`。这是 `docs/analysis` namespace 迁移任务的独立 `required_follow_up`，不因目录迁移完成而自动关闭。

## 目标

审阅并修复 phase-reference 与 phase-guided transport 诊断包中 `manifest.json` / `checksums.sha256` 与当前文件字节不一致的问题，恢复可验证的输出 metadata 合同。

## 已知范围

迁移前审阅（2026-08-18）记录：

- `phase_reference_current_state_freeze_v1`：`output_files` 8/8 不匹配，`checksums.sha256` 9/9 不匹配；
- `phase_reference_limited_evidence_audit_v1`：`output_files` 11/15 不匹配；
- `phase_reference_manual_overlay_promotion_audit_v1`：`output_files` 10/10 不匹配。

目录已迁移到 `docs/analysis/pnjl/phase_reference/`；上述计数是包内输出文件相对其自身 manifest 的 mismatch，不是路径迁移引入的变化。

本批迁移前追加审阅：

- `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_p128_xi001_analysis/manifest.json` 的 `outputs` 有 2/21 个条目与当前文件字节不一致（`README.md`、`figures/plot_manifest.json`）；该 mismatch 在目录迁移前已存在；
- v2 包 `manifest.json` 的 `outputs` 为 0/52 mismatch；不新增修复范围。

## Scope Lock

允许修改：

- phase-reference 包和 phase-guided transport v1 包内的 `manifest.json`、`checksums.sha256` 或明确属于 hash registry 的元数据；
- 为验证修复所需的只读脚本或审计记录。

明确不修改：

- CSV、JSON 结果表、PNG、README 的科学内容和历史 verdict；
- solver、Maxwell、C1/C2/reference、transport、正式图和外部 raw archive；
- 生成时 provenance、输入文件 hash 和 evidence 语义，除非有单独说明并经复核。

## 验收条件

- 修复前保存每个包的文件集合、当前字节 SHA-256、旧 manifest/checksum 和 mismatch 清单；
- 明确每个 mismatch 是 stale metadata、换行/序列化漂移还是实际内容变化；
- 修复后所有声明的 output 文件存在，manifest/checksum 与当前字节一致，文件集合无遗漏或额外文件；
- JSON 可解析，verdict、source provenance、输入 hash 和 claim ledger 语义未改变；
- 单独审阅、单独 `docs:` commit，不与后续目录迁移合并。

## 下一步

1. 生成不可变 pre-repair snapshot 和逐文件 mismatch report。
2. 只对可证明为 metadata stale 的条目重算 hash；无法解释的差异保持 `author_check`，不强行修复。
3. 运行 docs consistency、task ledger、JSON/hash audit 和 diff 检查后提交。
