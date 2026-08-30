---
title: Issue #130：RS old `prod_v1` path retirement implementation
archived: true
original: docs/dev/active/2026-08-28_Issue130_RS_old_reference_path_retirement.md
archived_date: 2026-08-30
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Issue #130：RS old `prod_v1` path retirement implementation

状态：accepted；PR #277 已通过 CI 并以 `main@74b53b47ebcca2b292cee72f70a70a84b0d2eea5`
合并。该 PR 只执行 canonical/path 退役；后续物理删除由独立任务和 PR 处理。

## 固定输入与目标

- audit 基线：PR #276，`main@ad47ad1482574faa4764af3913b8cfa3b6bae53e`
- implementation PR：#277，初始实现提交 `e5545bda`（当前 CI 待完成）
- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- source result commit：`05be2c05186f8e12baf3097b68f8619e53d19711`
- legacy case：`first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1`
- snapshot version：`legacy_prod_v1_snapshot_v1`
- current `prod_v2` canonical trees、数值内容和 runtime default 不变。

## 迁移合同

旧 case 的两个 mode result/figure tree 使用 `git mv` 原样迁入：

```text
data/outputs/results/relaxtime/transport/phase_guided/
  legacy_prod_v1_snapshot_v1/<mode>/<legacy_case>/
data/outputs/figures/relaxtime/transport/phase_guided/
  legacy_prod_v1_snapshot_v1/<mode>/<legacy_case>/
```

mode-specific canonical roots中不再出现旧 case。每棵 tree 的内部 manifest、CSV、PNG、
convergence sidecar 均保持原始字节；路径无关 tree SHA、manifest SHA、文件数和字节数写入
result snapshot 根目录的 `RETIREMENT_MANIFEST.json` 与 analysis 下的 figure snapshot manifest
记录上述信息；figure 输出根目录只保留图像和 `plot_manifest.json`，以满足 data-output-path
门禁。registry 的旧条目改为指向 snapshot，
并明确 `path_status=retired_to_versioned_legacy_snapshot`。

## Consumer 与回退边界

- 当前分析入口继续解析 `prod_v2`。
- 显式 legacy fallback/rollback 解析 versioned snapshot；不提供静默 canonical alias。
- 历史分析脚本可读取 snapshot；历史 evidence 包中的旧路径字符串保留为不可变 provenance，
  不回写或重生成。
- 物理删除不是本 PR 的目标，必须另立 deletion PR 并重新审核历史复现、fallback、rollback
  和 snapshot hash。

## 实现清单

- [x] 迁移 mode-A/mode-B 旧 result 与 figure tree，核对 canonical-root absence。
- [x] 新增 result snapshot README 与 retirement manifest；figure snapshot 的说明/manifest
  放在 `docs/analysis/relaxtime/issue130_rs_old_reference_path_retirement_v1/figure_snapshot/`，
  figure 输出目录仅保留允许的图像资产。
- [x] 更新 production registry、历史可执行分析/导入/审计工具和当前脚本指南的 snapshot 路径。
- [x] 增加 snapshot tree/hash、registry pointer、fallback/rollback 和旧 canonical absence 单元测试。
- [x] 运行 focused Python/Julia、consumer smoke、task-ledger、active-docs、docs consistency、
  data-output-path、script governance 与 `git diff --check`（本地 solver-free / governance checks）。
- [x] CI 全绿并经作者授权合并；PR #277 的 squash merge SHA 为
  `74b53b47ebcca2b292cee72f70a70a84b0d2eea5`。
- [x] path-retirement 分支及远端引用已清理；旧 snapshot 在本任务结束时仍由该 merge
  commit 保留，物理删除转交
  `2026-08-28_Issue130_RS_old_reference_physical_deletion.md`。

## 停止条件

若任一 snapshot tree hash 不匹配、旧 canonical 路径仍被可执行 consumer 使用、registry
指针无法解析、current/fallback 语义改变，立即停止；不重跑 solver、不改 production 数值、
不删除 snapshot。`prod_v1` snapshot 必须持续支持显式 rollback。
