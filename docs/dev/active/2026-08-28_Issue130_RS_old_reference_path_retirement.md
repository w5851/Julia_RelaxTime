# Issue #130：RS old `prod_v1` path retirement implementation

状态：active；作者已审核并接受 retirement audit，当前 PR 只执行 canonical/path
退役，不物理删除数值文件，不调用 solver，不触发 numerical Actions。

## 固定输入与目标

- audit 基线：PR #276，`main@ad47ad1482574faa4764af3913b8cfa3b6bae53e`
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
两个 snapshot 根目录的 `RETIREMENT_MANIFEST.json`。registry 的旧条目改为指向 snapshot，
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
- [x] 新增 result/figure snapshot README 与 retirement manifest。
- [x] 更新 production registry、历史可执行分析/导入/审计工具和当前脚本指南的 snapshot 路径。
- [x] 增加 snapshot tree/hash、registry pointer、fallback/rollback 和旧 canonical absence 单元测试。
- [x] 运行 focused Python/Julia、consumer smoke、task-ledger、active-docs、docs consistency、
  data-output-path、script governance 与 `git diff --check`（本地 solver-free / governance checks）。
- [ ] CI 全绿后等待作者审核；本 PR 合并后才将 audit item 记录为 path-retirement 完成。

## 停止条件

若任一 snapshot tree hash 不匹配、旧 canonical 路径仍被可执行 consumer 使用、registry
指针无法解析、current/fallback 语义改变，立即停止；不重跑 solver、不改 production 数值、
不删除 snapshot。`prod_v1` snapshot 必须持续支持显式 rollback。
