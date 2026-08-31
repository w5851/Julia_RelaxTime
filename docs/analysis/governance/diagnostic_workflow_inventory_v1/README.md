# Diagnostic workflow inventory v1

这是 2026-08-30 在 `main@cea8110cdc1e1bafbc9873aa79cb5b0954c3de76` 上生成的
workflow 生命周期盘点包。PR 基于随后快进到的
`main@6b2f4ffeb05001f6e39b6940d59e4cbc9cf9bfed`，其间没有 workflow 文件变化；盘点快照仍以
前一个 SHA 为唯一来源。它承接 Issue #130/RS 历史文档归档，但不改变任何 workflow、run、
artifact、solver 或数值结果。

## 包内文件

- `workflow_inventory.csv`：仓库内 47 个 workflow 定义、触发方式、调用角色、活动引用和
  建议 disposition。
- `run_artifact_snapshot.csv`：20 个一次性诊断 workflow 最近观察 run 的 ID、head SHA、
  artifact 数量、过期数量和 artifact 名称摘要。
- `reference_edges.csv`：可复用 workflow 的调用边和排除 archived docs 后的活动引用扫描结果。
- `retirement_proposal.md`：作者审核用的分组方案和动作边界；不是删除指令。
- `claim_ledger.csv`：每个盘点结论的证据、允许用途和限制。
- `manifest.json`：快照来源、命令、schema、文件清单和 hash policy。

## 读取口径

1. 先读 `manifest.json`，确认 source SHA、查询日期和扫描边界。
2. 再读 `workflow_inventory.csv`，判断 workflow 定义是否仍属于 core/reusable/production
   或 historical candidate。
3. 用 `run_artifact_snapshot.csv` 判断证据是否仍可从 GitHub 下载；`artifact_expired` 是
   当前观察，不是证据内容正确性的判断。
4. 用 `reference_edges.csv` 区分 workflow 调用、代码/文档引用和“当前扫描未发现活动引用”。
5. 只有作者批准 `retirement_proposal.md` 后，才可另立 retirement PR。

## 扫描边界和限制

- workflow 文件来自 `.github/workflows/*.yml`，数量按 source tree 统计。
- 活动引用扫描排除 `docs/dev/archived/**`，因为其中的旧路径是生成时 provenance；它们仍由
  各自 archive manifest 保留。
- GitHub run/artifact 通过 `gh workflow list`、`gh run list` 和 Actions REST API 在快照日期
  查询；网络/API 权限或 retention 变化可能使未来查询不同。
- “无 standalone run”只表示未观察到独立运行，不表示 workflow_call 或脚本没有用途。
- 本包 `solver_called=false`，没有调用 PNJL equilibrium solver。

## 结论边界

本包证明的是“当前 workflow 定义和历史证据可以被分组审核”，不是“任何候选都可以立即删除”。
历史 run/artifact、脚本、测试和文档必须在独立 retirement PR 中逐项闭合后才能改变生命周期。
