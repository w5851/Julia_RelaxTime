# Issue #130 task-document archive v1

状态：`completed`（治理归档；不包含新的数值计算）。

本包记录 2026-08-30 在 `main@bae8634f87df856eccfe333642b6a13e1abb6327` 上完成的
Issue #130/RS 历史任务文档归档。归档操作使用
`scripts/dev/archive_docs.jl`，每个目标文件均保留原始正文，并在头部加入
`title`、`archived`、`original` 和 `archived_date` 元信息。

## 归档范围

本次归档包含 17 个已经完成或被后续合同明确 superseded 的 active 文档。它们不再是
当前执行入口，但仍是 PR、Actions、solver-free replay、作者审核和 reference/transport
决策的历史记录。逐文件清单见 `manifest.json`。

归档不改变以下内容：

- `docs/analysis/**` 下的不可变 evidence、生成时路径、manifest 和 hash；
- `data/reference/**`、`data/outputs/**`、workflow 定义和数值结果；
- Issue #130 的 accepted reference/runtime 语义及 RS 的 archived 结论。

历史 evidence 中如果仍出现 `docs/dev/active/...`，那是生成时 provenance，不能回写。
当前治理台账中的 `task_file` 和 `evidence` 指针已改为 `docs/dev/archived/...`。

## 明确排除

以下文档没有移动：

- `docs/dev/active/2026-07-17_PNJL全温区相图热积分与参考生产.md`：仍含更广泛的
  production/迁移工作，不能由本次 Issue #130 归档替代；
- `docs/dev/active/2026-08-19_docs-analysis逻辑分组与物理迁移任务单.md`：metadata
  repair follow-up 仍独立存在；
- `docs/dev/active/2026-08-15_论文级绘图SOP与分层迁移任务单.md`：SOP 已 promotion，
  但逐图族迁移和 strict layout gate 仍是独立后续任务。

## 验证

- `julia --project=. scripts/dev/archive_docs.jl --check`：全部 archived 文件格式通过；
- `julia --project=. scripts/dev/check_task_ledger.jl`：台账 schema、引用和状态通过；
- `julia --project=. scripts/dev/check_task_ledger.jl --preflight`：当前分支/dirty
  状态可追溯；
- `julia --project=. scripts/dev/check_docs_consistency.jl`；
- `julia --project=. scripts/dev/check_active_docs_governance.jl`；
- `git diff --check`。

## 后续边界

归档完成后才启动 `diagnostic-workflow-lifecycle` 的只读 inventory。workflow 的归档或
删除必须另立 PR，先提交保留期、artifact/provenance 引用图和精确 allowlist，并经作者
批准；本包不删除任何 workflow。
