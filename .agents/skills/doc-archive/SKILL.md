---
name: doc-archive
description: 使用 scripts/dev/archive_docs.jl，将已完成、明确取消或被替代且已有终止依据的活动任务文档归档，保留正文并核对元数据与引用。
---

# 文档归档

## 契约

- 先核对任务已完成，或已明确取消/被替代；未完成的剩余工作不能通过归档消失。
- 使用 `scripts/dev/archive_docs.jl`，不要手动重写其移动和 frontmatter 处理逻辑。
- 保留原始内容，并记录 `title`、`archived`、`original` 和 `archived_date`。
- 遵循 `docs/dev/README.md` 中的命名和元数据规则。

## 工作流程

1. 核对终态与证据、文件范围以及仍有效的引用和 ledger 路径；取消或被替代必须记录原因。
2. 针对明确指定的活动文档文件名执行预演：

```powershell
julia --project=. scripts/dev/archive_docs.jl --dry-run <filename.md>
```

3. 用同一组参数执行。默认 `--status completed`；取消/被替代使用 `--status cancelled|superseded --reason <原因>`。可选 `--date YYYY-MM-DD`。脚本在预演和执行中均拒绝已有目标及范围外文件。
4. 确认源文件已移出 `docs/dev/active/`，目标文件存在于 `docs/dev/archived/`，并且所有必需的 frontmatter 字段齐全。
5. 优先检查本批目标文件，并更新仍有效的引用与 ledger 文件路径；全目录检查留给归档审计：

```powershell
julia --project=. scripts/dev/archive_docs.jl --check <archived-filename.md>
```

6. 汇报目标路径和验证结果。

## 边界

- 脚本不覆盖已有归档；遇到冲突先核对任务与目标，不自动改名绕过。
- 不为通过活动任务治理检查而归档只完成一部分的文档。
- 除非任务需要追溯来源，否则不要从归档历史中寻找当前事实。
