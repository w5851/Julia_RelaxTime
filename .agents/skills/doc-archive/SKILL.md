---
name: doc-archive
description: Archive docs/dev/active markdown tasks into docs/dev/archived with required metadata frontmatter, using scripts/dev/archive_docs.jl and validating format compliance.
---

# Skill: doc-archive

## Purpose

将 `docs/dev/active` 中“已完成”的开发任务文档归档到 `docs/dev/archived`，并确保符合项目归档规范（YAML 元信息、可追溯来源、归档日期）。

适用仓库约定：
- 规则来源：`docs/dev/README.md`
- 执行脚本：`scripts/dev/archive_docs.jl`

## When to Use

在用户出现以下意图时使用本 skill：
- “归档 active 文档”
- “把任务从 active 移到 archived”
- “按开发文档规范归档”
- “检查归档格式是否正确”

关键词（用于识别）：
- archive / 归档
- docs/dev/active
- docs/dev/archived
- archived metadata / frontmatter
- archive_docs.jl

## When NOT to Use

- 用户要“编辑任务内容”而非归档时。
- 文档尚未完成、仍在 active 迭代时。
- 用户希望手工自定义非标准归档结构（不走规范）时。

## Inputs

- 一个或多个 active 文档文件名（如 `2026-02-16_两味NJL模型实现.md`）
- 可选归档日期（`YYYY-MM-DD`）

归档文件命名约定：
- 推荐格式：`YYYY-MM-DD_描述.md`
- 日期部分使用连字符（ISO 8601），日期与描述之间使用下划线。

## Outputs

- 已移动到 `docs/dev/archived` 的归档文档
- 文档头部包含：`title`、`archived`、`original`、`archived_date`
- （可选）格式校验结果

## Standard Procedure

1. **定位目标文件**
	 - 确认文件存在于 `docs/dev/active`。

2. **执行归档脚本**
	 - 单文件归档：
		 - `julia --project=. scripts/dev/archive_docs.jl <filename.md>`
	 - 指定日期归档：
		 - `julia --project=. scripts/dev/archive_docs.jl -d 2026-02-17 <filename.md>`
	 - 批量归档：
		 - `julia --project=. scripts/dev/archive_docs.jl file1.md file2.md`

3. **验证结果**
	 - 确认 active 中目标文件已移除。
	 - 确认 archived 中新文件存在。
	 - 检查 frontmatter 是否含 4 个必需字段。
	 - 可运行全量校验：
		 - `julia --project=. scripts/dev/archive_docs.jl -c`

4. **回报用户**
	 - 给出归档文件路径与校验结论。
	 - 若发现命名/格式偏差，明确说明并给出修正建议。

## Quality Checklist

- `archived: true` 已写入
- `original` 指向 `docs/dev/active/<filename>`
- `archived_date` 为有效日期字符串
- 原始内容保留在元信息后（可追溯）

## Notes

- 脚本默认会移动文件并写入归档头，不是仅复制。
- 若目标归档文件已存在，脚本可能询问是否覆盖。
- 如需预览不落盘，使用：
	- `julia --project=. scripts/dev/archive_docs.jl --dry-run <filename.md>`

## Example

- 归档：`docs/dev/active/2026-02-16_两味NJL模型实现.md`
	1) `julia --project=. scripts/dev/archive_docs.jl 2026-02-16_两味NJL模型实现.md`
	2) 检查输出文件：`docs/dev/archived/2026-02-16_两味NJL模型实现.md`
	3) 确认头部字段完整且 active 中文件已移除
