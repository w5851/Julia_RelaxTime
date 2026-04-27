---
description: 开发文档全流程（设计/实现/台账/归档），管理 docs/dev/active 下的任务文档
---

# 开发文档管理

## 1. 从零设计任务单（doc-design-from-zero）

**输入要求：** 功能目标、范围边界、参考资料、时间/质量要求。

**输出格式（写入 `docs/dev/active/YYYY-MM-DD_描述.md`）：**
1. 背景与目标
2. 范围/非范围
3. 现状盘点
4. 任务分解（`- [ ]` 可勾选格式）
5. 测试与验收标准
6. 里程碑 / DoD
7. 风险与回退方案

**质量门禁：** 任务可执行、状态真实、范围清晰、可追踪。

## 2. 连续落地实现（doc-implementation）

**触发：** "按文档推进" / "执行下一步" / "继续推进"

1. 定位 `docs/dev/active/` 中匹配的活动文档。
2. 按任务顺序连续执行，默认不中途请示。
3. 每完成一项：主线自省（实现是否仍与文档/仓库结构一致）。
4. 每步验证后更新任务单状态。
5. 遇阻才停：文档与现状严重冲突、缺少关键依赖、需用户决策。
6. 全部完成后一次性汇总：完成项、验证结果、对齐总览、阻塞项、下一步建议。

## 3. 执行台账追加（doc-exec-log-append）

**Append-only 约束：**
- 只追加，不改写历史。
- 默认不回读台账正文。
- 每条记录自包含：目标、代码变更、验证命令、产物、结果、主线映射。

**脚本：** `julia --project=. scripts/dev/append_exec_log.jl`

## 4. 归档（doc-archive）

**脚本：** `julia --project=. scripts/dev/archive_docs.jl <filename.md>`
**校验：** `julia --project=. scripts/dev/archive_docs.jl -c`

归档后确认：active 中文件已移除、archived 中 frontmatter 含 4 个必需字段。
