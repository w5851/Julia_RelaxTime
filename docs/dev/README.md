# 开发文档总览

本目录用于记录开发过程中的计划、约定与变更说明。以下约定用于保持任务追踪与文档结构一致。

## active 与 archived 的区别

- active：进行中的任务与临时草案。
	- 只保留当前仍需推进的内容。
	- 任务完成后必须移出。
- archived：已完成任务的归档区。
	- 按统一格式记录“做了什么、验证了什么、对应了哪些文档”。
	- 保持可追溯，便于后续审计与复盘。

### 归档规则

完成的需求必须从 active 移出，放入 archived，并遵循以下格式：

#### 归档格式模板

归档文件需在头部使用统一的元信息格式，示例如下：

---
title: 任务名称
archived: true
original: docs/dev/active/原始任务文件.md
archived_date: 2026-01-19
---

归档内容规则：

- 头部元信息之后，直接原封不动复制开发文档内容，便于追溯。
- 当开发文档中包含多个需求时，应拆分为多个归档文档，并从 active 中移除对应需求点。

#### 归档示例

以下示例格式参考 [docs/dev/archived/2026-01-19_Dependency_Analysis_Recommendations.md](docs/dev/archived/2026-01-19_Dependency_Analysis_Recommendations.md)：

---
title: Dependency Analysis & Recommendations
archived: true
original: docs/dev/任务2.md
archived_date: 2026-01-19
---

以下为原始内容（保留，以便审阅与历史参考）：

***

（此处为原始开发文档内容原样粘贴）

## archived 文件命名规则

推荐命名格式：

- `[日期]_[描述]_[版本可选].扩展名`
- 日期格式：`YYYY-MM-DD`（ISO 8601），用于自然排序

示例：

- `2026-01-19_PNJL_Solver_Refactoring_Notes.md`
- `2026-01-19_Dependency_Analysis_Recommendations.md`

命名要点：

- 描述需明确，包含模块名/任务名/文档类型
- 可选追加版本号（如 `v1.2`）或状态（如 `deprecated`）
- 保留原扩展名

## 开发任务完成后的文档与测试补充

完成开发任务后，需同步补充以下内容：

- 必须：更新或新增 API 文档（docs/api/）
- 可选：更新或新增公式文档（docs/reference/formula/）
- 必须：补充单元测试（tests/unit/），或在归档说明中写明原因
- 可选：补充性能测试或分析（tests/perf/ 或 tests/analysis/）

## 项目结构约定（合并版）

### 顶层目录职责

- src/：核心可复用 Julia 源码（避免把一次性脚本放进来）。
- scripts/：可执行脚本、批处理与实验入口。
- tests/：测试与诊断。
	- tests/unit/：可自动化的单元测试。
	- tests/analysis/：分析/性能/调试脚本与报告。
- docs/：文档中心。
	- docs/api/：面向使用者的 API 文档。
	- docs/dev/：开发者文档。
	- docs/guides/：使用指南、排错与 Quickstart。
	- docs/reference/：公式、推导、理论背景材料。
- data/：数据与输出。
	- data/raw/：原始数据（只读/不做手改）。
	- data/processed/：处理后的数据。
	- data/outputs/：运行产物、结果与缓存输出。
- web/：前端静态资源与交互页面。
- config/：参数文件、扫描配置、模型参数集等。
- agent_rules/：代码风格、命名等规范（视作约束性文档）。

### 源码组织原则（src/）

- 以模块边界组织，而不是以功能碎片随意拆文件。
- 新增模块前先确认：是通用逻辑还是一次性实验。
	- 通用逻辑 → src/
	- 实验/临时对比 → scripts/ 或 tests/analysis/

### 单位与命名

- 单位约定以仓库根目录 README 为准。
- MeV 入参命名显式标注：T_MeV、mu_MeV、muB_MeV。
- fm⁻¹ 入参命名标注：T_fm、mu_fm。

### 文档与源码对齐策略

不强制 docs/api 与 src 逐文件 1:1 镜像：
- src 内部文件拆分会因重构而变化；逐文件镜像会带来路径震荡。
- API 文档应按稳定概念与入口组织。

推荐策略：
- docs/api/<domain>/ 与 src/<domain>/ 模块层级大致对齐。
- 每个稳定子模块/概念面提供一个页面。
- 当 src 重构但公共 API 不变时，仅更新文档链接即可。

### .github 规范

- .github/：协作相关文件统一管理。

### 允许/不允许提交的内容

- 建议：提交可复现脚本、小体积测试数据、关键结果摘要。
- 谨慎：大体积输出文件、重复生成的缓存、临时日志。
- 若必须提交大文件，请在 PR 说明来源、生成命令、用途与体积，并确认不会频繁变动。

### 变更流程建议

- 改代码 → 补测试（或说明为何不易测）→ 更新文档（docs/api 或 docs/dev）→ 在 PR 中写清验证方式。
