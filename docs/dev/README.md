# 开发文档总览

本目录用于记录开发过程中的计划、约定与变更说明。以下约定用于保持任务追踪与文档结构一致。

## 常用内部手册

- [Codex高阶使用手册](D:\Desktop\Julia_RelaxTime\docs\dev\Codex高阶使用手册.md)：面向本仓库的 Codex 协作约束、skill 设计与提示模板。

## active / backlog / archived 的区别

- active：进行中的任务与临时草案。
	- 只保留当前仍需推进的内容。
	- 应当是近期会继续执行的任务单，而不是长期路线图。
	- 任务完成后必须移出。
- backlog：长期路线图、候选需求、分阶段能力盘点。
	- 适用于“尚未进入当前执行批次、但仍需保留并滚动维护”的文档。
	- 不作为 active 治理检查对象。
	- 当 backlog 条目被正式拉起执行时，应拆成新的 active 任务单。
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

## active 命名与归档触发规则（新增）

- `docs/dev/active` 下任务文档命名统一为：`YYYY-MM-DD_描述.md`
- `docs/dev/backlog` 下文档建议沿用同一命名格式：`YYYY-MM-DD_描述.md`
- 归档触发条件（满足任一）：
	- 任务单 DoD 全部勾选完成；
	- 文档创建后超过 60 天仍停留在 active；
	- 任务被新任务单替代，旧文档仅保留历史价值。
- 长期路线图、候选能力盘点、未进入当前执行批次的规划文档，不应继续放在 `docs/dev/active`；应迁移到 `docs/dev/backlog`。
- 建议校验命令：

```powershell
julia --project=. scripts/dev/check_active_docs_governance.jl
julia --project=. scripts/dev/check_pnjl_migration_guard.jl
```

- 归档命令：

```powershell
julia --project=. scripts/dev/archive_docs.jl <filename.md>
```

## 开发任务完成后的文档与测试补充

完成开发任务后，需同步补充以下内容：

- 必须：更新或新增 API 文档（docs/api/）
- 可选：更新或新增公式文档（docs/reference/formula/）
- 必须：补充单元测试（tests/unit/），或在归档说明中写明原因
- 可选：补充性能测试或分析（benchmark/ 或 scripts/perf/ 或 scripts/analysis/）

## 项目结构约定（合并版）

### 顶层目录职责

- src/：核心可复用 Julia 源码（避免把一次性脚本放进来）。
- scripts/：可执行脚本、批处理与实验入口。
- tests/：测试与诊断。
	- tests/unit/：可自动化的单元测试。
	- tests/integration/：跨模块集成测试。
	- tests/regression/：数值回归测试。
	- tests/validation/：外部参考验证测试。
- scripts/analysis/：分析/诊断脚本与报告（非测试入口）。
- scripts/perf/：性能探针与 profiling 脚本（非测试入口）。
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
- `.github/copilot-instructions.md`：工作区级 Copilot 指令与项目级编码约束。
- `docs/dev/testing_governance.md`：测试组织、入口与 smoke/full 治理规范。
- `docs/architecture/models_solver_contract.md`：`Models` 统一求解接口契约（Problem-like + solve）。

### 源码组织原则（src/）

- 以模块边界组织，而不是以功能碎片随意拆文件。
- 新增模块前先确认：是通用逻辑还是一次性实验。
	- 通用逻辑 → src/
	- 实验/临时对比 → scripts/ 或 scripts/analysis/

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

## world-age 边界治理（迁移期）

- `src/models` 中 `Base.invokelatest(...)` 仅允许出现在已登记边界点，禁止新增散点调用。
- 当前登记基线为 7 处，详见 `docs/architecture/dependency_rules.md` 的 “world-age 动态调用边界” 小节。
- 边界登记的机读单一来源为 `config/ci/models_invokelatest_allowlist.toml`。
- 若确需新增，必须同时完成：
	- 更新 `scripts/dev/check_pnjl_migration_guard.jl` 白名单；
	- 在迁移任务单补充必要性与回归证据；
	- 提供 `models-invokelatest-audit` 输出（observed/baseline/allowlisted）。

## 运行时入口迁移现状（2026-02-24）

- `src/simulation/fullserver` 的扫描与单点求解主路径已切到 `Models` 统一入口。
- `src/pnjl/PNJL.jl` 当前定位为兼容层，保留对外导出以支持迁移过渡。
- 新增调用方默认应优先使用 `src/models/entrypoints.jl` 暴露的入口；仅在兼容需求下使用 `PNJL` 导出。
