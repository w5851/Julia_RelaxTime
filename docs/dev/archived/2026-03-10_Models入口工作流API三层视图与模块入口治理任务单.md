---
title: Models 入口工作流 API 三层视图与模块入口治理任务单
archived: true
original: docs/dev/active/2026-03-10_Models入口工作流API三层视图与模块入口治理任务单.md
archived_date: 2026-03-10
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Models 入口工作流 API 三层视图与模块入口治理任务单

更新日期：2026-03-10

> 目标：沿用已在 phase 主题验证过的 API 文档治理策略，对 `Models` 统一入口下的 transport workflow、meson workflow、模块访问器与参数适配层建立“面向用户入口 / 算法与职责核心 / 导出 API 全集”三层视图，并补齐当前导出索引中尚未被系统性吸收的公开入口。

---

## 1. 背景与问题定义

根据 [docs/api/generated/models/EntryPointsExportIndex.md](docs/api/generated/models/EntryPointsExportIndex.md)，`src/models/entrypoints.jl` 当前公开导出 15 个符号，其中相图主题已完成治理，但以下入口仍处于“有文档、但缺少统一主题治理”或“完全没有被人工文档吸收”的状态：

- `solve_gap_and_transport`
- `solve_transport_from_equilibrium`
- `solve_gap_and_meson_point`
- `transport_workflow_module`
- `meson_workflow_module`
- `workflow_param_adapters_module`
- `pnjl_module`

现状问题：

- transport 与 meson 各自已有页面，但仍分散在 `docs/api/relaxtime/workflow/` 与 `docs/api/pnjl/`，没有形成以 `Models` 统一入口为中心的主题导航；
- 模块访问器与参数适配层已被导出，但当前完全没有对应的人工说明；
- 使用者很难判断“应该从 `Models.*` 入口开始”，还是应该直接跳入底层 workflow 文件；
- 现有文档没有以导出索引为完整性基线，无法证明 `Models` 入口工作流主题已闭环。

因此，下一步仍应沿 API 治理方向推进，但目标不再是 phase，而是 `Models` 统一工作流入口主题。

---

## 2. 范围与非目标

### 2.1 范围

- `Models` 统一入口下的 transport workflow 文档治理；
- `Models` 统一入口下的 meson workflow 文档治理；
- 模块访问器函数的定位说明：
	- `transport_workflow_module`
	- `meson_workflow_module`
	- `workflow_param_adapters_module`
	- `pnjl_module`
- 对应导出 API 全集页的自动生成与覆盖检查；
- 旧路径页面与新主题目录之间的导航收口。

### 2.2 非目标

- 不在本任务中重写 transport 或 meson 底层实现；
- 不在本任务中统一重构整个 `docs/api/relaxtime/` 或 `docs/api/pnjl/` 目录；
- 不在本任务中引入 Documenter.jl；
- 不把所有 workflow 内部 helper 都提升为对外稳定入口。

---

## 3. 设计原则

- **仍以 `Models` 入口为第一视角**：用户入口页必须先讲 `Models.solve_gap_and_transport`、`Models.solve_transport_from_equilibrium`、`Models.solve_gap_and_meson_point`，而不是直接从底层 workflow 文件讲起。
- **模块访问器显式降级**：`*_workflow_module` 与 `workflow_param_adapters_module` 属于公开导出，但不是大多数用户的首选入口；文档应解释其职责边界，而不是把它们包装成常规业务入口。
- **导出完整性继续自动化**：完整性仍由导出索引页保证，不靠手工清单兜底。
- **保留原领域文档价值**：`docs/api/relaxtime/workflow/TransportWorkflow.md` 与 `docs/api/pnjl/MesonMassWorkflow.md` 可继续作为详细参考，但不应继续承担 `Models` 主题主入口职责。

---

## 4. 公开 API 基线

### 4.1 面向用户入口候选

- `Models.solve_gap_and_transport`
- `Models.solve_transport_from_equilibrium`
- `Models.solve_gap_and_meson_point`

### 4.2 模块访问器与适配层候选

- `Models.transport_workflow_module`
- `Models.meson_workflow_module`
- `Models.workflow_param_adapters_module`
- `Models.pnjl_module`

### 4.3 算法与职责核心候选

- transport workflow 的输入分层、参数优先级与输出合同
- meson workflow 的平衡态复用、通道扫描与 Mott 阈值合同
- workflow param adapters 的参数归一化职责边界
- module accessor 的 world-age / 聚合入口定位

说明：

- 用户入口层优先覆盖三个真正执行任务的 workflow 入口；
- 模块访问器与适配层进入“职责核心”或“进阶接口”层，而不作为首页主推荐；
- 公开导出完整性以 `src/models/entrypoints.jl` 与必要时 `src/models/Models.jl` 为准。

---

## 5. 建议目标结构

### 5.1 新结构

- `docs/api/models/workflows/README.md`
- `docs/api/models/workflows/Overview.md`
- `docs/api/models/workflows/CoreConcepts.md`
- `docs/api/models/workflows/generated/Exports.md`
- `docs/api/models/workflows/TransportWorkflow.md`
- `docs/api/models/workflows/MesonMassWorkflow.md`
- `docs/api/models/workflows/ModuleAccessors.md`

### 5.2 三层视图映射

| 视图 | 主页面 | 说明 |
| --- | --- | --- |
| 面向用户入口 | `workflows/Overview.md` | 统一说明 `Models` 入口下的 transport 与 meson 工作流怎么用 |
| 算法与职责核心 | `workflows/CoreConcepts.md` | 解释 workflow 输入分层、参数归一化、模块访问器职责边界 |
| 导出 API 全集 | `workflows/generated/Exports.md` | 自动生成，覆盖公开导出基线 |

### 5.3 旧路径保留策略

- `docs/api/relaxtime/workflow/TransportWorkflow.md` 保留为 transport 领域细节页；
- `docs/api/pnjl/MesonMassWorkflow.md` 保留为 meson 领域细节页；
- 新主题页负责从 `Models` 统一入口角度重新组织导航，而不是删除原页面。

---

## 6. 自动化策略

### 6.1 生成来源

- 首选来源：`src/models/entrypoints.jl`
- 必要时补充：`src/models/Models.jl`

### 6.2 初始过滤集合

- `solve_gap_and_transport`
- `solve_transport_from_equilibrium`
- `solve_gap_and_meson_point`
- `transport_workflow_module`
- `meson_workflow_module`
- `workflow_param_adapters_module`
- `pnjl_module`

### 6.3 生成目标

- 输出位置：`docs/api/models/workflows/generated/Exports.md`
- 仍使用 `scripts/dev/generate_api_export_index.jl`
- 页面必须显示：
	- 导出符号
	- 来源文件
	- 是否已被人工文档提及
	- 未被提及符号列表

---

## 7. 实施顺序

### 阶段 A：确定公开导出基线与主题边界

- [x] 确认 workflow 主题只覆盖 `Models` 统一入口，而不吞并整个 relaxtime/pnjl 领域文档
- [x] 固定导出过滤集合与生成命令
- [x] 确认模块访问器与适配层属于“进阶接口”，不放进首页首屏

阶段 A 落地说明（2026-03-10）：

- 主题边界已固定为 `Models` 统一入口层；
- transport 与 meson 原领域页继续保留，不被新主题吞并；
- 模块访问器与参数适配层已明确归入“进阶接口”层。

### 阶段 B：建立三层视图主结构

- [x] 新建 `docs/api/models/workflows/README.md`
- [x] 新建 `docs/api/models/workflows/Overview.md`
- [x] 新建 `docs/api/models/workflows/CoreConcepts.md`
- [x] 建立 `docs/api/models/workflows/generated/Exports.md`

阶段 B 落地说明（2026-03-10）：

- `README.md` 已建立为 workflow 主题主入口；
- `Overview.md` 已承接用户入口视图；
- `CoreConcepts.md` 已承接职责核心视图；
- 导出 API 全集页已由脚本生成，覆盖 7 个 workflow 主题导出。

### 阶段 C：迁入主题细节页

- [x] 新建 `docs/api/models/workflows/TransportWorkflow.md`
- [x] 新建 `docs/api/models/workflows/MesonMassWorkflow.md`
- [x] 新建 `docs/api/models/workflows/ModuleAccessors.md`
- [x] 在新页中明确链接回原领域细节页，而不是复制散落内容

阶段 C 落地说明（2026-03-10）：

- 已建立 transport、meson 与 module accessors 三个细节页；
- 新页已显式回链到 `docs/api/relaxtime/workflow/TransportWorkflow.md` 与 `docs/api/pnjl/MesonMassWorkflow.md`；
- 模块访问器与参数适配层已从“未被解释的导出”升级为有职责边界说明的进阶接口。

### 阶段 D：导航与覆盖校验

- [x] 更新 `docs/api/README.md`
- [x] 重新生成 `docs/api/generated/models/EntryPointsExportIndex.md`
- [x] 确认相关公开导出都能在人工文档或生成页中被找到

阶段 D 落地说明（2026-03-10）：

- `docs/api/README.md` 已纳入 `docs/api/models/workflows/` 主题入口；
- 已重新生成 `docs/api/generated/models/EntryPointsExportIndex.md`；
- workflow 主题的 7 个目标导出已在主题页或生成页中闭环，`meson_workflow_module`、`pnjl_module`、`transport_workflow_module`、`workflow_param_adapters_module` 已不再处于“未被提及”状态。

---

## 8. DoD

- [x] `docs/api/models/workflows/` 已形成三层视图
- [x] `Models` 入口下 transport 与 meson workflow 有统一主题入口
- [x] module accessor 与参数适配层已有职责说明
- [x] `docs/api/models/workflows/generated/Exports.md` 由脚本生成
- [x] `docs/api/generated/models/EntryPointsExportIndex.md` 中对应符号已不再处于“未被提及”状态
- [x] `docs/api/README.md` 已纳入新的 workflow 主题导航

---

## 9. 为什么下一步仍然是 API 方向

- phase 主题已经闭环，继续在同一主题内部追加收益很低；
- `Models` 统一入口仍然是仓库最重要的公开表面，当前恰好存在明确的导出覆盖缺口；
- 现有脚本与治理模板已经验证有效，把它扩到 workflow 主题的成本最低；
- 在转向新的代码治理任务之前，先把 `Models` 公开入口的文档面收拢，会显著提升整个仓库 API 的可发现性与一致性。
