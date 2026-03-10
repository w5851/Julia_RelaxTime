---
title: Models 扫描入口 API 三层视图与采样契约治理任务单
archived: true
original: docs/dev/active/2026-03-10_Models扫描入口API三层视图与采样契约治理任务单.md
archived_date: 2026-03-10
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Models 扫描入口 API 三层视图与采样契约治理任务单

更新日期：2026-03-10

> 目标：沿用 phase 与 workflow 主题已验证的 API 文档治理策略，对 `Models` 统一入口下的扫描能力建立“面向用户入口 / 职责核心 / 导出 API 全集”三层视图，并把当前仍主要挂在旧 `pnjl` 路径下的扫描说明，重组为以 `Models.run_tmu_scan`、`Models.run_trho_scan`、`Models.build_default_rho_grid` 为中心的主题导航。

---

## 1. 背景与问题定义

根据 [docs/api/generated/models/EntryPointsExportIndex.md](docs/api/generated/models/EntryPointsExportIndex.md)，`src/models/entrypoints.jl` 当前公开导出 15 个符号，其中 phase 主题与 workflow 主题已完成治理，剩余仍未形成独立 `Models` 主题入口的公开扫描能力为：

- `run_tmu_scan`
- `run_trho_scan`
- `build_default_rho_grid`

现状问题：

- 这三个入口虽然已有文档提及，但主要仍分散在 `docs/api/pnjl/TmuScan.md`、`docs/api/pnjl/TrhoScan.md`、`docs/api/pnjl/ScanSamplingGuidelines.md`、`docs/api/pnjl/ScanOutputContract.md` 等旧路径；
- 现有页面主要从 `PNJL.*` 兼容层视角出发，而不是从 `Models` 统一入口出发；
- 使用者难以快速判断：扫描主入口、采样网格生成、输出契约、相图流程之间的职责边界分别是什么；
- `build_default_rho_grid` 作为公开导出，当前只有零散说明，缺少在 `Models` 主题中的定位；
- 尽管导出索引已显示这些符号“被提及”，但尚未形成像 phase/workflows 那样的主题级三层闭环。

因此，下一步仍应沿 API 治理方向推进，但目标从 workflow 转向 `Models` 扫描入口主题。

---

## 2. 范围与非目标

### 2.1 范围

- `Models.run_tmu_scan` 的用户入口与职责说明；
- `Models.run_trho_scan` 的用户入口与职责说明；
- `Models.build_default_rho_grid` 的采样职责与推荐使用方式；
- 扫描输出契约、采样模板、相图流程衔接关系的主题化整理；
- 对应导出 API 全集页的自动生成与覆盖检查；
- 旧 `docs/api/pnjl/` 扫描页面与新 `docs/api/models/scans/` 主题之间的导航收口。

### 2.2 非目标

- 不在本任务中重写 `TmuScan.jl` 或 `TrhoScan.jl` 底层实现；
- 不在本任务中重构 CSV 输出契约本身；
- 不在本任务中新增 Documenter.jl；
- 不在本任务中把 phase 主题与 scans 主题合并成一个超大页面；
- 不在本任务中处理 `PNJL.solve` 等平衡单点求解入口。

---

## 3. 设计原则

- **仍以 `Models` 入口为第一视角**：用户入口页优先解释 `Models.run_tmu_scan`、`Models.run_trho_scan`、`Models.build_default_rho_grid`，而不是继续把 `PNJL.*` 兼容层作为主入口。
- **第二层按职责核心层组织，但本主题命名为 `Algorithms.md`**：扫描主题的第二层重点是采样策略、扫描顺序、连续性/phase-aware 策略、输出契约衔接，算法与流程判据占比更高，因此页面名采用 `Algorithms.md`；这与仓库统一规则“第二层语义为职责核心层”并不冲突。
- **保留旧领域页价值**：`docs/api/pnjl/TmuScan.md`、`docs/api/pnjl/TrhoScan.md`、`docs/api/pnjl/ScanSamplingGuidelines.md`、`docs/api/pnjl/ScanOutputContract.md` 继续保留为细节参考，但不再承担 `Models` 主主题入口职责。
- **导出完整性继续自动化**：完整性仍由 `scripts/dev/generate_api_export_index.jl` 生成页保证，而不是靠人工清单兜底。
- **把采样网格函数从“被顺带提及”升级为“有定位的公开接口”**：`build_default_rho_grid` 不应只埋在采样说明页里，而应成为 scan 主题的明确组成部分。

---

## 4. 公开 API 基线

### 4.1 面向用户入口候选

- `Models.run_tmu_scan`
- `Models.run_trho_scan`

### 4.2 进阶但公开的辅助入口

- `Models.build_default_rho_grid`

### 4.3 职责核心候选

- T-μ 与 T-ρ 扫描的输入口径与输出合同；
- `reverse_rho`、`phase_aware`、连续性种子/phase-aware 策略的职责边界；
- `build_default_rho_grid` 与 phase 产线、采样模板之间的衔接关系；
- 扫描输出契约与下游分析/作图/相图流程的接口边界。

说明：

- 本主题的第二层虽然命名为 `Algorithms.md`，但语义仍是“职责核心层”；
- `build_default_rho_grid` 属于公开导出，但不是首页首屏主推荐入口；
- 公开导出完整性以 `src/models/entrypoints.jl` 为主基线。

---

## 5. 建议目标结构

### 5.1 新结构

- `docs/api/models/scans/README.md`
- `docs/api/models/scans/Overview.md`
- `docs/api/models/scans/Algorithms.md`
- `docs/api/models/scans/generated/Exports.md`
- `docs/api/models/scans/TmuScan.md`
- `docs/api/models/scans/TrhoScan.md`
- `docs/api/models/scans/SamplingGrid.md`

### 5.2 三层视图映射

| 视图 | 主页面 | 说明 |
| --- | --- | --- |
| 面向用户入口 | `scans/Overview.md` | 统一说明 `Models` 入口下的扫描能力怎么选、怎么用 |
| 职责核心 | `scans/Algorithms.md` | 解释采样策略、扫描顺序、phase-aware/continuity、输出契约衔接 |
| 导出 API 全集 | `scans/generated/Exports.md` | 自动生成，覆盖公开导出基线 |

### 5.3 旧路径保留策略

- `docs/api/pnjl/TmuScan.md` 保留为 T-μ 扫描细节页；
- `docs/api/pnjl/TrhoScan.md` 保留为 T-ρ 扫描细节页；
- `docs/api/pnjl/ScanSamplingGuidelines.md` 与 `docs/api/pnjl/ScanOutputContract.md` 保留为支撑文档；
- 新主题页负责从 `Models` 统一入口角度重新组织导航，而不是立即删除旧页面。

---

## 6. 自动化策略

### 6.1 生成来源

- 首选来源：`src/models/entrypoints.jl`

### 6.2 初始过滤集合

- `run_tmu_scan`
- `run_trho_scan`
- `build_default_rho_grid`

### 6.3 生成目标

- 输出位置：`docs/api/models/scans/generated/Exports.md`
- 仍使用 `scripts/dev/generate_api_export_index.jl`
- 页面必须显示：
	- 导出符号
	- 来源文件
	- 是否已被人工文档提及
	- 未被提及符号列表

---

## 7. 实施顺序

### 阶段 A：确定扫描主题边界与公开导出基线

- [x] 确认 scan 主题只覆盖 `Models` 统一扫描入口，而不吞并整个 `pnjl` 扫描细节目录
- [x] 固定导出过滤集合与生成命令
- [x] 明确 `build_default_rho_grid` 属于公开辅助入口，而不是首页首屏主入口

阶段 A 落地说明（2026-03-10）：

- scan 主题边界已固定为 `Models` 统一扫描入口层；
- `docs/api/pnjl/` 下的扫描细节与契约页继续保留，不被新主题吞并；
- `run_tmu_scan`、`run_trho_scan`、`build_default_rho_grid` 已确定为主题导出基线，其中 `build_default_rho_grid` 明确归入公开辅助入口。

### 阶段 B：建立三层视图主结构

- [x] 新建 `docs/api/models/scans/README.md`
- [x] 新建 `docs/api/models/scans/Overview.md`
- [x] 新建 `docs/api/models/scans/Algorithms.md`
- [x] 建立 `docs/api/models/scans/generated/Exports.md`

阶段 B 落地说明（2026-03-10）：

- `README.md`、`Overview.md`、`Algorithms.md` 已建立；
- 第二层页面命名采用 `Algorithms.md`，但语义上仍对应统一的职责核心层；
- 导出全集页已由脚本生成，覆盖 `run_tmu_scan`、`run_trho_scan`、`build_default_rho_grid` 三个主题导出。

### 阶段 C：迁入主题细节页

- [x] 新建 `docs/api/models/scans/TmuScan.md`
- [x] 新建 `docs/api/models/scans/TrhoScan.md`
- [x] 新建 `docs/api/models/scans/SamplingGrid.md`
- [x] 在新页中明确链接回 `docs/api/pnjl/` 下的扫描细节与契约页，而不是复制散落内容

阶段 C 落地说明（2026-03-10）：

- 已建立 T-μ、T-ρ 与采样网格三个主题细节页；
- 新页均从 `Models` 统一入口角度重述定位，而不是沿用 `PNJL.*` 子模块视角；
- 旧路径细节页被保留为参考链接，避免重复维护同一批长篇实现说明。

### 阶段 D：导航与覆盖校验

- [x] 更新 `docs/api/README.md`
- [x] 重新生成 `docs/api/generated/models/EntryPointsExportIndex.md`
- [x] 确认 `run_tmu_scan`、`run_trho_scan`、`build_default_rho_grid` 在人工文档或生成页中形成主题级闭环

阶段 D 落地说明（2026-03-10）：

- `docs/api/README.md` 已纳入 `docs/api/models/scans/` 主题入口；
- 已重新生成 `docs/api/generated/models/EntryPointsExportIndex.md`；
- 总入口导出索引中，`run_tmu_scan`、`run_trho_scan`、`build_default_rho_grid` 已明确命中新建 scan 主题页面，不再只是旧 `pnjl` 路径下的零散提及。

---

## 8. DoD

- [x] `docs/api/models/scans/` 已形成三层视图
- [x] `Models` 入口下扫描能力已有统一主题入口
- [x] `build_default_rho_grid` 已有清晰职责说明
- [x] `docs/api/models/scans/generated/Exports.md` 由脚本生成
- [x] `docs/api/generated/models/EntryPointsExportIndex.md` 中对应 scan 符号已纳入新的主题导航语境
- [x] `docs/api/README.md` 已纳入新的 scans 主题导航

---

## 9. 为什么下一个主题选 scans

- phase 与 workflows 已形成 `Models` 主题主入口，扫描入口正好是下一个仍显著暴露在旧 `pnjl` 路径语境中的公开表面；
- `run_tmu_scan`、`run_trho_scan` 与 `build_default_rho_grid` 共同构成了相图产线与研究扫描的前置入口，价值足够集中；
- 现有支撑文档已经很多，但缺的是 `Models` 统一主题导航，而不是缺内容本身；
- 继续沿这一方向推进，能把 `src/models/entrypoints.jl` 的公开表面进一步整理成更清晰的主题集合。