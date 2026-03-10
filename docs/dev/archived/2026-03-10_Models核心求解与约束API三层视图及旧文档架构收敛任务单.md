---
title: Models 核心求解与约束 API 三层视图及旧文档架构收敛任务单
archived: true
original: docs/dev/active/2026-03-10_Models核心求解与约束API三层视图及旧文档架构收敛任务单.md
archived_date: 2026-03-10
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Models 核心求解与约束 API 三层视图及旧文档架构收敛任务单

更新日期：2026-03-10

> 目标：对 `Models` 统一入口下“核心求解、状态表示、约束模式、初值策略、隐式求解器”这一批高价值公开表面建立三层视图主题，并同步收敛旧 API 文档架构：把旧 `docs/api/pnjl/*.md` 中仍有价值的内容直接吸收入新主题，而不是继续让新主题依赖旧路径页面承载主体说明。

---

## 1. 背景与问题定义

根据 [docs/api/generated/models/ModelsExportIndex.md](docs/api/generated/models/ModelsExportIndex.md)，`src/models/Models.jl` 当前公开导出 144 个符号，其中 phase、workflows、scans 三个流程主题已建立新架构主题页，但以下高价值公共表面仍未形成独立 `Models` 主题：

- 模型创建与最小稳定入口：`create_model`
- 状态与合同：`MeanFieldState`、`meanfield_state`、`state_vector`、`normalize_mu_vec`
- 平衡求解入口：`solve_gap`、`solve`、`solve_multi`
- 约束模式求解：`solve_fixedmu_constraint`、`solve_fixedrho_constraint`、`solve_fixedentropy_constraint`、`solve_fixedsigma_constraint`、`solve_fixedasymrho_constraint`、`solve_constraint`
- 约束/残差构建：`ConstraintModes`、`GapParams`、`build_conditions`、`build_residual!`、`gap_state_dim`、`gap_residual`
- 初值策略：`SeedStrategy`、`DefaultSeed`、`MultiSeed`、`ContinuitySeed`、`HybridContinuitySeed`、`PhaseAwareSeed`、`PhaseAwareContinuitySeed`、`get_seed`、`get_all_seeds`、`update!`、`reset!`、`set_phase!`
- 求解器工厂：`create_implicit_gap_solver`、`create_flavor_mu_implicit_gap_solver`、`create_pnjl_implicit_solver`、`create_implicit_solver`、`solve_with_derivatives`

现状问题：

- 这些能力虽然部分已在 [docs/api/data_contracts.md](docs/api/data_contracts.md) 与 [docs/api/PARAMETER_TYPES_API.md](docs/api/PARAMETER_TYPES_API.md) 中被零散提及，但还没有形成以 `Models` 统一入口为中心的主题导航；
- 旧 `docs/api/pnjl/Solver.md`、`ImplicitSolver.md`、`ConstraintModes.md`、`Conditions.md`、`SeedStrategies.md`、`FluctuationADPath.md` 仍承载大量高价值说明；
- 当前新主题治理若继续采用“旧路径更细说明仍可参考”的方式，会把旧 API 文档架构长期保留下来，违背新架构主题收敛目标；
- `ModelsExportIndex.md` 中这一簇导出符号有大量仍是 0 命中或只有 1 次零散命中，说明“可发现性”仍不足。

因此，下一步最合理的方向不是继续做另一个流程主题，而是为 `Models` 的核心求解与约束公共表面建立主主题，并同步做旧文档架构吸收与降级。

---

## 2. 范围与非目标

### 2.1 范围

- `Models` 最小稳定入口与模型创建的主题化说明；
- 平衡求解、状态合同、约束模式、初值策略、隐式求解器工厂的三层视图整理；
- 将旧 `docs/api/pnjl/` 求解器相关页面中的高价值说明吸收入新主题；
- 新旧页面之间的导航收口与旧页降级策略；
- 对应导出 API 全集页的自动生成与覆盖检查。

### 2.2 非目标

- 不在本任务中重写 `src/models/solver/*` 或 `src/models/gap_solver.jl` 的实现；
- 不在本任务中处理 magnetic 主题的独立治理；
- 不在本任务中处理 conserved-charge / cumulant 全家桶的主题化；
- 不在本任务中立即删除所有旧 `docs/api/pnjl/*` 页面，但要明确哪些页在吸收完成后应降级为迁移页。

---

## 3. 设计原则

- **仍以 `Models` 入口为第一视角**：用户入口页必须优先解释 `Models.create_model`、`Models.solve_gap`、`Models.solve`、`Models.solve_multi` 以及状态/约束合同，而不是继续让 `PNJL.*` 顶层页承担主导航。
- **本主题的第二层采用 `CoreConcepts.md`**：因为这批公共表面不仅是算法判据，更是状态合同、约束模式、工厂函数和策略对象之间的职责边界。
- **新主题直接吸收旧页价值**：旧 `Solver.md`、`ImplicitSolver.md`、`ConstraintModes.md`、`Conditions.md`、`SeedStrategies.md`、`FluctuationADPath.md` 中仍有价值的说明应进入新主题页，而不是只在新页里放“详见旧页”的链接。
- **旧页完成吸收后降级**：吸收完成后，旧 `pnjl` 路径页面应尽量退化为迁移说明、兼容层说明或跳转页，避免新旧两套主架构长期并列。
- **导出完整性继续自动化**：完整性仍由导出索引页保证，不靠人工清单兜底。

---

## 4. 公开 API 基线

### 4.1 面向用户入口候选

- `Models.create_model`
- `Models.solve_gap`
- `Models.solve`
- `Models.solve_multi`

### 4.2 职责核心候选

- `MeanFieldState` / `meanfield_state` / `state_vector` / `normalize_mu_vec`
- `ConstraintModes`
- `GapParams` / `build_conditions` / `build_residual!` / `gap_state_dim` / `gap_residual`
- `SeedStrategy` 与各类 seed policy
- `create_implicit_gap_solver` / `create_flavor_mu_implicit_gap_solver` / `create_pnjl_implicit_solver` / `create_implicit_solver`

### 4.3 进阶公开接口候选

- `solve_fixedmu_constraint`
- `solve_fixedrho_constraint`
- `solve_fixedentropy_constraint`
- `solve_fixedsigma_constraint`
- `solve_fixedasymrho_constraint`
- `solve_constraint`
- `solve_with_derivatives`

说明：

- 本主题优先处理“最小稳定求解入口”和其直接支撑合同；
- conserved charge、magnetic、transport provider 等后续可再按独立主题治理；
- 若某些 seed 常量或 low-level 工厂最终不适合作为首页主入口，也必须在新主题中给出边界说明，而不是遗漏。

---

## 5. 建议目标结构

### 5.1 新结构

- `docs/api/models/solver/README.md`
- `docs/api/models/solver/Overview.md`
- `docs/api/models/solver/CoreConcepts.md`
- `docs/api/models/solver/generated/Exports.md`
- `docs/api/models/solver/StateContract.md`
- `docs/api/models/solver/ConstraintModes.md`
- `docs/api/models/solver/SeedStrategies.md`
- `docs/api/models/solver/ImplicitSolvers.md`

### 5.2 三层视图映射

| 视图 | 主页面 | 说明 |
| --- | --- | --- |
| 面向用户入口 | `solver/Overview.md` | 统一说明 `create_model`、`solve_gap`、`solve`、`solve_multi` 怎么选 |
| 职责核心 | `solver/CoreConcepts.md` | 解释状态合同、约束模式、seed policy、solver factory 边界 |
| 导出 API 全集 | `solver/generated/Exports.md` | 自动生成，覆盖 solver 主题导出基线 |

### 5.3 旧路径吸收与降级目标

待吸收的旧页面：

- `docs/api/pnjl/Solver.md`
- `docs/api/pnjl/ImplicitSolver.md`
- `docs/api/pnjl/ConstraintModes.md`
- `docs/api/pnjl/Conditions.md`
- `docs/api/pnjl/SeedStrategies.md`
- `docs/api/pnjl/FluctuationADPath.md`

吸收完成后的目标：

- 新 `docs/api/models/solver/*` 页面可独立承担主说明职责；
- 旧页面退化为迁移说明、兼容层说明或跳转页，不再作为新主题完成度前提。

---

## 6. 自动化策略

### 6.1 生成来源

- 首选来源：`src/models/Models.jl`

### 6.2 初始过滤集合

- `create_model`
- `MeanFieldState`
- `meanfield_state`
- `state_vector`
- `normalize_mu_vec`
- `solve_gap`
- `solve`
- `solve_multi`
- `ConstraintModes`
- `GapParams`
- `build_conditions`
- `build_residual!`
- `gap_state_dim`
- `gap_residual`
- `SeedStrategy`
- `DefaultSeed`
- `MultiSeed`
- `ContinuitySeed`
- `HybridContinuitySeed`
- `PhaseAwareSeed`
- `PhaseAwareContinuitySeed`
- `get_seed`
- `get_all_seeds`
- `update!`
- `reset!`
- `set_phase!`
- `create_implicit_gap_solver`
- `create_flavor_mu_implicit_gap_solver`
- `create_pnjl_implicit_solver`
- `create_implicit_solver`
- `solve_with_derivatives`
- `solve_fixedmu_constraint`
- `solve_fixedrho_constraint`
- `solve_fixedentropy_constraint`
- `solve_fixedsigma_constraint`
- `solve_fixedasymrho_constraint`
- `solve_constraint`

### 6.3 生成目标

- 输出位置：`docs/api/models/solver/generated/Exports.md`
- 仍使用 `scripts/dev/generate_api_export_index.jl`
- 页面必须显示：
	- 导出符号
	- 来源文件
	- 是否已被人工文档提及
	- 未被提及符号列表

---

## 7. 实施顺序

### 阶段 A：确定 solver 主题边界与导出基线

- [x] 确认本主题聚焦 `Models` 核心求解与约束公共表面，不吞并 magnetic / cumulant / transport provider 主题
- [x] 固定导出过滤集合与生成命令
- [x] 划分用户入口、职责核心与进阶接口层次

阶段 A 落地说明（2026-03-10）：

- solver 主题边界已固定为模型创建、状态合同、约束模式、seed strategy 与 implicit solver factory；
- magnetic、cumulant、transport provider 等能力仍保留为后续独立主题候选；
- 用户入口、职责核心与进阶接口层次已按任务单固定。

### 阶段 B：建立三层视图主结构

- [x] 新建 `docs/api/models/solver/README.md`
- [x] 新建 `docs/api/models/solver/Overview.md`
- [x] 新建 `docs/api/models/solver/CoreConcepts.md`
- [x] 建立 `docs/api/models/solver/generated/Exports.md`

阶段 B 落地说明（2026-03-10）：

- `README.md`、`Overview.md`、`CoreConcepts.md` 已建立；
- 第二层采用 `CoreConcepts.md`，用来承载合同、策略和工厂函数之间的职责边界；
- 导出全集页已由脚本生成，覆盖 solver 主题导出基线。

### 阶段 C：吸收旧页价值并建立细节页

- [x] 新建 `docs/api/models/solver/StateContract.md`
- [x] 新建 `docs/api/models/solver/ConstraintModes.md`
- [x] 新建 `docs/api/models/solver/SeedStrategies.md`
- [x] 新建 `docs/api/models/solver/ImplicitSolvers.md`
- [x] 将旧 `pnjl` 求解器相关页中的高价值说明直接吸收到新页

阶段 C 落地说明（2026-03-10）：

- 已建立状态合同、约束模式、seed strategy、implicit solver 四个细节页；
- 新页已直接吸收旧 `pnjl` 求解器相关页中的关键说明，而不是继续采用“详见旧页”的主说明方式；
- 旧页后续可进入迁移说明/跳转页定位。

### 阶段 D：旧文档架构收敛

- [x] 标出哪些旧 `pnjl` 页面在吸收完成后应降级为迁移说明/跳转页
- [x] 避免新 `models/solver` 页面继续以“详见旧页”作为主体补充方式
- [x] 必要时把旧页改写为薄跳转页而非主说明页

阶段 D 落地说明（2026-03-10）：

- `Solver.md`、`ImplicitSolver.md`、`ConstraintModes.md`、`Conditions.md`、`SeedStrategies.md`、`FluctuationADPath.md` 已降级为迁移说明页；
- 新 `models/solver` 页面已直接吸收旧页价值，不再以旧页作为主体说明依赖；
- 旧文档架构在 solver 主题上已从“并列主说明”降级为“兼容层路径说明”。

### 阶段 E：导航与覆盖校验

- [x] 更新 `docs/api/README.md`
- [x] 重新生成 `docs/api/generated/models/ModelsExportIndex.md`
- [x] 确认 solver 主题导出在新主题页中形成闭环

阶段 E 落地说明（2026-03-10）：

- `docs/api/README.md` 已纳入 `docs/api/models/solver/` 主题入口；
- 已重新生成 `docs/api/generated/models/ModelsExportIndex.md`；
- solver 主题导出在 `docs/api/models/solver/generated/Exports.md` 中已完成闭环，全部导出符号均已在至少一个非 generated API 页面中命中；
- `docs/api/generated/models/ModelsExportIndex.md` 已同步反映 solver 新主题的命中结果，`create_implicit_solver`、`solve_constraint`、`solve_with_derivatives` 等核心入口均已收敛到新路径。

---

## 8. DoD

- [x] `docs/api/models/solver/` 已形成三层视图
- [x] `Models` 的核心求解与约束已有统一主题入口
- [x] 旧 `pnjl` 求解器相关页中的高价值内容已被新主题吸收
- [x] 新 `models/solver` 页面不再以旧路径文档作为主体说明依赖
- [x] `docs/api/models/solver/generated/Exports.md` 由脚本生成
- [x] `docs/api/generated/models/ModelsExportIndex.md` 中对应 solver 符号已明显提升新主题命中
- [x] 旧 `pnjl` 求解器相关页已标出迁移后定位

---

## 9. 为什么下一步选 solver 主题

- 它是 `Models` 公开表面里尚未主题化的最大核心能力簇；
- 相关导出数量多、层次清晰，且有大量 0 命中或仅旧路径命中；
- 旧 `pnjl` 文档中已有丰富内容，正适合按照你刚澄清的目标做“新主题吸收旧文档价值”；
- 完成这一主题后，`Models` 的主入口将从“流程页为主”扩展到“核心求解合同也有主主题”，整体 API 架构会明显更完整。