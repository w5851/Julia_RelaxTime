---
title: Parameter Structs Phase C 任务单（可勾选）
archived: true
original: docs/dev/active/2026_02_18_parameter_structs_phaseC_plan.md
archived_date: 2026-02-18
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Parameter Structs Phase C 任务单（可勾选）

## 文档目的

在 Phase B（兼容层与入口统一）完成后，进入 Phase C：

- 收口参数契约（减少隐式兼容分支）
- 明确对外 API 的稳定输入形态
- 为后续移除历史接口做“可验证、可回滚”的准备

本阶段继续遵循策略：**对外渐进兼容、对内强约束统一**。

---

## 当前状态快照（2026-02-18）

- [x] Phase B 文档已归档（`docs/dev/archived/2026_02_17_parameter_structs_phaseB_plan.md`）
- [x] Workflow 参数适配层已集中到 `WorkflowParamAdapters`
- [x] Scan 已支持 config object（`TmuScanConfig` / `TrhoScanConfig`）
- [x] 关键等价性与边界诊断 smoke tests 已补齐
- [x] NamedTuple 路径已收口为兼容边界并带 depwarn（主路径为结构体）

---

## Phase C 目标与范围

### 目标

- 统一 PNJL 侧公共入口的参数输入规则（主路径以结构体为准）
- 保留最小必要兼容层，并给出明确弃用路线
- 建立“参数契约变更不会静默破坏”的测试与文档基线

### 本期范围（In Scope）

- `src/pnjl/workflows/*`
- `src/pnjl/scans/*`
- `docs/api/PARAMETER_TYPES_API.md`
- `docs/guides/PARAMETER_STRUCT_MIGRATION.md`
- `tests/unit/relaxtime/*` 与 `tests/unit/pnjl/*` 中参数契约相关测试

### 非目标（Out of Scope）

- 重写物理模型公式或求解器核心算法
- 大规模改动调用层业务脚本（仅做必要兼容）

---

## 任务清单（可勾选）

### C1 参数契约收口（入口签名与验证）

- [x] 梳理 workflows/scans 对外入口函数签名，形成“唯一推荐签名”列表
- [x] 在入口增加统一参数校验（类型、字段完整性、关键边界）
- [x] 对不满足契约的输入统一抛出 `ArgumentError`（信息包含字段名与期望类型）
- [x] 清理重复/分散的 `NamedTuple -> struct` 转换逻辑（保留单一边界）

推荐入口签名（Phase C 当前结论）：

- scans：优先使用 `run_tmu_scan(config::TmuScanConfig)` / `run_trho_scan(config::TrhoScanConfig)`；kwargs 入口保留兼容。
- workflows（transport）：`solve_gap_and_transport(T_fm::Real, mu_fm::Real; ...)` 与 `solve_transport_from_equilibrium(base, T_fm::Real, mu_fm::Real; ...)`。
- workflows（meson）：`solve_gap_and_meson_point(T_fm::Real, mu_fm::Real; ...)`。

### C2 弃用路径与兼容策略（软迁移）

- [x] 为历史调用形态增加 `Base.depwarn`（仅在兼容路径触发）
- [x] 为每个 depwarn 提供对应替代调用写法（文档示例可复制）
- [x] 确认不破坏当前测试与关键脚本入口

### C3 测试基线强化（防静默回归）

- [x] 新增“契约失败即报错”测试（缺字段、错类型、非法边界值）
- [x] 新增“弃用路径仍可运行但给出警告”测试
- [x] 保留并扩展“新旧路径数值等价”测试（在可比参数下）
- [x] 回归运行相关 smoke tests，确保核心路径稳定

### C4 文档同步（API + 迁移指南）

- [x] 更新 `docs/api/PARAMETER_TYPES_API.md`：标注推荐输入形态与弃用项
- [x] 更新 `docs/guides/PARAMETER_STRUCT_MIGRATION.md`：补充 Phase C 迁移步骤
- [x] 在文档中给出“最小改造示例”（旧写法 -> 新写法）

---

## 完成定义（Definition of Done）

- [x] workflows/scans 公共入口参数契约清晰且一致
- [x] 兼容路径仅保留单一边界并带可读 depwarn
- [x] 契约失败、弃用警告、等价性三类测试全部通过
- [x] API 与迁移文档完成同步且示例可运行
- [ ] 本任务单在完成后按规范归档到 `docs/dev/archived/`

---

## 执行记录

- 2026-02-18：创建 Phase C 任务单，作为 Phase B 完成后的下一阶段执行入口。
- 2026-02-18：完成 C1（入口参数契约收口）首轮落地：
	- `WorkflowParamAdapters` 增加字段级校验（`quark_params.m/μ` 与 `thermo_params.T/Φ/Φbar/ξ`）。
	- `run_tmu_scan` / `run_trho_scan` 增加入参向量与枚举值校验，统一 `ArgumentError` 语义。
	- 新增并通过契约失败测试与相关 smoke 回归。
- 2026-02-18：完成 C2 + C4 首轮落地：
	- Workflow 兼容路径（NamedTuple）增加 `Base.depwarn`。
	- 新增 depwarn smoke 测试，确认“兼容可运行 + 有警告”。
	- API 与迁移指南补充 Phase C 推荐签名、弃用说明与最小替代示例。
- 2026-02-18：完成 C3 等价性扩展与性能判断：
	- 扩展 `test_workflow_paramtypes_equivalence_smoke.jl`，增加 `as_legacy_inputs` 在 struct/NamedTuple 路径下的一致性断言。
	- 抽样性能结果显示：参数归一化调用约微秒级，而代表性求解/计算路径为百微秒及以上量级；当前不建议为“纯转换开销”进行全链路激进重写，优先保持兼容并做热点定向优化。
- 2026-02-18：完成 Phase C 收尾回归：
	- 通过 workflows 参数契约/弃用/等价性测试集。
	- 通过 scans config 等价、输入契约、边界错误与 T-μ / T-ρ smoke 测试集。
