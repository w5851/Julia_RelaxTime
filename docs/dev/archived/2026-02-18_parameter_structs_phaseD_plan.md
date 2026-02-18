---
title: Parameter Structs Phase D 任务单（可勾选）
archived: true
original: docs/dev/active/2026_02_18_parameter_structs_phaseD_plan.md
archived_date: 2026-02-18
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Parameter Structs Phase D 任务单（可勾选）

## 文档目的

在 Phase C 完成“契约收口 + 兼容边界 + depwarn”后，Phase D 聚焦：

- 用可复现 profiling 识别真实热点
- 仅对高收益路径做“结构体直通”定向特化
- 在不破坏兼容性的前提下获得可量化性能收益

---

## 当前状态快照（2026-02-18）

- [x] Phase C 任务单已归档
- [x] workflows/scans 入口契约已明确且有测试保护
- [x] NamedTuple 已收口为兼容边界并触发 depwarn
- [x] 已建立首版“参数归一化开销热点”基线报告

---

## 任务清单（可勾选）

### D1 建立 profiling 基线

- [x] 新增 `scripts/dev/profile_paramtypes_hotspots.jl`
- [x] 输出关键路径的单位调用耗时（结构体 vs NamedTuple）
- [x] 形成首版热点排序（Top N）

### D2 选择定向特化目标

- [x] 选取 1-2 个高频且收益明确的函数作为试点
- [x] 设计“结构体直通”实现，不破坏现有对外 API
- [x] 保留兼容边界与 depwarn 行为

当前试点候选：

- 候选 1：`TransportCoefficients` 中 `TransportRequest` 入口到 `transport_coefficients/shear/electric/bulk` 链路。
- 候选 2：`RelaxationTime` 链中重复 `_nt_quark/_nt_thermo` 归一化调用路径。

### D3 验证与回归

- [x] 为试点函数补充等价性测试（数值一致）
- [x] 增加最小性能回归脚本（前后对比）
- [x] 运行相关 smoke tests，确认无回归

### D4 文档同步

- [x] 更新 `docs/api/PARAMETER_TYPES_API.md`（新增性能/推荐说明）
- [x] 更新 `docs/guides/PARAMETER_STRUCT_MIGRATION.md`（Phase D 进展）

---

## 完成定义（Definition of Done）

- [x] 有可复现的 profiling 脚本与基线结果
- [x] 至少 1 个热点路径完成定向特化并通过等价性测试
- [x] 给出前后性能对比结论（可量化）
- [x] 文档同步更新

---

## 执行记录

- 2026-02-18：创建 Phase D 任务单，进入“热点驱动优化”阶段。
- 2026-02-18：完成 D1 基线 profiling（脚本：`scripts/dev/profile_paramtypes_hotspots.jl`）：
	- 代表性结果（修正后复测，us/call）：
		- `normalize(struct)` ≈ 1.86
		- `normalize(namedtuple)` ≈ 0.93
		- `A_from_equilibrium(struct)` ≈ 3.36
		- `A_from_equilibrium(namedtuple)` ≈ 720.76
		- `PNJL.solve` baseline ≈ 292.30
	- 结论：参数归一化本身仍是微秒级；优先做“热点函数定向特化”，不做全链路类型重写。
- 2026-02-18：完成 D2 试点优化（`TransportCoefficients` Request 入口快路径）：
	- 变更文件：`src/relaxtime/TransportCoefficients.jl`。
	- 去除 Request 入口对 `as_namedtuple` 的依赖，改为直接字段视图快路径。
	- 保持 NamedTuple 公共 API 与兼容边界行为不变。
- 2026-02-18：完成 D3 验证：
	- 通过 `tests/unit/relaxtime/test_transport_coefficients.jl`。
	- 通过 workflows/scans/contract 相关回归测试集。
	- profiling 显示：`transport_coeff(req struct)/transport_coeff(explicit nt) ≈ 0.9178`。
- 2026-02-18：完成 D4 文档同步：
	- 更新 `docs/api/PARAMETER_TYPES_API.md`（Phase D 增补）。
	- 更新 `docs/guides/PARAMETER_STRUCT_MIGRATION.md`（Phase D 进展与建议）。
