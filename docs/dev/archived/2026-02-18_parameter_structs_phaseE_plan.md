---
title: Parameter Structs Phase E 任务单（可勾选）
archived: true
original: docs/dev/active/2026_02_18_parameter_structs_phaseE_plan.md
archived_date: 2026-02-18
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Parameter Structs Phase E 任务单（可勾选）

## 文档目的

延续 Phase D 后的候选 2 优化：

- 收敛 RelaxationTime 链路中的重复参数归一化
- 保持对外兼容（NamedTuple + struct）
- 用最小改动实现可验证收益

---

## 当前状态快照（2026-02-18）

- [x] Phase D 已完成并归档
- [x] TransportRequest 入口快路径已落地并回归通过
- [ ] RelaxationTime 主入口仍存在重复归一化路径

---

## 任务清单（可勾选）

### E1 入口归一化收敛

- [x] 在 `RelaxationTime` 内引入单次归一化内部实现（NamedTuple 内核）
- [x] 移除 `relaxation_times` 的重复归一化调用
- [x] 保持 `compute_average_rates` 对外 API 不变

### E2 测试与性能验证

- [x] 增加 struct/NamedTuple 等价性测试（轻量 existing_rates 路径）
- [x] 更新 profiling 脚本增加 RelaxationTime 轻量路径对比
- [x] 运行回归测试确保无行为变化

### E3 文档同步与收尾

- [x] 更新 migration/API 文档中的 Phase E 进展
- [x] 完成任务单 DoD 并归档

---

## 完成定义（Definition of Done）

- [x] RelaxationTime 入口归一化路径单一且可读
- [x] struct 与 NamedTuple 输出等价测试通过
- [x] profiling 有可复现对比结果
- [x] 文档更新并完成归档

---

## 执行记录

- 2026-02-18：创建 Phase E 任务单，进入 RelaxationTime 链路重复归一化优化。
- 2026-02-18：完成 E1（入口归一化收敛）：
	- `RelaxationTime.jl` 增加 `_compute_average_rates_nt(...)` 内核。
	- `compute_average_rates(...)` 与 `relaxation_times(...)` 改为边界归一化一次后复用内核。
	- 移除 `relaxation_times(...)` 中重复归一化语句。
- 2026-02-18：完成 E2（测试与性能验证）：
	- `test_relaxation_time.jl` 增加 struct/NamedTuple existing-rates 等价性断言。
	- 通过 `test_relaxation_time.jl` 与 `test_transport_coefficients.jl` 回归。
	- profiling 更新后记录：`relaxation(struct)/relaxation(nt) ≈ 1.0642`。
- 2026-02-18：完成 E3（文档同步）：
	- 更新 `docs/api/PARAMETER_TYPES_API.md`（Phase E 增补）。
	- 更新 `docs/guides/PARAMETER_STRUCT_MIGRATION.md`（Phase E 进展与下一步建议）。
