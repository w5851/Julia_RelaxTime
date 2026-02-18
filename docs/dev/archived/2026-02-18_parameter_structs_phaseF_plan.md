---
title: Parameter Structs Phase F 任务单（可勾选）
archived: true
original: docs/dev/active/2026_02_18_parameter_structs_phaseF_plan.md
archived_date: 2026-02-18
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Parameter Structs Phase F 任务单（可勾选）

## 文档目的

在 Phase E 完成 RelaxationTime 入口归一化收敛后，Phase F 聚焦：

- 收敛 `AverageScatteringRate` / `TotalCrossSection` 中可避免的重复归一化
- 保持外部 API 兼容（struct + NamedTuple）
- 为高频路径减少边界转换噪声

---

## 当前状态快照（2026-02-18）

- [x] Phase E 已归档
- [ ] `AverageScatteringRate.get_sigma` 深层调用仍包含重复归一化
- [ ] `TotalCrossSection` 扫描/批量入口存在可收敛的重复归一化

---

## 任务清单（可勾选）

### F1 AverageScatteringRate 归一化收敛

- [x] 为 `get_sigma` 提供 NamedTuple 内核与兼容包装
- [x] 为 `design_w0cdf_s_grid` / `precompute_cross_section!` 提供内部 NamedTuple 内核
- [x] 保持 `average_scattering_rate` 对外签名不变

### F2 TotalCrossSection 归一化收敛

- [x] 为 `total_cross_section` 提供 NamedTuple 内核与兼容包装
- [x] 让 `calculate_all_total_cross_sections` / `scan_s_dependence` 复用内核，避免重复归一化

### F3 验证与文档

- [x] 增加/更新结构体与 NamedTuple 等价性 smoke 测试
- [x] 运行关键回归测试
- [x] 更新 API/迁移文档并归档任务单

---

## 完成定义（Definition of Done）

- [x] ASR/TCS 核心路径单边界归一化完成
- [x] 新旧参数路径行为一致（测试通过）
- [x] 关键回归通过
- [x] 文档更新并归档

---

## 执行记录

- 2026-02-18：创建 Phase F 任务单，进入 ASR/TCS 归一化收敛阶段。
- 2026-02-18：完成 F1（ASR 归一化收敛）：
	- `AverageScatteringRate.jl` 增加 `_precompute_cross_section_nt!`、`_design_w0cdf_s_grid_nt`、`_get_sigma_nt`。
	- 关键深层路径改为直接调用 NamedTuple 内核，避免重复 `_nt_*` 归一化。
- 2026-02-18：完成 F2（TCS 归一化收敛）：
	- `TotalCrossSection.jl` 增加 `_total_cross_section_nt` 与 NamedTuple 专用入口。
	- batch/scan 接口改为边界归一化一次后复用内核。
- 2026-02-18：完成 F3（验证与文档）：
	- `test_average_scattering_rate.jl` 新增 struct/NamedTuple 等价性（cache 路径）并通过。
	- `test_relaxation_time.jl`、`test_transport_coefficients.jl` 回归通过。
	- 更新 API/迁移文档；profiling 脚本保持可执行。
