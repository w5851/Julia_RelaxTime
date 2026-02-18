---
title: Parameter Structs Phase G 任务单（可勾选）
archived: true
original: docs/dev/active/2026_02_18_parameter_structs_phaseG_plan.md
archived_date: 2026-02-18
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Parameter Structs Phase G 任务单（可勾选）

## 文档目的

在 Phase F 完成 ASR/TCS 归一化收敛后，Phase G 聚焦：

- 建立 `TotalCrossSection` 的 model-ready 专项基准（含 `A` 与完整 `K_coeffs`）
- 避免通用 profiling 脚本中的模型依赖脆弱点
- 为后续 TCS 定向优化提供可复现基线

---

## 当前状态快照（2026-02-18）

- [x] Phase F 已归档
- [x] ASR/TCS 已完成单边界归一化收敛
- [ ] 缺少独立的 TCS model-ready 专项基准脚本

---

## 任务清单（可勾选）

### G1 新增 model-ready 基准脚本

- [x] 新增 `scripts/dev/profile_total_cross_section_model_ready.jl`
- [x] 自动构造 `quark_params.A` 与完整 `K_coeffs`
- [x] 输出 `total_cross_section` / `scan_s_dependence` / `calculate_all_total_cross_sections` 基线耗时

### G2 测试保障

- [x] 新增 TCS model-ready fixture smoke test
- [x] 运行相关 relaxtime 回归测试

### G3 文档与归档

- [x] 更新 API/migration 文档（Phase G 进展）
- [x] 完成任务单并归档

---

## 完成定义（Definition of Done）

- [x] TCS model-ready 基准脚本可复现运行
- [x] model-ready fixture smoke test 通过
- [x] 文档同步更新
- [x] 任务单归档完成

---

## 执行记录

- 2026-02-18：创建 Phase G 任务单，进入 TCS model-ready 专项基准阶段。
- 2026-02-18：新增 `scripts/dev/profile_total_cross_section_model_ready.jl`，并完成基线采集：
	- `total_cross_section(:uu_to_uu, s0)` = `0.0168 ms/call`
	- `scan_s_dependence(4 points)` = `0.0529 ms/call`
	- `calculate_all_total_cross_sections` = `1.0785 ms/call`
- 2026-02-18：新增并运行 `tests/unit/relaxtime/test_total_cross_section_model_ready_fixture_smoke.jl`，测试通过（6/6）。
- 2026-02-18：更新 API 与迁移文档，准备执行归档。
