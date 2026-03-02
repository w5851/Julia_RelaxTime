# Thermodynamics 收口未完成项任务单

更新日期：2026-03-02

> 来源说明：本任务单承接已归档文档
> - docs/dev/archived/2026-03-01_Thermodynamics原生化前置检查与实施任务单.md
> - docs/dev/archived/2026-03-01_Thermodynamics物理删除前迁移清单.md

---

## 1. 本轮目标

- 将 Thermodynamics 原生化收口到“可审计、可复现、可长期维护”的状态。
- 保持当前已通过的 C-4 smoke 基线，不在本任务单中引入新的物理功能。

---

## 2. 必做收口项（P0）

### 2.1 基线与留档（C0 缺口）

- [ ] 记录当前主线 commit 与工作区状态。
- [ ] 留档执行：`UNIT_PROFILE=smoke julia --project=. tests/unit/runtests.jl`。
- [ ] 留档执行：`julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base HEAD --head HEAD`。
- [ ] 留档执行：`julia --project=. scripts/dev/check_data_output_path_guard.jl --base HEAD --head HEAD`。
- [ ] 留档执行：`julia --project=. scripts/dev/run_prune_wave_gate.jl --base HEAD --head HEAD`。
- [ ] 固定点基线留档（至少 3 组 T-μ / T-ρ 代表点）。

### 2.2 接口与文档一致性

- [ ] 明确并文档化接口数据契约（输入/输出类型、单位、AD 约束）。
- [ ] 更新 API 文档与公式映射：`docs/reference/formula/models/pnjl` 与 `docs/api` 对应页面。
- [ ] 更新依赖图与审计文档，补充本轮结构收口后的映射关系。

### 2.3 门禁闭环

- [ ] 补跑并留档 dependency-audit / docs-consistency（或同等 CI 证据）。
- [ ] 在执行台账追加“命令-结果-结论”证据链。

---

## 3. 可选增强项（P1）

- [ ] 清理重复路径与同义函数，进一步收缩兼容层。
- [ ] 检查并治理 `thermo_backend` 在 `src/models/**` 的残留引用，形成“保留/移除”清单。
- [ ] 在不改变语义前提下，补充最小示例以覆盖关键入口（solver/scan/workflow）。

---

## 4. 验收标准（DoD）

- [ ] 数值一致性有固定点证据（满足既定 rtol/atol）。
- [ ] 迁移与输出路径门禁、prune-wave 门禁均有留档结果。
- [ ] 文档、依赖审计、执行台账三者一致且可追溯。
- [ ] 不引入新的 world-age/redefinition 问题。

---

## 5. 执行记录入口

- 建议将后续执行记录追加到：
  - `docs/dev/active/2026-02-25_新架构主线迁移执行台账.md`
  - 或新建同日台账文件并在本任务单中回链。
