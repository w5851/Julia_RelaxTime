# sysimage 产品化 Phase C：package 化 C2 第三批 phase/equilibrium 边界收敛任务单

更新日期：2026-05-04

当前状态：第三批已完成；`run_gap_transport_scan.jl` 中的 `phase tracker + equilibrium dispatch` 已抽离为 script-private helper。

> 目的：在不改动底层物理求解语义的前提下，继续压薄 `run_gap_transport_scan.jl`，把相对独立的相态跟踪与平衡态调度逻辑收敛成可复用 helper，并评估其是否已经具备上升为稳定 API 的条件。

---

## 1. 目标

- [x] C2-7 抽离 `phase tracker + equilibrium dispatch` helper
- [x] C2-8 保持原入口、输出契约与求解语义不变
- [x] C2-9 给出“是否升格为稳定 API”的阶段性判断

---

## 2. 本批范围

### 2.1 包含

- [x] 新增 `scripts/relaxtime/gap_transport_scan_phase_equilibrium.jl`
- [x] 更新 `scripts/relaxtime/run_gap_transport_scan.jl`
- [x] 更新 `tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl`

### 2.2 不包含

- [x] 不改动 `src/models/` / `src/relaxtime/` 中的底层求解实现
- [x] 不改动 mixed-meson governance、non-fixedmu unified joint-solve 语义
- [x] 不把该 helper 直接提升为 `src/` 稳定公共 API

---

## 3. 实现结果

### 3.1 新 helper 边界

- [x] 新增 `GapTransportScanPhaseEquilibrium` 模块
- [x] 收敛的职责包括：
  - `LocalPhaseTracker`
  - `build_phase_tracker`
  - `tracker_phase`
  - `phase_structure`
  - `describe_seed_source`
  - `classify_phase_from_solution`
  - `solve_models_equilibrium`
  - `solve_equilibrium_with_diagnostics`
  - `is_phase_transition`

收益：

- `run_gap_transport_scan.jl` 不再直接承载大段相态参考数据加载、seed 选择、equilibrium fallback 与诊断拼装逻辑
- 主脚本进一步向 thin-wrapper + orchestration body 收敛

### 3.2 当前边界判断

- [x] 已形成“脚本私有但清晰”的 helper 边界
- [x] 暂不建议直接升格到 `src/` 稳定 API

原因：

- 该层仍显著依赖脚本运行时上下文与 reference CSV 路径/缓存
- 复用场景目前主要仍限于 `run_gap_transport_scan.jl`
- 若过早提升到 `src/`，容易把尚未稳定的 orchestration 语义固化为公共契约

结论：

- 先保留在 `scripts/relaxtime/` 作为 script-private helper
- 待后续出现第二个稳定调用方、或 precompile/workload 明确需要绑定此层时，再评估上升

---

## 4. 验证

- [x] `julia --project=. -e 'include("tests/unit/relaxtime/test_run_gap_transport_scan_solver_entry.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_run_scan_subcommand_dispatch_smoke.jl"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`

---

## 5. 当前脚本边界状态

截至本批，`scripts/relaxtime/run_gap_transport_scan.jl` 的脚本层已拆出：

- [x] CLI 参数解析
- [x] IO / header / sidecar row 写入
- [x] provenance 汇总壳层
- [x] phase tracker / equilibrium dispatch helper

仍保留在主脚本中的主要职责：

- [ ] transport scan 主执行编排
- [ ] 渠道诊断行写入与结果拼装
- [ ] 具体扫描循环的 orchestration shell

---

## 6. 下一步建议

- [ ] C2-下一批：评估 `run_scan` 主执行编排是否还能继续压薄成 orchestration helper
- [ ] C3：让 precompile capability / workload 优先绑定已较稳定的 helper 边界，而不是脚本偶然路径
- [ ] 后续若出现第二个稳定调用方，再评估 phase/equilibrium helper 是否迁移到 `src/`

---

## 7. DoD

- [x] `phase tracker + equilibrium dispatch` 已抽离为独立 helper
- [x] 原入口、输出契约与求解语义保持不变
- [x] 已给出是否升格为稳定 API 的阶段性判断
