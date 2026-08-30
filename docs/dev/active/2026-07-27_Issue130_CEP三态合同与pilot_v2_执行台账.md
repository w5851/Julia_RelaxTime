# Issue130 CEP三态合同与pilot v2 执行台账

更新日期：2026-08-17

---

## 0. 台账定位

- 本文档仅记录“执行事实”：变更点、命令、产物、结果。
- 本阶段开发任务定义统一维护在：
  - `docs/dev/active/2026-07-27_Issue130_CEP三态合同与pilot_v2.md`

---

## 1. 记录规范（强约束）

- [x] 执行记录与开发任务分开保存：
  - 任务目标/范围/DoD 只写“开发任务单”；
  - 执行过程/命令/结果只写“执行台账”。
- [x] 追加记录时采用“直接追加”策略：
  - 不要求回读本台账历史上下文；
  - 仅按统一模板追加到文档末尾；
  - 每条记录必须自包含（含目标、命令、产物、结果）。
- [x] 每条记录必须可追溯到输出产物（`data/outputs/results/*`）。

---

## 2. 执行记录

- [ ] 批次号：phase-reference promotion audit v1
  - 目标：将 C2 自动证据与九点曲线、三点 CEP 人工 overlay 合并审计，判断是否允许 phase-reference promotion
  - 代码变更：新增 solver-free promotion audit 脚本与 versioned docs/analysis 审计包；不修改 C2、reference 或 transport
  - 验证命令：
    - `python -m py_compile scripts/analysis/pnjl/build_phase_reference_manual_overlay_promotion_audit.py`
    - `python scripts/analysis/pnjl/build_phase_reference_manual_overlay_promotion_audit.py --c1-root D:/Desktop/Julia_RelaxTime_issue130_artifacts/c1_complete_acceptance_31762201725_attempt2 --c2-root D:/Desktop/Julia_RelaxTime_issue130_artifacts/stagec_density_v2_c2_20260813_run31862752226 --c2-audit-root docs/analysis/pnjl/c2_blocking_audit_v2 --nine-review-root docs/analysis/pnjl/c2_targeted_manual_review_v1 --manual-overlay D:/Desktop/Julia_RelaxTime_issue130_artifacts/cep_manual_bisection_audit_v2_31999149922/author_decision_overlay_v2.json --output-root docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1`
    - `PowerShell manifest/hash/CSV uniqueness and finite-field validation`
    - `git diff --check`
    - `julia --project=. scripts/dev/check_docs_consistency.jl`
    - `julia --project=. scripts/dev/check_active_docs_governance.jl`
    - `julia --project=. scripts/dev/check_script_entrypoints.jl`
    - `julia --project=. scripts/dev/check_models_entry_contract.jl`
    - `julia --project=. scripts/dev/check_solver_contract_leakage.jl`
    - `julia --project=. scripts/dev/check_pnjl_migration_guard.jl`
    - `julia --project=. scripts/dev/check_data_output_path_guard.jl`
    - `python scripts/dev/check_python_script_syntax.py`
  - 产物：
    - `docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1/decision.json`
    - `docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1/manifest.json`
    - `docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1/tables/auto_gate_summary.csv`
    - `docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1/tables/grid_unresolved_by_reason.csv`
    - `docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1/tables/promotion_decision.csv`
  - 结果：blocked_manual_overlay_inconclusive; 5424 unresolved grid rows, 4661 geometry/interpolation-only and 763 non-geometry/unclassified or classification-transition rows; no reference promotion
  - 主线映射（N1/N2/N3/N4）：N3
  - 备注：C2 boundary/crossover/solver finite gates and C1-C2 physical stability pass; nine oracle curves and three CEP brackets are recorded as accepted diagnostic overlays only

- [ ] 批次号：phase-reference promotion audit governance closure
  - 目标：将 promotion audit 的阻塞 verdict、固定 run provenance 和下一步限制同步到 Issue #130 任务台账
  - 代码变更：更新 config/governance/task_tracks.toml 的 issue130-phase 与 full-hybrid-author-review evidence/next_action；不修改数值或 reference
  - 验证命令：
    - `julia --project=. scripts/dev/check_task_ledger.jl --preflight`
    - `julia --project=. scripts/dev/check_docs_consistency.jl`
    - `julia --project=. scripts/dev/check_active_docs_governance.jl`
    - `git diff --check`
  - 产物：
    - `config/governance/task_tracks.toml`
    - `docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1/decision.json`
    - `docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1/manifest.json`
  - 结果：记录完成；promotion 仍为 blocked_manual_overlay_inconclusive，RS transport/reference promotion 保持暂停
  - 主线映射（N1/N2/N3/N4）：N3
  - 备注：当前工作树含其他对话 dirty/untracked 路径，均未清理、未暂存、未修改
