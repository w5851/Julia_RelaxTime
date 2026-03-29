---
title: Phase-2 治理/可观测/分离部署任务单
archived: true
original: docs/dev/active/2026-03-29_phase2_governance_observability_split_deployment_tasklist.md
archived_date: 2026-03-29
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Phase-2 治理/可观测/分离部署任务单

## 1. 目标

- 对齐主任务单中的 Stage C / E / F，补齐 Phase-2 的工程化能力。
- 保持 localhost 默认可用的前提下，完成向分离部署迁移所需的最小治理与运维基础。

## 2. 关联文档

- Spec: `docs/superpowers/specs/2026-03-29-phase2-governance-observability-split-deployment-design.md`
- Plan: `docs/superpowers/plans/2026-03-29-phase2-governance-observability-split-deployment-implementation-plan.md`
- 主任务单: `docs/dev/active/2026-03-29_前后端分离中长期目标与前端短期待开发任务单.md`

## 3. Stage C（任务系统治理）

- [x] C1 状态机收敛：`queued/running/succeeded/failed/cancelled`。
- [x] C2 非法状态迁移拦截并可测试。
- [x] C3 创建任务幂等键支持（重放与冲突语义明确）。
- [x] C4 任务取消接口与行为语义明确（含终态处理）。
- [x] C5 超时策略可配置并在终态可见原因码。

## 4. Stage E（可观测与运维）

- [x] E1 生命周期结构化日志（含 `job_id`）。
- [x] E2 状态接口可见治理与队列关键指标。
- [x] E3 错误分类与前端错误映射对齐。
- [x] E4 提供最小运维排障说明（日志定位、常见故障）。

## 5. Stage F（分离部署落地准备）

- [x] F1 前后端分离部署配置矩阵（localhost/staging/remote）。
- [x] F2 反向代理与 CORS allowlist 基线策略文档化。
- [x] F3 前端 API 地址配置与后端环境变量策略可联动。
- [x] F4 保证 localhost 默认路径不回退。

## 6. 验收清单

- [x] V1 至少 1 条幂等重放场景通过。
- [x] V2 至少 1 条取消场景通过。
- [x] V3 至少 1 条超时场景通过。
- [x] V4 状态接口返回包含治理与诊断关键字段。
- [x] V5 分离部署 runbook 可按步骤执行（演练环境）。
- [x] V6 localhost 现有闭环能力保持可用。

## 7. 状态记录

- [x] 执行中持续回填每项完成证据（测试命令、日志片段、结果截图路径）。
- [x] 2026-03-29 Task1 证据：新增 `tests/unit/simulation/test_scan_job_state_machine.jl`；执行失败验证命令 `julia --project=. -e "ENV[\"UNIT_FILES\"]=\"simulation/test_scan_job_state_machine.jl\"; include(\"tests/unit/runtests.jl\")"`（初次失败，缺失状态机函数）；实现后复跑同命令通过（18/18）；回归 `julia --project=. -e "ENV[\"UNIT_FILES\"]=\"simulation/test_pnjl_scan_jobs.jl\"; include(\"tests/unit/runtests.jl\")"` 通过（69/69）。
- [x] 2026-03-29 Task2 证据：新增 `tests/integration/relaxtime/test_pnjl_scan_idempotency.jl`；先运行 `julia --project=. -e "include(\"tests/integration/relaxtime/test_pnjl_scan_idempotency.jl\")"` 观察失败（同 key 未重放/冲突未拦截）；实现幂等键提取、请求指纹与缓存后复跑通过（11/11）；并复跑 `simulation/test_scan_job_state_machine.jl`（18/18）与 `simulation/test_pnjl_scan_jobs.jl`（69/69）。
- [x] 2026-03-29 Task3 证据：新增 `tests/integration/relaxtime/test_pnjl_scan_cancel_timeout.jl`；先运行 `julia --project=. -e "include(\"tests/integration/relaxtime/test_pnjl_scan_cancel_timeout.jl\")"` 失败（缺失 cancel/timeout 处理）；实现 `/api/modules/pnjl-scan/jobs/{job_id}/cancel`、`handle_pnjl_scan_job_cancel` 与超时钩子 `_maybe_mark_job_timeout!` 后复跑通过（10/10）；回归 `simulation/test_scan_job_state_machine.jl`（18/18）与 `simulation/test_pnjl_scan_jobs.jl`（69/69）。
- [x] 2026-03-29 Task4 证据：新增 `tests/unit/simulation/test_scan_job_logging_contract.jl`；先运行 `julia --project=. -e "ENV[\"UNIT_FILES\"]=\"simulation/test_scan_job_logging_contract.jl\"; include(\"tests/unit/runtests.jl\")"` 失败（status 未暴露 events）；实现 `_new_job_event/_append_job_event!` 并在 create/start/progress/end 写入结构化事件后复跑通过（13/13）；回归 `test_pnjl_scan_idempotency.jl`（11/11）、`test_pnjl_scan_cancel_timeout.jl`（10/10）、`simulation/test_scan_job_state_machine.jl`（18/18）、`simulation/test_pnjl_scan_jobs.jl`（69/69）。
- [x] 2026-03-29 Task5 证据：新增 `tests/integration/relaxtime/test_pnjl_scan_metrics_diagnostics.jl`；先运行 `julia --project=. -e "include(\"tests/integration/relaxtime/test_pnjl_scan_metrics_diagnostics.jl\")"` 失败（缺少 runtime metrics 更新函数）；实现 `terminal` 与 `duration_buckets` 最小指标计数并在 status 暴露 `metrics` 后复跑通过（15/15）；回归 `test_pnjl_scan_idempotency.jl`（11/11）、`test_pnjl_scan_cancel_timeout.jl`（10/10）、`simulation/test_scan_job_logging_contract.jl`（13/13）、`simulation/test_pnjl_scan_jobs.jl`（69/69）。
- [x] 2026-03-29 Task6 证据：新增 `tests/integration/relaxtime/test_frontend_backend_config_contract.jl` 与 `web/js/runtime_config.contract.test.mjs`；先运行后端契约测试失败（缺失 `server_runtime_policy/runtime_policy_env`），实现 `ServerLauncher` 的 profile 策略矩阵与 env 导出后复跑通过（16/16）；前端配置契约测试通过（`runtime_config.contract.test.mjs: PASS`）；并复跑 `test_pnjl_scan_metrics_diagnostics.jl`（15/15）、`test_pnjl_scan_cancel_timeout.jl`（10/10）、`test_pnjl_scan_idempotency.jl`（11/11）确认 localhost 主流程未回退。
- [x] 2026-03-29 Task7 证据：新增 runbook `docs/reference/ops/phase2_split_deployment_runbook.md`，覆盖 localhost/staging/remote 启动、健康检查、回滚与排障；主任务单已同步 Stage C/E/F 与 Phase-2 对齐记录；执行 `julia --project=. scripts/dev/check_docs_consistency.jl` 输出 `[docs-consistency] OK`。
- [x] 完成后在主任务单中同步 Stage C/E/F 状态。
