# Phase-2 治理/可观测/分离部署任务单

## 1. 目标

- 对齐主任务单中的 Stage C / E / F，补齐 Phase-2 的工程化能力。
- 保持 localhost 默认可用的前提下，完成向分离部署迁移所需的最小治理与运维基础。

## 2. 关联文档

- Spec: `docs/superpowers/specs/2026-03-29-phase2-governance-observability-split-deployment-design.md`
- Plan: `docs/superpowers/plans/2026-03-29-phase2-governance-observability-split-deployment-implementation-plan.md`
- 主任务单: `docs/dev/active/2026-03-29_前后端分离中长期目标与前端短期待开发任务单.md`

## 3. Stage C（任务系统治理）

- [ ] C1 状态机收敛：`queued/running/succeeded/failed/cancelled`。
- [ ] C2 非法状态迁移拦截并可测试。
- [ ] C3 创建任务幂等键支持（重放与冲突语义明确）。
- [ ] C4 任务取消接口与行为语义明确（含终态处理）。
- [ ] C5 超时策略可配置并在终态可见原因码。

## 4. Stage E（可观测与运维）

- [ ] E1 生命周期结构化日志（含 `job_id`）。
- [ ] E2 状态接口可见治理与队列关键指标。
- [ ] E3 错误分类与前端错误映射对齐。
- [ ] E4 提供最小运维排障说明（日志定位、常见故障）。

## 5. Stage F（分离部署落地准备）

- [ ] F1 前后端分离部署配置矩阵（localhost/staging/remote）。
- [ ] F2 反向代理与 CORS allowlist 基线策略文档化。
- [ ] F3 前端 API 地址配置与后端环境变量策略可联动。
- [ ] F4 保证 localhost 默认路径不回退。

## 6. 验收清单

- [ ] V1 至少 1 条幂等重放场景通过。
- [ ] V2 至少 1 条取消场景通过。
- [ ] V3 至少 1 条超时场景通过。
- [ ] V4 状态接口返回包含治理与诊断关键字段。
- [ ] V5 分离部署 runbook 可按步骤执行（演练环境）。
- [ ] V6 localhost 现有闭环能力保持可用。

## 7. 状态记录

- [ ] 执行中持续回填每项完成证据（测试命令、日志片段、结果截图路径）。
- [ ] 完成后在主任务单中同步 Stage C/E/F 状态。
