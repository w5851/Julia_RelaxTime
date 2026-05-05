---
title: sysimage 产品化 Phase E：E3 CLI 与服务职责边界任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseE_E3_CLI与服务职责边界任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase E：E3 CLI 与服务职责边界任务单

更新日期：2026-05-05

当前状态：已完成首批收口；已把 point/job 的服务边界固化为机器可读 discovery 契约，并同步写回数据契约文档。

> 目的：避免 CLI 入口和服务入口长期漂移，把“谁是权威执行入口、谁是默认在线入口、哪些能力只能走 job”收敛成可读、可测、可发现的固定边界。

---

## 1. 前置结论

- [x] E1 已完成 service 候选能力盘点
- [x] E2 已完成 server warmup profile / hook 首批收口
- [x] 当前已有明确 service 路由：
  - `pnjl-gap/run`
  - `pnjl-scan/jobs`
- [x] 当前已有明确 CLI 白名单：
  - `phase`
  - `unified-scan`
  - `transport-scan`
  - `relaxtime-orchestrator`
  - `suscept`
  - `server`

---

## 2. 本批目标

- [x] E3-1 明确 CLI 和 service 的职责矩阵
- [x] E3-2 固化 point 与 job 的服务表面边界
- [x] E3-3 让 `/api/modules` discovery 返回可机器消费的边界元数据
- [x] E3-4 把边界规则同步进 `docs/api/data_contracts.md`

---

## 3. 范围与非范围

### 3.1 本批范围

- [x] service module discovery contract
- [x] CLI/service 边界文档
- [x] integration contract test
- [x] backlog 状态同步

### 3.2 非范围

- [ ] 本批不新增 `transport-point` endpoint
- [ ] 本批不让 `phase / scan / orchestrator` 变成同步接口
- [ ] 本批不改写现有 CLI 参数语法
- [ ] 本批不碰活跃物理求解主线

---

## 4. 边界结论

### 4.1 service 默认入口

- [x] 面向前端、在线 agent、交互式调用
- [x] 只暴露两类表面：
  - `point`
  - `job`

### 4.2 CLI 默认入口

- [x] 面向离线批处理、研究扫描、产物落盘、人工调参
- [x] 对 `phase / transport-scan / orchestrator / suscept / unified-scan` 保持权威脚本入口地位

### 4.3 固定规则

- [x] `point-family`
  - 可同时拥有 CLI/API
  - service 为默认在线入口
- [x] `scan / phase / orchestrator`
  - CLI 为权威执行入口
  - service 仅可提供 job 壳
  - 不提供同步长任务接口

---

## 5. 已落地实现

- [x] 更新 [src/simulation/fullserver/shared.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/simulation/fullserver/shared.jl)
  - `MODULE_REGISTRY` 新增：
    - `invocation_style`
    - `service_surface`
    - `default_client_surface`
    - `stable_entrypoint`
    - `http`
- [x] 新增 [tests/integration/relaxtime/test_service_module_registry_contract.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/tests/integration/relaxtime/test_service_module_registry_contract.jl)
  - 校验 `pnjl-gap` 为 `sync/point`
  - 校验 `pnjl-scan` 为 `async/job`
  - 校验 discovery 返回的 HTTP 路径与稳定入口标记
- [x] 更新 [docs/api/data_contracts.md](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/docs/api/data_contracts.md)
  - 增补 E3 固化后的 CLI/service 边界规则
  - 增补 `GET /api/modules` discovery contract

---

## 6. 验收标准

- [x] CLI/service 的默认职责已写清
- [x] `/api/modules` 返回的元数据足以让前端/agent 判断该走 point 还是 job
- [x] 数据契约文档与实际 discovery 返回不再分离
- [x] 后续 E4 若继续扩 service，只需在这套边界内补新 module

---

## 7. DoD

- [x] E3 已不再停留在口头约定
- [x] 服务边界已具备机器可读 discovery 契约
- [x] 后续可继续进入下一阶段 service 扩展，而不必重新争论入口职责

---

## 8. 验证

- [x] `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_service_module_registry_contract.jl"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`
