---
title: sysimage 产品化 Phase E：E2 服务预热策略任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseE_E2服务预热策略任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase E：E2 服务预热策略任务单

更新日期：2026-05-05

当前状态：已完成首批收口；已形成 service warmup profile 契约，并将最小 server warmup hook 接入现有 `server_full` 壳层。

> 目的：在不扩新 endpoint、不回改物理主线求解逻辑的前提下，把 Phase E 的“长驻服务进程可复用预热态”落实为明确 profile、默认映射和最小启动 hook。

---

## 1. 前置结论

- [x] E1 已确认 service 主体沿用现有 `server_full` / `FullServerApp`
- [x] E1 已确认 point-family 是第一批优先服务化白名单
- [x] E1 已确认 `phase / scan / orchestrator` 不宜先走同步 endpoint
- [x] 已有 sysimage / precompile registry 能力可复用，无需发明第二套预热路径

---

## 2. 本批目标

- [x] E2-1 定义 service warmup profile 契约
- [x] E2-2 明确 deploy profile 到 warmup profile 的默认映射
- [x] E2-3 把 warmup hook 以最小侵入方式接入 `ServerLauncher`
- [x] E2-4 为后续 E3 endpoint 扩展保留稳定的 step/profile 边界

---

## 3. 范围与非范围

### 3.1 本批范围

- [x] `ServerWarmup` 模块
- [x] `ServerLauncher` runtime policy / env contract 扩展
- [x] integration contract test 补充
- [x] E2 文档与 backlog 状态同步

### 3.2 非范围

- [ ] 本批不新增 transport-point HTTP endpoint
- [ ] 本批不重写 precompile registry
- [ ] 本批不让 scan / phase / orchestrator 变成同步 endpoint
- [ ] 本批不追求完整启动耗时 benchmark 收口

---

## 4. 设计结论

### 4.1 warmup profile 契约

- [x] `none`
  - 不做服务启动预热
- [x] `point`
  - 预热 `pnjl_point`
  - 预热 `transport_point`
- [x] `service_core`
  - 包含 `point`
  - 额外预热 `phase_pipeline`

对应 step 集合：

- [x] `:none -> []`
- [x] `:point -> [:pnjl_point, :transport_point]`
- [x] `:service_core -> [:pnjl_point, :transport_point, :phase_pipeline]`

### 4.2 deploy profile 默认映射

- [x] `localhost -> warmup_profile=none`
- [x] `staging -> warmup_profile=point`
- [x] `remote -> warmup_profile=service_core`

设计理由：

- [x] 本地开发默认不强行增加启动等待
- [x] staging 优先收敛 point-family 的交互延迟
- [x] remote 才承担更完整的 phase pipeline 预热成本

### 4.3 hook 接入方式

- [x] 继续沿用 `ServerLauncher.run_full_server(...)` 作为唯一接入点
- [x] 通过 env 契约暴露：
  - `JRT_SERVER_WARMUP_PROFILE`
  - `JRT_SERVER_WARMUP_STRICT`
- [x] 仅在 profile 非 `none` 时执行 warmup
- [x] 默认 `strict=false`，避免服务因单步 warmup 失败而无法启动

---

## 5. 已落地实现

- [x] 新增 [src/simulation/ServerWarmup.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/simulation/ServerWarmup.jl)
  - `resolve_server_warmup_profile`
  - `list_server_warmup_steps`
  - `run_server_warmup`
  - `run_server_warmup_from_env`
- [x] 更新 [src/simulation/ServerLauncher.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/simulation/ServerLauncher.jl)
  - runtime policy 新增 `warmup_profile / warmup_strict`
  - env map 新增 `JRT_SERVER_WARMUP_PROFILE / JRT_SERVER_WARMUP_STRICT`
  - `run_full_server(...)` 在启动前执行最小 warmup hook
- [x] 更新 [tests/integration/relaxtime/test_frontend_backend_config_contract.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/tests/integration/relaxtime/test_frontend_backend_config_contract.jl)
  - 校验 warmup profile 默认映射
  - 校验 env contract
  - 校验 warmup step 集合稳定性

---

## 6. 与 precompile/sysimage 的关系

- [x] server warmup 不替代 sysimage
- [x] sysimage 负责跨进程/跨重启复用编译产物
- [x] warmup 负责单进程启动后的热点路径实跑与缓存填充
- [x] `transport_point` 直接复用现有 `Models.run_precompile_capability(:transport_point_api)`，避免再造第三套“服务专用预热代码”

---

## 7. 验收标准

- [x] 已有明确 warmup profile 契约
- [x] deploy profile 与 warmup profile 的默认映射已固定
- [x] `server_full` 已具备最小 warmup hook
- [x] 后续 E3 可以直接围绕 `point` / `service_core` profile 扩 endpoint 或补 benchmark

---

## 8. DoD

- [x] E2 不再停留在“口头预热策略”
- [x] 现有 server 壳已具备最小可执行 warmup 入口
- [x] 后续可进入 E3：CLI 入口与服务入口职责边界的进一步落地

---

## 9. 验证

- [x] `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_frontend_backend_config_contract.jl"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`
