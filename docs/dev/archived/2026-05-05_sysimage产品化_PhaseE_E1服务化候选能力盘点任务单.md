---
title: sysimage 产品化 Phase E：E1 服务化候选能力盘点任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseE_E1服务化候选能力盘点任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase E：E1 服务化候选能力盘点任务单

更新日期：2026-05-05

当前状态：已完成首轮盘点收口；已形成 Phase E 第一批服务化白名单、暂缓名单与同步/异步边界，可直接进入 E2 服务预热策略设计。

> 目的：在不回改当前物理主线与 D3 launcher PoC 结论的前提下，先盘点哪些稳定能力值得走长驻服务化，以及这些能力分别适合怎样的服务边界与预热方式。

---

## 1. 前置结论

- [x] Phase A 已完成 wrapper 默认入口收口
- [x] Phase B 已完成 metadata / bootstrap / release workflow 基础
- [x] Phase C 已完成 package 化主路径收敛
- [x] Phase D 已完成最小 launcher app PoC 主收口
- [x] 当前不再把 D3 的 launcher sysimage 长时编译结果作为主线阻塞项

---

## 2. 本批目标

- [x] E1-1 盘点现有稳定能力中哪些最适合优先服务化
- [x] E1-2 区分 `scan / point / phase / orchestrator / server` 的调用形态、生命周期与预热收益
- [x] E1-3 给出第一批服务化候选白名单与暂缓名单
- [x] E1-4 明确进入 E2 预热策略设计前所需的最小边界结论

---

## 3. 范围与非范围

### 3.1 本批范围

- [x] 稳定 CLI / helper / launcher 能力盘点
- [x] 服务化收益与风险对照
- [x] 长驻进程价值判断
- [x] Phase E 文档与 backlog 状态同步

### 3.2 非范围

- [ ] 本批不直接实现 HTTP / daemon / worker 服务进程
- [ ] 本批不修改活跃物理求解主线
- [ ] 本批不重新打开 D3 的 launcher 编译验证
- [ ] 本批不把所有稳定脚本一次性提升为服务入口

---

## 4. 盘点维度

### 4.1 能力分类

- [x] 单点求解能力
- [x] 扫描型批处理能力
- [x] phase 结构能力
- [x] orchestrator/工作流能力
- [x] 现有 server 能力

### 4.2 服务化判断口径

- [x] 是否有明显冷启动 / 预热收益
- [x] 是否适合重复调用、共享进程态缓存或共享 sysimage
- [x] 是否已有较稳定 API / 输入输出契约
- [x] 是否会引入状态污染、并发或资源占用风险
- [x] 是否需要与 CLI 入口长期并存

---

## 5. 验收标准

- [x] 已给出第一批优先服务化候选
- [x] 已给出不建议优先服务化的能力类别及原因
- [x] 已明确 E2 预热策略设计围绕哪些入口展开
- [x] backlog 已同步 Phase E 首批状态

---

## 6. DoD

- [x] Phase E 已从“概念阶段”进入“有白名单/边界的执行阶段”
- [x] 后续可直接进入 E2 服务预热策略设计
- [x] 不与当前活跃物理开发主线冲突

---

## 7. 盘点输入与现状证据

- [x] 稳定 API 主要已收敛在 [src/models/entrypoints.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/models/entrypoints.jl)
  - `solve_pnjl_point`
  - `solve_gap_and_transport`
  - `solve_transport_from_equilibrium`
  - `run_phase_pipeline`
  - `run_scan_pipeline`
  - `run_relaxtime_orchestrator_pipeline`
- [x] 稳定 CLI 已基本形成薄壳：
  - [scripts/pnjl/calculate_phase_structure.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/pnjl/calculate_phase_structure.jl)
  - [scripts/models/run_unified_scan.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/models/run_unified_scan.jl)
  - [scripts/relaxtime/run_relaxtime_orchestrator.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/relaxtime/run_relaxtime_orchestrator.jl)
- [x] transport scan 仍是脚本装配层较厚，但 C2 已把扫描计划、orchestration、phase/equilibrium、provenance 拆成 helper：
  - [scripts/relaxtime/run_gap_transport_scan.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/relaxtime/run_gap_transport_scan.jl)
- [x] 现有 service 壳已经存在：
  - [scripts/server/server_full.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/scripts/server/server_full.jl)
  - [src/simulation/ServerLauncher.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/simulation/ServerLauncher.jl)
  - [src/simulation/FullServerApp.jl](C:/Users/Wmzx/.codex/worktrees/c75d/Julia_RelaxTime/src/simulation/FullServerApp.jl)
- [x] 已有现成服务边界证据：
  - 同步 point：`POST /api/modules/pnjl-gap/run`
  - 异步 scan job：`POST /api/modules/pnjl-scan/jobs` + status/result/cancel

---

## 8. 分类盘点结论

| 能力类别 | 当前代表入口 | 调用形态 | 预热收益 | 契约稳定度 | 并发/状态风险 | E1 结论 |
|---|---|---|---|---|---|---|
| `server` | `scripts/server/server_full.jl` | 长驻进程 | 高 | 中高 | 中 | `优先作为 service 壳层继续演进` |
| `point` | `Models.solve_pnjl_point`、`Models.solve_transport_from_equilibrium` | 同步 request/response | 高 | 高 | 中低 | `第一批优先服务化白名单` |
| `phase` | `Models.run_phase_pipeline`、`phase` CLI | 任务型、产物型 | 高 | 中高 | 中 | `第二批候选，宜走异步 job` |
| `scan` | `Models.run_scan_pipeline`、`run_gap_transport_scan.jl` | 长时批处理 | 高 | 中高 | 高 | `不做同步服务；保留 CLI + job 化` |
| `orchestrator` | `Models.run_relaxtime_orchestrator_pipeline` | 配置驱动批处理 | 中高 | 中高 | 高 | `暂缓为直接服务入口；先视作 worker/job` |

### 8.1 单点求解能力

- [x] `solve_pnjl_point` 已经具备最典型的服务化形态：
  - 输入是有限参数集
  - 输出是单个结构化结果
  - 现有 HTTP handler 已经存在
  - 冷启动/JIT 摊销收益明确
- [x] `solve_transport_from_equilibrium` 与 `solve_gap_and_transport` 虽尚未作为现成 HTTP endpoint 暴露，但从 API 稳定性和 sysimage 受益角度看，属于最值得在 E2/E3 继续展开的 point-family 候选。

结论：

- [x] point-family 是 Phase E 第一优先级
- [x] 其服务目标应是“低延迟同步调用 + 共享预热态”，而不是再包一层批处理语义

### 8.2 扫描型批处理能力

- [x] `run_scan_pipeline` / `run_tmu_scan` / `run_trho_scan` 已经更像作业系统，而不是同步 API
- [x] 现有 `FullServerApp` 已明确采用：
  - 创建任务
  - 查询状态
  - 拉取结果
  - 取消任务
  的 job 模式
- [x] `run_gap_transport_scan.jl` 仍具有明显的：
  - 长执行时间
  - 文件输出依赖
  - provenance/sidecar 产物
  - 并发资源占用

结论：

- [x] scan-family 不应作为第一批同步服务入口
- [x] 更合理路径是“CLI 继续作为权威批处理入口，service 侧只承接 job 提交/队列/状态管理”

### 8.3 phase 能力

- [x] `phase` 已经有薄 CLI 和稳定 pipeline 边界
- [x] phase 产物天然带有：
  - 较长运行时
  - 文件/目录产物
  - reference/promote 等流程语义
- [x] 因此 phase 更接近“任务型服务”，而不是 point API

结论：

- [x] phase 属于第二批候选
- [x] 若进入服务化，优先采用异步 job / artifact endpoint，而不是先做同步接口

### 8.4 orchestrator 能力

- [x] `run_relaxtime_orchestrator.jl` 当前以 `config_path / outdir / resume / overwrite / processes` 为主
- [x] 这类入口本质上是“批处理工作流编排”，更像 worker，而不是面向交互式调用的细粒度服务
- [x] 若过早直接服务化，容易把 CLI 批处理契约、路径约束、队列治理和在线请求语义耦在一起

结论：

- [x] orchestrator 暂不进入第一批直接服务化白名单
- [x] 后续若需要上线，应优先挂到 job/worker 层，而不是直接扩成同步 endpoint

### 8.5 server 能力

- [x] 当前 `ServerLauncher + FullServerApp` 已经是最自然的 service 壳
- [x] 其中已经内置了：
  - deploy profile / runtime policy
  - scan queue / timeout / idempotency / cancel / metrics
  - point 与 scan 的初步分流

结论：

- [x] Phase E 不需要另造第二套 service 容器
- [x] 后续 service 化主线应优先扩展现有 server 壳，而不是新开并行 daemon 体系

---

## 9. 第一批白名单与暂缓名单

### 9.1 第一批优先服务化白名单

- [x] `server` 作为统一长驻 service 壳层
- [x] `pnjl-gap` / `solve_pnjl_point` 作为同步 point endpoint
- [x] `transport-point` family 作为下一批最值得补进 server 的同步 point endpoint

### 9.2 第二批候选

- [x] `phase` pipeline job 化
- [x] `suscept` 单点能力，可作为 point-family 扩展候选

### 9.3 暂缓名单

- [x] `unified-scan` 直接同步服务化
- [x] `transport-scan` 直接同步服务化
- [x] `relaxtime-orchestrator` 直接顶层 endpoint 化

暂缓原因：

- [x] 长时运行、文件输出和队列治理需求明显强于同步请求语义
- [x] 更适合作为 job/worker 挂到现有 server 壳之后
- [x] 若过早服务化，会把 CLI 权威入口、文件产物治理和在线调用契约混在一起

---

## 10. 进入 E2 的最小边界结论

- [x] Phase E 的 service 主体沿用现有 `server_full` / `FullServerApp`
- [x] 请求形态明确分两类：
  - 同步：point-family
  - 异步：phase / scan / orchestrator
- [x] E2 预热策略应优先围绕：
  - `Models.solve_pnjl_point`
  - `Models.solve_transport_from_equilibrium`
  - `Models.run_phase_pipeline`
- [x] CLI 入口继续保留，service 不替代 CLI 的批处理/研究工作流职责
- [x] 后续若要扩大 endpoint 集，应先判断其是否属于 point-family；若不是，优先考虑 job 化而不是同步化

---

## 11. 治理校验

- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`
