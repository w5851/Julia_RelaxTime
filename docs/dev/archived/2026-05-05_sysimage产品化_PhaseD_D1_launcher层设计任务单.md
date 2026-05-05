---
title: sysimage 产品化 Phase D：D1 launcher 层设计任务单
archived: true
original: docs/dev/active/2026-05-05_sysimage产品化_PhaseD_D1_launcher层设计任务单.md
archived_date: 2026-05-05
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# sysimage 产品化 Phase D：D1 launcher 层设计任务单

更新日期：2026-05-05

当前状态：已完成首轮设计收口；已形成 launcher 职责边界、关系梳理与首批白名单初稿，可直接进入 D2。

> 目的：把当前已经验证有效的 wrapper/sysimage/stable CLI 能力，收敛为一个对用户与零知识 agent 都更稳定、更统一的 launcher 层，为后续 `create_app(...)` 评估与最小 app 分发 PoC 建立固定入口。

---

## 1. 背景

- [x] Phase A 已完成 wrapper 默认化首轮收口
- [x] Phase B 已完成 metadata / bootstrap / release workflow 基础
- [x] Phase C 已完成 stable CLI 与 stable API 的主路径收敛
- [x] 当前 C3 已进入 perf refinement 区间，不再适合作为主线继续深挖

当前缺口：

- [x] 用户侧仍然需要理解“脚本入口 / wrapper / sysimage / bootstrap”之间的关系
- [x] 零知识 agent 仍缺少一个更强约束的“默认稳定启动入口”
- [x] 后续若要 app 化，先需要定义 launcher 的职责边界，避免直接把当前脚本集合粗暴打包

---

## 2. 目标

- [x] D1-1 明确 launcher 层的职责边界
- [x] D1-2 定义 launcher 与现有 wrapper / bootstrap / stable CLI 的关系
- [x] D1-3 选出进入 launcher 首批白名单的稳定能力
- [x] D1-4 形成可供 D2 / D3 继续执行的 launcher 设计结论

---

## 3. 范围与非范围

### 3.1 本批范围

- [x] launcher 概念与职责定义
- [x] launcher 命令面 / 参数面 / fallback 面设计
- [x] 与 `run_with_sysimage.*`、`bootstrap_sysimage.*` 的关系梳理
- [x] 对首批稳定能力的 launcher 白名单筛选

### 3.2 非范围

- [ ] 本批不直接实现完整 app
- [ ] 本批不直接切到 service 化
- [ ] 本批不修改活跃物理主线求解逻辑
- [ ] 本批不把尚不稳定的脚本强行纳入 launcher 白名单

---

## 4. 关键设计问题

### 4.1 launcher 是什么

- [x] launcher 是“统一稳定入口层”，而不是新的求解实现层
- [x] launcher 负责：
  - 选择稳定能力
  - 调度 wrapper / sysimage / bootstrap
  - 提供一致的错误提示与回退策略
- [x] launcher 不负责：
  - 重新实现业务求解逻辑
  - 替代 `src/models` 的 stable API

设计结论：

- [x] launcher 应定义为“用户/agent 看到的唯一默认入口壳层”
- [x] launcher 下层继续复用：
  - `scripts/dev/run_with_sysimage.ps1/.sh`
  - `scripts/dev/bootstrap_sysimage.ps1/.sh`
  - 既有稳定 CLI 脚本
- [x] launcher 自身不承载求解参数解释之外的物理业务逻辑
- [x] 后续若进入 `create_app(...)`，app 应优先打包 launcher，而不是把全部脚本直接暴露为 app 顶层命令

### 4.2 launcher 与现有组件关系

- [x] `run_with_sysimage.*`
  - 作为 launcher 的执行基础设施，负责优先走 sysimage
- [x] `bootstrap_sysimage.*`
  - 作为 launcher 的准备/获取路径，负责匹配或下载产物
- [x] stable CLI
  - 作为 launcher 的被调度目标，而不再直接面向用户暴露为“推荐默认入口”

关系结论：

- [x] `bootstrap_sysimage.*` 负责“准备可用 sysimage”
- [x] `run_with_sysimage.*` 负责“在本机决定是否使用 sysimage 执行 Julia”
- [x] launcher 负责“把用户意图映射到稳定能力，并选择调用哪个稳定 CLI”
- [x] stable CLI 保持为可独立调用的底层入口，但文档与 agent 默认策略应优先引用 launcher

推荐分层：

1. 用户/agent：调用 launcher
2. launcher：解析能力名与高层子命令
3. wrapper：根据 mismatch policy / 本地兼容性选择 sysimage 或 plain julia
4. stable CLI：调用 `Models` / `src/models/entrypoints.jl` 暴露的稳定能力

### 4.2.1 launcher 命令面建议

- [x] 首轮建议采用“能力名 + 原脚本参数透传”的薄命令面，而不是重新发明全套参数语法
- [x] 推荐形态：
  - `launcher phase ...`
  - `launcher unified-scan ...`
  - `launcher transport-scan ...`
  - `launcher relaxtime-orchestrator ...`
  - `launcher suscept ...`
  - `launcher server ...`
- [x] launcher 应保留高层控制参数：
  - `--sysimage-policy=<fallback|strict|rebuild>`
  - `--bootstrap-if-missing`
  - `--dry-run`
  - `--show-target`
- [x] 业务参数默认透传给底层稳定 CLI，避免在 D1 阶段制造第二套参数契约

### 4.2.2 fallback 策略建议

- [x] 默认策略沿用 wrapper 当前契约：`fallback`
- [x] 对零知识 agent，建议默认：
  - 本机已有 sysimage 时走 launcher + wrapper 默认策略
  - fresh clone / 新机器先执行 bootstrap，再走 launcher
- [x] launcher 应统一错误语义：
  - 找不到匹配能力：列出白名单能力
  - sysimage 不兼容：回显 wrapper 的 mismatch reason
  - 目标 CLI 不在白名单：明确提示“可直接脚本调用，但不是 launcher 稳定入口”

### 4.3 launcher 首批白名单

- [x] 至少评估以下能力是否进入首批：
  - transport point
  - gap transport scan
  - phase / pipeline 相关稳定入口
- [x] 每个候选项都要标注：
  - 是否已有稳定 API
  - 是否已有 sysimage 收益证据
  - 是否已有清晰输入/输出契约

#### 首批 launcher 白名单建议

| launcher 能力名 | 底层稳定入口 | 稳定 API / 边界 | sysimage 证据 | 契约清晰度 | 首批结论 |
|---|---|---|---|---|---|
| `phase` | `scripts/pnjl/calculate_phase_structure.jl` | 已对齐 `Models.run_phase_pipeline` 薄 CLI | 有；README / Quickstart 已默认经 wrapper 启动 | 高 | `纳入` |
| `unified-scan` | `scripts/models/run_unified_scan.jl` | 已属 `Models` 统一扫描主链 | 有；已列入稳定白名单 | 高 | `纳入` |
| `transport-scan` | `scripts/relaxtime/run_gap_transport_scan.jl` | 已完成 C2 主路径收敛，脚本基本退回装配职责 | 有；C3 多批 residual / cold-start 证据 | 中高 | `纳入` |
| `relaxtime-orchestrator` | `scripts/relaxtime/run_relaxtime_orchestrator.jl` | 属稳定工作流编排入口 | 有；已列入稳定白名单 | 中高 | `纳入` |
| `suscept` | `scripts/pnjl/run_conserved_charge_susceptibilities.jl` | 已列稳定白名单 | 有；当前以 wrapper 默认化治理为主 | 中高 | `纳入` |
| `server` | `scripts/server/server_full.jl` | 入口稳定，但更接近服务/联调场景 | 有；Quickstart 已统一入口 | 中高 | `纳入，但单列为 service-like 能力` |
| `transport-point` | 无单独稳定 CLI；当前主要通过 API/workload/probe 证据存在 | 稳定 API 已有，但用户入口尚不独立 | 有；C3 证据最充分 | 中 | `不单独纳入；并入 transport-scan / 后续服务化` |

白名单结论：

- [x] launcher 首批优先承载“已经在 README / scripts guide 白名单中存在的稳定 CLI”
- [x] `transport-point` 目前更适合作为：
  - `transport-scan` 的一个高层模式
  - 或后续 service 化时的 point endpoint
- [x] 未进入现有稳定白名单的研究型/分析型脚本，当前不纳入 launcher 首批范围

---

## 5. 拟交付物

- [x] launcher 设计说明
- [x] launcher 首批白名单表
- [x] launcher / wrapper / bootstrap / stable CLI 关系图或结构化说明
- [x] 对 D2 `create_app(...)` 评估所需的前置结论

### 5.1 D2 前置结论

- [x] 若评估 `PackageCompiler.create_app(...)`，目标不应是“打包整个仓库脚本集合”
- [x] 更合理的对象是：
  - 一个 launcher 顶层命令
  - 若干白名单子命令
  - 对 wrapper/bootstrap 语义的内嵌或替代实现
- [x] D2 需要重点回答：
  - app 内是否仍保留 sysimage mismatch policy 概念
  - app 打包后如何处理 fresh clone 与 release 资产分发边界
  - app 是否替代 bootstrap，还是仅替代 launcher + wrapper

---

## 6. 验收标准

- [x] 能明确回答“用户默认该怎么启动稳定能力”
- [x] 能明确回答“零知识 agent 默认该走哪个入口”
- [x] 能明确回答“哪些能力现在能进 launcher，哪些还不能”
- [x] 结论足以支撑 D2 / D3，而不需要重新回到 Phase C 重做边界定义

当前答案：

- [x] 用户默认入口：launcher
- [x] launcher 底层默认执行：`run_with_sysimage.*`
- [x] fresh clone / 新机器：先 bootstrap，再 launcher
- [x] 首批仅纳入现有稳定白名单 CLI，不扩大到研究/分析型脚本

---

## 7. 风险

- [x] 风险 R1：把 launcher 做成又一层脚本别名，没形成真正的职责边界
- [x] 风险 R2：把尚未稳定的脚本过早纳入白名单
- [x] 风险 R3：launcher 设计若脱离当前 wrapper / bootstrap 契约，会制造第二套入口语义

应对：

- [x] launcher 只新增“能力选择与统一错误语义”，不复制业务实现
- [x] 白名单严格继承现有稳定入口治理
- [x] wrapper / bootstrap 契约保持下沉复用，而不是在 launcher 中重写第二份

---

## 8. DoD

- [x] launcher 的职责边界已清晰写明
- [x] launcher 与 wrapper / bootstrap / stable CLI 的关系已固定
- [x] 首批 launcher 白名单已形成
- [x] 已能直接进入 D2：`PackageCompiler.create_app(...)` 适配性评估
