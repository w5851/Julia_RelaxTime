---
title: Issue #130：phase-reference adapter contract 与引用迁移
archived: true
original: docs/dev/active/2026-08-23_Issue130_phase-reference-adapter-contract.md
archived_date: 2026-08-30
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Issue #130：phase-reference adapter contract 与引用迁移

状态：accepted。PR253 的 solver-free compatibility audit、PR259 consumer migration 与 PR260
runtime switch 已完成并接受；legacy 数值现以 byte-preserving versioned snapshot 保留，未被删除。
本任务把 versioned candidate 的 schema 语义收束为显式 adapter 合同，并在独立 runtime-switch PR
中验证默认切换与回滚。本任务不改 phase 数值、solver 或容差；RS numerical production 由独立任务单承接。

## 已知边界

- candidate：`data/reference/pnjl/issue130_phase_reference_v1/`
- layer：`strict`、`derived`；`render` 仅可用于绘图，不可作为 runtime 输入。
- `mu_MeV`/`mu_CEP_proxy_MeV` 是 `mu_q`；adapter 同时派生 `muB=3*mu_q`，不得由调用方猜测。
- Maxwell 的 `grid_unresolved`、`geometry_converged`、`finite_and_converged` 和状态字符串
  共同决定 `certified`；unresolved 行必须保留为诊断，但默认不能进入 strict candidate
  runtime view。当前 runtime adapter 的逐键顺序为 strict certified、作者接受且在
  common support 内的 accepted fallback，最后才是 legacy snapshot。
- derived 的 `interpolated_noncertified` 只能显式允许用于派生图或诊断，不能伪装成 strict。
- CEP 同时保留 `[T_low_MeV,T_high_MeV]`、宽度和 `T_midpoint_MeV`；使用中点必须显式选择
  `cep_mode=estimated_midpoint`，内部审计可以选择 `cep_mode=bracket`。
- crossover 仅在 `physical_region=crossover_below_CEP` 时赋予 crossover 语义；Maxwell 区的
  导数峰不能被 adapter 自动暴露为 crossover。

## 任务项

- [x] 新增 solver-free `scripts/pnjl/phase_reference_adapter.py`，校验 import/layer manifest、
  CSV 字段、finite、重复键和 runtime_consumption gate。
- [x] 归一化 boundary/crossover/CEP/spinodal，明确 `mu_q`/`mu_B` 和 strict/derived 语义。
- [x] 提供只在内存中生成的 legacy view；默认拒绝 unresolved/non-certified 行，不写旧 reference。
- [x] 增加 fixture 单元测试，覆盖三态、CEP 双模式、单位换算、重复键、NaN/Inf 和 runtime gate。
- [x] 将 plotting consumer 改为显式 candidate root/layer 入口；candidate
  入口只生成诊断图，不改变 runtime 语义。
- [x] 为 `run_gap_transport_scan.jl` / phase-guided plan 定义 Julia 侧等价 adapter，并完成
  legacy-vs-candidate parity fixtures；runtime switch 仍由独立 PR 控制。
- [x] 审核 Paper P1 的 tagged-file contract，并提供显式 candidate overlay 与 parity fixture；
  candidate 仍只用于显式诊断/图层输入，不改变 P1 默认 legacy 路径。
- [x] 完成 transport、phase-guided、paper 和 legacy plot 的显式引用入口迁移及 solver-free
  回归；runtime resolver 已支持 certified-only candidate、逐键 legacy fallback 和显式 rollback。
- [x] versioned runtime-switch PR #260 默认使用 `strict` candidate 的 certified-only view，
  对缺失键保留 legacy fallback/rollback；作者随后接受限定 numerical pilot 并授权 canonical-path
  retirement。legacy 数值迁入 versioned snapshot，不在本任务中物理删除。

## 验收与非目标

本阶段验收是 adapter contract 可读、可拒绝不安全输入、消费者 parity 可追溯，并能以
strict candidate + 受控 accepted fallback + legacy fallback 运行或显式回滚；不是
phase-reference promotion，也不是 RS transport production。验证只运行 Julia/Python fixture、
candidate schema replay、CLI/dry-run、`git diff --check`、task-ledger/docs governance；本机不调用新的
equilibrium solver 扫描。runtime-switch 合并并经作者审核后才能标记 candidate 为默认
runtime reference；accepted fallback 的实际覆盖和旧 reference retirement 仍需独立审计。

## Contract revision (2026-08-29)

早期合同曾明确禁止 accepted 进入 runtime。该禁止已被新的活动任务
`2026-08-29_Issue130_phase-reference-accepted-runtime-fallback.md` 有限地修订：
accepted 不能 primary，但可作为作者接受、非外推、common-support 且状态有效的逐键
fallback。进入 runtime view 的 accepted 行保持 `certified=false` 并标记
`runtime_eligible=true`；直接 `allow_runtime=true, layer=accepted` 仍拒绝，避免把
下游派生图层伪装成 strict 数值 reference。
