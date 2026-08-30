---
title: Issue #130：accepted-primary runtime 与 PNJL legacy retirement
archived: true
original: docs/dev/active/2026-08-29_Issue130_phase-reference-accepted-primary-runtime.md
archived_date: 2026-08-30
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# Issue #130：accepted-primary runtime 与 PNJL legacy retirement

状态：completed。作者已明确 phase-guided transport 属于 phase-reference 的下游消费者，
`accepted` 是默认 runtime source；`strict` 仅在显式开启时作为 certified-only source。
旧 PNJL legacy snapshot 不再作为 runtime fallback 或 rollback，先完成 path retirement，
再以精确 allowlist 做物理删除。

## 1. 固定输入与目标

- accepted package：`data/reference/pnjl/issue130_phase_reference_v2/`；
- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`；
- source run：`32354095831`；solver-free replay：`32451053476`；
- 现有 strict/render/accepted 数据和旧 legacy snapshot 均保持 byte-stable；
- 不调用 PNJL solver，不改 Maxwell、CEP、equilibrium solver 或物理容差。

目标是让实际 phase-guided/gap transport runtime、phase-map plotting 和后续 consumer
共享同一个明确的 v2 accepted source，同时把严格证书和作者接受语义分开记录，避免较旧
legacy 数值因缺少 exact key 而被自动优先使用。

## 2. 三层与 runtime 合同

### strict

strict 是数值证据层。显式命令
`--phase-reference-mode strict --phase-reference-layer strict` 只暴露
`certified=true` 行；未认证行被省略，不自动补到 legacy。

### render

render 是结构化展示/诊断层，不允许作为 runtime 输入。

### accepted

accepted 是作者接受的下游主输入，包含完整 render 支持上的 native 和
`interpolated_noncertified` 行。每行仍保留 `certified=false`、`source_status`、
`coverage_status` 和 `extrapolation`，因此“可供当前下游消费”不等于“通过 strict 证书”。
loader 只接受 `author_accepted_for_downstream`、非外推、共同 support 内且不含
`unresolved/ambiguous/not_converged` 的行；通过后标记
`runtime_source_layer=accepted_primary`。

### legacy

`data/reference/pnjl/legacy_phase_reference_v1/` 只用于不可变历史审计和 Git 恢复边界。
正常 runtime 不加载它；`load_phase_reference_runtime_with_fallback`、
`--phase-reference-mode legacy` 和 `source=:legacy` 不再是可用的 runtime 路径。
历史审计脚本可以继续读取 snapshot，直到物理删除 PR 完成。

## 3. exact key 与“覆盖”语义

accepted 与 legacy 的网格并非逐键相同，但这不再构成自动 legacy fallback：

| 表 | legacy 行 | accepted 行 | legacy 未 exact 命中 | 同 ξ 最近 accepted 距离 |
|---|---:|---:|---:|---:|
| boundary/Maxwell | 48 | 12,537 | 34 | 1 MeV（温度） |
| crossover（CEP 以下） | 336 | 3,135 | 190 | 中位约 5.62 MeV、最大约 12.63 MeV（μq） |
| CEP | 11 | 161 | 0 | 0 |
| spinodals | 57 | 11,989 | 31 | 0.31–2 MeV（温度） |

这些距离来自 solver-free 邻近检查，不是新的数值结果。boundary/crossover consumer
按 accepted 曲线执行已有的温度/μ 插值和 ξ 最近层选择；不允许把邻近值伪装成 exact
legacy key。对超出 accepted support 的请求返回缺失/错误，由调用方显式处理。

## 4. 实现与验证状态

- [x] Python adapter：accepted-primary default、strict explicit、legacy fallback 拒绝；
- [x] Julia adapter：accepted-primary default、strict explicit、legacy fallback API 进入
  retired error；
- [x] phase-guided/gap transport CLI：默认 `accepted/runtime`，支持 `strict` 和
  `diagnostic`，不再接受 legacy mode；
- [x] solver-free v2 accepted loader 与实际 package smoke；
- [x] Python adapter tests：16 passed；Julia adapter：全部通过；run-gap/phase-guided
  focused tests：通过；
- [x] accepted-primary legacy coverage audit v3 已生成（历史基线）：
  `docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v3/`；
  accepted/strict 输入、legacy snapshot bytes/hash、exact-key 与同 ξ 邻近支持、
  runtime fallback=0 均已验证。其 verdict
  `accepted_primary_runtime_ready_path_retirement_pending` 仅表示当时仍有 21 个
  active legacy path contract；这不是 runtime fallback 依赖。
- [x] path-retirement implementation：6 个 active 默认路径契约已迁移到 accepted；旧
  CSV 仅保留显式 fixture/历史输入能力，snapshot 仍保留；solver-free audit v4 已重跑。
- [x] physical-deletion PR：PR #285 在 path-retirement CI、静态扫描、恢复 manifest 和作者
  授权后删除 snapshot，并 squash 合并为
  `dd943779a9ddc170491a61e69e2e41198532084d`。

## 5. 非目标与停止条件

- 不把 accepted 插值行改写为 strict certified；
- 不因 exact-key 缺失重跑 C0/C1/C2 或 PNJL solver；
- 不删除历史 `docs/analysis` evidence、source run 或 v2 candidate；
- 若 accepted manifest/hash/support 校验失败、active consumer 仍直接依赖 legacy，或
  发现未登记的 runtime 入口，停止 path retirement 和物理删除，先修复合同。

## 6. 可追溯证据

本任务的 accepted-primary 实现已在 PR #280 合并；本次路径迁移位于独立
`codex/issue130-pnjl-phase-reference-path-retirement` 分支，合并前不把未合并分支
写成 main provenance。迁移后的 solver-free evidence 为
`docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v4/`，其
`active_path_contract_count=0`、`physical_deletion_eligible=true` 曾表示可以提出
独立删除 PR；PR #285 已完成物理删除。
v1/v2 legacy coverage audit 保留为历史基线；v3 audit 仍作为历史基线，v4 audit 保存摘要、邻近覆盖、
consumer matrix、manifest 和 claim ledger，不复制原始全量曲线。legacy 现在只能从 Git
recovery ref 在临时分支恢复，不能作为 runtime fallback/rollback。
