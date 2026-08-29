# Issue #130：accepted phase-reference runtime fallback 合同修订

状态：active。该任务修订此前“accepted 不得进入 runtime”的明确禁止，原因是
legacy snapshot 的部分数值来自较旧、较稀疏的 reference，而当前 v2 `accepted`
层已经过作者物理审核并覆盖了大量 strict certificate 缺口。本任务只调整
adapter 的来源优先级和 provenance，不修改任何 phase 数值、equilibrium solver、
Maxwell、容差或 C0/C1/C2 artifact。

## 1. 当前合同

runtime view 的固定顺序为：

1. `strict` candidate 中通过现有 certificate predicate 的行；
2. v2 `accepted` 中 `author_accepted_for_downstream`、非外推、处于
   `native_support` 或 `interpolated_common_support` 的行；
3. versioned `legacy_phase_reference_v1` snapshot 的最后回退行。

`accepted` 行可以是 `interpolated_noncertified`。进入 runtime view 后必须保留
`certified=false`，同时显式标记 `runtime_eligible=true`、
`source_layer=accepted_fallback` 和原始 `source_status`。这表示“允许当前下游
消费”，不表示通过 strict 数值证书，也不等于 accepted 可以作为主 runtime
source。

直接以 `accepted` 或 `render` 作为 primary runtime source 仍然拒绝。显式
`--phase-reference-mode legacy` 仍是 rollback；传入自定义 candidate root 时，
除非同时传入 accepted root，否则保持旧的 strict→legacy 行为。

## 2. 为什么需要修订

legacy coverage audit v2 记录了 382 个历史 fallback key，而 accepted v2 与
strict candidate 的完整 key 集一致地覆盖了 strict 的 4,071 个未认证 boundary
行，并在相同 support 内覆盖了大多数 legacy 请求。继续把所有 strict 缺口直接
交给 legacy，会让更旧的 snapshot 优先于作者审核过的 accepted 非证书派生值。
修订后的顺序解决这个优先级问题，同时不把 accepted 的插值状态抹成 certified。

当前 solver-free 默认加载的诊断字段包括：

- `fallback_order=\"strict_candidate>accepted_downstream>legacy_snapshot\"`；
- `accepted_fallback_row_counts` 与 `legacy_fallback_row_counts`；
- `accepted_manifest_sha256`、candidate manifest SHA 和显式 rollback 状态；
- 每个 fallback 行的 `source_layer`、`accepted_source_status`、`certified` 与
  `runtime_eligible`。

## 3. 资格边界

- 未标记 `author_accepted_for_downstream` 的行不能进入 accepted fallback；
- `extrapolation=true`、support 外行、`unresolved`、`ambiguous` 或
  `not_converged` 状态不能进入 accepted fallback；
- crossover 只有 `physical_region=crossover_below_cep`（或空的历史兼容值）可
  进入；Maxwell/CEP/spinodal 保留各自字段和 bracket 语义；
- accepted 行与 strict/legacy 使用同一 adapter key 去重，不能按行号配对；
- accepted 资格失败或 key 仍缺失时才使用 legacy，并记录最终来源。

## 4. 实现与验证

- [x] Python adapter 支持显式 `accepted_bundle`/accepted fallback，保留旧
  strict→legacy API 和 direct accepted runtime rejection；
- [x] Julia adapter 支持默认 v2 accepted root、显式 accepted fallback 和
  legacy rollback；`run_gap_transport_scan.jl` sidecar 写入顺序、hash 和逐表计数；
- [x] Python/Julia fixture 覆盖 accepted 非证书行优先于 legacy，以及 strict/
  render/accepted runtime 边界；
- [ ] 增加 solver-free consumer smoke/replay v3，确认实际 active consumers 的
  fallback 矩阵与 accepted 资格统计；
- [ ] 作者审核 fallback 矩阵后，才决定是否可以继续 legacy path retirement；
  在此之前不创建物理删除 PR。

## 5. 非目标与回退

本任务不把 accepted 提升为 strict，不修改 `data/reference/**` 的数值，不删除
legacy snapshot，也不启动新的 PNJL/RS numerical run。若 accepted manifest、
support 或状态字段校验失败，adapter 回退到已验证的 strict→legacy 路径，并在
diagnostics 中保留错误；显式 legacy source 始终可用。
