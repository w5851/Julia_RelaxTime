# Issue #130：accepted phase-reference runtime fallback 合同修订

状态：superseded（2026-08-29）。该任务曾修订此前“accepted 不得进入 runtime”的明确禁止，原因是
legacy snapshot 的部分数值来自较旧、较稀疏的 reference，而当前 v2 `accepted`
层已经过作者物理审核并覆盖了大量 strict certificate 缺口。本任务只调整
adapter 的来源优先级和 provenance，不修改任何 phase 数值、equilibrium solver、
Maxwell、容差或 C0/C1/C2 artifact。当前合同已进一步收束为 accepted-primary；
strict 仅显式启用，legacy 不再作为 runtime fallback/rollback。当前执行入口是
`2026-08-29_Issue130_phase-reference-accepted-primary-runtime.md`。

## 1. 已废止的 v1 合同（历史记录）

PR #280 之前曾拟定的 runtime view 顺序为：

1. `strict` candidate 中通过现有 certificate predicate 的行；
2. v2 `accepted` 中 `author_accepted_for_downstream`、非外推、处于
   `native_support` 或 `interpolated_common_support` 的行；
3. versioned `legacy_phase_reference_v1` snapshot 的最后回退行。

当时 `accepted` 行可以是 `interpolated_noncertified`。进入 runtime view 后必须保留
`certified=false`，同时显式标记 `runtime_eligible=true`、
`source_layer=accepted_fallback` 和原始 `source_status`。这表示“允许当前下游
消费”，不表示通过 strict 数值证书，也不等于 accepted 可以作为主 runtime
source。

当时直接以 `accepted` 或 `render` 作为 primary runtime source 仍然拒绝。显式
`--phase-reference-mode legacy` 仍是 rollback；传入自定义 candidate root 时，
除非同时传入 accepted root，否则保持旧的 strict→legacy 行为。

## 2. 为什么当时需要修订

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
- crossover 只有 `physical_region=crossover_below_cep` 可进入；历史空字段只在
  同行存在有限 `mu_CEP_proxy_MeV` 且 `muq_MeV <= mu_CEP_proxy_MeV` 时兼容，
  以确保 CEP 以上的导数峰不被赋予 crossover 物理含义。Maxwell/CEP/spinodal
  保留各自字段和 bracket 语义；
- accepted 行与 strict/legacy 使用同一 adapter key 去重，不能按行号配对；
- accepted 资格失败或 key 仍缺失时才使用 legacy，并记录最终来源。

## 4. 实现与验证

- [x] Python adapter 支持显式 `accepted_bundle`/accepted fallback，保留旧
  strict→legacy API 和 direct accepted runtime rejection；
- [x] Julia adapter 支持默认 v2 accepted root、显式 accepted fallback 和
  legacy rollback；`run_gap_transport_scan.jl` sidecar 写入顺序、hash 和逐表计数；
- [x] Python/Julia fixture 覆盖 accepted 非证书行优先于 legacy，以及 strict/
  render/accepted runtime 边界；
- [x] fixture 覆盖 CEP 以上历史 crossover 导数峰：该类行保留在 legacy snapshot，
  但从 runtime crossover fallback 中排除并单独计数；
- [ ] 增加 solver-free consumer smoke/replay v3，确认实际 active consumers 的
  fallback 矩阵与 accepted 资格统计；
- [ ] 作者审核 fallback 矩阵后，才决定是否可以继续 legacy path retirement；
  在此之前不创建物理删除 PR。

## 5. 当前 solver-free 默认 materialization

在当前 candidate/accepted/legacy 文件树上只读加载得到的来源计数为：

| 表 | strict candidate | accepted fallback | legacy retained | legacy excluded |
|---|---:|---:|---:|---:|
| boundary/Maxwell | 3,091 | 9,446 | 34 | 0 |
| crossover | 1,343 | 1,792 | 190 | 125 |
| CEP | 90 | 71 | 0 | 0 |
| spinodals | 6,886 | 5,103 | 31 | 0 |

`legacy retained` 是当前 runtime view 实际仍会读取的逐键 fallback；
`legacy excluded=125` 是历史 crossover 表中 `muq_MeV > mu_CEP_proxy_MeV` 的
导数峰，它们不是待补齐的物理 crossover 行。原始 snapshot 仍保持
`boundary 48 / crossover 336 / CEP 11 / spinodals 57` 的字节级完整性。

这里要区分“参数范围覆盖”和“runtime key 覆盖”：accepted package 在其声明的共同
support 内提供派生网格，但 adapter 当前按各表的 exact key 去重（boundary/spinodals
为 `(xi,T_MeV)`、crossover 为 `(xi,muq_MeV)`、CEP 为 `(xi)`），不会在 runtime
入口隐式插值。因此 accepted 目前不能逐键完全替代 legacy；它已经优先接管了
strict 的大量 non-certified gap，但剩余的 34/190/0/31 行仍需 legacy 或另一个
明确的插值适配合同。

## 6. 非目标与回退

本任务不把 accepted 提升为 strict，不修改 `data/reference/**` 的数值，不删除
legacy snapshot，也不启动新的 PNJL/RS numerical run。若 accepted manifest、
support 或状态字段校验失败，adapter 回退到已验证的 strict→legacy 路径，并在
diagnostics 中保留错误；显式 legacy source 始终可用。

## 7. 历史保留说明

本文件中的 strict→accepted→legacy 表格、fallback 计数和 rollback 语义是
superseded 的历史合同，仅用于解释 PR #280 之前的决策。不得据此配置当前 runtime；
legacy snapshot 的读取只允许留在 retirement audit 或 byte-preserving 历史恢复工具中。

## 8. 当前合同指针

当前不再启用本文件第 1 节的 fallback 顺序。请使用
`2026-08-29_Issue130_phase-reference-accepted-primary-runtime.md`：
accepted 是默认 runtime source，strict 是显式 certified-only 选项，legacy 既不
fallback 也不 rollback；v3 覆盖审计位于
`docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v3/`。
