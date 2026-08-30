# Issue #130：PNJL legacy phase-reference retirement 与物理清理

状态：active；这是 RS `prod_v1` 物理删除（PR #278）之后的独立 PNJL
phase-reference follow-up。作者已确认 v2 `accepted` 是 phase-reference 下游的默认
runtime source，strict 只在显式模式下使用；legacy 不再作为 runtime fallback/rollback。
当前已完成 path retirement，并在独立分支准备 PNJL legacy snapshot 物理删除 PR；
该 PR 尚待 CI 与作者合并授权。

## 1. 背景与目标

Issue #130 已将
`data/reference/pnjl/issue130_phase_reference_v2/accepted` 设为 runtime 的默认
下游输入，并保留 v2 strict 作为显式 certified-only 选项。历史
`data/reference/pnjl/legacy_phase_reference_v1/` 不再承担逐键 fallback 或 rollback，
但在 path retirement 完成前仍作为不可变历史审计输入保留。PR #264 只完成了旧
canonical 根路径退役，并不等于 snapshot 可以立即物理删除。

本任务的目标是按顺序完成：

1. 证明 candidate strict、adapter 和所有实际消费者已经不再需要 legacy fallback；
2. 在独立 PR 中完成 fallback/path retirement（如审计发现仍有路径需要迁移）；
3. 只有在审计证据和作者授权均满足后，另开精确 allowlist 的物理删除 PR，并保留
   可验证的 Git 恢复边界。

## 2. 范围与非目标

### 范围

- 审计 candidate strict 的 certified/uncertified 行、legacy 的逐表 key 覆盖和
  实际 fallback 计数；
- 盘点 Julia/Python/workflow 的 active consumer、显式 legacy rollback 和历史
  evidence 引用；
- 固化 candidate/legacy manifest 与文件 SHA、schema、单位、缺失/失败语义；
- 生成不可变的 retirement evidence package，并据此决定是否可以创建后续 PR。

### 非目标

- 不调用 PNJL equilibrium solver，不重跑 C0/C1/C2、CEP、RS numerical 或 transport；
- 不修改 Maxwell、crossover、CEP、spinodal 数值或容差；当前 adapter 合同为
  accepted primary、strict 显式 opt-in，legacy 仅允许历史审计/恢复边界读取；
- 不删除 `data/reference/pnjl/crossover.csv`（它是另一组历史 fixed-point 输入）；
- 不删除 `docs/analysis/` 历史 evidence、raw rho-mu 外部归档指针或当前 candidate；
- 不把 candidate 的 unresolved/non-certified 行强行标为 certified，也不以图像完整性
  代替 runtime key 覆盖证明。

## 3. 当前能力与缺口（2026-08-28 solver-free 基线）

权威输入为：

- candidate：`data/reference/pnjl/issue130_phase_reference_v1/`，calculation SHA
  `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`；
- legacy：物理删除前的 `data/reference/pnjl/legacy_phase_reference_v1/` snapshot，
  其来源、字节数和 SHA-256 由
  `docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v4/tables/legacy_snapshot_inventory.csv`
  及恢复引用固定保存；合并后工作树中不再保留该目录。
- adapter：`scripts/pnjl/phase_reference_adapter.py` 与
  `scripts/relaxtime/phase_reference_adapter.jl`。

按 adapter 的 certified predicates 对当前 CSV 做的只读预检为：

| 表 | candidate 行 / certified / uncertified | legacy 行 | legacy 中未被 certified candidate 精确 key 覆盖 |
|---|---:|---:|---:|
| boundary/Maxwell | 7162 / 3091 / 4071 | 48 | 35 |
| crossover | 1343 / 1343 / 0 | 336 | 315 |
| CEP | 93 / 90 / 3 | 11 | 1 |
| spinodals | 6886 / 6886 / 0 | 57 | 31 |

这些数字是审计起点而不是最终 gate；它们已经足以说明当前不能物理删除 snapshot：
runtime view 仍可能从 legacy 提供缺键或未认证键。尤其 candidate 与 legacy 的网格
并不相同，不能只用行数或文件是否存在判断“已覆盖”。

## 4. 可勾选任务分解

### 阶段 A：solver-free retirement audit

- [x] 固定 candidate/legacy 文件树、manifest、calculation/postprocess SHA 和
  tree hash；检查 NaN/Inf、重复键、schema/单位及 snapshot 字节完整性。
- [x] 以 adapter 的真实 key 语义（而非 CSV 行号）重建每表 certified、uncertified、
  fallback、overlap 和 consumer-requested coverage 矩阵。
- [x] 静态扫描 Julia、Python、workflow、配置和文档，区分 runtime consumer、显式
  `legacy` rollback、历史复现和已失效路径；检查 fallback 是否仍可达。
- [x] 输出 versioned evidence package：`manifest.json`、`coverage.csv`、
  `fallback_matrix.csv`、`consumer_matrix.csv`、`claim_ledger.json` 和停止原因；
  package 不复制原始全量曲线。

#### 阶段 A 实际结果（2026-08-28）

审计包位于
`docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v1/`，
固定 candidate calculation SHA 为
`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`，source run 为 `32354095831`，
replay run 为 `32451053476`。审计本身 `solver_called=false`、
`reference_write=false`，candidate/legacy tree 与
`RETIREMENT_MANIFEST.json` 的 bytes/SHA 均通过。

当前 verdict 为 `retirement_inconclusive`：legacy fallback key 共 `382`（
boundary `35`、crossover `315`、CEP `1`、spinodals `31`），candidate
uncertified 行 `4074`；静态 tracked-source 扫描记录 `24` 个仍需迁移或保留
rollback 的 active contract，未知 active 引用为 `0`。因此 snapshot 继续保留，
不得创建物理删除 PR。candidate-only 迁移不能把这些 unresolved 行自动提升为
certified；阶段 B 必须先定义请求键/适配层合同并验证显式 legacy rollback。

#### 阶段 A post-acceptance 重审（2026-08-29）

作者已确认 v2 package 的 `accepted` 层作为 phase-map/分析下游默认，但这不改变
strict runtime 的 certified-key 合同。基于同一 candidate calculation SHA
`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`、source run `32354095831` 和
solver-free replay `32451053476`，重新生成的审计包位于
`docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v2/`。

v2 结果为 `retirement_inconclusive`：legacy fallback key `382`
（boundary `35`、crossover `315`、CEP `1`、spinodals `31`），candidate
uncertified row `4074`，active consumer/rollback blocker `25`，unknown active
reference `0`。新增的旧审计契约测试已明确归类为 snapshot contract test；它仍需在
物理删除前迁移或随删除 allowlist 一并处置，但不再构成未登记引用。`accepted`
插值行仍为 non-certified，故不能用来宣称 strict fallback 已消失。

v2 审计保持 `solver_called=false`、`reference_write=false`、`runtime_consumption=false`，
不修改 candidate/legacy tree，也不创建物理删除 PR。下一步仍是迁移 strict runtime
consumer 的请求键/适配层，并审计 accepted fallback 的资格与实际来源；只有新的
fallback 矩阵确认无 active legacy 依赖、消费者/恢复边界完整后，才可提出独立
physical-deletion PR。旧的 strict→legacy audit 数字继续保留为历史基线。

#### 阶段 A.1：accepted fallback 合同修订（2026-08-29，历史 superseded）

作者指出 strict certificate 缺口直接落到 legacy 会优先使用较旧 snapshot，
而 v2 `accepted` 中的作者审核非证书派生值在当前研究用途下更可信。该历史任务曾
提出 `strict candidate → accepted_downstream → legacy snapshot` 的逐键 fallback，
但随后被作者的 accepted-primary 决策取代：accepted 现在直接作为下游 runtime
primary，strict 只在显式模式下启用，legacy 不再作为 fallback/rollback。原任务单、
旧 audit 数字和 snapshot 仍保留作历史 provenance，不再指导当前实现。

### 阶段 B：fallback/path retirement（以阶段 A 与 accepted-primary 证据为输入）

- [x] 将 6 个活动默认路径契约迁移到 v2 `accepted`：Julia gap/phase-equilibrium、
  Python 相图/验证器、ξ 平滑采样和 CEP baseline exporter；旧 CSV 仅保留显式
  fixture/历史输入能力。
- [x] 为 accepted schema 增加必要的 production-parity 字段映射、路径参数和 focused
  unit coverage；不改变数值 solver、Maxwell 或物理容差。
- [x] 在独立 path-retirement 分支生成 solver-free audit v4，验证活动 legacy 路径为
  `0`、runtime fallback/rollback 为 `0`、accepted/strict/snapshot provenance 完整。
  此阶段仍不物理删除 snapshot。

### 阶段 C：物理清理（再次单独授权）

- [x] 生成精确 deletion allowlist、删除前 tree/hash/bytes manifest 和 Git restore
  指令；allowlist 只包含已确认不再被任何 consumer 使用的 PNJL legacy snapshot。
- [x] 在物理删除分支应用 allowlist，运行 deletion validator、data-path/active-doc/task-ledger governance 和
  focused adapter tests；确认 candidate、历史 evidence 与 `crossover.csv` 不受影响。
- [x] 创建独立 physical-deletion PR；未经作者再次授权不合并、不删除本地或远端
  其他分支/文件。

#### 阶段 C 实际结果（2026-08-30，待合并）

物理删除提案位于
`docs/analysis/pnjl/phase_reference/issue130_phase_reference_physical_deletion_v1/`，
allowlist 精确包含 8 个 legacy snapshot 文件，共 50,749 bytes。删除前 tree hash 为
`c128ef6358a5813533fc5a9726047585a750a4421f162293564bea8e363764e6`，恢复引用固定为
`9aa4c313901ca0c91e851f58514e3df9aa124df4`。分支中的 validator 报告
`physical_deletion_proposal_valid`，`solver_called=false`、`production_write=false`，
candidate/accepted/strict、历史 evidence 和独立 `crossover.csv` 均保持存在。
合并后 legacy 不再作为 runtime fallback/rollback；如需复原只能从 Git 恢复引用在
临时分支恢复，不承诺运行时回退。

## 5. 测试与验收标准

- 阶段 A 必须 `solver_called=false`、`reference_write=false`，双 SHA/hash 与
  `RETIREMENT_MANIFEST.json` 一致；无 NaN/Inf、重复键或路径漂移。
- consumer matrix 必须逐项给出 active/fallback/rollback/history 分类；任何未解释
  的 fallback key、未知消费者或 hash mismatch 都使 verdict 为
  `retirement_inconclusive`，不得进入物理删除。
- 阶段 B 只有在 runtime candidate-only smoke、legacy rollback 和 schema/units
  parity 全部通过后才可接受。
- 阶段 C 的 deletion PR 必须证明 allowlist 外零变更，并保留可达的 Git 恢复引用；
  合并后不再声称 PNJL legacy fallback 可用。

## 6. 里程碑与停止条件

1. `legacy_audit_v1`：只读覆盖/消费者审计完成并由作者审核；当前结果已归档为
   `retirement_inconclusive`，故保持现状并进入阶段 B 设计。
2. `legacy_path_retirement_v1`：若需要，迁移 active consumer 并保留 rollback。
3. `legacy_physical_deletion_v1`：仅在 fallback=0、consumer migration 完成、
   evidence/hash 完整且作者明确授权后执行。

以下任一情况立即停止：candidate 仍有 active key 缺口、fallback 仍被运行时使用、
snapshot hash 不一致、发现未登记 consumer、恢复引用不可达，或任何数值/solver 变化。

## 7. DoD、风险与回退

DoD 是阶段性而非自动连跳：每一阶段都有独立 evidence、PR 和作者审核；在阶段 C
合并前 legacy snapshot 必须保持可读、可复现，并作为 Git 恢复边界保留。风险主要是
candidate 的细网格与旧粗网格不共键、隐藏 consumer 依赖以及误把诊断层 unresolved
当作可运行数据；缓解方式是按真实 adapter key 审计、静态引用矩阵和精确 allowlist。
任何阶段失败都保留当前 snapshot 和历史 evidence，不重算数值、不修改物理容差。

## 8. 当前 accepted-primary 重审结果（2026-08-29）

作者已明确接受 `accepted` 作为所有 phase-reference 下游（包括 phase-guided transport）
的默认输入，`strict` 只在显式模式下使用。当前 adapter 已移除可用的
`strict→accepted→legacy` fallback/rollback 入口；legacy 仅由历史审计和恢复边界读取。

solver-free v3 审计包（历史基线）：
`docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v3/`。
审计固定 calculation SHA `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`，验证 accepted
manifest promotion、strict certified-only、legacy snapshot bytes/SHA、同 ξ 邻近支持以及
运行时 `legacy_fallback_rows=0`。legacy 与 accepted 的 exact key 仍不是同一网格，
因此 coverage 表仅用于解释邻近支持，不能把 nearest row 伪装成 exact row。

当前 v3 verdict `accepted_primary_runtime_ready_path_retirement_pending` 仅为历史基线。
路径迁移后的 solver-free v4 审计包位于
`docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v4/`，其
`active_path_contract_count=0`、`runtime_legacy_fallback_rows=0`、
`runtime_legacy_rollback_enabled=false`，当前 verdict 为
`legacy_physical_deletion_candidate`。这只表示可以提出独立 physical-deletion PR；删除
PR 仍只允许精确 allowlist、恢复 manifest 和作者授权，且保留全部历史 evidence。
