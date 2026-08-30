# Issue #130：phase reference 三层语义与 accepted 下游数据层

状态：accepted；作者已确认 `accepted` 作为所有 phase-reference 下游（含 runtime）默认层，
当前转入 PNJL legacy path-retirement。三层数据仍按 `strict → render → accepted` 的
构建/展示关系保留，但这不是 runtime 优先级；strict 仅显式开启，legacy 不再 fallback/rollback。
现有 v1 package 和历史 evidence 保持不可变；`derived` 降为 render 的内部构建
输入与 provenance，不再作为下游公共语义层。

## 1. 背景与目标

当前 `data/reference/pnjl/issue130_phase_reference_v1/` 同时保存 strict、derived
和 render 三个历史层。作者已经审核当前 phase-surface 图像，认为其整体物理形状
和覆盖效果足以支持下游研究；但 render package 仍缺少实际 spinodal 坐标表，且
当前 derived 的插值/未认证状态没有一个明确的下游接受语义。

本任务的目标是建立一个不重跑 PNJL solver 的、可追溯的三层 v2 package：

1. `strict` 保留原生 C2/v6 与真实 Maxwell 补点及 unresolved 语义；
2. `render` 形成完整结构化展示数据（含 spinodal）和已审核图像；
3. `accepted` 从结构化 render 数据生成，供允许派生/插值的下游研究消费；
4. 原 derived 构建输入继续作为内部 provenance 保存，但不再冒充第四个公共层。

`accepted` 的“接受”表示在当前研究用途、覆盖范围和作者物理审核下可使用，
不表示所有行通过 strict 数值证书。

## 2. 范围与非目标

### 范围

- 建立 versioned v2 schema、manifest、README、claim ledger 和 review record；
- 从现有 v1 strict/derived/render 做 solver-free render 完整化；
- 将 derived 的 spinodal 结构表纳入 render，保持值、单位、key 和 hash 可追溯；
- 生成 accepted 下游表，保留 native/unresolved/interpolated 状态、来源端点和
  coverage；
- 为 Python/Julia adapter 增加显式 accepted 读取和 consumer smoke contract；
- 更新 Issue #130 任务台账和 legacy fallback 审计依赖关系。

### 非目标

- 不调用 PNJL equilibrium solver，不重跑 C0/C1/C2、CEP、Maxwell 或 RS numerical；
- 不修改 Maxwell、crossover、CEP、spinodal 数值或任何物理容差；
- 不把 interpolated/unresolved 行改写为 strict certified；
- 不从 PNG 反向提取数据，不用 render 图像替代结构化数据；
- 不在本任务中物理删除 PNJL legacy snapshot；runtime 默认切换与物理清理由后续独立任务承接；
- 不修改现有 v1 package、v1 figure 或历史 manifest。

## 3. 三层数据合同

### 3.1 strict

`strict` 是数值证据层。它保留原生 C2/v6 行、276 条真实 Maxwell 补点和已有
unresolved 状态。`certificate_status` 只能反映既有严格门禁，不因作者审核图像而
改变。

### 3.2 render

`render` 是完整结构化展示层，而不只是 PNG：

- crossover、Maxwell/first-order、CEP bracket/boundary、spinodal；
- ordered xi-slice 数据和 coverage mask；
- 图像、plot manifest 与作者审核 provenance。

现有 render 缺失的 spinodal 坐标必须从已有 derived 表补齐；不得跨缺失 support
做无界外推。

### 3.3 accepted

`accepted` 是 render 之后的下游派生层。每行至少保留：

- `source_status`：`strict_certified`、`unresolved` 或 `interpolated_noncertified`；
- `acceptance_status`：`author_accepted_for_downstream` 或待审核状态；
- `source_layer`、来源 xi/axis 端点、插值方法和 `coverage_status`；
- `extrapolation=false` 以及 CEP/phase-exclusive 规则。

accepted 可以满足完整相结构下游查询，但不得被表述为 strict convergence evidence。

## 4. 可勾选任务分解

### 阶段 A：合同与冻结

- [x] 从最新 `origin/main` 创建隔离 branch，记录 repo SHA、现有 v1 tree/hash 和
  legacy-retirement active 依赖。
- [x] 定义 v2 package/schema、层字段、状态机、runtime/consumer 边界和回退语义。
- [x] 为所有历史 v1 文件建立 byte/hash unchanged 检查。

### 阶段 B：render 完整化

- [x] 新增 solver-free materializer，把现有三类 render 表与 derived spinodal 合并
  为完整 render 表。
- [x] 验证 key、单位、行数、NaN/Inf、重复键、coverage 和来源 hash；不修改 strict。
- [x] 生成 render manifest、plot manifest、claim ledger 和完整性报告。

### 阶段 C：accepted 派生层

- [x] 从结构化 render 表生成 accepted 四类表，不从图像提取。
- [x] 保存 native/unresolved/interpolated 状态、来源端点、局部插值和 phase-exclusive
  规则；禁止无界外推。
- [x] 生成 accepted manifest、coverage/fallback 对照表和作者审核记录入口。

### 阶段 D：adapter 与 consumer smoke

- [x] Python/Julia adapter 支持显式 `accepted` layer；render/accepted 仍拒绝作为
  runtime input。
- [x] 对 phase-map/Paper P1 display path 做 solver-free accepted consumer smoke；
  phase-guided/RS 的 runtime smoke 保持阻断，直到 accepted 通过作者审核并另立
  runtime-switch 变更。
- [x] 保留 strict candidate 的历史 provenance；旧 legacy fallback/rollback 合同已由
  accepted-primary 任务取代，不切换 solver 数值路径。

### 阶段 E：作者审核与后续 retirement

- [x] 审核完整 spinodal、accepted 覆盖、插值状态和代表图。
- [x] 将 accepted 标为 `accepted_for_downstream`，并设为下游分析/派生相图默认层。
- [x] 重新运行 legacy fallback coverage audit；只有 fallback=0 且恢复证据完整时，才另立
  PNJL legacy physical-deletion PR。

`accepted` 的默认范围仅为 phase-map、Paper P1 和其他明确允许派生/插值的分析消费者。
Julia transport/runtime 当前使用 accepted-primary view；strict certified-only view 仅在
显式模式开启，legacy snapshot 只供 retirement audit/history，promotion 不改变 solver 数值。

审计 v2 已生成于
`docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v2/`，
verdict 为 `retirement_inconclusive`：legacy fallback `382`、active
consumer/rollback blocker `25`、unknown active reference `0`。因此 accepted
默认已经完成下游语义切换，但 strict runtime 的 consumer migration 和 legacy
fallback coverage 仍未完成，不能进入物理删除。

## 5. 测试与验收标准

- 所有 materializer/replay 必须 `solver_called=false`，且不写入现有 v1 文件；
- strict/render/accepted 的 schema、单位、key、NaN/Inf、重复键和状态语义通过 focused tests；
- render 必须包含四类结构化表及 coverage，accepted 必须保留来源和插值状态；
- 无超出共同 support 的外推；CEP 与 crossover/first-order 排他规则保持有效；
- v1 文件 SHA、legacy snapshot SHA 和现有 figure SHA 不变；
- adapter 的 render runtime rejection、accepted explicit opt-in、strict fallback/rollback 均通过；
- 文档、task ledger、data-path 和 active-doc governance 通过。

## 6. 风险与回退

- 若 spinodal 与 render key 无法无歧义配对，停止并保留 v1；
- 若 accepted 需要跨 support 外推，停止，不用插值掩盖缺口；
- 若 consumer 依赖 strict certification 语义，accepted 只能作为显式研究视图；
- 任一阶段失败均停止并保留现有 accepted/strict evidence，不删除数据、不调用 solver；legacy
  snapshot 只作为物理删除前的历史恢复边界。

## 7. DoD

本任务完成的最低条件是：v2 三层 package、完整 render、accepted 数据、manifest/claim
ledger、adapter contract 和 solver-free consumer smoke 均有可追溯证据，并经作者审核。
它不自动等价于 strict 数值 convergence，也不自动完成 PNJL legacy 物理删除。

## 8. 本次执行证据（2026-08-29）

- 隔离分支：`codex/issue130-phase-reference-three-layer-v2`；基线为
  `origin/main@ead55ece4a799948e54799288177f45908ff9d88`。
- source manifest：`data/reference/pnjl/issue130_phase_reference_v1/manifest.json`，
  SHA-256 `7d37aedd8ec7a1c4715f5f0b9462460ed5619d26586b27e5326e32a777b06519`。
- v2 输出：`data/reference/pnjl/issue130_phase_reference_v2/`；root/layer manifest
  均声明 `solver_called=false`、`runtime_consumption=false`、`reference_write=false`。
- 行数：strict=`boundary 7162 / crossover 1343 / cep 93 / spinodals 6886`；
  render=`12537 / 3135 / 161 / 11989`；accepted 与 render 相同。
- accepted 状态计数：boundary=`strict_certified 7162 + interpolated_noncertified 5375`；
  crossover=`1343 + 1792`；cep=`93 + 68`；spinodals=`6886 + 5103`。状态未被升级，
  且不包含无界外推。
- focused evidence：Python adapter/materializer/P1 smoke `12 passed`；相关
  runtime/layer tests `22 passed`；Julia adapter `45 passed`（含 v2 三层 fixture），
  phase-guided plan `44 passed`。本地未调用 PNJL equilibrium solver。
- 兼容性修正：import candidate 的 legacy byte snapshot 已改为读取
  `data/reference/pnjl/legacy_phase_reference_v1/`，与已合并的物理 retirement
  路径一致；原顶层 dense 路径继续保持不存在。
- accepted 已完成作者审核并记录为 `accepted_for_downstream`；行级插值状态仍为
  `interpolated_noncertified`，strict certificate 未被升级。旧的“strict primary +
  accepted/legacy fallback”文字仅是历史快照；当前 runtime 已改为 accepted-primary，
  strict 仅显式启用，legacy 只保留给 retirement audit，不能据此自动回退或 rollback。

## 9. accepted promotion 与 legacy audit v2（2026-08-29）

- promotion 脚本：`scripts/analysis/pnjl/promote_issue130_phase_reference_accepted.py`；
  只更新 accepted 行级决策字段、layer/root manifest、claim ledger 和说明文字。
- v2 accepted package：`data/reference/pnjl/issue130_phase_reference_v2/`；
  `promotion_status=accepted_for_downstream`、`downstream_default_layer=accepted`、
  `runtime_consumption=false`、`reference_write=false`、`solver_called=false`。
- strict 四表、render 四表、源 manifest 和 legacy snapshot 均未由 promotion 脚本写入；
  accepted 的 12,537/3,135/161/11,989 行分别保留原 source status 和非外推约束。
- solver-free legacy coverage audit v2 已生成于
  `docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v2/`；
  其 `retirement_inconclusive` 结果只进入 consumer/path retirement；accepted 插值行不
  自动提升为 certified。accepted-primary runtime 合同、strict 显式模式和 legacy
  path/physical retirement 由 `2026-08-29_Issue130_phase-reference-accepted-primary-runtime.md`
  承接；v2 audit 数字仍作为历史基线，不覆盖新合同。
