# Issue #130：phase reference 三层语义与 accepted 下游数据层

状态：active；本任务承接作者已确认的三层方案：`strict → render → accepted`。
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
- 不在本任务中物理删除 PNJL legacy snapshot，不自动改变 strict runtime 默认；
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
- [x] 保留 strict candidate、legacy fallback 和显式 rollback；不自动切换 solver runtime。

### 阶段 E：作者审核与后续 retirement

- [ ] 审核完整 spinodal、accepted 覆盖、插值状态和代表图。
- [ ] 审核通过后将 accepted 标为 `accepted_for_downstream`，并决定是否作为下游分析默认。
- [ ] 重新运行 legacy fallback coverage audit；只有 fallback=0 且恢复证据完整时，才另立
  PNJL legacy physical-deletion PR。

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
- 任一阶段失败均回退到现有 strict candidate + legacy fallback，不删除数据、不调用 solver。

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
- 当前仍待作者审核 accepted 覆盖与语义；在此之前不切换 runtime 默认、不删除
  legacy、不启动新的数值重验。
