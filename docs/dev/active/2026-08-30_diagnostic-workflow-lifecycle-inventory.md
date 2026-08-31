# Diagnostic workflow lifecycle inventory v1

创建日期：2026-08-30

状态：`wave2_author_review`

这是 Issue #130/RS 历史任务文档归档后的第二阶段任务单。归档阶段已经在
`main@cea8110cdc1e1bafbc9873aa79cb5b0954c3de76` 完成；本 PR 基于
`main@6b2f4ffeb05001f6e39b6940d59e4cbc9cf9bfed`，中间提交没有修改 workflow 文件。本任务只盘点 workflow 生命周期，
不执行 workflow 删除、停用、改触发条件或 artifact 迁移。

## 1. 本阶段目标

- 固定当前仓库内 47 个 `.github/workflows/*.yml` 的静态 inventory。
- 记录一次性诊断 workflow 的最近 run、artifact 数量、过期状态和可追溯入口。
- 绘制 workflow-to-workflow 调用边和仓库内活动引用边。
- 将 workflow 分为 `retain`、`preserve_historical_replay`、`retain_archive_recovery` 和
  `retirement_candidate_review`，形成作者可审核的精确候选清单。
- 在独立 retirement PR 之前冻结保留期、provenance 迁移和删除/归档边界。

详细静态结果见
`docs/analysis/governance/diagnostic_workflow_inventory_v1/`。

## 2. 明确非目标

- 不删除或禁用任何 workflow，不修改 `.github/workflows/` 文件。
- 不删除 GitHub Actions run/artifact，不假定当前 retention 等于永久保存。
- 不删除或迁移诊断脚本、测试、分析表和已归档 Issue #130 文档。
- 不重跑 PNJL、Maxwell、CEP、RS 或任何 solver；本阶段 `solver_called=false`。
- 不把“近期没有 standalone run”解释为“没有下游依赖”；`workflow_call`、脚本和文档引用
  分别核对。

## 3. 快照边界

| 项目 | 固定值 |
| --- | --- |
| repository | `w5851/Julia_RelaxTime` |
| source branch | `main` |
| source tree SHA | `cea8110cdc1e1bafbc9873aa79cb5b0954c3de76` |
| PR base SHA | `6b2f4ffeb05001f6e39b6940d59e4cbc9cf9bfed`（中间无 workflow 文件变化） |
| inventory date | `2026-08-30` |
| local workflow files | 47 |
| one-shot diagnostic run rows | 20 |
| solver called by inventory | `false` |

run/artifact 数量是 2026-08-30 的 GitHub API 观察快照，不是未来保留承诺。artifact 可能
在 GitHub retention 到期后消失；需要长期复现的数值证据必须依赖仓库内或外部不可变 evidence
包及其 hash，而不是依赖 workflow 定义本身。

## 4. 分类结果

### 4.1 立即保留

- CI、治理、模型入口、依赖、数据路径、benchmark、nightly 和 precompile workflow。
- Dense reference 的主 workflow、resume/replay 入口及其 `workflow_call` 内部 shards。
- raw rho-mu archive、Zenodo draft/restore 及 archive shard；它们承担恢复和外部归档职责。
- relaxtime manual/meson/phase-guided production、sysimage 等仍可能被使用的生产入口。
- `pnjl-phase-diagram.yml` 虽然最近 standalone run 较早，但仍被 legacy audit 脚本引用，暂不
  作为删除候选。

### 4.2 盘点时的保留定义和证据

以下 workflow 有已完成的 Issue #130/CEP/Maxwell 诊断证据，或仍有测试/活动任务引用。即便
未来归档定义，也必须先把可复现入口、输入 hash、run ID、artifact provenance 固化到 analysis 包：

- `pnjl-c2-cep-limited-feasibility.yml`
- `pnjl-c2-cep-three-midpoint.yml`
- `pnjl-c2-cep-xi05-high-side-extension.yml`
- `pnjl-c2-limited-feasibility.yml`
- `pnjl-c2-targeted-closure-v1.yml`
- `pnjl-cep-deep-oracle.yml`
- `pnjl-cep-hybrid-production-shadow.yml`
- `pnjl-cep-narrow-pilot.yml`
- `pnjl-cep-production-shadow.yml`
- `pnjl-issue130-crossover-mu-endpoint-pilot.yml`
- `pnjl-issue130-maxwell-cep-local-expansion.yml`
- `pnjl-issue130-maxwell-cep-local-pilot.yml`
- `pnjl-maxwell-endpoint-local-production-shadow-v4.yml`
- `pnjl-maxwell-endpoint-production-shadow.yml`
- `pnjl-maxwell-endpoint-refinement.yml`
- `pnjl-stagec-density-certificate-feasibility-v2.yml`
- `pnjl-stagec-density-certificate-feasibility.yml`（盘点时 v1 workflow-specific test 仍引用；
  wave 1 已迁移到 retired definition）

盘点时 `relaxtime-issue130-rs-numerical-pilot-v1.yml` 的最新 run 失败且仍有不可变诊断
artifact，因此先列为候选审核。wave 1 已保留该失败 evidence 和原字节定义，只移除 Actions
触发入口。

### 4.3 首批退役候选（已获作者批准）

首轮提出以下四个候选；作者已于 2026-08-31 批准移出 Actions 目录并保留 versioned definition：

1. `pnjl-c2-cep-manual-bisection.yml`：无仓库活动引用，aggregate artifact 仍可下载；
   先保存 provenance，再决定 archive definition 或 delete。
2. `pnjl-issue130-crossover-mu-endpoint-expansion.yml`：无仓库活动引用，aggregate replay
   artifact 仍可下载；同样先保存 provenance。
3. `relaxtime-issue130-rs-numerical-pilot-v1.yml`：属于已归档 RS pilot，需保留失败/成本
   诊断记录；只有确认下游不再需要 workflow 重放后才可归档定义。
4. `pnjl-stagec-density-certificate-feasibility.yml`：v2 已替代，但历史测试和证据仍需
   明确迁移；在测试引用迁移前不得删除。

执行 PR 使用路径、当前文件 hash、Git blob SHA、替代 evidence 和引用迁移结果逐项记录；
inventory PR 本身没有执行这些动作。

## 5. 作者审核问题

首批四个候选已经完成作者决定。第二轮需审核
`diagnostic_workflow_retirement_wave1_v1/manual_dispatch_inventory_v2.csv` 和
`consolidation_groups.csv`：

1. 已闭环 Issue #130 专用入口是否进入 wave 2 直接退役。
2. 三个高结构相似组是整体退役，还是保留一个 versioned 参数化入口。
3. 两个 CEP oracle/pilot 是否仍有明确的未来复用需求。
4. 是否把未来新增的一次性 workflow 强制登记到本 inventory/ledger 规则。

在这些问题有明确答案前，不能改变第二轮候选的 workflow 定义。

## 6. 验收与后续入口

- [x] 47 个 workflow 文件已盘点并分类。
- [x] 20 个一次性诊断 run/artifact 已记录快照。
- [x] workflow call graph 和活动引用边已记录。
- [x] 精确 retirement candidate 清单已提出。
- [x] 未执行删除、停用、solver 或数值重跑。
- [x] 作者于 2026-08-31 批准首批四个候选退役，并要求继续审计其余手动 workflow。

首批执行包位于
`docs/analysis/governance/diagnostic_workflow_retirement_wave1_v1/`：四个 YAML 移出
`.github/workflows/` 并原字节保存在 `definitions/`，不删除历史 run/artifact。剩余 25 个纯
手动入口已形成第二轮分类和结构相似组；本 PR 不执行第二批退役或合并。

## 7. 第二波审阅包（2026-08-31）

wave 1 已在 PR #296 合并，合并 SHA 为
`15e02b8e8f9b23629908b6ee164d1e858f7ac463`。在不修改剩余 YAML 的前提下，第二波审阅包
`docs/analysis/governance/diagnostic_workflow_retirement_wave2_review_v1/` 对 25 个纯手动入口
中的 17 个给出实际功能、触发/输入、后端脚本、历史 evidence、可退役理由和必须保留的
provenance 边界。

- 8 个：已闭环 Issue #130/legacy 专用入口，进入直接退役审核；
- 7 个：phase shadow/Maxwell local 结构重复，进入参数化合并或整体退役审核；
- 2 个：CEP oracle/pilot，进入历史重放保留审核；
- 8 个：dense reference、archive/recovery 和 relaxtime production，继续保留。

本轮只生成 solver-free 审阅 evidence；不退役第二波 YAML。作者按组决定后，另立精确
allowlist 的 retirement 或 parameterization PR。
