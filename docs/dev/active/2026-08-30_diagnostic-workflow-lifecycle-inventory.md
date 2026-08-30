# Diagnostic workflow lifecycle inventory v1

创建日期：2026-08-30

状态：`review_pending_author_retirement_approval`

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

### 4.2 保留定义和证据，先不退役

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
- `pnjl-stagec-density-certificate-feasibility.yml`（v1 已 superseded，但测试仍引用）

`relaxtime-issue130-rs-numerical-pilot-v1.yml` 的最新 run 失败且仍有不可变诊断 artifact；在
没有完成证据保留确认前也只列为候选审核，不执行删除。

### 4.3 退役候选（仅供作者批准）

当前只提出四个候选，不等于允许删除：

1. `pnjl-c2-cep-manual-bisection.yml`：无仓库活动引用，aggregate artifact 仍可下载；
   先保存 provenance，再决定 archive definition 或 delete。
2. `pnjl-issue130-crossover-mu-endpoint-expansion.yml`：无仓库活动引用，aggregate replay
   artifact 仍可下载；同样先保存 provenance。
3. `relaxtime-issue130-rs-numerical-pilot-v1.yml`：属于已归档 RS pilot，需保留失败/成本
   诊断记录；只有确认下游不再需要 workflow 重放后才可归档定义。
4. `pnjl-stagec-density-certificate-feasibility.yml`：v2 已替代，但历史测试和证据仍需
   明确迁移；在测试引用迁移前不得删除。

候选动作必须在单独 PR 中用路径、当前文件 hash、替代 evidence、引用迁移结果和作者批准逐项
列出。inventory PR 本身不执行这些动作。

## 5. 作者审核问题

作者只需对 `retirement_proposal.md` 中的候选组作决定：

1. 四个候选是否允许归档 workflow 定义、仅停用触发入口，或完全删除。
2. 哪些 run/artifact 必须先复制到仓库 analysis 包或外部不可变 archive。
3. 统一 artifact 保留期和复现入口；过期 artifact 是否已经有等价 hash-bound evidence。
4. 是否把未来新增的一次性 workflow 强制登记到本 inventory/ledger 规则。

在这些问题有明确答案前，下一步只能是补充引用审计或制作 dry-run allowlist，不能改变 workflow
定义。

## 6. 验收与后续入口

- [x] 47 个 workflow 文件已盘点并分类。
- [x] 20 个一次性诊断 run/artifact 已记录快照。
- [x] workflow call graph 和活动引用边已记录。
- [x] 精确 retirement candidate 清单已提出。
- [x] 未执行删除、停用、solver 或数值重跑。
- [ ] 作者审核 proposal 并批准下一独立 retirement PR。

作者批准后，另立 `diagnostic-workflow-retirement` PR，先执行 dry-run 路径/hash 校验，再按
allowlist 归档或删除；同时更新 `docs/guides/scripts/run_script_catalog.md` 中受影响脚本状态，
并重新运行 workflow/schema、docs、ledger、active-docs 和 data-path governance。
