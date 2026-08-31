# Diagnostic workflow retirement proposal v1

状态：`wave1_author_approved_2026-08-31`

本提案只定义后续动作边界，不是删除或停用指令。任何实际改变 `.github/workflows/` 的操作
必须另立 PR，使用精确路径和 SHA allowlist，并在合并前重新运行引用扫描。

## A. Wave 1 不进入退役清单

以下组在当前阶段继续保留：

- `retain_governance`：CI、ledger、docs、模型入口、依赖、数据路径和 smoke/benchmark gate。
- `retain_scheduled`：nightly regression、precompile trace budget。
- `retain_production`：dense reference、manual/meson/phase-guided production、sysimage。
- `retain_reusable`：dense/reference/raw archive 的 `workflow_call` 内部 job。
- `retain_archive_recovery`：raw rho-mu archive、Zenodo draft/restore 和 archive shard。
- `retain_legacy_audit`：本轮不执行 `pnjl-phase-diagram.yml` 的退役，因为 legacy audit 脚本仍
  保留其历史入口名。它不属于 wave 1 四文件 allowlist；在 wave 2 清单中仍可重新审核，若要
  退役必须先迁移该审计入口的历史引用。

这些分类不代表所有 artifact 永久保留；artifact retention 仍由各 workflow 的配置和 GitHub
实际到期时间决定。

## B. 历史重放组

CEP、Maxwell、Stage-C 及 crossover pilot/expansion workflow 目前保留定义和已归档证据。保留
理由是测试、活动 task、ledger pointer 或后续审计仍引用它们。若未来归档定义，必须先完成：

1. 将输入文件、calculation/postprocess SHA、run ID、artifact 名称和结果 hash 写入不可变
   analysis evidence。
2. 把仍有测试的 workflow 引用迁移到 versioned schema fixture 或明确的 replay entry。
3. 运行 `rg` 活动引用扫描和 workflow YAML schema 检查，确认没有隐式调用。
4. 在 retirement PR 中逐文件列出替代入口和 rollback 说明。

## C. 当前四个候选

| workflow | 当前观察 | 建议动作（待批准） | 前置条件 |
| --- | --- | --- | --- |
| `pnjl-c2-cep-manual-bisection.yml` | run `31999149922` success，1 个 aggregate artifact，未发现活动引用 | archive definition 或 delete | 复制 manifest/provenance，确认无作者复现需求 |
| `pnjl-issue130-crossover-mu-endpoint-expansion.yml` | run `32255786553` success，1 个 aggregate artifact，未发现活动引用 | archive definition 或 delete | 保留 crossover expansion evidence 和输入 hash |
| `relaxtime-issue130-rs-numerical-pilot-v1.yml` | run `32684074876` failure，1 个 artifact，历史 audit 有引用 | archive definition；不删除 evidence | 确认失败诊断仍可从 analysis 包重建，明确无 numerical rerun 需求 |
| `pnjl-stagec-density-certificate-feasibility.yml` | v1 run `31313661297` failure，无 artifact，仍有 v1 tests | archive definition after test migration | 迁移/删除 v1 tests，保留 v1 failure provenance |

“archive definition”表示让 YAML 不再成为可触发 workflow，同时在 Git 历史或明确的 analysis
快照中保存内容；具体采用移动到文档快照、禁用触发器还是删除文件，必须由作者逐项选择。

## D. Artifact/provenance 保留合同

- GitHub artifact 的 `expires_at` 只记录当前观察，不作为长期证据。
- 需要长期复现的 aggregate/manifest/claim ledger 应保存在仓库 analysis 包或已有外部不可变
  archive，并记录 SHA-256、source run、calculation SHA 和 schema。
- 不下载或重打包原始全量 solver 曲线，除非已有 task 明确授权；本 inventory 只记录名称和
  计数。
- 失败 run 也必须保留其失败状态和无 artifact 事实，不能用删除 workflow 掩盖失败。

## E. 后续独立 PR 的最小验证

1. `git diff --check` 和精确 allowlist/hash preflight。
2. `rg` workflow、脚本、测试、active docs 和 ledger 引用扫描。
3. workflow YAML/schema、docs consistency、task ledger、active-docs 和 data-output governance。
4. 对受影响脚本运行 `check_script_entrypoints.jl`，先将 `deprecate-candidate` 状态写回
   `docs/guides/scripts/run_script_catalog.md`。
5. 确认现有 run/artifact/evidence URL、manifest 和恢复说明仍可追溯。

## 作者批准栏

- [x] 批准四个候选进入 retirement PR。
- [x] 批准将四个 YAML 移出 `.github/workflows/`，原字节定义保存在 versioned analysis 包。
- [x] 批准沿用现有仓库 evidence/manifest 和历史 run URL，不删除 GitHub artifact。
- [x] 批准历史失败 run 继续以 diagnostic-only 语义保留。

执行包：`docs/analysis/governance/diagnostic_workflow_retirement_wave1_v1/`。
剩余手动入口的合并或第二批退役不属于本次批准范围，必须另行审核。
