# Diagnostic Workflow Retirement Wave 2 v1

本包执行作者于 2026-08-31 对第二波手动诊断 workflow 的决定：将 16 个已完成、一次性或历史专用入口移出 `.github/workflows/`，保留原始 YAML 定义、测试合同和历史 run/artifact provenance；当时额外保留了一个带 target-list 参数化能力的 Maxwell-local expansion 入口，该保留决定已于 2026-11-04 撤销并删除该入口。

## 执行语义

- “退役”只表示 YAML 不再位于 `.github/workflows/`，因此不再出现在 Actions 手动入口；不删除 GitHub run/artifact，不改 solver、reference、production 数据或数值结果。
- 每个退役 YAML 在 `definitions/` 中按原字节保存，并记录 SHA-256、字节数、来源提交和历史 evidence。
- workflow-specific Python/Julia tests 改为读取 versioned definition，继续验证原有输入、SHA、scope、failed-only、artifact 和 `solver_called` 合同。
- `pnjl-issue130-maxwell-cep-local-expansion.yml` 当时作为 Maxwell-local family 唯一的 target-list/failed-target-only 参数化入口保留；该入口已于 2026-11-04 与其 runner、collector 一并删除，YAML 可从删除前的 Git 提交取回。

## 本轮范围

退役 16 个入口：

- 8 个 Issue #130/C2/phase-diagram 一次性入口；
- 4 个 CEP/Maxwell phase-shadow 入口；
- 2 个 CEP oracle/pilot 历史入口；
- 2 个 Maxwell-local 专用入口（pilot、endpoint refinement）。

迁移 `pnjl-phase-diagram.yml` 前已将 legacy-retirement audit 的活动路径契约改为版本化历史 definition；旧 audit/evidence 中出现的原路径继续作为历史记录保留。

## 计数与边界

基于 wave 1 合并后的状态，本轮预期 active workflow 从 43 降至 27，含 `workflow_dispatch` 的入口从 39 降至 23，纯手动入口从 25 降至 9。计数不包含已保存的 versioned definitions。

本包不宣称任何数值 workflow 成功或失败被改变；它只记录入口生命周期变化。未来若需要 CEP oracle 或历史 pilot，必须从其 versioned definition 复制出新的、带明确 scope 和 provenance 的 workflow，而不是恢复旧入口。

## 证据文件

- `retired_workflows.csv`：16 个退役入口的功能、理由、历史 run 和 hash；
- `reference_migrations.csv`：测试、活动任务单和 legacy audit 的路径迁移；
- `claim_ledger.csv`：本轮声明与可复核证据；
- `definitions/`：原始 YAML 的字节保持副本；
- `manifest.json`：来源提交、计数、执行语义和文件 hash。

本轮 `solver_called=false`、`numerical_results_changed=false`、`github_runs_deleted=false`、`github_artifacts_deleted=false`。
