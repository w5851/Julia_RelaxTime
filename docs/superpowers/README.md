# AI 过程产物（已退役，只读）

**状态：`retired_readonly`（已退役，只读）。不再新增文件，也不作为现行权威规范引用。** 本目录保留 56 份早期 AI 规划工作流产物（`plans/` 29 份、`specs/` 27 份）。因 [`asset_registry.json`](../analysis/governance/figure_asset_registry_v1/asset_registry.json) 等历史注册表直接引用了目录内路径作为图件 provenance，所有文件原地只读保留，不再移动或并入 `docs/dev/archived/`；禁止新增由 [`check_active_docs_governance.jl`](../../scripts/dev/check_active_docs_governance.jl) 与 [`retired_docs_allowlist.toml`](../../config/governance/retired_docs_allowlist.toml) 强制门禁。

## 现行做法

- 设计草案与任务合同放 `docs/dev/`：进行中放 `docs/dev/active/`，延后放 `docs/dev/backlog/`，完成后按 [归档脚本](../../scripts/dev/archive_docs.jl) 归档。
- 一次性计划优先留在对话、issue 或 PR，不落盘进 `docs/`。
- 落点判断规则见 [AGENTS.md](../../AGENTS.md)。

## 历史处置说明

这些文件是早期 AI 规划工作流的产物，未接入 `docs/dev/` 的生命周期管理；保留原路径用于追溯历史 commit / PR 与图件 provenance，不再迁并。
