# AI 过程产物（已冻结）

**状态：`frozen_historical`。不再新增，也不作为现行权威规范引用。**

- `plans/`：29 份历史实施计划。
- `specs/`：27 份历史规格说明。

这些文件是早期 AI 规划工作流的产物，未接入 `docs/dev/` 的生命周期管理。保留原路径只为追溯历史 commit / PR 上下文。

## 现行做法

- 设计草案与任务合同放 `docs/dev/`：进行中放 `docs/dev/active/`，延后放 `docs/dev/backlog/`，完成后按 [归档脚本](../../scripts/dev/archive_docs.jl) 归档。
- 一次性计划优先留在对话、issue 或 PR，不落盘进 `docs/`。
- 落点判断规则见 [AGENTS.md](../../AGENTS.md)。

## 后续处置

若确认无活跃引用，可在 `docs/dev/README.md` 中登记为只读历史区，或整体并入 `docs/dev/archived/`；在此之前保持原地不动，避免破坏历史链接。
