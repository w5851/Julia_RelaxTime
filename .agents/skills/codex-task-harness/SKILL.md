---
name: codex-task-harness
description: 将宽泛的 Julia_RelaxTime 任务收束为可执行 harness：边界、风险、验证、交付物与推荐 skill。适用于跨后端、前端契约、文档、论文、研究分析的混合任务。关键词：task harness, scope, validation, skill routing, execution plan
---

# codex-task-harness

## Purpose

把“大而散”的需求快速整理成一个可以直接执行的任务框架，减少中途漂移。

本 Skill 不追求长计划，而是尽快给出范围、约束、风险、验证和推荐 skill，然后进入执行。

## Apply When

- 用户需求跨多个方向：后端实现、文档、前端契约、论文、实验分析。
- 任务描述宽泛，容易在实施中扩边界。
- 需要先判断该调用哪些 skill，再落地。

## Do Not Use

- 任务已经非常窄且边界明确时；直接执行即可。
- 用户明确只要 brainstorming，不要实现时；保留为轻量分析。

## Hard Rules

- Harness 只保留执行真正需要的信息，避免大段泛化分析。
- 先界定不该改什么，再界定要改什么。
- 必须给验证方式，不能只给实现思路。
- 如果已有合适 skill，优先路由而不是重复发明流程。

## Output Frame

每次使用时，至少给出：

- 目标：这次要完成什么
- 边界：允许修改 / 不允许修改
- 契约：必须保持不变的行为或入口
- 风险：最可能出问题的 1 到 3 点
- 验证：最小必要测试、文档更新或结果核对
- 路由：建议同时使用的其他 skill

## Repository Routing Heuristics

- 数值回归与 baseline：`transport-regression-keeper` 或 `baseline-regression-governance`
- API 文档与稳定入口：`api-doc-authoring`
- 前端字段与 payload：`frontend-contract-author`
- 论文与研究写作：`paper-traceable-coauthor` 或 `doc-coauthoring`
- 实验日志：`experiment-logbook-append` 或 `doc-exec-log-append`
- Julia 数值实现细节：`julia-pro`

## Standard Workflow

1. Collapse the request
   - 用一句话写主目标。
2. Freeze boundaries
   - 文件边界、契约边界、非目标。
3. Pick the verification floor
   - 至少要跑什么、检查什么、更新什么。
4. Route to skills
   - 只选最少但足够的 skill。
5. Execute
   - 当边界足够清楚时，直接进入实现而不是停留在 plan。

## Prompt Skeleton

可将需求整理成：

```text
目标：
- ...

允许修改：
- ...

不要修改：
- ...

必须保持不变：
- ...

验证要求：
- ...

建议同时使用：
- ...
```
