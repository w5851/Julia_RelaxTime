---
name: codex-task-harness
description: 将宽泛、跨模块的 Julia_RelaxTime 需求收束为可执行 harness：目标、边界、契约、风险、验证、交付物与最小 skill 路由。仅用于容易扩 scope 的混合任务；边界清楚的窄任务应直接执行。
---

# Codex Task Harness

## Purpose

把“大而散”的请求压缩成可以立即执行的边界与验证框架，避免计划本身成为交付物。

## Hard rules

- 先写不允许改变的范围和契约，再写允许修改的内容。
- 只保留会影响执行的风险和假设。
- 给出最小必要验证，不以“以后补测试”代替验收。
- 只路由最少但足够的 skill；边界清楚后立即执行。

## Harness contract

- 目标：本轮必须完成的结果
- 允许修改：文件、模块、数据或文档范围
- 不要修改：显式非目标和外部影响边界
- 必须保持：稳定入口、求解语义、单位或数据契约
- 风险：最可能失败的 1–3 点
- 验证：测试层、治理检查、结果核对和文档义务
- 路由：本轮实际需要的其他 skill

## Repository routing

- 数值 baseline 制度：`baseline-regression-governance`
- transport/relaxtime 漂移：`transport-regression-keeper`
- API 文档：`api-doc-authoring`
- 前端数据契约：`frontend-contract-author`
- 论文证据写作：`paper-traceable-coauthor`
- 技术文档：`doc-coauthoring`
- 实验日志：`experiment-logbook-append`
- 执行台账：`doc-exec-log-append`
- Julia 实现：`julia-pro`

## Workflow

1. 用一句话锁定主目标。
2. 冻结文件、语义和外部影响边界。
3. 选择验证下限和文档义务。
4. 选择最小 skill 组合。
5. 执行并在范围或风险实质变化时重新锁定边界。
