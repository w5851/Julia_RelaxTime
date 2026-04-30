---
name: frontend-contract-author
description: 为 Julia_RelaxTime 提炼前端可消费的后端数据契约。适用于单点求解、扫描结果、脚本输出、HTTP 接口、JSON/CSV 字段表、单位说明、错误语义、mock payload。关键词：frontend contract, API schema, payload, units, docs/api, docs/dev
---

# frontend-contract-author

## Purpose

在不急着做 UI 的前提下，把后端入口整理成前端可直接消费的契约文档。

本 Skill 的目标是交付字段、单位、错误语义和稳定性边界，而不是直接设计界面。

## Apply When

- 需要为未来前端页面准备后端契约。
- 需要从 `Models` 入口、脚本入口或 server 路径抽取 request/response 结构。
- 需要整理 JSON/CSV 输出字段、单位、排序、缺失值语义。
- 需要生成 mock payload 或前端集成说明。

## Do Not Use

- 用户明确要求先实现 UI：先完成契约后再转前端实现。
- 后端能力尚未存在且任务目标是新增接口实现：此 Skill 只负责契约整理与缺口标注。

## Hard Rules

- 不虚构后端已存在的字段或接口。
- 每个数值字段默认写清单位。
- 区分稳定契约与内部实现细节。
- 若发现契约缺口，标注为“建议新增”，不要伪装成现有能力。

## Repository Defaults

- 优先入口：
  - `src/models/entrypoints.jl`
  - `src/simulation/`
  - `scripts/server/server_full.jl`
  - `scripts/`
- 优先参考文档：
  - `docs/api/README.md`
  - `docs/api/data_contracts.md`
  - `docs/guides/FRONTEND_DEBUG.md`
- 推荐落点：
  - 稳定公开契约：`docs/api/`
  - 过渡性前后端协作文档：`docs/dev/`

## Standard Workflow

1. Identify the consumer flow
   - 明确页面或交互：单点求解、扫描提交、结果列表、图表展示。
2. Locate the authoritative backend surface
   - 优先统一入口与现有 server/script，而不是深层内部 helper。
3. Extract the contract
   - 请求参数：名称、类型、是否必填、单位、默认值、枚举口径
   - 响应字段：名称、类型、单位、排序/索引语义、缺失值语义
   - 错误模型：参数错误、求解失败、空结果、降级结果
4. Classify stability
   - 稳定契约：可推荐前端直接依赖
   - 内部细节：仅实现态存在，不建议绑定
   - 建议项：尚未实现但前端需要
5. Produce examples
   - 给最小 request/response payload
   - 给边界状态样例，如空数据、失败、部分结果

## Output Contract

最终文档至少包含：

- 场景说明
- 后端权威入口
- 请求参数表
- 响应字段表
- 单位与取值语义
- 错误与边界状态
- 最小 payload 示例
- 稳定契约 / 内部细节 / 建议新增 三分标注
