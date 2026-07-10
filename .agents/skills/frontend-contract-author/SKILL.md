---
name: frontend-contract-author
description: 从 Julia_RelaxTime 的 Models、脚本或 server 权威入口提炼前端可消费的数据契约，包括字段、类型、单位、顺序、缺失值、错误语义、稳定性和 mock payload。用于契约文档与缺口标注，不直接实现 UI，也不虚构尚不存在的后端接口。
---

# Frontend Contract Author

## Hard rules

- 只记录源码、脚本或服务中实际存在的接口和字段。
- 每个数值字段写明单位、值域和缺失值语义。
- 区分稳定契约、内部实现细节和建议新增项。
- 不把建议项写成已经可用的接口。

## Authoritative sources

- `src/models/entrypoints.jl`
- `src/simulation/`
- `scripts/server/server_full.jl`
- 相关稳定脚本和现有 `docs/api/`

稳定契约写入 `docs/api/`；过渡性协作契约写入 `docs/dev/`。

## Workflow

1. 定义前端消费场景和交互边界。
2. 定位统一后端入口，而不是深层 helper。
3. 提取请求参数：名称、类型、必填性、默认值、单位和枚举。
4. 提取响应字段：类型、单位、顺序、索引和缺失值。
5. 描述参数错误、求解失败、空结果、部分结果和降级结果。
6. 标注稳定性并生成最小成功与边界 payload。
7. 列出前端需要但后端尚未实现的明确缺口。

## Output contract

- 场景与权威入口
- 请求参数表
- 响应字段表
- 单位、值域、顺序与缺失语义
- 错误和边界状态
- 最小成功、空结果、失败和部分结果 payload
- 稳定契约 / 内部细节 / 建议新增分类
