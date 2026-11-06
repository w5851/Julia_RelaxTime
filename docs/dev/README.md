# 开发文档总览

本目录保存需要跨对话共享的设计、开发约定与历史依据。普通工作使用对话计划或已有 issue/PR，不默认创建任务单、执行台账或归档副本。

## 维护入口

| 需要解决的问题 | 权威入口 |
| --- | --- |
| 项目约束、目录职责与编码约定 | [AGENTS.md](../../AGENTS.md) |
| 计划、任务恢复与跨主线状态 | [任务跟踪治理](task_tracking_governance.md) |
| Codex 与 skill 的使用边界 | [Codex 协作手册](Codex高阶使用手册.md) |
| 环境、测试、治理与 benchmark 命令 | [命令参考](agent_command_reference.md) |
| 验证层级与测试选择 | [测试治理](testing_governance.md) |
| 产物位置与来源记录 | [交付物管理](项目交付物定位与管理规范.md) |
| 科研分析、实验事实与结论留存 | [科研证据治理](../analysis/governance/research_evidence.md) |
| Julia 模块加载与依赖 | [模块规范](Julia模块引用规范准则.md)、[架构规则](../architecture/dependency_rules.md) |

## 临时文档的使用

仅在用户要求交付文档、需要共享设计/验收依据，或现有 ledger 合同要求文件时使用任务文档。是否保留文档取决于信息的用途，不取决于任务是否复杂或执行轮数。

- `active/`：当前仍需维护的任务合同、设计草案或恢复说明。保持当前范围、未完成项、验收条件与证据链接清楚。
- `backlog/`：尚未执行的候选需求和研究方向。启动时可直接用对话计划或 issue；只有存在上述文件需求时才新建 active 文档。
- `archived/`：已完成、取消或被替代工作的历史记录，明确其终止原因。普通续接不读取该目录；需要历史依据时再定向查找。

当前状态摘要可原位更新，详细实验/验收事实放在相应证据位置并链接，避免在任务单末尾逐轮追加。不要为了缩短摘要重写原始实验记录、作者决策或哈希绑定文档；具体保留边界见[任务跟踪治理](task_tracking_governance.md)。

active 和 backlog 文档使用 `YYYY-MM-DD_描述.md`。超过 60 天只触发复核提醒，不自动表示任务完成或失效。未完成工作退出当前批次时转入 backlog，并保留剩余工作与依赖；不要以完成归档的方式清掉待办。

## 已有任务文档的归档

已完成、取消或被替代且仍有历史用途的任务文档，使用
[`archive_docs.jl`](../../scripts/dev/archive_docs.jl) 迁移，保留正文和终止原因。
未完成事项继续由当前任务或 backlog 承接。没有任务文档的工作无需补建归档。

```powershell
julia --project=. scripts/dev/archive_docs.jl --dry-run <filename.md>
julia --project=. scripts/dev/archive_docs.jl <filename.md>
julia --project=. scripts/dev/archive_docs.jl --check <filename.md>
```

取消或被替代使用 `--status cancelled|superseded --reason "终止原因"`；
脚本帮助和[归档技能](../../.agents/skills/doc-archive/SKILL.md)说明完整参数。
迁移后更新当前引用和 ledger，检查本批文件即可。

## 文档和验证随变更维护

稳定入口或数据契约变化时更新 `docs/api/` 或 `docs/guides/`；公式、方法和单位变化时维护 `docs/reference/`。内部拆分不要求 API 文档逐文件镜像，也不要求为每轮工作新增过程文档。

按[测试治理](testing_governance.md)选择能验证本次行为的覆盖。数值漂移、性能主张和生产晋升分别遵循其专项证据要求。
