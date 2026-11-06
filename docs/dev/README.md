# 开发文档总览

本目录保存需要跨对话共享的设计、开发约定与历史依据。普通工作使用对话计划或已有 issue/PR，不默认创建任务单、执行台账或归档副本。

## 工作区收尾状态（2026-11-06）

当前计算项目登记两个 worktree：本会话 `D:/w/transport-v15-paper` 和主工作区
`D:/Desktop/Julia_RelaxTime`。远程仅有 `main`，无待合并 PR（本次整理 PR 另计）。
已合入主干的本地 `codex/workflow-guidance-trim` 已回收。

桌面主工作区仍为 `codex/docs-architecture-cleanup`，tip
`5946aaf73a4d26f0f392bc7d43b49419f0ee6f26`。相对整理前主干
`703c23d234bfebd107a9f1acc7179a8290bd1640` 有 6 个未合入提交、2 个已跟踪脏文件和
1783 个未跟踪文件；暂存区为空。它不是可直接删除的旧隔离目录。

| 范围 | 收尾处理 |
| --- | --- |
| `feed2c99` 的已合并任务说明 | 冻结线 PR #327 已合并的事实已同步到当前入口；审阅 pending 继续保留。旧图清理任务的归档可随桌面分支收尾处理 |
| `e5c21ebf` 的 Web/API 导航 | 独立文档改动，需按当前入口核对后选入收尾 PR |
| `4b91ee73` 的早期冻结线文献与交接 | 与已保留的 PR #327 五个归档查重，提取仍独有的研究解释；不恢复重复数值产物 |
| `c498d85e` 的旧 FIG4 留存 | v15 已替代；保留提交恢复来源，不把旧输运图重新带入当前目录 |
| `9d1c637d` 的早期 APS 包归档 | 先核对现行相图交付是否仍需要，不因分支收尾重建整套备份 |
| `5946aaf7` 的 Issue130 工具退役 | 独立工程改动；重新检查当前脚本／workflow 消费者并验证后再合入 |
| 两处脏文档及早期输运研究（约 22 MiB） | 含旧 Mott 推导、k=0 完整 RS 比较、Julia–Fortran 分析及源码，尚无完整已核验的持久替代；保留原处逐项收口。旧台账里的“待整批重算”已由 PR #331–334 完成，当前状态看输运专题 |

当前[任务台账](../../config/governance/task_tracks.toml)继续以冻结线研究为主线；
输运 v15 已交付，GBU 与冻结线的作者审阅、完整输运正则性和 A 降维等研究问题分别保留。
桌面分支后续按上述内容选择合入，不能整体 merge 后重新引入已退役资产。

`D:/w` 另外三个目录均已失去 Git worktree 身份：`charged-gbu-contours` 只剩空标记，
其 GBU 产物已移至 `D:/DevData/Julia_RelaxTime/artifacts/charged_gbu_pr325`；
`project-skills-wrapup` 只剩空标记；`meson-freezeout-review-closeout-20261022` 只剩
PR #327 的临时检查 ZIP 和收据。三处残留的删除被自动审批以 `blocked by policy` 拒绝，未执行。
本会话 worktree 的旧临时缓存继续按作者意见留待回收时处理。

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

已完成或已明确取消/被替代、且需要保留历史的文档通过[归档脚本](../../scripts/dev/archive_docs.jl)迁移。没有创建任务文档的普通工作无需补建文档再归档。

归档使用显式文件名，先预览，再执行并检查：

```powershell
julia --project=. scripts/dev/archive_docs.jl --dry-run <filename.md>
julia --project=. scripts/dev/archive_docs.jl <filename.md>
julia --project=. scripts/dev/archive_docs.jl --check <archived-filename.md>
```

默认终态是 `completed`。取消或被替代使用 `--status cancelled` 或 `--status superseded`，并传入 `--reason "终止原因"`；预览和执行使用同一组参数。`--date YYYY-MM-DD` 可指定归档日期。

元信息由脚本生成，正文按原始字节保留：

```yaml
---
title: "任务名称"
archived: true
original: "docs/dev/active/原始任务文件.md"
archived_date: "2026-10-07"
task_status: "completed"
---
```

脚本以自身位置定位仓库，只接收 `active/` 的直接 Markdown 子文件，拒绝链接重定向和已有目标。帮助、检查和预览不写文件；批量操作先完成整批参数与路径检查，但执行期间失败不自动回滚已经归档的文件。

写入时先校验同目录临时文件，通过硬链接发布完整目标，再核对原文未变后移除源文件。目标文件系统必须支持硬链接；发布失败会保留源文件，目标已发布后的校验失败会保留两份文件供核对。

`--check <文件...>` 只验证本批目标；省略文件名才检查全部归档。检查器支持平面的标量 frontmatter，核对必需字段、类型、日期及终止原因，拒绝重复字段；不解析任意嵌套 YAML。错误返回非零退出码。旧文件没有 `task_status` 仍可检查，不因此批量重写历史。

迁移前核对状态与验收证据，迁移后由调用者更新仍有效的引用和 ledger 路径，脚本不推断这些关系。取消或被替代不能标为验收完成；不要破坏冻结来源。已有历史可通过 Git 和归档路径恢复，不因减少上下文而批量删除。

## 文档和验证随变更维护

稳定入口或数据契约变化时更新 `docs/api/` 或 `docs/guides/`；公式、方法和单位变化时维护 `docs/reference/`。内部拆分不要求 API 文档逐文件镜像，也不要求为每轮工作新增过程文档。

按[测试治理](testing_governance.md)选择能验证本次行为的覆盖。数值漂移、性能主张和生产晋升分别遵循其专项证据要求。
