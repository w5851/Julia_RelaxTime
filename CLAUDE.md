# CLAUDE.md

本文件是 Claude Code 的薄适配层，只保留 Claude Code 特有说明。跨工具的权威规则在 [AGENTS.md](AGENTS.md)；可执行命令在 [命令目录](docs/dev/agent_command_reference.md)；模块与子系统结构在 [架构总览](docs/architecture/agent_overview.md)。

## 必须遵守

- 与用户沟通默认使用中文。
- 规则冲突时以 `AGENTS.md` 为准；本文件不重复其中的编码规范、单位约定、测试分层与提交信息规则。
- 命令一律查 `docs/dev/agent_command_reference.md`，不要在文档或示例里另写一套并行版本。

## Claude Code 专属

- 技能：`.agents/skills/` 是权威定义，`.claude/skills/` 是由 `scripts/dev/sync_claude_skills.jl` 生成的适配副本；改技能只改 `.agents/skills/`，再重新生成。
- 命令宏：`.claude/commands/` 存放参数化动作宏。
- Copilot 的对应入口是 `.github/copilot-instructions.md`，同样指向 `AGENTS.md`。

## 需要细节时

- 架构与子系统：`docs/architecture/agent_overview.md`
- 测试分层与 profile：`docs/dev/testing_governance.md`
- 任务台账与跨主线分诊：`docs/dev/task_tracking_governance.md`
- 科研证据留存：`docs/analysis/governance/research_evidence.md`
