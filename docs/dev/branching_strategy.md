# Git 分支与提交工作流

当前项目围绕 `main` 和短期任务分支开发。原先的 `feature → develop → main` 模板已退出当前规范；历史使用方式可通过 Git 记录追溯。

## 开始工作

- 先检查 branch、HEAD、dirty paths 和 worktree；保留用户已有改动。
- 基准以当前任务、PR 或已明确的 `main` 提交为准。
- Codex 新建任务分支默认使用 `codex/<topic>`；用户指定的名称优先。
- 当前目录被另一项工作使用或有无关改动时，用独立 worktree；普通小改动不强制新建 worktree。
- 聚焦问题可以直接由 issue/PR 追踪；跨主线依赖按 [任务治理](task_tracking_governance.md) 记录。

## 实现与提交

- 提交只包含本任务的文件，选择明确路径并审阅 staged diff。
- 提交主题准确描述变化；使用 `docs:`、`fix:`、`refactor:`、`ci:`、`feat:` 等已有前缀即可。
- 为任务选择 [受影响的验证](testing_governance.md)；验证足够时结束，不按文件数量追加测试。
- 数值语义变化需要对应回归或误差证据；性能主张需要相应测量。
- 公开稳定入口、脚本合同或数据字段变化时同步相应说明。

## PR 与集成

默认将任务分支合入 `main`。PR 说明范围、用户影响和相关验证；拆分依据独立可审阅的变化，不按固定行数或时间限额拆分。

实际 required checks、审批人数和分支保护以 GitHub 设置为准。仓库文档不预设当前远端已启用哪些保护；额外审查按数值、架构或外部影响选择。

文档、小范围修复和内部重构不附带全套性能基准、覆盖率配额或额外 Release。GitHub Actions 的触发范围见 `.github/workflows/`；sysimage 发布由对应显式入口管理。

## 收尾

核对已提交范围、验证结果及仍需用户操作的步骤。分支/worktree 清理前确认没有未保存工作或其他任务占用；按已授权范围清理。

相关入口：[代码审查](code_review_guidelines.md)、[命令参考](agent_command_reference.md)。