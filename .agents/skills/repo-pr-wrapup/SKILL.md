---
name: repo-pr-wrapup
description: 为 Julia_RelaxTime 的“提交改动 + 开 PR + 收尾验证”提供标准化流程。适用于用户要求“提交”“开 PR”“收尾”“准备合并”“整理本次改动并发 PR”“处理 review comments 后更新 PR”等场景。该 skill 会先判断是否应优先处理现有 review comments，再按仓库约束完成分支命名、选择性 stage、历史提交风格匹配、分层测试与治理检查、PR 标题与固定正文整理，以及 `docs/dev/active/` 归档检查。
---

# Skill: repo-pr-wrapup

## 核心规则

1. 先判定当前工作是“更新已有 PR”还是“新开 PR”。若存在 review comments 或用户明确要求处理评论，优先处理评论，不重复走全新的 PR 创建流程。
2. 默认分支名使用 `codex/{topic}`。若当前已在合适分支上，沿用现有分支；不要为了形式切换到无关分支。
3. 提交前先运行 `git status --short` 审视改动范围，再按 scope 选择性 stage；默认不要使用 `git add -A`。
4. 提交信息必须先读取最近 10 条历史提交并匹配已有前缀风格，再拟定本次 commit message。
5. 开 PR 前必须执行与改动层级匹配的 tests / governance checks；无法运行时要明确写出原因与风险。
6. PR title 按 squash merge 的最终提交标题来写，简洁、可直接进入主干历史。
7. PR body 固定覆盖：变更范围、用户影响、实现方式、验证项、已知非目标。
8. 若改动涉及 `docs/dev/active/`，必须检查是否应归档、更新状态，或保持 active 并说明原因。

## 工作流

### 1. 判断流程入口

先做以下判断：
- 是否已有对应分支和未提交改动
- 是否已有打开中的 PR
- 是否存在待处理的 review comments / review threads
- 用户是要“提交并开 PR”，还是只做其中一部分

若存在 review comments，先进入“评论优先”路径：
- 先读取评论与关联文件
- 只实现与评论相关的必要修正
- 重新运行受影响范围的验证
- 更新已有分支 / PR，而不是重复构造新一轮 PR 文案

### 2. 整理分支与改动范围

1. 查看当前分支。若尚未进入任务分支且用户没有指定命名，使用 `codex/{topic}`。
2. 运行 `git status --short`，按文件与任务目标划定本次提交 scope。
3. 识别无关改动、用户已有改动、生成物或临时文件；除非明确属于本次 scope，不要 stage。
4. stage 后再次检查 `git diff --cached --stat` 或等价信息，确认暂存内容与本次任务一致。

### 3. 匹配提交信息风格

1. 运行 `git log -10 --oneline`。
2. 提取最近历史里已存在的前缀模式，优先匹配与本次改动类别最接近的样式。
3. 仅在无法可靠判断时，使用仓库允许的 fallback：`docs:`, `fix:`, `refactor:`, `ci:`。
4. 保持标题单行、意图明确，避免空泛表述。

### 4. 选择验证集

按改动层级选择最小但足够的验证：
- `tests/unit/` 变更：优先 unit smoke/core 或受影响单文件
- `tests/integration/` 触达的工作流改动：补 integration smoke/core
- 数值路径、solver、scan、`src/relaxtime/`、`src/models/`、`src/simulation/` 改动：评估 regression / validation / 相应 governance checks
- 文档与治理改动：运行对应 `scripts/dev/check_*.jl`
- 涉及 `docs/dev/active/`：额外检查是否需要归档，必要时转交或联动 `doc-archive`

不要为了省事跳过高风险验证，也不要为了过测试而放宽容差。若跳过验证，必须在最终汇报和 PR body 中说明。

### 5. 生成 PR 标题与正文

PR title:
- 按 squash merge 目标来写
- 与最终期望提交标题同风格
- 不写冗余前缀如 “WIP” 或聊天式说明，除非用户明确要求 draft

PR body 固定包含以下五段：
1. `变更范围`
2. `用户影响`
3. `实现方式`
4. `验证项`
5. `已知非目标`

若已有 PR，则更新已有正文时保持这五段结构，不随意丢字段。

### 6. 检查活动文档归档义务

当改动涉及 `docs/dev/active/` 时：
- 判断该文档对应任务是否已完成
- 若已完成，检查是否应转 `doc-archive` 流程
- 若未完成但本次更新了状态，确保文档状态、代码状态、验证结果一致
- 不要把尚未验证完成的事项写成已完成

## 输出要求

完成该 skill 所对应任务时，最终汇报至少包含：
- 当前分支与是否新建分支
- 选择性 stage 的文件范围
- 使用的提交信息及其历史风格依据
- 已运行验证与结果
- PR title
- PR body 五段摘要
- 是否检查了 `docs/dev/active/` 归档义务
- 若未实际执行 commit / push / 开 PR，明确停在了哪一步

## 边界

- 本 skill 负责仓库内提交与 PR 收尾流程，不负责大范围功能设计。
- 若核心任务是“处理 GitHub 评论并回填修复”，优先结合现有 comment 处理流程推进。
- 若核心任务是归档活动开发文档，转 `doc-archive`。
- 若用户只要求审查而非执行，不要擅自提交、推送或创建 PR。

## 触发示例

- “用这个仓库的规范把这次改动提交并开 PR”
- “帮我收尾一下，先看该 stage 什么，再准备 PR”
- “处理 review comments 后更新 PR”
- “按仓库约束完成 commit 和 PR 文案”
