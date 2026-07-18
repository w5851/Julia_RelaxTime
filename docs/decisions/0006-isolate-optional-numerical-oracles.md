# ADR-0006: 隔离可选数值 oracle 与精简 agent 指令

## 状态

已接受

## 背景

QuadGK 自提交 `f08256cc` 起作为根环境直接依赖存在，最初用于 B0 测试和对照积分。当前 production 源码已改用仓库内 Gauss--Legendre、奇点减法和显式误差控制；PR #137 又从 PNJL phase thermodynamics 核心路径移除了 QuadGK。继续把它留在根 `Project.toml` 会扩大主环境解析与预编译面，并模糊 production 算法和外部对照 oracle 的边界。

根 `AGENTS.md` 同时累积了大量易变命令清单，达到 456 行。Agent 每次读取这些命令会挤占用于任务事实和代码上下文的空间，而命令本身更适合由可链接的开发指南维护。

## 决策

1. 根 `Project.toml` 不直接依赖 QuadGK；`src/`、`scripts/` 和 `tests/` 不导入或调用 QuadGK。
2. QuadGK 只允许作为 `benchmark/Project.toml` 中的隔离对照 oracle。先用 `--project=benchmark` 实例化；需要同时 include 根源码时，必须显式把 `benchmark` 环境叠加到仅该 benchmark 进程的 `LOAD_PATH`。benchmark 结果不得自动成为 production 正确性证据。
3. 根环境测试和分析脚本使用仓库内规则或根环境已有的 `FastGaussQuadrature` 构造独立固定/收敛对照。
4. 根 `AGENTS.md` 只保存高优先级、跨仓库适用的强约束，行数上限为 220。完整命令清单迁到 `docs/dev/agent_command_reference.md`。
5. 只有子树存在真正局部且不适用于全仓的约束时才增加 scoped `AGENTS.md`；不得用 scoped 文件复制根规则或命令目录。
6. `scripts/dev/check_dependency_policy.jl` 与 `scripts/dev/check_agent_instructions.jl` 在 dependency-audit CI 中执行，并由 unit 测试覆盖有效/无效 fixture。

## 理由

- 将 production 依赖与验证 oracle 分离，可降低冷启动和预编译风险，也让 manifest provenance 更容易解释。
- 保留 benchmark 独立环境，使高精度外部对照仍可用，而不把可选工具传播到核心运行时。
- 精简 agent 指令遵循 progressive disclosure：先加载稳定规则，需要命令时再读取权威指南。
- 机读门禁避免该约定依赖人工记忆。

## 备选方案

### 保留 QuadGK 根依赖但禁止核心导入

实现成本最低，但仍扩大根环境并让“可用依赖”与“允许依赖”相互矛盾，因此不采用。

### 完全禁止 QuadGK，包括 benchmark

边界最简单，但会失去一个成熟的独立对照工具。预编译问题来自根环境传播而不是隔离 benchmark 本身，因此不采用。

### 把全部 AGENTS 内容拆成多个 scoped AGENTS

会引入继承和重复规则的审计成本。当前主要问题是命令目录过长，不是子树规则差异，因此先采用一个精简根文件和一个权威命令指南。

## 后果

正面影响：

- 根 manifest 不再包含无其他依赖需要的 QuadGK；核心 `Models` 加载不会触达该包。
- unit、analysis 与 production 的数值依赖边界明确。
- 根 agent 指令显著缩短，命令仍可追溯且由 CI 验证。

负面影响与风险：

- benchmark 使用者必须先实例化 `benchmark/` 环境。
- 固定节点 oracle 需要显式检查节点收敛，不能把单次高节点结果当成自带误差估计。
- 若未来核心算法确实需要 QuadGK，必须通过新 ADR 或替代本决策，并提供预编译、数值和 AD 兼容证据。

## 相关决策与实现

- [依赖规则](../architecture/dependency_rules.md)
- [Agent command reference](../dev/agent_command_reference.md)
- `config/ci/dependency_policy.toml`
- `scripts/dev/check_dependency_policy.jl`
- `scripts/dev/check_agent_instructions.jl`
- PR #137：phase thermodynamics 核心 QuadGK 移除

---

**日期**：2026-07-18
**作者/审查者**：项目作者与 Codex 协作审查
