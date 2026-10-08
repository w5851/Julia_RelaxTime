# 依赖规则（目录级）

本页约束当前源码的文件加载方向。`scripts/dev/analyze_deps.jl` 直接扫描 `src/`，生成的 Mermaid 图只用于阅读，不作为审计输入。

## 分组与方向

`base` 合并 `src/` 根部兼容文件及 `config/`、`constants/`、`types/`；这些目录存在受控的常量/类型转发。其余分组保留实际目录名。

| 来源 | 允许依赖 |
| --- | --- |
| `base` | `base` |
| `utils` | `base`、`utils` |
| `integration` | `base`、`utils`、`integration` |
| `models` | `base`、`utils`、`integration`、`models` |
| `relaxtime` | `base`、`utils`、`integration`、`relaxtime` |
| `simulation` | `base`、`utils`、`integration`、`simulation` |

两项按职责限定的跨层入口：

- `src/models/workflow_apps/` 可依赖 `relaxtime`，用于工作流编排。
- `src/simulation/` 的服务启动可加载 `src/models/Models.jl`；业务调用经统一入口，不直接加载 solver 内部文件。

现存五条 `src → scripts` 桥接在审计器 `SCRIPT_BRIDGES` 中逐文件列明，报告单列为待迁移依赖：charged GBU 内核与绘图两条、orchestrator 配置和流程三条。它们不构成对其他反向依赖的许可；后续内核迁移需要独立验证数值和产物合同。

## 静态检查的范围

- 支持字面量 include、路径变量、`joinpath`/`normpath`/`dirname`、`@__DIR__`/`@__FILE__` 及绝对路径的 `abspath`。
- 条件分支中的 include 合并展示；无法唯一确定的路径留作未解析项，不执行待分析源码。
- 相对 `using`/`import` 只显示名称，没有解析完整模块身份；宏展开、调用链及被加载脚本内部的依赖仍需人工核查。
- `DEPS_STRICT=1` 拒绝已定位的违规边、缺失文件和解析错误；动态 include 明列并提示覆盖不完整。通过检查不代表完整运行时依赖无环。

规则变更同步修改本页与审计器，并运行 `tests/unit/config/test_gen_deps.jl`、`tests/unit/config/test_analyze_deps.jl`。更新阅读图用 `scripts/dev/gen_deps.jl`；检查当前源码只需 `scripts/dev/analyze_deps.jl`。

## 第三方数值 oracle 的环境边界

- 根 `Project.toml` 只声明 production、稳定 CLI 和常规测试实际需要的依赖。
- QuadGK 不属于根运行时/测试依赖；`src/`、`scripts/`、`tests/` 不得导入或调用它。
- QuadGK 可在 `benchmark/Project.toml` 中作为隔离对照 oracle；必须先独立实例化，并在需要根源码依赖时通过显式 `LOAD_PATH` 叠加，仅对该 benchmark 进程可见。
- 外部 oracle 只提供交叉验证证据，不能替代节点/容差自收敛和 production provenance。
- 机读规则位于 `config/ci/dependency_policy.toml`，门禁为 `scripts/dev/check_dependency_policy.jl`。
- 决策背景和重新评估条件见 [ADR-0006](../decisions/0006-isolate-optional-numerical-oracles.md)。

## Models 入口契约联动

- `Models` 统一求解接口契约见：`docs/architecture/models_solver_contract.md`。
- 新增公开求解入口时，除依赖规则外，还应运行：
  - `julia --project=. scripts/dev/check_models_entry_contract.jl`

## world-age 动态调用边界（迁移期规则）

迁移期允许在 `src/models` 保留少量 `Base.invokelatest(...)` 作为 world-age 安全边界，但必须集中、可审计。

当前必要边界（Phase5-8 后）：
- `src/models/gap_solver.jl`：1 处（legacy adapter 的 `omega` 边界）
- `src/models/factory.jl`：1 处（legacy 构造器边界）
- `src/models/entrypoints.jl`：workflow bridge 边界

Phase5-8 结论（2026-03-03）：
- 旧 PNJL 物理实现已迁移；`src/models/pnjl/` 仅保留 capabilities、API 和 adapter 锚点。
- 物理实现迁移到 `src/models/pnjl_physics/`，不再保留 `module PNJL` 运行时入口。

单一来源：
- 机读清单位于 `config/ci/models_invokelatest_allowlist.toml`。
- 门禁脚本与审计输出以该清单为准，文档仅做解释性说明。

准入规则：
- 不允许新增分散 `invokelatest` 点位；新增必须先更新门禁白名单并附迁移任务证据。
- 相关门禁：`scripts/dev/check_pnjl_migration_guard.jl`（规则5 + `models-invokelatest-audit`）。
- 建议在 PR 中附审计输出：`observed`、`allowlist_baseline`、`allowlisted`。
