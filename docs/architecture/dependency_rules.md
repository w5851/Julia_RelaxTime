# 依赖规则（目录级）

本规则用于限制 `src/` 目录内模块之间的依赖方向，避免出现“底层依赖上层”。

## 规则形式

- **文档位置**：`docs/architecture/dependency_rules.md`
- **规则粒度**：目录级（`src/<group>/`）
- **核心原则**：允许跨目录，但必须**单向**；违反时需要重构或加入明确的例外说明。

## 分层与允许依赖矩阵

分组说明：
- `root`：`src/` 根目录下的文件（如 `Constants_PNJL.jl`、`QuarkDistribution*.jl`）
- `utils`：通用工具与常量
- `integration`：数值积分相关
- `simulation`：运动学与服务接口
- `models`：QCD 模型实现与求解/扫描/工作流入口（含 `pnjl_physics`）
- `relaxtime`：弛豫时间与散射链路

允许依赖（✅ 允许 / ❌ 不允许）：

| From \ To | root | utils | integration | simulation | models | relaxtime |
|---|---|---|---|---|---|---|
| root | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| utils | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| integration | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ |
| simulation | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| models | ✅ | ✅ | ✅ | ❌ | ✅ | ⚠️ 仅 workflows |
| relaxtime | ✅ | ✅ | ✅ | ❌ | ❌ | ✅ |

**例外约定**：
- `src/models/workflows/` 允许依赖 `src/relaxtime/`（用于输运流程编排）。

迁移期补充约束（2026-02-24）：
- `src/simulation/fullserver` 新增调用应优先经 `src/models/entrypoints.jl` 进入扫描/工作流，不再新增对 `PNJL.run_*` 的运行时依赖。
- `src/models/pnjl_module()` 仅作兼容别名入口，允许存量调用逐步迁移，但不作为新功能默认入口。

## 变更流程

- 如需新增例外或调整矩阵，请在本文件记录理由与影响范围。
- 依赖图更新后，请运行 `scripts/dev/analyze_deps.jl` 生成依赖审计报告。

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
- `src/models/pnjl/` 已下线删除。
- 物理实现迁移到 `src/models/pnjl_physics/`，不再保留 `module PNJL` 运行时入口。

单一来源：
- 机读清单位于 `config/ci/models_invokelatest_allowlist.toml`。
- 门禁脚本与审计输出以该清单为准，文档仅做解释性说明。

准入规则：
- 不允许新增分散 `invokelatest` 点位；新增必须先更新门禁白名单并附迁移任务证据。
- 相关门禁：`scripts/dev/check_pnjl_migration_guard.jl`（规则5 + `models-invokelatest-audit`）。
- 建议在 PR 中附审计输出：`observed`、`allowlist_baseline`、`allowlisted`。
