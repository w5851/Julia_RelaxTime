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
- `pnjl`：PNJL 求解与扫描
- `relaxtime`：弛豫时间与散射链路

允许依赖（✅ 允许 / ❌ 不允许）：

| From \ To | root | utils | integration | simulation | pnjl | relaxtime |
|---|---|---|---|---|---|---|
| root | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| utils | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ |
| integration | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ |
| simulation | ✅ | ✅ | ✅ | ✅ | ❌ | ❌ |
| pnjl | ✅ | ✅ | ✅ | ❌ | ✅ | ⚠️ 仅 workflows |
| relaxtime | ✅ | ✅ | ✅ | ❌ | ❌ | ✅ |

**例外约定**：
- `src/pnjl/workflows/` 允许依赖 `src/relaxtime/`（用于输运流程编排）。

迁移期补充约束（2026-02-24）：
- `src/simulation/fullserver` 新增调用应优先经 `src/models/entrypoints.jl` 进入扫描/工作流，不再新增对 `PNJL.run_*` 的运行时依赖。
- `src/pnjl/PNJL.jl` 视为兼容层入口，允许存量调用逐步迁移，但不作为新功能默认入口。

## 变更流程

- 如需新增例外或调整矩阵，请在本文件记录理由与影响范围。
- 依赖图更新后，请运行 `scripts/dev/analyze_deps.jl` 生成依赖审计报告。

## world-age 动态调用边界（迁移期规则）

迁移期允许在 `src/models` 保留少量 `Base.invokelatest(...)` 作为 world-age 安全边界，但必须集中、可审计。

当前必要边界（基线=7）：
- `src/models/pnjl/PNJLModel.jl`：2 处（`FixedMu` 构造、`solve` 调用）
- `src/models/legacy/LegacyPNJLModel.jl`：2 处（`FixedMu` 构造、`solve` 调用）
- `src/models/gap_solver.jl`：1 处（legacy adapter 的 `omega` 边界）
- `src/models/factory.jl`：1 处（legacy 构造器边界）
- `src/models/entrypoints.jl`：1 处（统一流程入口边界）

单一来源：
- 机读清单位于 `config/ci/models_invokelatest_allowlist.toml`。
- 门禁脚本与审计输出以该清单为准，文档仅做解释性说明。

准入规则：
- 不允许新增分散 `invokelatest` 点位；新增必须先更新门禁白名单并附迁移任务证据。
- 相关门禁：`scripts/dev/check_pnjl_migration_guard.jl`（规则5 + `models-invokelatest-audit`）。
- 建议在 PR 中附审计输出：`observed`、`allowlist_baseline`、`allowlisted`。
