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

## 变更流程

- 如需新增例外或调整矩阵，请在本文件记录理由与影响范围。
- 依赖图更新后，请运行 `scripts/dev/analyze_deps.jl` 生成依赖审计报告。
