# 初值策略与分支跟踪

本页吸收旧 `SeedStrategies` 与相关求解器页面中的核心说明，统一从 `Models` 公共表面解释 seed policy。

## `SeedStrategy` 家族

当前对外可见的核心策略包括：

- `SeedStrategy`
- `DefaultSeed`
- `MultiSeed`
- `ContinuitySeed`
- `HybridContinuitySeed`
- `PhaseAwareSeed`
- `PhaseAwareContinuitySeed`

配套接口包括：

- `get_seed`
- `get_all_seeds`
- `update!`
- `reset!`
- `set_phase!`

## 各策略的职责

### `DefaultSeed`

提供固定经验初值，适合单点求解或作为其它策略的回退基线。

### `MultiSeed`

显式尝试多个候选初值，适合多值解、分支竞争或需要更稳健选解的区域。它与 `solve_multi` 的关系最紧密。

### `ContinuitySeed`

把“上一个收敛解”作为“下一个点”的初值，适合参数扫描。

### `HybridContinuitySeed`

连续性优先，但在必要时允许更强的回退策略，适合比 `ContinuitySeed` 更复杂的扫描环境。

### `PhaseAwareSeed` 与 `PhaseAwareContinuitySeed`

这两个策略用于相变敏感区：

- `PhaseAwareSeed`：更偏单点或无状态初值选择
- `PhaseAwareContinuitySeed`：更偏扫描，跨相变线时自动切换更合理的默认相种子

## `solve`、`solve_multi` 与 seed 的关系

这部分是旧文档里最值得保留的说明之一：

- 若 `seed_strategy isa MultiSeed`，`solve` 实际上会转向 `solve_multi`
- 若 `seed_strategy isa PhaseAwareContinuitySeed` 且 `bootstrap_multiseed=true` 且当前是首点，求解器会在首点使用多初值自举
- 其它策略则通过 `get_seed` 生成单一初值进入常规求解

这说明 seed policy 不只是“给一个向量”，它会直接改变入口行为。

## `update!` 与 `reset!`

连续性类策略是有状态对象：

- `update!`：在收敛后推进跟踪状态
- `reset!`：切换扫描线、切换温度切片或重新开始时重置状态

如果不理解这一点，就会误把连续性策略当作无状态配置对象使用。

## Seed 常量的定位

诸如 `HADRON_SEED_5`、`QUARK_SEED_5`、`HADRON_SEED_8`、`QUARK_SEED_8` 这类导出常量属于公开但进阶的策略支持对象。它们不应作为首页主入口，但应在导出全集中被完整收录，并在需要时在本主题中说明其定位。