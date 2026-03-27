# Rotation 主题职责核心

本页说明 rotation 为什么应作为 `Models` 的模型变体主题，以及模型适配层和核心热力学层的职责边界。

## 1. 变体定位

rotation 不是独立的通用求解框架，而是特定物理背景下的模型变体分支：

- 通过 `RotationModel` 对接 `Models` 统一接口
- 通过 `RotationWorkflow` 提供单点 workflow 封装
- 核心数值细节收敛在 `rotation/core/RotationThermo.jl`

因此 rotation 在目录上归入 `docs/api/models/variants/`，而不是与 `solver` 或 `workflows` 并列为一级主题。

## 2. 两层职责边界

- 适配层：`RotationModel`
  - 对齐 `Models` 抽象接口（`solve_gap`、`omega_components`、`number_densities` 等）
  - 管理状态与统一调用口径
- 核心层：`RotationThermo`
  - 承载旋转相关热力学计算与残差公式
  - 供适配层和 workflow 复用

## 3. workflow 的语义边界

`solve_rotation_point` 负责把“模型求解 + 可观测量整理”打包成可直接消费的结果结构，用于脚本与服务层快速接入。

它是面向调用方的稳定入口，不应在新调用方中被核心公式函数替代。

## 4. 维护建议

- 新增旋转能力时，优先先判断是否属于适配层 API，还是核心层实现细节。
- 若是对外入口，先更新本主题的用户入口层；若是核心公式，优先更新职责核心层说明。
