---
title: 重构顺序建议
archived: true
original: docs/dev/重构顺序.md
archived_date: 2026-01-19
---

以下为原始内容（保留，以便审阅与历史参考）：

***

重构顺序建议：

先创建目录结构和基础模块（core/、solver/、derivatives/）
提取 Thermodynamics.jl - 从 AnisoGapSolver 中分离热力学计算
实现 ConstraintModes.jl - 定义求解模式类型
实现 Conditions.jl - 重构条件函数
实现 SeedStrategies.jl - 初值策略
实现 ImplicitSolver.jl - 整合 ImplicitDifferentiation
迁移 ThermoDerivatives.jl - 适配新接口
更新 PNJL.jl - 统一导出
测试验证 - 确保行为一致