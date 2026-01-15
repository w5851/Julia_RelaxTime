# 条件构建模块 `PNJL.Conditions`

代码位置：`src/pnjl/solver/Conditions.jl`

本模块负责构建能隙方程/约束方程的残差函数，并为 `ImplicitSolver` 提供可用于 `NLsolve` 的 `build_residual!` 等接口。

## 相关类型

- `GapParams`：把积分节点、温度、各向异性参数等打包，便于在残差中复用。

## 相关函数（摘要）

- `gap_conditions(...)`
- `build_conditions(...)`
- `build_residual!(...)`

建议直接从 `PNJL.solve(...)` 入口使用；需要自定义求解流程时再直接调用本模块。
