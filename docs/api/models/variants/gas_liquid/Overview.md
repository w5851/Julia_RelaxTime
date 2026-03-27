# GasLiquid 主题总览

本页优先回答“如何在 `Models` 入口下使用 gas_liquid 变体”，而不是先展开 RMF 方程细节。

## 何时使用该主题

当你需要以下任一能力时，应从本主题开始：

- 构造气液变体模型对象
- 在给定 `(T, mu)` 下完成气液单点求解并获取热力学量
- 在统一 `Models` 工作流下接入 RMF/Walecka 变体分支

## 首选公开入口

优先关注：

- `GasLiquidModel`
- `solve_gas_liquid_point`

完整导出基线见 [generated/Exports.md](generated/Exports.md)。

## 最短使用路径

```julia
using .Models

model = Models.create_model(:GasLiquid)
result = Models.solve_gas_liquid_point(0.60, 0.20)
```

## 输入与单位口径

- `T_fm`、`mu_fm` 使用 `fm^-1`
- workflow 输出中 `pressure`、`rho`、`entropy`、`energy` 为可直接消费的单点结果

## 非首页首选入口

`EquationSet.jl` 与 `Thermodynamics.jl` 中的底层函数属于职责核心层，不建议作为新调用方首选入口。
