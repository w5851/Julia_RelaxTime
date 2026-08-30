# GasLiquid 主题总览

本页优先回答“如何在 `Models` 入口下使用 gas_liquid 变体”，而不是先展开 RMF 方程细节。

## 何时使用该主题

当你需要以下任一能力时，应从本主题开始：

- 构造气液变体模型对象
- 在给定 `(T, mu)` 下完成气液单点求解并获取热力学量
- 使用 RMF 的固定化学势/固定密度单点和 T-mu/T-rho diagnostic scan
- 在统一 `Models` 工作流下接入 RMF/Walecka 变体分支

## 首选公开入口

优先关注：

- `GasLiquidModel`
- `solve_gas_liquid_point`
- `solve_gas_liquid_rmf_point`
- `run_gas_liquid_tmu_scan` / `run_gas_liquid_trho_scan`

完整导出基线见 [generated/Exports.md](generated/Exports.md)。

## 最短使用路径

```julia
using .Models

model = Models.create_model(:GasLiquid)
result = Models.solve_gas_liquid_point(0.60, 0.20)

rmf = Models.solve_gas_liquid_rmf_point(
    T_MeV=15.0,
    mode=:fixed_rho,
    rho_B_fm3=0.16,
    alpha=0.2,
    profile="DiToro_NLrhoDelta",
)
```

## 输入与单位口径

- `T_fm`、`mu_fm` 使用 `fm^-1`
- 兼容 workflow 输出中 `pressure`、`rho`、`entropy`、`energy` 使用自然单位
- RMF workflow 的 `T_MeV`、`mu_B_MeV`、`mu_3_MeV` 和 row 字段单位见 [RMF 核心公式](../../../../reference/formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md)
- RMF row 默认 `formal_status=:diagnostic_only`，不能直接作为生产结果

## 非首页首选入口

`EquationSet.jl` 与 `Thermodynamics.jl` 中的底层函数属于职责核心层，不建议作为新调用方首选入口。
