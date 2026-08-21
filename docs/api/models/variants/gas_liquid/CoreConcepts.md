# GasLiquid 主题职责核心

本页说明 gas_liquid 为什么应作为 `Models` 的模型变体主题，以及适配层与核心方程层的职责边界。

## 1. 变体定位

gas_liquid 代表 RMF/Walecka 路径下的一条模型变体分支：

- 通过 `GasLiquidModel` 对接 `Models` 抽象接口
- 通过 `GasLiquidWorkflow` 提供单点 workflow 封装
- 核心方程与热力学实现在 `gas_liquid/core/`

因此它应作为 `variants` 主题治理，而不是独立一级流程主题。

## 2. 两层职责边界

- 适配层：`GasLiquidModel`
  - 将气液平衡态映射到 `MeanFieldState` 语义
  - 对齐统一接口（`solve_gap`、`omega_components`、`number_densities`）
- 核心层：`EquationSet` / `Thermodynamics`
  - 承载平衡求解、有效质量与热力学量计算
  - 供适配层和 workflow 复用

RMF 核心状态保存 `(S,D,mu_p_star,mu_n_star)`，其中 `S/D` 是场贡献；
`W/R` 由 `f_omega*rho_B` 和 `f_rho*rho_3` 生成。`f_i` 始终表示
`g_i^2/m_i^2`，旧 `g_i=f_i*m_i` 解释只在历史文档中保留，不是当前 API 语义。

## 3. workflow 的语义边界

`solve_gas_liquid_point` 是对外稳定入口，负责把模型求解和可观测量组织成统一结果结构。

`solve_gas_liquid_rmf_point`、`run_gas_liquid_tmu_scan` 和
`run_gas_liquid_trho_scan` 是 RMF diagnostic workflow；每行包含 solver
status、失败原因、场贡献、核子密度和 quadrature 设置。

新调用方优先使用该入口，而不是直接耦合 core 层函数。

## 4. 维护建议

- 若新增能力面向调用方，应优先在 `Models` 聚合层暴露并更新用户入口文档。
- 若属于方程与数值细节演进，应优先更新职责核心层说明并保持接口语义稳定。
