# GasLiquid RMF/RMFT SOP

## Scope

该 SOP 覆盖 `Models.solve_gas_liquid_rmf_point`、`run_gas_liquid_tmu_scan` 和
`run_gas_liquid_trho_scan`。它不吸收 Fortran 文件、固定宽度输出、高阶有限差分
涨落或 freeze-out 参数化。Fortran 只作为 provenance 和 diagnostic 设计参考。

## Inputs and profiles

- `T_MeV`、`mu_B_MeV`、`mu_3_MeV` 使用 MeV；密度使用 `fm^-3`。
- `f_i` 定义为 `g_i^2/m_i^2`，单位 `fm^2`。
- `alpha=(rho_n-rho_p)/rho_B`，因此 `rho_3=-alpha*rho_B`。
- `DiToro_NLrho`、`DiToro_NLrhoDelta` 和 `Thesis_NLrho` 分开配置；`default`
  仍是 legacy/diagnostic，不自动晋升 formal。

## Point flow

1. 读取 profile，记录 effective config hash 和 source matrix IDs。
2. 使用 `GaussLegendre` 计算粒子/反粒子数密度、标量密度、动能和熵积分。
3. 联合求解 `S,D,mu_p_star,mu_n_star`；线性 omega/rho 方程代数消元为
   `W=f_omega*rho_B`、`R=f_rho*rho_3`，并在 row 中保留独立场贡献和残差证据。
4. 由统一 `Omega` 计算 `P=-Omega`、`s` 和
   `epsilon=Ts+mu_p*rho_p+mu_n*rho_n-P`。
5. 返回 `gas_liquid_rmf_row_v1`。失败点保留 `solver_status`、
   `failure_reason`、`iterations`、`residual_norm` 和 `attempts`。

## T-mu and T-rho

`T-mu` 网格输入温度、重子化学势和可选同位旋化学势；支持 continuation seed。
`T-rho` 网格输入温度、`rho_B` 和 `alpha` 或 `rho_3`，反解 `mu_B/mu_3` 与四场。
每个 scan 返回 rows 与 `gas_liquid_rmf_manifest_v1`，manifest 记录网格、积分截断/节点、
solver 设置、失败点、source IDs 和 output hash 占位字段。

## Validation gates

- 对称基线：`rho_p=rho_n`、`rho_3=rho_s3=D=R=0`、`M_p^*=M_n^*`。
- 非对称 delta：`rho_s3 != 0`、`D=f_delta*rho_s3`、质量分裂和 `R` 符号闭环；
  与 `f_delta=0` 的 NLrho 对照。
- 热力学：`rho_i=partial P/partial mu_i`、`s=partial P/partial T` 和
  `epsilon+P=Ts+sum(mu_i*rho_i)`，并做积分节点/截断 refinement。
- schema：字段、单位、缺失值、失败状态和 manifest/hash 完整。

在这些 gate 与人工 profile 审阅完成前，`formal_status` 必须是
`diagnostic_only`。未收敛、非有限、负有效质量、残差超限或 provenance 不完整的
结果不能进入正式 production，也不能把 Fortran 输出写成 regression truth。
