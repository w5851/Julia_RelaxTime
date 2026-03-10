# 模型与配置合同

本页集中说明 magnetic 主题的模型对象、配置对象与参数模板。

## `PNJLMagneticModel`

定义位于 [src/models/pnjl_physics/PNJLMagneticModel.jl](../../../../src/models/pnjl_physics/PNJLMagneticModel.jl#L1)。

结构：

- `base::PNJLModel`
- `magnetic::MagneticConfig`

它的职责不是重写整套求解链，而是：

- 复用零磁场 PNJL 的稳定 gap 解路径
- 把 magnetic 配置与磁场热力学能力接入 `Models` 聚合表面

## `MagneticIMCParams`

定义位于 [src/models/pnjl_physics/core/MagneticThermodynamics.jl](../../../../src/models/pnjl_physics/core/MagneticThermodynamics.jl#L45)。

字段：

- `a`
- `b`
- `c`
- `d`
- `Lambda_QCD_MeV`

默认构造器：

- `default_imc_params()`

## `MagneticConfig`

定义位于 [src/models/pnjl_physics/core/MagneticThermodynamics.jl](../../../../src/models/pnjl_physics/core/MagneticThermodynamics.jl#L53)。

关键字段：

- `eB_fm2`
- `n_max`
- `p_num`
- `pz_max`
- `cutoff_N`
- `imc`

默认构造器：

- `default_magnetic_config(; eB_fm2=0.0)`

## 单位与向量合同

- `x_state`：5 维，`(phi_u, phi_d, phi_s, Phi, PhiBar)`
- `mu_vec`：3 维，`(mu_u, mu_d, mu_s)`
- `T_fm`、`mu_vec`、质量、动量：`fm^-1`
- `eB_fm2`：`fm^-2`

## 配置模板

推荐配置模板： [config/models/pnjl/magnetic_default.toml](../../../../config/models/pnjl/magnetic_default.toml)

关键段落：

- `[magnetic]`
- `[magnetic.imc]`

典型使用顺序：

1. 用 `default_magnetic_config` 构造最小配置
2. 必要时从 TOML 模板同步 `eB_fm2`、`n_max`、`p_num`、`pz_max`
3. 运行后用 `magnetic_nmax_convergence_report` 检查截断质量