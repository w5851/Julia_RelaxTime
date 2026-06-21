# 介子数密度：FixedAsymmetricRho 路径

本文固定 `FixedAsymmetricRho` 接入介子数密度 workflow 的公式口径。

## 1. 路径定位

`FixedAsymmetricRho` 是 density-constrained equilibrium source。它不是新的
介子求解器，也不是新的 branch selection 策略。

组合式介子数密度扫描中的 `--path trho_asymmetric` 执行：

```text
(T, rho_target, asym_ud_ratio_target, asym_s_target)
-> Models.solve(model, FixedAsymmetricRho(...), T)
-> Models.solve_meson_point_from_equilibrium(...)
-> meson density regimes
```

其中第二步使用已经求出的 equilibrium result，不重新求 gap。

## 2. 约束方程

令 `rho_i` 为通过模型压强对 flavor 化学势求导得到的净味密度：

```math
\rho_i = \frac{\partial P}{\partial \mu_i}.
```

当前约束为：

```math
\frac{\rho_u + \rho_d + \rho_s}{3\rho_0} = \rho_{\mathrm{target}},
```

```math
\frac{\rho_u}{\rho_d} = \mathrm{asym\_ud\_ratio\_target},
```

```math
\rho_s = \mathrm{asym\_s\_target}.
```

实现层的约束审计字段必须复用 `Models.model_rho`，即 `model_rho = dP/dmu_i`。
不得用介子数密度后处理中的 `number_densities(...).quark` 或任意 meson
density 结果反推这些约束字段。

### 2.1 `0.876` 默认值的来源

项目内把 `asym_ud_ratio_target = 0.876` 解释为 `rho_Q/rho_B = 0.4`
同位旋不对称口径在 flavor density 约束中的历史默认值，而不是新的拟合参数。

在 quark density 记号下：

```math
\rho_B = \frac{\rho_u + \rho_d + \rho_s}{3},
```

```math
\rho_Q = \frac{2}{3}\rho_u - \frac{1}{3}\rho_d - \frac{1}{3}\rho_s.
```

因此：

```math
\frac{\rho_Q}{\rho_B}
= \frac{2\rho_u - \rho_d - \rho_s}{\rho_u+\rho_d+\rho_s}.
```

当默认同时取 `asym_s_target = 0` 时，令
`r = rho_u/rho_d`，有：

```math
\frac{\rho_Q}{\rho_B} = \frac{2r - 1}{r + 1}.
```

若严格令 `rho_Q/rho_B = 0.4`，则 `r = 0.875`。当前默认值
`0.876` 是沿用组内老 Fortran 路径的历史近似口径，对应
`rho_Q/rho_B \approx 0.40085`。因此文档、脚本和产物中应把
`0.876` 说明为 `rho_Q/rho_B = 0.4` 口径的默认实现，而不要把它
误解为独立可调的物理拟合常数。

## 3. 接入介子数密度

平衡态求解完成后，workflow adapter 将 `SolverResult` 转换为：

```text
quark_params  = (m_u,m_d,m_s, mu_u,mu_d,mu_s)
thermo_params = (T, Phi, Phibar, xi)
meson_results = masses / widths / Mott diagnostics
```

随后 stable、strict BW、current BU、generalized BU 等 density regime 只消费
同一个 meson point。这样可以让 fixed-`mu` 路径与 density-constrained 路径
共享同一套介子后处理，而不在脚本层重新拼装质量、宽度或相移积分。

## 4. 输出诊断字段

`trho_asymmetric` smoke 输出应记录：

- `constraint_mode = FixedAsymmetricRho`
- `rho_target`
- `rho_norm`
- `rho_u_fm3`, `rho_d_fm3`, `rho_s_fm3`
- `rho_u_over_rho_d`
- `asym_ud_ratio_target`, `asym_s_target`
- `constraint_residual_norm`
- `mu_u_MeV`, `mu_d_MeV`, `mu_s_MeV`
- `muB_MeV`, `muQ_MeV`, `muS_MeV`

守恒化学势组合沿用 `TrhoScan` 现有输出口径：

```math
\mu_B = \mu_u + 2\mu_d,\qquad
\mu_Q = \mu_u - \mu_d,\qquad
\mu_S = \mu_d - \mu_s.
```

## 5. 当前边界

当前 `trho_asymmetric` 仅作为 smoke / diagnostic path strategy 接入 combined
meson density scan。它不生产正式高精度产物，不更新 numerical regression
baseline，也不引入 `muB/muQ/muS` fixed-`mu` 新路径。
