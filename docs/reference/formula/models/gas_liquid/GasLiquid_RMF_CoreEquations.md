# Gas-Liquid (RMF/Walecka) 核心公式（面向 Julia_RelaxTime 内化）

## 1. 适用范围与定位

本文档提炼 legacy `Function_Gas_Liquid` 中对当前主线最有迁移价值的核心：

- 质子-中子双组分相对论平均场（RMF）
- 自洽场方程 + 约束方程联立求解
- 压强与高阶导数（用于涨落量）

不保留 legacy 中“脚本即长时运行”的形态，统一转为“模型函数 + workflow 调用”。

## 2. 基本自由度与单位

- 费米子：`p, n`（质子、 中子）
- 介子平均场：`sigma, omega, rho, delta`
- 内部单位：自然单位（`fm^-1`）
- 外部接口若使用 MeV，需显式转换：`x_inv_fm = x_MeV / hbarc_MeV_fm`

## 3. 有效质量与单粒子能量

设质子/中子有效质量为：

```math
M_p^* = M_N - g_\sigma\sigma - g_\delta\delta,
\qquad
M_n^* = M_N - g_\sigma\sigma + g_\delta\delta.
```

单粒子能量：

```math
E_i^*(p) = \sqrt{p^2 + (M_i^*)^2}, \quad i\in\{p,n\}.
```

有效化学势（平均场近似）：

```math
\mu_p^* = \mu_p - g_\omega\omega - g_\rho\rho,
\qquad
\mu_n^* = \mu_n - g_\omega\omega + g_\rho\rho.
```

## 4. 分布函数与密度

```math
f_i = \frac{1}{\exp((E_i^* - \mu_i^*)/T)+1},
\qquad
\bar f_i = \frac{1}{\exp((E_i^* + \mu_i^*)/T)+1}.
```

数密度与标量密度：

```math
\rho_i = 2\int\frac{d^3p}{(2\pi)^3}(f_i-\bar f_i),
\qquad
\rho_{s,i} = 2\int\frac{d^3p}{(2\pi)^3}\frac{M_i^*}{E_i^*}(f_i+\bar f_i).
```

## 5. 平均场方程（最小闭环）

推荐以“残差联立”形式实现，不在文档里绑定某个求解器：

```math
R_\sigma = m_\sigma^2\sigma + b\,\sigma^2 + c\,\sigma^3 - g_\sigma(\rho_{s,p}+\rho_{s,n}) = 0,
```

```math
R_\delta = m_\delta^2\delta - g_\delta(\rho_{s,p}-\rho_{s,n}) = 0,
```

```math
R_\omega = m_\omega^2\omega - g_\omega(\rho_p+\rho_n) = 0,
```

```math
R_\rho = m_\rho^2\rho - g_\rho(\rho_p-\rho_n) = 0.
```

## 6. 约束模式（建议与 Models 统一）

为与当前 `Models.solve_constraint` 生态兼容，建议优先支持：

- `FixedMu`：给定 `(T, mu_p, mu_n)`
- `FixedRho`：给定 `(T, rho_B)`，其中 `rho_B=(rho_p+rho_n)/2`
- `FixedAsymRho`：给定 `(T, rho_B, alpha)`，`alpha=(rho_n-rho_p)/(rho_n+rho_p)`

等价残差可写为：

```math
R_{\rho_B}=\rho_B^{calc}-\rho_B^{target}=0,
\qquad
R_\alpha=\alpha^{calc}-\alpha^{target}=0.
```

## 7. 压强与热力学导数（用于涨落）

压强定义保持主线口径：

```math
P = -\Omega.
```

迁移阶段建议先保留“原始导数定义”，再决定 AD 或有限差分实现：

```math
\chi_1 = \frac{\partial P}{\partial \mu_B},
\chi_2 = \frac{\partial^2 P}{\partial \mu_B^2},
\chi_3 = \frac{\partial^3 P}{\partial \mu_B^3},
\chi_4 = \frac{\partial^4 P}{\partial \mu_B^4}.
```

并可构造观测比值：

```math
\kappa\sigma^2 = \chi_4/\chi_2,
\qquad
S\sigma = \chi_3/\chi_2.
```

## 8. 与仓库实现的映射建议

- 公式层：本页
- 代码层（计划）：`src/models/variants/gas_liquid/`
- 测试层（计划）：
  - `tests/unit/models/test_gas_liquid_model.jl`
  - `tests/integration/models/test_gas_liquid_workflow_smoke.jl`

## 9. 迁移注意事项

- 先实现“最小可运行核”，不要一次迁入全部高阶扫描脚本。
- 明确 `rho_B` 归一化口径（不同文献常见 `1/2` 或 `1/3` 差异）。
- 保留对初值敏感性的诊断字段，避免把收敛失败静默吞掉。
