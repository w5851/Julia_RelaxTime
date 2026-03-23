# Gas-Liquid (RMF/Walecka) 核心公式（文献对齐内化版）

## 1. 适用范围与定位

本文档提炼 legacy `Function_Gas_Liquid` 与学位论文第 2.2 节（Walecka 模型）中对当前主线最有迁移价值的核心：

- 质子-中子双组分相对论平均场（RMF）
- 自洽场方程 + 约束方程联立求解
- 压强与高阶导数（用于涨落量）

不保留 legacy 中“脚本即长时运行”的形态，统一转为“模型函数 + workflow 调用”。

## 2. 基本自由度与单位

- 费米子：`p, n`（质子、 中子）
- 介子平均场：`sigma, omega, rho, delta`
- 内部单位：自然单位（`fm^-1`）
- 外部接口若使用 MeV，需显式转换：`x_inv_fm = x_MeV / hbarc_MeV_fm`

注：论文 2.2 节的基线模型是 `sigma-omega-rho`（NLrho）；本文保留 `delta` 作为与仓库现有 Gas-Liquid 壳兼容的扩展通道。

## 3. 有效质量、有效化学势与单粒子能量

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

有效化学势（平均场近似，含等标量与同位旋矢量通道）：

```math
\mu_p^* = \mu_p - g_\omega\omega_0 - g_\rho\rho_{03},
\qquad
\mu_n^* = \mu_n - g_\omega\omega_0 + g_\rho\rho_{03}.
```

在常见符号里也可写成统一形式

```math
\mu_i^*=\mu_i-g_\omega\omega_0-g_\rho\tau_{3,i}\rho_{03},\quad \tau_{3,p}=+1,\;\tau_{3,n}=-1.
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

与论文第 2.2 节一致的非线性写法常写作

```math
U(\sigma)=\frac{1}{3}b\,m_N\,(g_\sigma\sigma)^3+\frac{1}{4}c\,(g_\sigma\sigma)^4,
```

并把该势导数并入 `R_\sigma`。

## 6. 巨热力学势与压强

在 RMF 近似下，可写为

```math
\Omega_{RMF}=U(\sigma)-\frac{1}{2}m_\omega^2\omega_0^2-\frac{1}{2}m_\rho^2\rho_{03}^2-\frac{1}{2}m_\delta^2\delta^2-2T\sum_{i=p,n}\int\frac{d^3p}{(2\pi)^3}\Bigl[\ln(1+e^{-(E_i^*-\mu_i^*)/T})+\ln(1+e^{-(E_i^*+\mu_i^*)/T})\Bigr],
```

```math
P=-\Omega_{RMF}.
```

对应能量密度常写为

```math
\epsilon=U(\sigma)+\frac{1}{2}m_\omega^2\omega_0^2+\frac{1}{2}m_\rho^2\rho_{03}^2+\frac{1}{2}m_\delta^2\delta^2+2\sum_{i=p,n}\int\frac{d^3p}{(2\pi)^3}E_i^*(f_i+\bar f_i).
```

## 7. 约束模式（建议与 Models 统一）

为与当前 `Models.solve_constraint` 生态兼容，建议优先支持：

- `FixedMu`：给定 `(T, mu_p, mu_n)`
- `FixedRho`：给定 `(T, rho_B)`
- `FixedAsymRho`：给定 `(T, rho_B, alpha)`，`alpha=(rho_n-rho_p)/(rho_n+rho_p)`

文献（论文 2.2 节）采用的核物质口径为：

```math
\rho_B=\rho_p+\rho_n,
\qquad
\alpha=\frac{\rho_n-\rho_p}{\rho_B}.
```

兼容说明：若某些历史脚本使用 `\tilde\rho_B=(\rho_p+\rho_n)/2`，应在接口层显式换算，避免把“变量定义差异”混入模型方程。

等价残差可写为：

```math
R_{\rho_B}=\rho_B^{calc}-\rho_B^{target}=0,
\qquad
R_\alpha=\alpha^{calc}-\alpha^{target}=0.
```

## 8. 压强与热力学导数（用于涨落）

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

## 9. 论文 2.2 节逐式对照（Eq. 2.40-2.70）

由于 `thesis.pdf` 的中文字体编码（`Adobe-GB1`）在文本提取时有损，下表采用“方程编号 + 可确认结构 + 仓库映射”的保守对齐方式。

| 论文编号 | 可确认物理结构 | 本仓库文档映射 |
|---|---|---|
| Eq. (2.40) | `L_NLWM = L_N + L_M + L_int - U(sigma)` 总拉氏量分解 | 本文 `## 5` 与 `## 6` |
| Eq. (2.41) | 核子 Dirac 项 | 本文 `M_i^*` 与 `E_i^*` 定义 |
| Eq. (2.42) | `sigma/omega/rho` 介子自由项 | 本文 `R_sigma,R_omega,R_rho` 与 `Omega_RMF` |
| Eq. (2.43) | `g_sigma,g_omega,g_rho` 耦合项（含同位旋） | 本文 `mu_i^*` 统一式 |
| Eq. (2.44) | 非线性标量势 `U(sigma)` 三次+四次项 | 本文 `U(sigma)` 定义 |
| Eq. (2.45)-(2.48) | Euler-Lagrange 场方程组 | 本文残差方程组 `R_x=0` |
| Eq. (2.49)-(2.52) | 平均场替换（均匀静态） | 本文“平均场近似”语义 |
| Eq. (2.53)-(2.58) | 平均场代数方程（时间分量/同位旋分量） | 本文 `omega_0,rho_03` 记号 |
| Eq. (2.59)-(2.61) | 进一步化简后的闭合场方程 | 本文最小闭环实现建议 |
| Eq. (2.63) | 标量密度 `rho_s` | 本文 `rho_{s,i}` |
| Eq. (2.64) | 重子数密度 `rho_B` | 本文 `rho_i, rho_B` |
| Eq. (2.65)-(2.67) | `f, fbar, mu_i^*` 分布与有效化学势 | 本文 `## 3`, `## 4` |
| Eq. (2.68) | `Omega_NLWM`（含势场项 + 热对数积分） | 本文 `Omega_RMF` |
| Eq. (2.69) | 压强表达（等价于 `P=-Omega`） | 本文 `P=-Omega_RMF` |
| Eq. (2.70) | 能量密度表达 | 本文 `epsilon` |

## 10. 2.2 节符号映射（实现口径）

为避免后续实现出现“同式异名”，约定如下映射：

- 论文 `omega_0` -> 文档/代码 `omega`（若需区分分量，显式写 `omega_0`）
- 论文 `rho_03` -> 文档/代码 `rho`（若需区分分量，显式写 `rho_03`）
- 论文同位旋标记 `tau_3` -> 代码中采用 `tau3_p=+1`, `tau3_n=-1`
- 论文 `rho_B`（重子密度）-> 约束层统一 `rho_B`
- 论文 `mu_i^*` -> 代码层保持 `mu_p`, `mu_n` 输入，内部计算 `mu_i^*`

结合 `tmp/公式文档temp.md`，额外补充两条可直接用于实现边界的说明：

- 在标准 NLrho 口径下，`m_p^*=m_n^*=m_N-g_sigma sigma`；若启用 `delta` 通道，则恢复 `p/n` 质量分裂写法。
- 若当前阶段不实现反核子热项裁剪，建议保持 `f_i` 与 `fbar_i` 成对保留，以便在中高温区维持热力学一致性。

## 11. 文献参数对齐（第 2.2 节 / 表 2.3）

按论文文本中可辨识参数，给出一组可用于复现实验曲线量级的基准值：

- `m_N = 939 MeV`
- 非线性势系数：`b = 0.00692`, `c = -0.0048`
- 耦合比参数（按原文常用 `f_i` 记号）：
  - `f_sigma = g_sigma / m_sigma = 10.329`
  - `f_omega = g_omega / m_omega = 5.423`
  - `f_rho = g_rho / m_rho = 3.15`

实现建议：

- 在 `config/models/gas_liquid/` 中保存“原表参数 + 映射后参数（显式 `g_x`, `m_x`）”。
- 在文档或配置中明确这 3 个无量纲项对应 `g^2/m^2` 的哪一组（`sigma/omega/rho`），避免口径漂移。

备注：当前 `tmp/公式文档temp.md` 给出的结构与论文提取结果相互印证，支持将三组参数按 `(f_sigma, f_omega, f_rho)` 顺序定稿。若后续拿到原始 TeX/高清表格，可再做一次符号层复核（不影响当前文档默认口径）。

## 12. 与仓库实现的映射建议

- 公式层：本页
- 代码层（计划）：`src/models/variants/gas_liquid/`
- 测试层（计划）：
  - `tests/unit/models/test_gas_liquid_model.jl`
  - `tests/integration/models/test_gas_liquid_workflow_smoke.jl`

## 13. 迁移注意事项

- 先实现“最小可运行核”，不要一次迁入全部高阶扫描脚本。
- 明确 `rho_B` 归一化口径（不同文献常见 `1/2` 或 `1/3` 差异）。
- 保留对初值敏感性的诊断字段，避免把收敛失败静默吞掉。
- 对论文方程编号保持可追溯（建议在代码注释或 docs 中保留 `Eq.(2.xx)` 标签）。

## 14. 基于当前实现反推的不确定项（必须在参数化阶段收敛）

以下结论来自现有实现：

- `src/models/variants/gas_liquid/core/EquationSet.jl`
- `src/models/variants/gas_liquid/core/Thermodynamics.jl`
- `src/models/variants/gas_liquid/GasLiquidModel.jl`

### 14.1 已确认是“占位实现”而非论文定式的部分

- 当前仅有 `g_sigma`, `g_delta`，尚未引入 `g_omega`, `g_rho`, `m_sigma`, `m_omega`, `m_rho`, `m_delta` 全参数集。
- `solve_equilibrium` 里 `sigma` 由 `_sigma_target` 人工目标函数驱动（`tanh` 形式），不是由 Eq. (2.59)-(2.61) 自洽求解。
- `pressure_density_entropy_energy` 目前是代数近似表达，不是由 `Omega_RMF` 的积分表达（Eq. (2.68)-(2.70)）直接导出。
- `mu_p`, `mu_n` 目前在模型壳里默认取同一个 `muB`，同位旋分裂尚未启用。

### 14.2 可用代码反推的参数映射约束

1) **非线性势参数口径**

- 文档写法：`U(sigma)=1/3 b (g_sigma sigma)^3 + 1/4 c (g_sigma sigma)^4`
- 代码现状：尚未显式实现 `U(sigma)`，因此 `b,c` 的“是否含 `g_sigma` 幂次吸收”仍未定稿。
- 约束建议：落地时二选一并写进配置注释：
  - 方案 A：`b,c` 直接对应论文定义（含 `g_sigma` 幂）
  - 方案 B：重定义为 `b_tilde,c_tilde`（把 `g_sigma` 幂吸收）

2) **三组无量纲耦合比（10.329, 5.423, 3.15）**

- 结合 `tmp/公式文档temp.md` 的整理结果，可按 `(f_sigma, f_omega, f_rho)` 锁定顺序。
- 这三项不应直接写成 `g_x`，应先按 `f_x = g_x/m_x` 入配置，再由 `m_x` 还原 `g_x`。

3) **`rho_03` 与 `delta` 通道并存关系**

- 文档已保留 `delta`（质量分裂）与 `rho_03`（化学势分裂）双通道。
- 代码现状只有 `delta` 参与有效质量，尚未把 `rho_03` 引入 `mu_i^*`。
- 约束建议：后续实现中应保持“质量分裂由 `delta`、化学势分裂由 `rho_03`”的正交职责，不要混并。

### 14.3 参数配置落地前的最小核验清单

- 在 `config/models/gas_liquid/` 新增参数文件时，显式写出：
  - 原文符号
  - 配置键名
  - 单位（`MeV` 或 `fm^-1`）
  - 是否吸收 `g_sigma` 幂次
- 用一个最小回归点验证：`T=50 MeV` 下 `P(mu_B)` 曲线应出现论文描述的 S 形区间（若尚未出现，说明参数映射仍有误）。
