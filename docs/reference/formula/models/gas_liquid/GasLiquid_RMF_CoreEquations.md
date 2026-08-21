# GasLiquid RMF/RMFT 核心公式与来源矩阵

本页是当前 Julia RMF 核心的公式、符号、单位和证据边界。实现使用自然单位：能量/动量/质量为 `fm^-1`，密度为 `fm^-3`，`Omega/pressure/energy` 为 `fm^-4`；稳定入口同时提供 MeV-facing workflow。

## 1. Source matrix

| source id | 权威范围 | 当前映射 | 限制 |
|---|---|---|---|
| `thesis_2.2_Eq2.40-2.70` | NLWM/NLrho 拉氏量、非线性 sigma、有限温度分布、Omega/P/epsilon | `GasLiquidEquationSet`、`GasLiquidThermodynamics` | Eq.2.69 漏反粒子熵项；Eq.2.70 的 delta 写法与 NLrho-delta 拉氏量不一致 |
| `DiToro_2006_Eq1-13_Eq21-22` | Di Toro et al., NPA 775 (2006), DOI [10.1016/j.nuclphysa.2006.04.007](https://doi.org/10.1016/j.nuclphysa.2006.04.007)：NLrho-delta 拉氏量、delta 场、对称能和约束 | `DiToro_NLrho`、`DiToro_NLrhoDelta`、非对称 validation | 混合相 Gibbs 条件仅作未来 adapter 参考 |
| `QHDI_nucl-th_9311005_Eq1-4` | [arXiv:nucl-th/9311005](https://arxiv.org/abs/nucl-th/9311005)：线性 Walecka/QHD-I 交叉核验 | sigma/omega 线性极限 | 不含 rho/delta 的完整有限温度 EOS |
| `fortran_archive_gas_liquid` | 实现路线、continuation、诊断字段和失败输出 | provenance/diagnostic 设计 | 不是公式权威或 regression truth |

Fortran 价值矩阵：四场残差结构、`T-mu/T-rho` 分离、continuation seed、字段级
诊断和失败点记录属于“可吸收/可借鉴”，但全部以 Julia 重写；旧文件、固定宽度输出、
Visual Fortran 工程、未审计 `.dat`、高阶有限差分涨落和 freeze-out 参数化属于“不可采用”；
`.obj/.exe/.vs`、停滞/overflow/NaN 输出和旧 mixed/quark phase 脚本仅作历史 provenance。

| 分类 | Fortran 内容 | Julia 处理 |
|---|---|---|
| 可吸收 | 四场残差组织、`RMFT-T-mu`/`RMFT-T-rho` 分离、continuation seed、参数拟合输入 | 重写为 `EquationSet`、`solve_equilibrium` 和 scan workflow |
| 可借鉴 | `S/D/W/R` 场贡献、EOS 分项、失败点记录、低阶物理诊断 | 进入 row/manifest 字段和 validation 设计 |
| 不可采用 | 旧 Fortran 文件/变量名/固定宽度输出、Visual Fortran 工程、未审计 `.dat`、高阶有限差分涨落、重复 freeze-out | 不进入 `src`、不生成 regression truth |
| 仅历史参考 | `.obj/.exe/.vs`、旧 mixed/quark phase 脚本、overflow/NaN/停滞输出 | 只保留外部目录与审计记录指针 |

外部证据指针：`D:\w\fortran-archive-validation\gas-liquid`；历史审计记录：
`D:\Desktop\Julia_RelaxTime-record\fortran_gas_liquid_zip_analysis_2026-08-19.md`。

## 2. 耦合比和场贡献

所有 profile 的参数语义冻结为

```math
f_i \equiv \frac{g_i^2}{m_i^2},\qquad [f_i]={\rm fm^2},
```

因此 `g_i = sqrt(f_i) m_i`，不是 `f_i m_i`。代码入口是 `GasLiquidCoreParams` 与 `coupling_ratio`。`b,c` 采用不吸收额外 `g_sigma` 幂次的 `S=g_sigma sigma` 形式：

```math
U_\sigma(S)=\frac13 b m_N S^3+\frac14 c S^4.
```

场贡献变量为

```math
S=g_\sigma\sigma,\quad D=g_\delta\delta_3,\quad
W=g_\omega\omega_0,\quad R=g_\rho b_0.
```

`GasLiquidState` 中保存 `(S,D,mu_p^*,mu_n^*)`；裸 meson field 不进入结果 contract。代码映射：`src/models/variants/gas_liquid/core/EquationSet.jl`。

求解器将 `W/R` 的线性场方程代数消元为密度函数，因此数值未知量是
`(S,D,mu_p^*,mu_n^*)`；这不是省略 omega/rho 方程，row 仍显式返回
`W/R` 和其密度残差。

## 3. 有效质量、化学势和密度

采用 `tau3_p=+1`、`tau3_n=-1`：

```math
M_p^*=M_N-S-D,\qquad M_n^*=M_N-S+D,
```

```math
\mu_p^*=\mu_p-W-R,\qquad
\mu_n^*=\mu_n-W+R.
```

```math
\rho_i=2\int\frac{d^3k}{(2\pi)^3}(n_i-\bar n_i),
\qquad
\rho_{s,i}=2\int\frac{d^3k}{(2\pi)^3}\frac{M_i^*}{E_i^*}(n_i+\bar n_i),
```

```math
\rho_B=\rho_p+\rho_n,\quad \rho_3=\rho_p-\rho_n,
\quad \rho_s=\rho_{s,p}+\rho_{s,n},\quad
\rho_{s3}=\rho_{s,p}-\rho_{s,n}.
```

动量积分复用 `src/integration/GaussLegendre.jl`；节点、截断和 `p_num` 被写入结果 row。

## 4. 四场方程与 delta 自洽

矢量场的密度方程为

```math
W=f_\omega\rho_B,\qquad R=f_\rho\rho_3.
```

非线性 sigma 和 delta 方程为

```math
S=f_\sigma[\rho_s-bm_NS^2-cS^3],\qquad
D=f_\delta\rho_{s3}.
```

delta 是同位旋标量介子。Julia 中 `g_delta` 只在 legacy 配置边界转换为 `f_delta=(g_delta/m_delta)^2`；方程内部只使用 `D` 和 `f_delta`。因此 `gdelta` 场贡献、裸耦合和 `D` 不可混写。

对 `alpha=(rho_n-rho_p)/rho_B`，`rho_3=-alpha*rho_B`。自洽闭环必须同时满足 `rho_s,p/n`、`M_p/n^*`、`D` 和密度积分；当前最小验证位于 `tests/unit/models/test_gas_liquid_rmf_core.jl`。

## 5. Omega、压强、熵和能量

当前 on-shell 约定为

```math
\Omega=U_\sigma(S)-\frac{S^2}{2f_\sigma}
-\frac{W^2}{2f_\omega}-\frac{R^2}{2f_\rho}-\frac{D^2}{2f_\delta}
-2T\sum_i\int\frac{d^3k}{(2\pi)^3}
[\log(1+e^{-(E_i^*-\mu_i^*)/T})
+\log(1+e^{-(E_i^*+\mu_i^*)/T})].
```

当 `f_delta=0` 时 delta 项为零。统一定义为

```math
P=-\Omega,\qquad s=-\left(\frac{\partial\Omega}{\partial T}\right)_\mu,
\qquad \epsilon=Ts+\mu_p\rho_p+\mu_n\rho_n-P.
```

熵的实现必须含粒子和反粒子两项；对 `x=(E-\mu^*)/T`、`\bar x=(E+\mu^*)/T`，积分 integrand 使用 `log(1+e^{-x})+x n + log(1+e^{-\bar x})+\bar x\bar n`。这修正 thesis Eq.2.69 的遗漏。入口为 `pressure_density_entropy_energy` 与 `omega_components`。

## 6. 约束和 workflow

固定化学势：`T_MeV, mu_B_MeV, mu_3_MeV`，联合求解四个场/有效化学势，输出物理 `mu_p/n`、`rho_B/3` 和 EOS。

固定密度：`T_MeV, rho_B_fm3, alpha` 或 `rho_3_fm3`，联合反解 `mu_B,mu_3` 和四场。入口为 `solve_gas_liquid_rmf_point`、`run_gas_liquid_tmu_scan`、`run_gas_liquid_trho_scan`；旧 `solve_gas_liquid_point(T_fm,mu_fm)` 保持兼容。

## 7. 结果 contract 和 formal gate

row schema 为 `gas_liquid_rmf_row_v1`，至少包含：

`schema_version`, `run_id`, `point_id`, `mode`, `profile`, `T_MeV`, `mu_B_MeV`, `mu_3_MeV`, `mu_p_MeV`, `mu_n_MeV`, `mu_p_star_MeV`, `mu_n_star_MeV`, `rho_B_fm3`, `rho_3_fm3`, `rho_p_fm3`, `rho_n_fm3`, `rho_s_fm3`, `rho_s_p_fm3`, `rho_s_n_fm3`, `M_p_MeV`, `M_n_MeV`, `S_inv_fm`, `D_inv_fm`, `W_inv_fm`, `R_inv_fm`, `omega_fm4`, `pressure_fm4`, `entropy_fm3`, `energy_fm4`, `converged`, `solver_status`, `iterations`, `residual_norm`, `failure_reason`, `quadrature_settings`, `formal_status`。

manifest schema 为 `gas_liquid_rmf_manifest_v1`，记录 `git_sha`、effective config hash、profile、source matrix IDs、grid、solver/quadrature 设置、row count、failed points 和 output hashes。当前所有 workflow 默认 `formal_status=diagnostic_only`。

只有 profile 语义、对称 baseline、非对称 delta、热力学恒等式、全网格收敛、quadrature/solver refinement、schema/manifest/hash 和人工审阅全部通过后，结果才可晋升 `formal`；Fortran 历史输出不能替代任何 gate。

## 8. profile 语义

- `DiToro_NLrho`: `f_rho=0.95 fm^2`, `f_delta=0`。
- `DiToro_NLrhoDelta`: `f_rho=3.150 fm^2`, `f_delta=2.500 fm^2`。
- `Thesis_NLrho`: `f_rho=3.150 fm^2`, `f_delta=0`。
- `default`: legacy/diagnostic profile，仍保留 `f_rho=3.15` 与零 delta，不作为 formal 语义。

## 9. 当前源码与测试证据

| 层 | 入口/证据 | 状态 |
|---|---|---|
| 参数/密度/场 | `core/EquationSet.jl` | 已实现，profile 语义需继续审阅 |
| 热力学 | `core/Thermodynamics.jl` | 已实现统一 Omega/P/s/epsilon；导数 refinement 仍需扩大 |
| 兼容 API | `GasLiquidModel.jl`, `entrypoints.jl` | 已实现 |
| T-mu/T-rho/schema | `workflows/GasLiquidWorkflow.jl` | 已实现 diagnostic workflow |
| unit/integration | `tests/unit/models/test_gas_liquid_rmf_core.jl`, `tests/integration/models/test_gas_liquid_rmf_workflows.jl` | 已覆盖最小闭环 |
| regression | `tests/regression/models/test_gas_liquid_rmf_symmetric_regression.jl` | diagnostic fixed point，非 formal truth |
| validation | `tests/validation/models/test_gas_liquid_rmf_formula_mapping.jl` | 公式/单位映射 smoke |

未完成：外部文献数值点的独立复现、热力学偏导 refinement、完整生产网格和人工 formal promotion。未完成前不得将结果用于正式论文数值结论。
