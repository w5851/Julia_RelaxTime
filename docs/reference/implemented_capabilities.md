# 已实现计算能力与方法追踪清单

> 快照日期：2026-08-23
> 审计基线：`origin/main@1ccf29310fb20c30bcd154f0b4966e25a7565225`；磁场五维 FixedMu
> solver/production adapter 已在该主线，外部参数与数值验证仍按下文的 diagnostic gate 管理。
> 权威聚合入口：[`src/models/Models.jl`](../../src/models/Models.jl) 与 [`src/models/entrypoints.jl`](../../src/models/entrypoints.jl)
> 用途：论文/毕业答辩的方法索引、开发者路线图、数值产物审计入口
> 边界：本页记录“仓库已经实现什么、如何计算、证据到哪一层”，不把 smoke、内部回归或已有产物自动解释为外部物理验证。

## 1. 如何使用本清单

本页是当前实现的路线级总账，不替代各 API、公式和 SOP 页面。阅读时应同时区分四类事实：

1. **实现事实**：源码中存在可调用计算核或 workflow。
2. **接口事实**：能力是否经 `Models` 公共 API 或稳定 CLI 暴露。
3. **数值事实**：是否有 unit、integration、regression 或 validation 证据。
4. **物理结论**：是否已有针对具体参数、网格和文献目标的收敛与外部验证。

本文使用以下状态标签：

| 标签 | 含义 |
|---|---|
| **公共 API** | 由 `Models` 或对应领域模块导出；调用合同有 API 文档/测试。 |
| **稳定 CLI** | 位于 [`docs/guides/scripts/README.md`](../guides/scripts/README.md) 的 Phase A 白名单。 |
| **domain-candidate** | 有专题脚本、SOP 和测试，但尚未进入核心稳定 CLI 白名单。 |
| **diagnostic/research** | 只用于比较、审计、局部可行性或后处理，不能自动替代生产主线。 |
| **internal regression-grade** | 有固定点/路径回归；证明当前仓库口径未漂移，不等于外部正确性。 |
| **external/legacy validated** | 有文献数字化点、外部 oracle 或 legacy 对照；结论只覆盖对应目标和容差。 |
| **formal case** | 某一命名 case 已通过该 case 的收敛、manifest、回归和人工审阅；不能外推为全参数域保证。 |

### 1.1 全局计算骨架

多数路线共享以下上游：

```text
配置/物理常数
-> 构造模型
-> 定义约束模式
-> 求平衡态 x*(T, mu, ...)
-> 计算 Omega、P、rho、s、epsilon
-> 进入相结构 / 导数与涨落 / 介子 / 输运 / 服务编排
-> CSV/JSON/manifest/diagnostics
```

核心热力学合同为

$$
\frac{\partial\Omega}{\partial X_a}=0,
\qquad P=-\Omega,
\qquad \rho_i=-\frac{\partial\Omega}{\partial\mu_i},
\qquad s=-\frac{\partial\Omega}{\partial T},
\qquad \epsilon=Ts+\sum_i\mu_i\rho_i-P.
$$

公式来源与推导骨架见 [`Lagrangian_to_GrandPotential_MFA_Workflow.md`](formula/models/shared/Lagrangian_to_GrandPotential_MFA_Workflow.md)，有效实现由各模型的 `omega`、`gap_residual`、`model_thermo` 和 `number_densities` 方法决定。

## 2. 总览：目前有哪些独立计算路线

| 路线 | 主要输出 | 推荐入口 | 当前证据边界 |
|---|---|---|---|
| A. 模型平衡态单点 | 序参量、质量、`Omega/P/rho/s/epsilon` | `Models.solve_gap`、`Models.solve`、`Models.model_thermo` | PNJL/NJL/RPNJL 主线较完整；变体需看第 13-15 节。 |
| B. FixedMu / T-mu 扫描 | 网格化平衡态 CSV | `Models.run_tmu_scan`；`run_unified_scan.jl` | 公共 API + 稳定 CLI + integration/regression。 |
| C. 非 FixedMu / T-rho 扫描 | 联合求解的 flavor `mu`、密度和热力学量 | `Models.solve_constraint`、`Models.run_trho_scan` | 公共 API + 稳定 CLI；必须保持 8 维联合求解语义。 |
| D. 相结构 | S-shape、spinodal、Maxwell、crossover、CEP、manifest | `Models.run_phase_pipeline` / `run_production_phase_pipeline` | 稳定 CLI；有 regression/validation，但 formal 结论按 case 审计。 |
| E. 热力学/质量导数 | 一、二阶导数，bulk 组合量 | `thermo_derivatives`、`mass_derivatives` | Taylor-series 隐函数主线；固定点回归。 |
| F. 守恒荷涨落 | `chi_BQS`、cumulant、`Ssigma`、`kappa_sigma2` | `Models.chi_BQS` 等；susceptibility CLI | 公共 API + 稳定 CLI；单轴与 mixed 使用不同导数后端。 |
| G. 介子谱与 Mott | 质量、宽度、阈值、Mott gap/line | `solve_gap_and_meson_point`；Mott 专题脚本 | 公共 API；专题 CLI 为 domain-candidate；含外部/legacy validation 子集。 |
| H. 介子数密度 | stable/BW/phase-shift `n_pi,n_K,K/pi` | `Models.solve_*_meson_density_*` | 公共 API + domain-candidate SOP；formal verdict 按 path/regime/case。 |
| I. 介子热力学 | `P_M`、QP/LD、`P_total,s,epsilon,I` | `Models.solve_*_meson_thermo_*` | 公共 API + domain-candidate；当前无 external validation gate。 |
| J. 完整输运 | `sigma(s)`、rate、`tau`、`eta/zeta/sigma` | `solve_gap_and_transport`；transport CLI | 公共 API + 稳定 CLI + 回归/部分外部对照。 |
| K. phase-guided transport | 相线附近两种采样路线的输运产物 | `run_phase_guided_transport_scan.jl` | 专题生产流程；formal case 仅覆盖命名参数集。 |
| L. 外磁场 PNJL | Landau 能级热力学、完整五维 FixedMu 的 `(T,mu,eB)` 扫描 | `Models.run_magnetic_scan` / `scan magnetic` | 公共 API + 稳定 CLI；保留 selected/candidates；FixedRho/phase 明确非目标。数值资格仍按代表点收敛审计。 |
| M. Rotation-PNJL | `(T,mu,omega)` 单点热力学 | `Models.solve_rotation_point` | 可运行最小核 + unit/integration smoke。 |
| N. GasLiquid/RMF | 兼容 `(T,mu)` 单点；diagnostic RMF `T-mu/T-rho` EOS rows | `Models.solve_gas_liquid_point`、`solve_gas_liquid_rmf_point`、`run_gas_liquid_tmu_scan`、`run_gas_liquid_trho_scan` | 核心四场/热力学/契约已实现并有分层 smoke；完整外部 validation/refinement/formal production 仍未完成。 |
| O. 工作流、manifest、HTTP | 声明式 pipeline、后台 job、API JSON | `run_*_pipeline`、`server_full.jl` | 稳定服务入口；适合联调，不代表全部物理域验收。 |

## 3. 路线 A：公共模型与平衡态底座

### 3.1 从模型到平衡态

三味 NJL/PNJL 家族的共同质量关系为

$$
M_i=m_i-4G\phi_i+2K\phi_j\phi_k,
\qquad i\ne j\ne k,
$$

并以平均场驻点条件求 `x*`。PNJL 在 NJL 凝聚、真空和热项之外加入 Polyakov 势：

$$
\Omega=2G\sum_i\phi_i^2-4K\phi_u\phi_d\phi_s
-2N_c\sum_i\int_0^\Lambda\frac{d^3p}{(2\pi)^3}E_i
-2T\sum_i\int\frac{d^3p}{(2\pi)^3}(\mathcal Q_1+\mathcal Q_2)
+\mathcal U(\Phi,\bar\Phi;T).
$$

其中

$$
\mathcal Q_1=\ln\left[1+3\Phi e^{-x_+}+3\bar\Phi e^{-2x_+}+e^{-3x_+}\right],
\quad x_+=\frac{E_i-\mu_i}{T},
$$

$$
\mathcal Q_2=\ln\left[1+3\bar\Phi e^{-x_-}+3\Phi e^{-2x_-}+e^{-3x_-}\right],
\quad x_-=\frac{E_i+\mu_i}{T}.
$$

内部单位为自然单位；`T`、`mu`、质量和动量为 `fm^-1`。对外 MeV 字段通过 `hbar*c` 显式转换。公式主索引：

- NJL：[`NJL_core.md`](formula/models/njl/NJL_core.md)，源为高雪艳博士论文第 2.1 节。
- PNJL：[`Omega_各向同性.md`](formula/models/pnjl/Omega_各向同性.md)。
- RS-PNJL：[`Omega_RS各向异性.md`](formula/models/pnjl/Omega_RS各向异性.md)。
- rPNJL：[`rPNJL_core.md`](formula/models/rpnjl/rPNJL_core.md)，源为贺伟博博士论文第 3.3.1 节。

### 3.2 RS 各向异性进入位置

当前标量 PNJL 热力学采用 Romatschke-Strickland 变形：

$$
E_{\mathrm{dist}}=E_\xi
=\sqrt{p^2+M_i^2+\xi p^2\cos^2\theta}.
$$

它只替换热分布/热对数中的自变量；真空项仍使用

$$
E_{\mathrm{kin}}=\sqrt{p^2+M_i^2}.
$$

这是一项模型约定，不应描述为各向异性色散关系。相图的 `rs_reduced_adaptive` 也只适用于该标量热核，不能外推到 magnetic 或 transport。来源：[`Omega_RS各向异性.md`](formula/models/pnjl/Omega_RS各向异性.md) 与 [`pnjl_phase_structure.md`](../guides/sop/workflows/pnjl_phase_structure.md)。

### 3.3 已注册模型与实际成熟度

模型工厂当前注册 `NJL`、`NJL2`、`PNJL`、`RPNJL`、`PNJLMagnetic`、`Rotation`、`GasLiquid`。同构接口测试要求每个模型暴露 `solve_gap`、`number_densities`、`calculate_mass_vec` 和状态维数。

| 模型 | 状态维数 | 已实现的共同表面 | 特有能力与边界 | 代表证据 |
|---|---:|---|---|---|
| NJL | 3 | 模型构造、gap、Omega、热力学、密度 | 三味 NJL；无 Polyakov 两变量 | `test_njl_core.jl`、`test_njl_omega.jl`、implicit smoke |
| NJL2 | 2 | 同上 | 两味模型；第三味 `mu_s` 仅兼容占位，不应作三味物理解读 | `test_njl2_core.jl`、`test_njl2_omega.jl` |
| PNJL | 5 | 完整公共主链 | 当前相图、导数、介子、输运的 canonical 背景 | unit + integration + regression + validation |
| RPNJL | 5 | 同构 gap/Omega/热力学 | 八夸克项 + Vandermonde 项；有 fixed-point regression | `test_rpnjl_model.jl`、`test_rpnjl_gap_fixedpoint_regression.jl` |
| PNJLMagnetic | 5 | 五维 magnetic equilibrium + magnetic API | Landau 能级、`eB`、`n_max` 收敛；非零 `eB` 的密度为专用净密度，不满足普通独立 `quark/antiquark` 输运合同；不等于普通 PNJL 的 RS 路线 | magnetic unit/regression；见第 14 节 |
| Rotation | 3 | 单点 gap/热力学 workflow | 最小旋转 PNJL 核；见第 15 节 | unit + integration smoke |
| GasLiquid | 4 | 单点 workflow、同构接口、diagnostic RMF core | `f_i`/四场/统一 EOS/T-mu/T-rho 已实现；尚未达到 formal production；见第 16 节 | GasLiquid unit/integration + RMF unit/integration/regression/validation |

代码事实源：[`factory.jl`](../../src/models/factory.jl)、[`test_model_interface_homomorphism.jl`](../../tests/unit/models/test_model_interface_homomorphism.jl)。稳定全功能路线仍以 PNJL 为主；“工厂可构造”不能外推成“所有下游 workflow 对所有模型均物理完备”。

## 4. 路线 B：FixedMu 单点与 T-mu 扫描

### 4.1 单点计算

`FixedMu` 给定 `T` 与 flavor 化学势，未知量仅为模型序参量。对 PNJL：

$$
F_a(x;T,\mu_u,\mu_d,\mu_s,\xi)
=\frac{\partial\Omega}{\partial x_a}=0,
\quad x=(\phi_u,\phi_d,\phi_s,\Phi,\bar\Phi).
$$

求解完成后在同一状态上读取 `Omega/P/rho/s/epsilon`。调用层次：

```text
create_model(kind)
-> solve(model, FixedMu(), T; mu_vec=...)
-> SolverResult + diagnostics
-> model_thermo / number_densities
```

当前实现支持单初值、continuation、phase-aware 和多初值候选治理；需要热力学稳定相时不能只把上一点 continuation 当唯一候选。公共求解合同见 [`ConstraintModes.md`](../api/models/solver/ConstraintModes.md) 与 [`ResultDiagnosticErrorContracts.md`](../api/models/solver/ResultDiagnosticErrorContracts.md)。

### 4.2 T-mu 网格

`Models.run_tmu_scan` 对 `(T,mu,xi)` 网格重复上述单点求解，负责扫描顺序、分支 seed、失败行、resume 和 CSV 合同。其稳定输出包含序参量、质量、压力、密度、熵、能量、迭代数和 residual。

- **公共 API**：`Models.run_tmu_scan`。
- **稳定 CLI**：`scripts/models/run_unified_scan.jl scan tmu`。
- **算法/合同**：[`TmuScan.md`](../api/models/scans/TmuScan.md)。
- **测试**：`tests/unit/models/test_tmu_scan.jl`、`tests/integration/models/test_models_unified_entrypoints_smoke.jl`、`tests/regression/pnjl/test_scan_fixedpoint_regression.jl`。
- **边界**：continuation 路径依赖遍历顺序；失败点写出不等于可用于物理拟合。

## 5. 路线 C：FixedRho、FixedAsymmetricRho、FixedEntropy、FixedSigma 与 T-rho

### 5.1 联合求解语义

非 `FixedMu` 模式不是“先猜 `mu`，再独立解 gap”的分层默认算法。当前主线将 5 个 PNJL 状态量与 3 个 flavor 化学势展平成 8 维未知量：

$$
y=(\phi_u,\phi_d,\phi_s,\Phi,\bar\Phi,\mu_u,\mu_d,\mu_s),
$$

并同时闭合 5 个驻点方程与 3 个约束方程。各模式的物理残差为：

- `FixedRho(rho_target)`：固定重子密度，并施加 equal-mu 约束；核心密度残差可写为
  $$
  R_\rho=\rho_B^{calc}/\rho_0-\rho_{target}=0.
  $$
- `FixedAsymmetricRho(rho_target,ud_ratio_target,s_target)`：
  $$
  R_B=\rho_B^{calc}/\rho_0-\rho_{target}=0,
  \quad R_{ud}=\rho_u/\rho_d-r_{ud}^{target}=0,
  \quad R_s=\rho_s-\rho_s^{target}=0.
  $$
- `FixedEntropy(s_target)`：固定 `s=s_target`，其余化学势关系由 mode 组件定义。
- `FixedSigma(sigma_target)`：固定比熵
  $$
  \sigma=s/n_B=\sigma_{target}.
  $$

实际方程装配以 `build_constraint_components`、`gap_core_residual!` 和 `ProblemSpec` 为准；上式用于说明物理目标，不替代实现中的归一化和 equal-mu 组件。来源：[`ConstraintModes.md`](../api/models/solver/ConstraintModes.md)、[`ConstraintModes.jl`](../../src/models/solver/spec/ConstraintModes.jl)、[`ConstraintComponents.jl`](../../src/models/solver/spec/ConstraintComponents.jl)。

### 5.2 T-rho 扫描

```text
T/rho/xi grid
-> reverse-rho / seed policy
-> 8D joint solve at every point
-> same-config thermo and residual postprocess
-> trho CSV
```

- **公共 API**：`Models.run_trho_scan`、`Models.build_default_rho_grid`。
- **稳定 CLI**：`scripts/models/run_unified_scan.jl scan trho`。
- **算法/合同**：[`TrhoScan.md`](../api/models/scans/TrhoScan.md)。
- **回归**：`test_constraint_fixedpoint_regression.jl`、`test_constraint_selection_regression.jl`、`test_wavee_unified_scan_stability.jl`。
- **边界**：`rho` 分辨率直接影响 S-shape/Maxwell；默认 `reverse_rho=true` 是稳定性合同。`rs_reduced_adaptive` 目前只允许 PNJL + fixed-rho。

## 6. 路线 D：相结构、Maxwell 与 CEP

### 6.1 计算流程

```text
T-rho joint scan
-> group mu(rho) by T
-> S-shape detection
-> spinodal extraction
-> Maxwell equal-area construction
-> first-order boundary
-> temperature-side CEP bracket/refinement
-> optional crossover line
-> artifacts + manifest + optional reviewed promotion
```

### 6.2 S-shape 与 spinodal

离散曲线对 `dmu/drho` 做符号分析。出现

$$
+\ \longrightarrow\ -\ \longrightarrow\ +
$$

时判定存在 S-shape；两个导数换号点分别给出 hadron-side 与 quark-side spinodal。实现还用邻近三点二次插值细化极值，失败时回退到离散点。来源：[`PhaseTransition.md`](../api/models/phase/PhaseTransition.md)。

### 6.3 Maxwell 构造

在唯一三交点拓扑下，寻找共存化学势 `mu*` 使

$$
A(\mu_*)=\int_{\rho_h}^{\rho_q}[\mu(\rho)-\mu_*]d\rho=0.
$$

输出 `mu_transition`、`rho_hadron`、`rho_quark`、`area_residual` 与 `converged`。production 路线只接受带非零宽度面积变号 bracket 的唯一候选；单个近零采样点只记录诊断，不制造候选根。来源：[`Algorithms.md`](../api/models/phase/Algorithms.md) 与 [`PhaseCore.jl`](../../src/models/phase/PhaseCore.jl)。

### 6.4 CEP 与 crossover

CEP 以温度切片状态从 `confirmed_first_order` 到 `confirmed_monotone` 的证据区间收口。公开状态仅为：

- `resolved`：有独立同点 `(T,mu)` 临界估计；
- `ambiguous`：保留 `[T_last_first_order,T_first_monotone]`；
- `not_found`：扫描范围内未建立证据，不能自动解释为物理上不存在 CEP。

普通 bracket 的 `temperature_resolution_target` 是数值搜索目标，不是物理误差条。高温支路可用 `peak` 或 `inflection` 判据寻找 crossover；它不替代 Maxwell 一阶判据。

### 6.5 通用、production 与诊断路线

| 子路线 | 入口 | 用途 | 不能声称什么 |
|---|---|---|---|
| 通用研究 | `Models.run_phase_pipeline` | 插值/direct CEP、crossover、策略比较 | 低网格输出不是正式相图。 |
| production | `Models.run_production_phase_pipeline` | 三态分类、rho/T/xi 几何收敛、工件治理 | 生成工件不等于自动晋升 reference。 |
| PM 双分支 | `Models.analyze_pm_branch_competition` | 同配置 hadron-like/quark-like 分支与 `Delta Omega` 对照 | compare-only，不改变主 pipeline CEP 选择。 |
| cascade/hybrid/endpoint-local | production 内显式策略 | CEP 附近 support 诊断和证书可行性 | 内部证书不能无审阅替代 canonical production。 |

- **稳定 CLI**：`scripts/pnjl/calculate_phase_structure.jl`。
- **SOP**：[`pnjl_phase_structure.md`](../guides/sop/workflows/pnjl_phase_structure.md)。
- **Integration**：`test_phase_pipeline_smoke.jl`、`test_phase_cli_smoke.jl`。
- **Regression**：`tests/regression/phase/` 与 `tests/regression/pnjl/`。
- **Validation**：`tests/validation/pnjl/`，包括 literature targets、reference targets 与 thermal quadrature oracle；每个 target 的适用范围单独解释。

## 7. 路线 E：Taylor-series 热力学与质量导数

### 7.1 方法

当前 PNJL 导数主线不是在每个扰动点重新做全局 gap 求解，而是：

1. 求目标点平衡态 `x0`，或复用 workflow 提供的同一分支状态；
2. 构造 `J0=dF/dx` 并分解；
3. 令 `T(delta)`、`mu(delta)`、`x(delta)` 携带 Taylor series；
4. 逐阶求解
   $$
   F(x(\delta),T(\delta),\mu(\delta))=0;
   $$
5. 从压力、质量或状态量 series 读出导数并检查 series residual。

该路线锁定零阶分支；bulk 所需的 `T`、`mu`、`T+mu` 方向共享零阶状态与线性化。旧 `ForwardDiff + ImplicitDifferentiation` production fallback 已下线。

### 7.2 已实现输出

- `mass_derivatives`：`M`、`dM/dT`、`dM/dmu` 及二阶导数。
- `thermo_derivatives`：`P/epsilon/rho/s` 与一阶导数。
- `dP_dT`、`dP_dmu`：单方向高阶压力导数。
- `bulk_derivative_coeffs`：`(dP/depsilon)_n`、`(dP/dn)_epsilon` 等组合。
- `bulk_viscosity_coefficients`：`v_n^2`、`dmuB/dT|sigma`、质量导数、`c_p` 和分支诊断。

例如

$$
v_n^2=
\frac{s\chi_{\mu\mu}-n_B\chi_{\mu T}}
{T(\chi_{TT}\chi_{\mu\mu}-\chi_{\mu T}^2)}.
$$

公式与 API：[`CoreConcepts.md`](../api/models/derived/derivatives/CoreConcepts.md)、[`ThermoDerivatives.md`](../api/models/derived/derivatives/ThermoDerivatives.md)、[`BulkViscosityDerivatives.md`](../api/models/derived/derivatives/BulkViscosityDerivatives.md)。测试事实源：`tests/unit/pnjl/test_thermo_derivatives.jl`、`tests/regression/pnjl/test_pnjl_thermo_derivatives_regression.jl`。

## 8. 路线 F：守恒荷 susceptibility 与 cumulant

### 8.1 化学势映射

压力以 flavor 化学势为底层变量；守恒荷方向使用

$$
\mu_u=\frac13\mu_B+\frac23\mu_Q,
\quad
\mu_d=\frac13\mu_B-\frac13\mu_Q,
\quad
\mu_s=\frac13\mu_B-\frac13\mu_Q-\mu_S.
$$

总阶数 `n=i+j+k` 的广义 susceptibility 采用主题合同

$$
\chi_{ijk}^{BQS}
=T^{n-4}
\frac{\partial^n P}
{\partial\mu_B^i\partial\mu_Q^j\partial\mu_S^k}.
$$

### 8.2 两条导数路线

- **纯 B/Q/S 单轴**：单变量 TaylorDiff fast path。
- **mixed BQS**：内部 multivariate Taylor jet；变量数为非零轴数量，阶数为 `sum(orders)`。

两者都先求 primal gap、构造 `J0`、逐阶闭合 gap series，再从 pressure series 取系数。显式 `derivative_backend=:forwarddiff` 已退休。

### 8.3 已实现读出

$$
C_{ijk}^{BQS}=VT^3\chi_{ijk}^{BQS},
\qquad S\sigma=\frac{\chi_3^B}{\chi_2^B},
\qquad \kappa\sigma^2=\frac{\chi_4^B}{\chi_2^B}.
$$

公共入口包括 `chi_BQS`、`chi1_B...chi4_B`、Q/S 同构入口、`chi11_BQ/BS/QS`、`cumulant_BQS`、`baryon_Ssigma` 与 `baryon_kappa_sigma2`。

- **稳定 CLI**：`scripts/pnjl/run_conserved_charge_susceptibilities.jl`。
- **API**：[`ConservedChargeSusceptibilities.md`](../api/models/derived/susceptibility/ConservedChargeSusceptibilities.md)、[`Cumulants.md`](../api/models/derived/susceptibility/Cumulants.md)。
- **Regression**：`test_pnjl_chi_b_taylordiff_regression.jl`、`test_pnjl_chi_bqs_mixedjet_regression.jl`、`test_pnjl_conserved_charge_regression.jl`。
- **边界**：`chi2` 接近零时比值需额外稳定性检查；高阶 mixed 成本随 jet 维数和总阶数增长。

## 9. 路线 G：介子质量、宽度与 Mott

### 9.1 一圈积分到极化函数

平衡态先给出 `(m_f,mu_f,T,Phi,PhiBar,xi)`。单线积分在各向同性记号下为

$$
A_f=-4\int_0^\Lambda dp\,\frac{p^2}{E_f}
[1-f_f(E_f-\mu_f)-f_{\bar f}(E_f+\mu_f)].
$$

两线 `B0_ff'` 是有序 flavor 对象，核心能量组合为

$$
\lambda_{f f^{\prime}}=k_0+\mu_f-\mu_{f^{\prime}}.
$$

赝标量/标量极化函数为

$$
\Pi_{f f^{\prime}}^{P,S}
=-\frac{N_c}{8\pi^2}\left[
A_f+A_{f^{\prime}}+\bigl((m_f\mp m_{f^{\prime}})^2-\lambda_{f f^{\prime}}^2+k^2\bigr)B_{0,f f^{\prime}}
\right].
$$

来源：[`OneLoopIntegral_A.md`](formula/relaxtime/integrals/OneLoopIntegral_A.md)、[`OneLoopIntegral_B0.md`](formula/relaxtime/integrals/OneLoopIntegral_B0.md)、[`Polarization_极化函数byB0.md`](formula/relaxtime/polarization/Polarization_极化函数byB0.md)。

### 9.2 RPA 传播子与 pole

简单通道采用

$$
\mathcal D_M(k_0,k)=\frac{2K_M}{1-4K_M\Pi_M(k_0,k)},
$$

零动量 pole 条件为

$$
1-4K_M\Pi_M(k_0=M_M,k=0)=0.
$$

`pi/K/sigma_pi/sigma_K` 走简单通道；`eta/eta_prime` 与标量 partner 走 `0-8` 混合矩阵和 determinant 条件。workflow 会复用同一平衡态，按通道返回 `mass/gamma/converged/residual/threshold/gap`，并维护 mixed-pair identity/continuation，而不是独立按最小根排序。

当前底层极化主接口的 `k0` 仍以实轴数值路径为基础；复 pole 通过显式实/虚方程或上层封装求解。不能仅因结果字段含 `gamma` 就声称所有通道都使用了完整复分析延拓。

### 9.3 Mott 条件

$$
M_M(T,\mu)=M_{q_1}(T,\mu)+M_{q_2}(T,\mu),
\quad
\Delta_M=M_M-(M_{q_1}+M_{q_2})=0.
$$

阈值以上传播子出现虚部，介子由束缚态转为共振态。来源：[`MesonMass_RPA_Pole.md`](formula/relaxtime/propagator/MesonMass_RPA_Pole.md)、[`MottTransition.md`](formula/relaxtime/propagator/MottTransition.md)。

- **公共 API**：`Models.solve_gap_and_meson_point`、`solve_meson_point_from_equilibrium`。
- **专题入口**：`run_gap_meson_mass_scan.jl`、`run_mott_phase_scan.jl`，均非核心稳定 CLI。
- **测试**：meson mass/Mott unit、integration、fixed-point regression 和 literature digitized validation。
- **文献信源**：NJL/RPA 背景 `Hatsuda:1994pi`；SU(3) pole 方程 `Rehberg:1996da`。正式论文引用前仍需从主 BibTeX 核对完整元数据。

## 10. 路线 H：介子数密度的四种物理口径

所有 density 路线共享：

```text
equilibrium
-> meson mass/pole workflow
-> select density regime
-> n_pi, n_K, K/pi + status/metadata
-> path scan and provenance
```

### 10.1 Stable-particle

$$
n_M=d_M\int_0^\infty\frac{dq\,q^2}{2\pi^2}
\frac{1}{\exp[(E_M-\mu_M)/T]-1},
\qquad E_M=\sqrt{q^2+m_M^2}.
$$

聚合通道默认 `d_pi=3`、`d_K=4`；电荷分辨通道 `d=1`。直接来源：`Blaschke:2020bzh` Eq. (23)、`Blaschke:2021yml` Eq. (7)，内部映射见 [`MesonDensity_稳定粒子与KPi比值.md`](formula/relaxtime/meson_density/MesonDensity_稳定粒子与KPi比值.md)。

### 10.2 Reduced strict BW (Stage 1)

$$
n_M^{BW,red}=d_M\int\frac{dq\,q^2}{2\pi^2}
\int d\omega\,g_M(\omega)
\frac{\Gamma_M}{2\pi[(\omega-E_M(q))^2+\Gamma_M^2/4]}.
$$

Stage 1 读取 `q=0` pole 的 `mass/gamma`，用 `E(q)=sqrt(q^2+m^2)` 与常宽度；内层用 `atan` 变量变换稳定解析窄峰。正能量窗口必须满足 `omega_min>mu_M`，跨 Bose pole 的全实轴积分不是默认生产处方。来源：[`MesonDensity_BW过渡验证.md`](formula/relaxtime/meson_density/MesonDensity_BW过渡验证.md)，文献锚点 `Blaschke:2020bzh` Eqs. (26)-(27)。

### 10.3 q-dependent strict BW (Stage 2)

Stage 2 在每个 `q` 点重解

$$
z_p(q)=E_M(q)-i\Gamma_M(q)/2,
\qquad 1-4K_M\Pi_M(z_p,q)=0,
$$

用前一 `q` 点做 continuation，再把 `E_M(q),Gamma_M(q)` 代回 BW 双积分。输出同时保留每个 `q` 的 residual、convergence 与 acceptance，不把失败点静默平滑掉。

### 10.4 Phase-shift / BU

BU 原式为

$$
n_M=d_M\int\frac{d^3q}{(2\pi)^3}
\int\frac{d\omega}{2\pi}g_M(\omega)
\frac{d\delta_M}{d\omega}.
$$

当前数值主线优先使用分部积分后的相移本体形式：

$$
n_M=\frac{d_M}{T}\int\frac{dq\,q^2}{2\pi^2}
\int_0^\infty\frac{d\omega}{2\pi}
g_M(\omega)[1+g_M(\omega)]F(\delta_M).
$$

- `scheme=:current`：`F(delta)=delta`，当前项目默认 production branch。
- `scheme=:gbu_reference`：`F(delta)=delta-sin(2delta)/2`，stricter reference branch。

相移可取 `arg(propagator)` 或审计用 `arg(inverse propagator)`；实轴后端区分 finite-eta 与 `pv_b0_eta0`。`strict_normal_domain` 遇到 `omega<=mu_M` 返回 `NaN/status=:unsafe_bose_domain`；`x_min_cut`、`excitation_only_E_gt_mu` 和 no-anomalous subtraction 都是显式诊断/延拓策略，不能写成唯一文献处方。公式与治理见 [`MesonDensity_BU相移公式.md`](formula/relaxtime/meson_density/MesonDensity_BU相移公式.md) 和 [`MesonDensity.md`](../api/relaxtime/meson_density/MesonDensity.md)。

### 10.5 扫描路径必须与 regime 分开

当前组合 workflow 支持：

- `tmu`：固定 `mu_q` 的温度路径；
- `trho_asymmetric`：复用 `FixedAsymmetricRho` 的联合平衡态；
- 独立 API：freezeout、crossover、external path；
- chemical profiles：flavor 与 charged meson 化学势配置。

freeze-out baseline 使用

$$
T_{fo}(\mu_B)=a-b\mu_B^2-c\mu_B^4,
\qquad
\mu_B(\sqrt{s_{NN}})=\frac{d}{1+e\sqrt{s_{NN}}},
$$

但 pseudo-critical、fixed `mu/T` 或 stitched path 是额外路径规则，不能与 baseline 参数化混称。来源：[`MesonDensity_Freezeout路径参数化.md`](formula/relaxtime/meson_density/MesonDensity_Freezeout路径参数化.md)，直接参数化信源为 `Cleymans:2006qe`。

- **公共 API**：六组 stable/BW/phase-shift point/from-equilibrium 入口。
- **domain-candidate CLI**：`run_combined_meson_density_scan.jl`。
- **SOP**：[`meson_density.md`](../guides/sop/workflows/meson_density.md)。
- **验证**：unit、integration、regime regression、plot-review regression、legacy meson targets、literature digitized K/pi smoke。
- **formal 边界**：现有 `trho_asymmetric` production-grade verdict 只适用于对应 grid/profile/policy/branch repair；不能外推到别的路径或 Bose policy。

## 11. 路线 I：介子热力学与总 EOS

### 11.1 三条 pressure 路线

稳定粒子压强：

$$
P_M=d_M\int_0^\infty\frac{dq\,q^2}{2\pi^2}
T\ln\left(1-e^{-(E_M-\mu_M)/T}\right)^{-1}.
$$

Reduced strict BW 将同一玻色 pressure 核与 Stage-1 Lorentzian 谱权卷积。Phase-shift/GBU 路线为

$$
P_M=d_M\int\frac{d^3q}{(2\pi)^3}
\int_0^\infty\frac{d\omega}{\pi}
\frac{F(\delta_M)}{e^{\omega/T}-1}.
$$

### 11.2 QP 与 Landau damping

令 `s=omega^2-q^2`：

$$
P_M^{LD}=d_M\int\frac{d^3q}{(2\pi)^3}
\int_0^q\frac{d\omega}{\pi}\frac{F(\delta_M)}{e^{\omega/T}-1},
$$

$$
P_M^{QP}=d_M\int\frac{d^3q}{(2\pi)^3}
\int_q^{\omega_{max}}\frac{d\omega}{\pi}\frac{F(\delta_M)}{e^{\omega/T}-1},
\qquad P_M=P_M^{QP}+P_M^{LD}.
$$

LD 的低频热权重很强，因此 `ld_cutoff`、`ld_cutoff_mode` 与 `ld_threshold_mode` 必须进入输出和 manifest。来源：[`MesonThermo_BU_EOS_OffShell_LD.md`](formula/relaxtime/meson_thermo/MesonThermo_BU_EOS_OffShell_LD.md)、[`MesonThermo_QP_LD_Cutoff_Governance.md`](formula/relaxtime/meson_thermo/MesonThermo_QP_LD_Cutoff_Governance.md)，外部主线锚点 `Maslov:2023wul`。

### 11.3 并入总热力学势

$$
P_{total}=P_{MF}+\sum_M P_M,
\qquad
\Omega_{total}=\Omega_{MF}-\sum_M P_M.
$$

`s/epsilon/trace_anomaly` 必须从 `Omega_total -> Models.model_thermo -> AD` 统一读出，不能在脚本层手工拼差分。来源：[`OmegaTotal_并入介子压强后的统一AD热力学流程.md`](formula/models/shared/OmegaTotal_并入介子压强后的统一AD热力学流程.md)。

- **公共 API**：stable、strict BW、phase-shift 的 point/from-meson-point 入口及 `build_meson_thermo_contract_row`。
- **domain-candidate CLI**：`run_phase_shift_meson_thermo_scan.jl`。
- **canonical 当前口径**：`phase_shift_current + pi/sigma_pi + mu_B=0 + xi=0`。
- **内部证据**：point、fixed-point、canonical-path、plot-review regression。
- **明确缺口**：尚无 external literature validation gate；`pi/sigma_pi` 也不是 full nonet。

## 12. 路线 J：从平衡态到完整输运系数

### 12.1 端到端数据流

```text
equilibrium
-> quark_params / thermo_params
-> A_f, B0_ff'
-> Pi_ff'^(P/S)
-> meson propagators
-> scattering |M|^2
-> d sigma/dt
-> sigma_process(s)
-> average rates wbar_ij
-> tau_i
-> eta / zeta / electric sigma
```

逐步输入、单位与对象合同见 [`Transport_EndToEnd_Pipeline.md`](formula/relaxtime/transport/Transport_EndToEnd_Pipeline.md)。canonical workflow 默认把单个 `mu_fm` 扩展为 `mu_u=mu_d=mu_s`；底层 A/B0/polarization 可消费不等 flavor 化学势，但主 transport workflow 不是 charge-resolved kaon workflow。

### 12.2 振幅与截面

同一过程内按允许的 `s/t/u` 道组合标量与赝标量传播子并做色自旋平均。不同散射过程互不相加。微分截面为

$$
\frac{d\sigma}{dt}
=\frac{\overline{|\mathcal M|^2}}
{16\pi\lambda(s,m_i^2,m_j^2)},
$$

总截面为

$$
\sigma_{ij\to cd}(s)
=\int_{t_-}^{t_+}dt\,\frac{d\sigma}{dt}
[1-f_c][1-f_d].
$$

`t_±` 由 Kallen 函数确定；相同末态使用对称因子。来源：[`ScatteringAmplitude_FromTotalPropagator.md`](formula/relaxtime/scattering/ScatteringAmplitude_FromTotalPropagator.md)、[`DifferentialCrossSection_FromScatteringAmplitude.md`](formula/relaxtime/scattering/DifferentialCrossSection_FromScatteringAmplitude.md)、[`TotalCrossSection_FromDifferentialCrossSection.md`](formula/relaxtime/scattering/TotalCrossSection_FromDifferentialCrossSection.md)。

### 12.3 平均散射率与 tau

轴对称 RS 分布下，平均率保留两个动量、两个极角和相对方位角积分：

$$
\bar w_{ij}=\frac{d_q^2}{\rho_i\rho_j(2\pi)^5}
\int dp_i\,p_i^2dp_j\,p_j^2d\cos\theta_i d\cos\theta_j d\phi\,
f_i f_j v_{rel}\sigma_{ij}(s).
$$

当前 cache 主策略是 `w0cdf` 取点 + PCHIP，默认 60 个 `sqrt(s)` 点；这是实现策略，不是新物理公式。弛豫率为

$$
\tau_i^{-1}=\sum_j\rho_j\bar w_{ij}.
$$

实现利用 `u<->d` 同位旋对称返回六分量 `(u,d,s,ubar,dbar,sbar)`。来源：[`AverageScatteringRate_FromCrossSection.md`](formula/relaxtime/scattering/AverageScatteringRate_FromCrossSection.md)、[`RelaxationTime_FromAverageScatteringRate.md`](formula/relaxtime/transport/RelaxationTime_FromAverageScatteringRate.md)。

### 12.4 eta、electric sigma 与 zeta

各向同性 RTA 代表公式：

$$
\eta=\frac{1}{15T}\sum_a\int\frac{d^3p}{(2\pi)^3}
\frac{p^4}{E_{kin,a}^2}\tau_a f_a(1-f_a),
$$

$$
\sigma_{el}=\frac{1}{3T}\sum_a\int\frac{d^3p}{(2\pi)^3}
\frac{p^2q_a^2}{E_{kin,a}^2}\tau_a f_a(1-f_a).
$$

体粘滞采用 Albright-Kapusta 等熵平方核：

$$
\zeta=\frac{1}{9T}\sum_a\int d\Gamma_a
\frac{\tau_a}{E_{kin,a}^2}
\left[p^2+3v_n^2T^2E_{kin,a}
\frac{\partial}{\partial T}\left(\frac{E_{kin,a}-\mu_a}{T}\right)_\sigma
\right]^2 f_a(1+d_af_a).
$$

各向异性时只把 `f -> f_RS(E_dist)` 并恢复显式 `cos(theta)` 积分；运动学分母、`x=E∓mu` 和 bulk `B` 核继续使用 `E_kin=sqrt(p^2+m^2)`。`prefer_energy_aniso` 只改变分布调用路由，不改变运动学核。

公式外部映射：

- eta/sigma 与 RS 能量职责：Thakur et al. (2017), PRD 95, 096009, Eqs. (18), (21), (25)-(26), (37)-(38), DOI `10.1103/PhysRevD.95.096009`。
- zeta：Albright & Kapusta (2016), PRC 93, 014903, Eq. (138)，经典统计 Eq. (111)，DOI `10.1103/PhysRevC.93.014903`。
- `mu_B=0` 交叉核验：Mykhaylova & Sasaki (2021), PRD 103, 014007, Eq. (10)，DOI `10.1103/PhysRevD.103.014007`。
- 备用热力学导数形式：Sasaki & Redlich (2009), PRC 79, 055207, Eq. (2.19)，DOI `10.1103/PhysRevC.79.055207`。

“A26”只是无法映射到外部文献的 legacy 别名，不得在论文或答辩中当作公式信源。完整逐项映射见 [`TransportCoefficients_FromRelaxationTime.md`](formula/relaxtime/transport/TransportCoefficients_FromRelaxationTime.md)。

### 12.5 守恒荷扩散矩阵与热导率

守恒荷输运以 `B/Q/S` 三个驱动方向组织扩散矩阵：

$$
\Delta J_X^i=\sum_{Y\in\{B,Q,S\}}\kappa_{XY}D^i\alpha_Y,
\qquad \alpha_Y=\frac{\mu_Y}{T}.
$$

当前公开 `kappa_BB`、`kappa_BQ`、`kappa_BS`、`kappa_QQ`、`kappa_QS`、`kappa_SS`、`diffusion_matrix` 与由 `kappa_BB` 派生的

$$
\lambda=\kappa_{BB}\left(\frac{\epsilon+P}{n_BT}\right)^2.
$$

公式和接口来源：[`KappaSeries_AndLambda_FromLegacyFortran.md`](formula/relaxtime/transport/KappaSeries_AndLambda_FromLegacyFortran.md)、[`TransportCoefficients.md`](../api/relaxtime/transport/TransportCoefficients.md)。外部理论锚点是 Das, Mishra & Mohapatra (2022), PRD 106, 014013；该文在 HRG 框架推导扩散矩阵，本仓库落地为夸克侧 RTA，因此只能声称扩散矩阵结构与投影口径有来源映射，不能声称模型细节完全等同。

### 12.6 接口与证据

- **公共 API**：`Models.solve_gap_and_transport`、`solve_transport_from_equilibrium`；底层 `transport_coefficients`、`diffusion_coefficient`、`diffusion_matrix` 等。
- **稳定 CLI**：`run_gap_transport_scan.jl`、`run_relaxtime_orchestrator.jl`。
- **SOP**：[`relaxtime_transport.md`](../guides/sop/workflows/relaxtime_transport.md)。
- **Regression**：cross section、tau、transport fixed-point、diffusion matrix、workflow-vs-direct consistency。
- **Validation**：legacy sigma/tau/ratios、扩散/热导率 legacy 对照与 digitized eta/s/sigma/tau targets；只覆盖目标文件所列点和容差。
- **formal 边界**：低节点 pilot、非默认 propagator/cache policy、未解释 failed points 均为 diagnostic-only。

## 13. 路线 K：phase-guided transport

该路线不是新的输运公式，而是把相结构证据转成输运采样计划：

```text
phase reference / direct coexistence anchor
-> build sampling plan near crossover/first-order regions
-> solve stable branch-specific equilibrium
-> reuse the full transport chain
-> convergence + channel diagnostics + plot manifest
```

两条路线必须分开解释：

1. **mode A: fixed-muB-phase-scaled**
   固定 `mu_B`，以 `alpha_T=T/T_phase` 采样；一阶区域默认直接求两分支等热力学势锚点。
2. **mode B: fixed-T-sparse-muB**
   固定 `T`、离散 `mu_B`，沿 `xi` 采样。

一阶共存处严格 `xi=0` 的单一输运值可能不唯一，当前 formal case 用经过热力学与更高节点认证的正/负近零 `xi` 分别表示两侧极限。phase reference 负责指导采样，不自动晋升新的相结构 reference。

- **专题入口**：`scripts/relaxtime/run_phase_guided_transport_scan.jl`。
- **当前 formal case**：`first_canonical_v1_p128_xi001_validated_anchored_prod_v1`；低 `xi` 分辨率收敛锚为 `...xi005...`。具体路径以 [`scripts/README.md`](../guides/scripts/README.md) 的当前登记为准。
- **治理边界**：GitHub workflow 默认产出 diagnostic candidate；formal 资格来自命名 case 的 convergence/import gate，而不是 workflow 成功状态。

## 14. 路线 L：外磁场 PNJL

### 14.1 流程与公式

```text
eB + magnetic config
-> Landau levels and smooth cutoff
-> Omega_vac(B) + Omega_T(B)
-> 5D magnetic stationarity (required by the source model)
-> pressure/densities + n_max convergence report
-> selected state + all converged candidate branches
-> (T,mu,eB) FixedMu production scan
```

Landau 能谱为

$$
E_{f,n}=\sqrt{p_z^2+M_f^2+2n|q_f|B},
\qquad \alpha_n=2-\delta_{n0}.
$$

巨势 flavor 项为

$$
\Omega_f^0=-N_c\frac{|q_f|eB}{2\pi}\sum_n\alpha_n
\int\frac{dp_z}{2\pi}f_\Lambda^2(p_{f,n})E_{f,n},
$$

$$
\Omega_f^T=-T\frac{|q_f|eB}{2\pi}\sum_n\alpha_n
\int\frac{dp_z}{2\pi}(Z_f^++Z_f^-).
$$

可选 IMC 耦合参数化：

$$
G(B)=G_0\frac{1+a\zeta^2+b\zeta^3}{1+c\zeta^2+d\zeta^4},
\qquad \zeta=eB/\Lambda_{QCD}^2.
$$

当前配置的 `a=0.108805` 来自高雪艳博士论文表 5-1；Ferreira 2014/2018 以及
旧 Fortran、`pnjl_mag` 使用 `a=0.0108805`。这不是单位转换差异：内部
`eB_fm2` 与 `Lambda_QCD/hbarc` 在 `zeta` 中使用同一单位，转换因子会抵消。
在 profile 决策冻结前，磁场外部数值验证保持 `diagnostic_only`，不生成
acceptance target，也不把两个参数值静默混用。

磁场巨势还包含

$$
\Omega=\sum_{f=u,d,s}(\Omega_f^0+\Omega_f^T)
+2G(eB)\sum_f\phi_f^2-4K\phi_u\phi_d\phi_s+U(\Phi,\bar\Phi,T),
$$

磁场路线沿用第 2 章的对数 Polyakov 势和动力学质量关系：

$$
\frac{U(\Phi,\bar\Phi;T)}{T^4}
=-\frac{a(T)}{2}\bar\Phi\Phi
+b(T)\ln\left[1-6\bar\Phi\Phi+4(\bar\Phi^3+\Phi^3)-3(\bar\Phi\Phi)^2\right],
$$

$$
a(T)=a_0+a_1\frac{T_0}{T}+a_2\left(\frac{T_0}{T}\right)^2,
\qquad b(T)=b_3\left(\frac{T_0}{T}\right)^3,
\qquad M_i=m_i-4G(eB)\phi_i+2K\phi_j\phi_k.
$$

一般温度的逐味密度按式(5-10)--(5-12)由 Landau level 上的夸克/反夸克 Polyakov 分布积分得到；当 `Phi≈PhiBar≈0` 且温度足够低时，才可使用式(5-13)的三夸克近似。五个平均场变量的候选平衡条件是

$$
\partial_{(\phi_u,\phi_d,\phi_s,\Phi,\bar\Phi)}\Omega=0.
$$

完整密度路线来自论文式(5-10)--(5-12)，低温禁闭近似来自式(5-13)；五变量平衡条件和对数 Polyakov 势来自第 2 章式(2-50)--(2-63)。固定化学势约定为 `mu_u=mu_d=mu_s`。完整公式和逐式来源见 [`PNJL_magnetic_core.md`](formula/models/pnjl_magnetic/PNJL_magnetic_core.md)；开发者逐项审核表见该文档的“开发者审核表：公式到实现的对应关系”一节。

来源核验范围：高雪艳博士论文《强相互作用物质相变与重子数涨落的研究》第 2.2 节印刷页 21--24、第五章第 5.1 节印刷页 65--68。原文有符号/排版不一致（式(2-63)重复 `phi_u`，式(5-10)导数符号与式(2-65)冲突，式(5-11)使用带符号 `q_f B`），所以本页不把这些冲突字符升级为无条件代码合同。IMC 参数还存在来源分歧：当前 Julia 配置为 `a=0.108805`，而 Ferreira 2014/2018、旧 Fortran 和 `pnjl_mag` 为 `a=0.0108805`；Julia 的 `MeV^2 -> fm^-2` 转换在 `zeta=eB/Lambda_QCD^2` 中抵消，未发现隐藏十倍补偿。该 profile 决策尚未冻结，详见 [`magnetic_imc_parameter_provenance_v1.md`](../analysis/historical/legacy/legacy_extraction_v1/magnetic_imc_parameter_provenance_v1.md)。

- **公共 API**：`PNJLMagneticModel`、`solve_magnetic_gap`、`Models.run_magnetic_scan`、`MagneticConfig`、magnetic Omega/pressure/density 与 `magnetic_nmax_convergence_report`。磁场 API 只接受正标量 `eB`：外部生产单位为 `MeV^2`，要求 `eB >= MAGNETIC_EB_MIN_MEV2 = 100`；核心内部单位为 `fm^-2`，要求 `eB_fm2 >= MAGNETIC_EB_MIN_FM2 = 100/hbarc^2`。低于门槛、零值和负值统一抛出 `ArgumentError`，不会路由到普通 PNJL；模型 capability 恒定只承诺磁场专用净密度能力，调用方应使用 `calculate_magnetic_number_densities` 的 `net` 字段。
- **稳定 CLI**：`scripts/models/run_unified_scan.jl scan magnetic` 是 `(T,mu,eB)` 完整五维 FixedMu equilibrium 产线，写出 selected CSV 和 candidates CSV；`mu` 表示共同的 `mu_u=mu_d=mu_s`，外部单位为 MeV/MeV^2。
- **固定态诊断 CLI**：`run_magnetic_point.jl`、`run_magnetic_eb_scan.jl`、`run_magnetic_stability_scan.jl` 仍使用固定 `x_state`，只负责内核、`n_max` 或稳定性诊断，不应解释为 equilibrium 扫描。
- **证据**：magnetic unit、thermodynamics unit、fixed-point regression，以及低节点 `solve_magnetic_gap` stationarity/branch probe；固定点证据覆盖固定 `x_state` 的内核/回归，不等于默认高节点磁场 equilibrium 或全分支全集已验收。Fortran 与 `pnjl_mag` 的外部交叉核对台账和代表点计划见 [`magnetic_external_crosscheck_v1.md`](../analysis/historical/legacy/legacy_extraction_v1/magnetic_external_crosscheck_v1.md) 与同目录 CSV；当前仍等待参数 profile、ensemble、截断和 branch gate。
- **当前实现边界**：受支持正 `eB` 的 `PNJLMagneticModel.solve_gap` 通过 `solve_magnetic_gap` 对磁场 `Omega` 的五维驻点做多 seed 求解，并以候选集合保留分支、残差、`n_max` 和可选 Hessian 稳定性标签；普通 `solve_gap` 在未启用稳定性分类时按已找到候选中的最低 `Omega` 选择一个 convenience state，但 branch-aware API 仍保留全部可行候选。`classify_stability=true` 只启用有限差分 Hessian 诊断/显式研究策略，不是 PNJL 系列模型的默认生产过滤条件；`saddle_or_maximum` 标签不能单独否定一个已收敛驻点。`T_fm` 必须为正，磁场模型只接受 `xi=0`。磁场 `calculate_magnetic_rho` 与 `calculate_magnetic_number_densities` 共用含 `Phi/PhiBar` 的净密度语义，后者的 `net` 和历史 `quark` 字段均表示 `q-qbar`，`antiquark` 明确为 `nothing`；该结果不是普通 PNJL 的独立夸克/反夸克输运输入，且模型 capability 会将通用 `number_densities` 标为不支持；需要磁场密度时应调用专用 API。共享 `solve_constraint`/ProblemSpec 约束链目前显式拒绝 `PNJLMagneticModel`，避免把普通 PNJL residual 误用于磁场；磁场完整平衡态入口是 `solve_magnetic_gap`。
- **统一入口拒绝边界**：普通 `run_tmu_scan`、`run_trho_scan` 以及 phase/Maxwell/CEP pipeline 遇到 `model_kind=:PNJLMagnetic` 会显式报错；当前只承诺磁场 `FixedMu`，不把普通 PNJL 的 `FixedRho` 或 phase pipeline 伪装成磁场实现。
- **数值边界**：生产路径使用预先验证的保守 `n_max`、`p_z_max`、平滑截断和 `eB` 分辨率配置；收敛报告用于代表性极端点和发布前数值审计，不要求每个生产点重复验证。论文的低温占据 `n_max` 估计不等同于真空项的 cutoff-based 求和上限；`p_num/pz_max/rtol` 等为实现参数，不能从论文公式推断。不能套用 RS 标量热核的 quadrature 结论。
- **外场能量边界**：当前 `omega/pressure` 是固定外部磁场背景下的物质巨势/压力，只组装式(5-2)中的夸克 Landau 项、手征项和 Polyakov 势；拉格朗日量中的外部 Maxwell 项不作为独立的 `B^2/2` 能量加入输出。因此这些量不能直接解释为包含电磁场自能的总 EOS，若未来需要该口径必须新增显式、版本化的 Maxwell 选项并重新验证压力、能量和相平衡。
- **压力方向边界**：`calculate_magnetic_pressure=-Omega` 当前只定义固定外部 `B` 背景下的标量物质压力；代码没有磁化强度导数、横向/纵向压力拆分或磁场介质自洽 EOS。不能把该标量直接解释为包含磁化各向异性的完整压力张量。

## 15. 路线 M：Rotation-PNJL

旋转准粒子能量与模权为

$$
\epsilon_n=\sqrt{M^2+p_t^2+p_z^2}-(n+1/2)\omega,
\qquad
\mathcal W_n=J_n^2(p_tr)+J_{n+1}^2(p_tr).
$$

巨势将该能谱和模权代入 PNJL 双对数热核，并加多项式 Polyakov 势。平衡态由

$$
\partial_\phi\Omega_{rot}=
\partial_\Phi\Omega_{rot}=
\partial_{\bar\Phi}\Omega_{rot}=0
$$

确定，且输入应满足因果约束 `omega*r<1`。完整来源：Sun, Xu & Huang, *Quarkyonic phase induced by Rotation*, arXiv:2307.14402v1, Eqs. (12)-(15)；仓库映射见 [`Rotation_PNJL_CoreEquations.md`](formula/models/rotation/Rotation_PNJL_CoreEquations.md)。

- **公共 API**：`Models.create_model(:Rotation)`、`Models.solve_rotation_point`。
- **证据**：`test_rotation_model.jl`、`test_rotation_workflow.jl`、integration smoke。
- **边界**：当前是三变量最小核与单点 workflow；高阶 rotation 导数和系统性外部 validation 尚未列为完成能力。

## 16. 路线 N：GasLiquid / RMF

### 16.1 目标理论路线

完整 RMF/Walecka 目标应计算

$$
M_{p,n}^*=M_N-g_\sigma\sigma\mp g_\delta\delta,
\qquad
\mu_{p,n}^*=\mu_{p,n}-g_\omega\omega_0\mp g_\rho\rho_{03},
$$

并联合求解 `sigma/omega/rho/delta` 场方程，再由 `Omega_RMF` 得到 `P` 与 `epsilon`。公式和论文 Eq. (2.40)-(2.70) 映射见 [`GasLiquid_RMF_CoreEquations.md`](formula/models/gas_liquid/GasLiquid_RMF_CoreEquations.md)。

### 16.2 当前真实实现边界

文档声明（历史状态）曾将本路线标为接口壳；当前源码和测试证据已增加一个
独立的 diagnostic RMF 核心，但这不等价于 formal production：

- `GasLiquidModel` 可构造；
- `solve_gas_liquid_point` 可运行并返回结构化 `pressure/rho/entropy/energy`；
- `GasLiquidCoreParams` 明确 `f_i=g_i^2/m_i^2`，并提供 `DiToro_NLrho`、
  `DiToro_NLrhoDelta`、`Thesis_NLrho` 三个 diagnostic profile；
- `GasLiquidEquationSet` 已实现稳定粒子/反粒子积分、sigma/omega/rho/delta
  场闭环、固定 `T-mu` 与固定 `T-rho` 求解诊断；
- `GasLiquidThermodynamics` 已从统一 Omega 读出 `P/s/epsilon`，熵含粒子和反粒子项；
- `solve_gas_liquid_rmf_point`、`run_gas_liquid_tmu_scan`、
  `run_gas_liquid_trho_scan` 已返回版本化 row/manifest，默认 `diagnostic_only`；
- unit/integration/regression/validation 文件覆盖耦合比、对称基线、最小非对称 delta
  闭环和公式映射。

尚未验证、因此不能声称 formal RMF 生产结果的部分：

- 外部文献数值点的独立复现、完整网格收敛与 refinement 仍未完成；
- `rho_i=partial P/partial mu_i`、`s=partial P/partial T` 的偏导检查尚未形成生产 gate；
- profile 物理范围、饱和性质和 compression 尚未人工审阅并晋升 baseline；
- Fortran 输出没有被转成 regression truth，也不能替代上述验证。

因此本路线当前状态为 **可运行、可审计的 diagnostic RMF core**；所有新结果仍不可用于论文正式数值结论。

## 17. 路线 O：声明式工作流、manifest 与 HTTP/API

### 17.1 声明式 pipeline

- `Models.run_workflow_pipeline(:transport; ...)`：transport workflow facade。
- `Models.run_scan_pipeline(:tmu|:trho; ...)`：扫描 facade。
- `Models.run_relaxtime_orchestrator_pipeline(...)`：relaxtime 子命令编排。
- phase、meson density/thermo、rotation、gas-liquid 也有各自的 `Models` 入口或模块访问器。

pipeline 层负责参数适配、统计、manifest 和工件定位；物理公式仍由底层模型/领域模块定义。职责说明见 [`workflows/CoreConcepts.md`](../api/models/workflows/CoreConcepts.md)。

### 17.2 可复现性产物

稳定/专题 workflow 通常记录：

- `git_commit`；
- Julia/环境信息；
- argv 与 effective config；
- physics/model profile 与 config hash；
- 行级状态、failed points、diagnostics；
- artifacts 列表和 plot manifest。

这些字段使结果可定位到计算实现，但 manifest 中存在一个 SHA 仍不自动证明当前 worktree、外部数据或物理结论有效。

### 17.3 HTTP 服务

`scripts/server/server_full.jl` 当前提供：

- `POST /compute`：2->2 动量/椭球计算；
- `POST /api/modules/pnjl-gap/run`：PNJL 单点；
- `POST /api/modules/transport-point/run`：transport 单点；
- `/api/modules/pnjl-scan/jobs`：T-mu/T-rho 后台 job 的 create/status/result/cancel；
- `/api/modules/script-tasks`：受控脚本目录与后台 job。

路由事实源：[`routing.jl`](../../src/simulation/fullserver/routing.jl)、[`shared.jl`](../../src/simulation/fullserver/shared.jl)。服务有 handler、状态机、路径保护和 job contract 单测；当前定位仍是演示与联调，不代表全部参数域、长任务资源和前端交互已完成产品级验收。

## 18. 验证证据如何解读

| 证据层 | 能证明什么 | 不能证明什么 |
|---|---|---|
| Unit | 单一公式核、参数校验、结构合同未破坏 | 完整计算链在物理目标上正确 |
| Integration smoke | 模块能串联、CLI/API/产物合同可运行 | 网格/积分/求解器已收敛 |
| Regression | 固定点或路径相对仓库 baseline 未漂移 | baseline 本身等于文献或实验真值 |
| Validation | 对指定文献、legacy 或 oracle target 满足声明容差 | 未覆盖参数区也同样有效 |
| Formal case audit | 指定 SHA、配置、网格、策略和产物通过审阅 | 同算法任意新 case 自动获得 formal 资格 |

代表性验证目录：

- `tests/unit/`：确定性计算核和合同；
- `tests/integration/`：workflow/CLI/API 串联；
- `tests/regression/` + `tests/baselines/`：内部基线；
- `tests/validation/pnjl/`：相结构与 quadrature 外部/参考目标；
- `tests/validation/relaxtime/`：介子、截面、tau、输运比值和 RS 能量公式映射；
- `benchmark/`：性能与可选外部数值 oracle，QuadGK 只能在该环境显式启用。

## 19. 公式信源与引用登记

### 19.1 可直接按 DOI/arXiv 复核的核心信源

| 主题 | 信源 | 本仓库用途 |
|---|---|---|
| PNJL Polyakov 势 | Ratti et al., PRD 73, 014019 (2006), DOI `10.1103/PhysRevD.73.014019` | PNJL/rPNJL Polyakov 热核背景；具体参数仍以 config 为准。 |
| Rotation | Sun, Xu, Huang, arXiv:2307.14402v1 | Rotation-PNJL Eqs. (12)-(15)。 |
| transport eta/sigma RS | Thakur et al., PRD 95, 096009 (2017), DOI `10.1103/PhysRevD.95.096009`, arXiv:1703.03142 | RTA eta/sigma 与 RS 分布职责。 |
| transport zeta | Albright & Kapusta, PRC 93, 014903 (2016), DOI `10.1103/PhysRevC.93.014903`, arXiv:1508.02696 | 等熵平方核 Eq. (138)。 |
| zeta 交叉核验 | Mykhaylova & Sasaki, PRD 103, 014007 (2021), DOI `10.1103/PhysRevD.103.014007`, arXiv:2007.06846 | `mu_B=0` Eq. (10)。 |
| zeta 备用形式 | Sasaki & Redlich, PRC 79, 055207 (2009), DOI `10.1103/PhysRevC.79.055207`, arXiv:0806.4745 | 备用热力学导数形式 Eq. (2.19)。 |
| meson RPA | Hatsuda & Kunihiro, Phys. Rept. 247 (1994) 221; `Hatsuda:1994pi` | NJL/RPA 传播子和极化背景。 |
| SU(3) meson pole | Rehberg, Klevansky & Huefner, PRC 53 (1996) 410; `Rehberg:1996da` | `1-4K_M Pi_M=0` 上游信源。 |
| off-shell/LD thermo | Maslov & Blaschke, PRD 107, 094010 (2023), DOI `10.1103/PhysRevD.107.094010`; `Maslov:2023wul` | phase-shift pressure、QP/LD 与 cutoff 方法背景。 |
| freeze-out | `Cleymans:2006qe` | `T_fo(mu_B)`、`mu_B(sqrt(s_NN))` 参数化 baseline。 |

### 19.2 论文引用前必须补齐的书目信息

以下内部公式页有清楚的章节/方程号，但当前仓库没有足以直接生成规范 bibliography 的完整元数据：

- 高雪艳博士论文：NJL、PNJL、magnetic PNJL；
- 贺伟博博士论文：rPNJL、transport 备用公式与散射章节；
- GasLiquid/Walecka 学位论文第 2.2 节；
- `Blaschke:2020bzh`、`Blaschke:2021yml` 等 citekey 的完整题名、卷页和 DOI 应从 `D:/Desktop/paper/bib` 主文献库核对。

因此本页可以作为方法定位和答辩追踪索引，但在提交论文前应把 citekey 替换/连接到主 BibTeX 的已核验条目；不要在本仓库另建一份竞争性的主 bibliography。

## 20. 已知文档/实现债务

1. 部分早期公式页仍保留旧绝对 worktree 链接或“待实现”标题，即使当前模块已经存在；判断实现状态应以 `Models` 入口、当前源码和测试为准。
2. PNJL/NJL/磁场的学位论文信源缺完整书目信息；论文写作前需要补齐。
3. `MesonPropagator` 的旧 API 页对质量/宽度提取的描述落后于独立 `MesonMass` workflow；本页按当前 workflow 事实登记。
4. meson thermo 尚无 external validation gate。
5. GasLiquid 已超出接口壳并具备可审计 diagnostic RMF core，但尚未与 Rotation/magnetic 的成熟度等同；外部数值 validation、refinement 和 formal gate 仍缺失。
6. phase/transport 的 formal 资格与具体 case、SHA、配置和收敛证据绑定，不能从目录名或一次成功运行推断。

## 21. 后续维护规则

每次新增或改变稳定入口、公式语义、数值策略、输出合同或验证等级时，按以下顺序更新：

1. 更新对应源码与 unit/integration/regression/validation。
2. 更新专题 API 或 `docs/reference/formula/` 权威页。
3. 更新稳定 CLI 白名单或 `run_script_catalog.md` 的状态。
4. 在本页修改对应路线的“流程、公式、入口、证据、边界”，不要只在总览表改一个词。
5. 更新页首快照日期与审计 SHA。
6. 运行文档、Models 入口、脚本入口和 SOP 治理检查。

复核清单：

- [ ] 新能力是否区分了公共 API、稳定 CLI 与内部 helper？
- [ ] 每一步是否有代码入口和公式/方法信源？
- [ ] `mu_q`、`mu_B`、flavor `mu` 与单位是否明确？
- [ ] 数值方法变更是否说明非变更约束和回归层？
- [ ] internal regression 是否被误写成 external validation？
- [ ] diagnostic candidate 是否被误写成 formal production？
- [ ] 文献 citekey 是否能在主 bibliography 中解析？
- [ ] manifest 与工件是否绑定计算 SHA、effective config 和状态诊断？
