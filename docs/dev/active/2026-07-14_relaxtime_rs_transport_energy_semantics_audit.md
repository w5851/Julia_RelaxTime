# RS 输运能量语义审计与实施任务单

创建日期：2026-07-14

状态：PR 1 源码、测试、稳定文档、registry 与经作者批准的 regression baseline 刷新已完成；等待 Draft PR CI/review

基线提交：`ea706548e9167db61e0cb7537bab2d2d4daf4cad`

关联 issue：[Issue #130](https://github.com/w5851/Julia_RelaxTime/issues/130)

## 1. 背景与目标

当前 `src/relaxtime/TransportCoefficients.jl` 在 $\xi\ne0$ 时通过 `_energy_for_kernel` 得到

$$
E_{\mathrm{dist}}
=E_\xi
=\sqrt{p^2+M^2+\xi(p\cos\theta)^2},
$$

并把同一个值同时用于 RS 分布和部分输运积分核。本文采用的目标物理约定则是“各向异性只进入分布函数，不重新定义准粒子色散关系”，因此必须明确拆分

$$
E_{\mathrm{kin}}=\sqrt{p^2+M^2},
\qquad
E_{\mathrm{dist}}=\sqrt{p^2+M^2+\xi(p\cos\theta)^2}.
$$

其中：

- $E_{\mathrm{kin}}$ 是普通在壳准粒子能量，用于 $\eta$、$\sigma$、$\zeta$、$\kappa_{XY}$ 的运动学分母，并用于 $\zeta$ 的色散响应/$B$ 核和 $\kappa_{XY}$ 的 Landau-Lifshitz 投影；
- $E_{\mathrm{dist}}$ 只作为 RS 分布 $f_\xi$ 的变形自变量，并由此进入各公式采用的分布占据核；
- $M=M(T,\mu_B,\xi)$ 仍可来自当前各向异性热力学背景。这里排除的是在普通色散关系上额外加入 $\xi(p\cos\theta)^2$，不是把质量强制换成各向同性热力学解。

作者已于 2026-07-15 确认 5.6 节的完整映射，PR 1 已据此完成源码、测试、稳定文档和 registry 实施。当前 gate 是 review 已量化的 $\xi\ne0$ 数值漂移；在作者批准前不刷新 regression baseline。新的正式 production 仍须等待代码 PR 合并，并通过独立数据/图像 PR 导入。

## 2. 范围与非目标

### 2.1 范围

- 审计 `TransportCoefficients.jl` 中 $\eta$、$\sigma$、等熵平方核 $\zeta$、守恒荷扩散矩阵 $\kappa_{XY}$ 及派生热导率 $\lambda$ 的能量角色。
- 把当前 `E` 的实际含义、目标能量、仓库公式和外部文献逐项对应。
- 设计 $E_{\mathrm{kin}}$ 与 $E_{\mathrm{dist}}$ 的最小拆分路径，并让所有共享 helper 调用方显式选择正确能量角色。
- 规定 unit、integration、regression、validation 各层验证责任。
- 规定旧 production 的保留、外部 registry 标记和后续新 case 晋升顺序。
- 在代码实施阶段修正 `docs/api/integrals/OneLoopIntegralsAniso.md` 中 `1/E_{\rm aniso}` 与实际普通能量分母不一致的问题。

### 2.2 非目标

- 不修改散射过程的 $s$、$v_{\rm rel}$、$t_\pm$ 或普通质量壳外线运动学。
- 不改变初末态 RS 分布、Pauli blocking 或 `propagator_xi_policy = match_thermo`。
- 不改变弛豫时间和截面的物理语义。
- 不在本 issue 中引入纵向、横向或更完整的空间各向异性输运张量分解。本次修正现有角积分后的 $\eta$、$\sigma$、$\zeta$、守恒荷空间扩散矩阵 $\kappa_{XY}$ 及其派生 $\lambda$；守恒荷标签矩阵不等同于空间纵/横向张量响应。
- 不通过 `prefer_energy_aniso=false` 伪装新的物理策略；该字段当前只是分布函数调用路由偏好。
- 本次人工确认同步只更新文档与 Issue/PR 描述，不修改 `src/`、测试基线、production registry、正式 CSV 或图像。
- 不修改 `D:\Desktop\paper` 下的主 bibliography；本仓库只记录已经核验的 citation key 和公式对应关系。

## 3. 当前能力与缺口

### 3.1 已有分层能力

- `AverageScatteringRate.jl` 使用普通 $E_i=\sqrt{p_i^2+m_i^2}$、$E_j=\sqrt{p_j^2+m_j^2}$ 构造 $s$ 和 $v_{\rm rel}$。
- `TotalCrossSection.jl` 的 $t_\pm$ 由普通二体运动学决定；末态普通动量再进入 RS 分布计算 Pauli blocking。
- `TransportCoefficients.jl:127-140` 已同时提供普通 `energy_from_p` 和 RS `energy_from_p_aniso`。
- `distribution_for_species_from_E` 支持把预计算的 $E_{\mathrm{dist}}$ 传给普通分布接口，也支持调用 `*_distribution_aniso`；两条路径对默认 PNJL provider 表达同一个 RS 分布。

### 3.2 当前混用

- `TransportCoefficients.jl:579-585` 的 `_energy_for_kernel` 在 $\xi\ne0$ 且 provider 提供 `energy_from_p_aniso` 时返回 $E_{\mathrm{dist}}$。
- `TransportCoefficients.jl:723-745` 的 `_species_transport_state` 只返回一个 `E`，同时把它交给分布函数和积分核。
- `TransportCoefficients.jl:849-853`、`903-907` 因而分别在 $p^6/E^2$、$p^4q_f^2/E^2$ 中使用 $E_{\mathrm{dist}}$。
- `TransportCoefficients.jl:988-1038` 还把同一个 $E_{\mathrm{dist}}$ 用于 $\zeta$ 的色散导数、$x=E\mp\mu$、$B$ 核和 $p^2/E^2$。

### 3.3 `prefer_energy_aniso` 不能解决该问题

`_energy_for_kernel` 是否返回 $E_{\mathrm{dist}}$ 只取决于 $\xi$ 和 provider 是否具有 `energy_from_p_aniso`，不取决于 `prefer_energy_aniso`。后者只在 `distribution_for_species_from_E` 中选择：

- 把已经得到的 $E_{\mathrm{dist}}$ 传给 `quark_distribution` / `antiquark_distribution`；或
- 调用内部自行计算同一 $E_{\mathrm{dist}}$ 的 `quark_distribution_aniso` / `antiquark_distribution_aniso`。

因此，只修改配置并重新运行不会把输运核运动学能量切换成 $E_{\mathrm{kin}}$。

### 3.4 共享 helper 延伸到 $\kappa_{XY}$ 与 $\lambda$

`_integrate_species_sum` 不只服务 $\eta$ 和 $\sigma$。`TransportCoefficients.jl:1128-1133` 的 `diffusion_coefficient` 也消费当前单一 `E`，并把它用于：

- $p^4/E^2$；
- `_kappa_projection` 中的 $q_X-(n_X/h)E$；
- 分布函数。

仓库公式文档 `KappaSeries_AndLambda_FromLegacyFortran.md` 已将扩散矩阵锚定到 Das et al. (2022) Eq. (55)，并明确使用普通

$$
E_f(p)=\sqrt{p^2+M_f^2}
$$

作为 $p^4/E_f^2$ 和 Landau-Lifshitz 投影 $q_X-n_XE_f/(\epsilon+P)$ 中的能量。因此，作者于 2026-07-15 确认把 `diffusion_coefficient` 和全部 $\kappa_{XY}$ 纳入本 issue：上述运动学与投影能量改用 $E_{\mathrm{kin}}$，分布继续使用 $E_{\mathrm{dist}}$。

`lambda_from_kappa_BB` 不直接进行单粒子动量积分，而是计算

$$
\lambda=\kappa_{BB}\left(\frac{\epsilon+P}{n_BT}\right)^2.
$$

其中函数参数 `energy` 是热力学能量密度 $\epsilon$，不是单粒子 $E_{\mathrm{kin}}$ 或 $E_{\mathrm{dist}}$。$\lambda$ 无需新增单独的 energy policy，但必须由修正后的 $\kappa_{BB}$ 重新计算；Lorenz number、legacy Lorentz 和 Prandtl number 等下游量随之重新派生。

## 4. 公式证据与编号治理

### 4.1 权威公式来源

| 项目公式 | 仓库权威位置 | 外部来源 | 使用边界 |
| --- | --- | --- | --- |
| $\eta$ 与 RS 分布职责分离 | `docs/reference/formula/relaxtime/transport/TransportCoefficients_FromRelaxationTime.md:78-99` | Thakur et al. (2017), Eqs. (18), (21), (25)-(26) | 外部文献的显式各向异性结果是小 $\xi$ 展开；它支持普通 $E_f^2=p^2+M_f^2$ 与 RS 分布的职责分离，但不等同于仓库的完整全角 RS 积分。 |
| $\sigma$ 与 RS 分布职责分离 | 同上 `:104-124` | Thakur et al. (2017), Eqs. (37)-(38) | 运动学分母保持普通 $E_f^2$，RS 变形进入分布。 |
| 等熵平方核 $\zeta$ | 同上 `:129-180` | Albright and Kapusta (2016), Eq. (138)；经典统计版本 Eq. (111) | Eq. (138) 是当前量子统计主公式；Eqs. (54)-(55)、(59)-(61) 给出等熵导数与准粒子色散响应所需定义。 |
| $\mu_B=0$ 的 $\zeta$ 交叉核验 | 同上 `:129-180` | Mykhaylova and Sasaki (2021), Eq. (10) | 支持加号形式及其 $\mu_B=0$ 极限。 |
| 备用体黏滞公式 | 同上备用公式小节 | Sasaki and Redlich (2009), Eq. (2.19) | 对应贺伟博博士论文 Eq. (5.5)，不是本次采用的等熵平方核；博士论文 Eq. (5.6) 是电导率。 |
| $\kappa_{XY}$ 的分母与 Landau-Lifshitz 投影 | `docs/reference/formula/relaxtime/transport/KappaSeries_AndLambda_FromLegacyFortran.md` §3-§4 | Das et al. (2022), Eq. (55) | 外部公式给出普通准粒子能量结构；完整 RS 角积分是仓库采用的分布替换约定，不声称由 Das 2022 直接推导。 |
| $\lambda$ | 同上热导率小节 | legacy Fortran + Landau-Lifshitz 框架关系 | $\lambda$ 由 $\kappa_{BB}$ 与热力学量派生，不含独立的单粒子能量积分。 |

仓库公式文档的符号表明确规定

$$
E_a=\sqrt{p^2+m_a^2}.
$$

### 4.2 “A26”与正负号

- “A26”无法在已核对的外部文献或贺伟博博士论文中找到可验证来源，只能保留为 legacy Fortran/内部别名。
- 公式本身可追溯到 Albright and Kapusta Eq. (138)，经典统计版本为 Eq. (111)。代码、审计和论文不得把“A26”继续当成外部方程号。
- 本地稿件 *Bulk viscosity of quark matter across the QCD phase transitions* 的 Eq. (1) 曾把 $B$ 核第二项印成减号；作者已确认这是文稿错误。目标公式采用

$$
B=p^2+3v_n^2T^2E\,\partial_T[(E\mp\mu)/T]_\sigma,
$$

即第二项前为加号。这与 Albright and Kapusta Eq. (138)、Mykhaylova and Sasaki Eq. (10) 的 $\mu_B=0$ 等价形式及当前 Julia/legacy Fortran 实现一致。

## 5. 公式—实现逐项审计

下表中的“当前能量”描述 `main@ea706548...` 的实际执行语义；“目标能量”描述待作者确认后才允许实施的语义。

### 5.1 公共状态与分布路由

| 代码位置/表达式 | 数学角色 | 当前能量 | 目标能量 | 公式依据 | 结论与理由 |
| --- | --- | --- | --- | --- | --- |
| `energy_from_p`, `:127` | 普通在壳色散 | $E_{\mathrm{kin}}$ | $E_{\mathrm{kin}}$ | 仓库符号表；Thakur et al. 的普通 $E_f^2$ | 保留。质量可来自当前 $\xi$ 热力学背景。 |
| `energy_from_p_aniso`, `:129-131` | RS 分布自变量 | $E_{\mathrm{dist}}$ | $E_{\mathrm{dist}}$ | RS 分布定义 | 保留，但名称和调用方不得暗示它是物理色散关系。 |
| `_energy_for_kernel`, `:579-585` | 当前共用能量 helper | $E_{\mathrm{dist}}$（$\xi\ne0$） | 拆成 $E_{\mathrm{kin}}$ 与 $E_{\mathrm{dist}}$ | 全部目标公式 | 当前名称掩盖了两种职责，是本次最小重构点。 |
| `_species_transport_state`, `:723-745` | 物种积分状态 | 单一 `E` 同时服务分布与核 | 同时返回/提供 $E_{\mathrm{kin}}$、$E_{\mathrm{dist}}$ | 全部目标公式 | 调用方必须按数学角色显式取值。 |
| `distribution_for_species_from_E`, `:588-625` | $f_\xi$ | $E_{\mathrm{dist}}$ | $E_{\mathrm{dist}}$ | RS 分布；Thakur et al. Eqs. (21), (25)-(26), (37) | 保留两条 provider 路由的等价分布语义。 |

### 5.2 剪切黏度 $\eta$

| 代码位置/表达式 | 数学角色 | 当前能量 | 目标能量 | 仓库公式 | 外部文献 | 结论与理由 |
| --- | --- | --- | --- | --- | --- | --- |
| `p6 / (E * E)`, `:849-853` | 速度/运动学核 $p^6/E^2$ | $E_{\mathrm{dist}}$ | $E_{\mathrm{kin}}$ | $\eta$ 各向异性公式 `:93-99` | Thakur et al. Eqs. (18), (25)-(26) | 各向异性不重定义色散；分母必须是普通在壳能量。 |
| `ff = f(1-f)`, `:741-742` | 占据与 Pauli 因子 | 由 $E_{\mathrm{dist}}$ 得到 | 由 $E_{\mathrm{dist}}$ 得到 | 同上 | Thakur et al. Eqs. (21), (25)-(26) | RS 变形保留在分布中。 |
| `τ`, `degeneracy`, 角测度 | 弛豫时间、简并和全角积分 | 不涉及本次能量拆分 | 不变 | 同上 | 同上 | 本 issue 不改变这些因子。 |

### 5.3 电导率 $\sigma$

| 代码位置/表达式 | 数学角色 | 当前能量 | 目标能量 | 仓库公式 | 外部文献 | 结论与理由 |
| --- | --- | --- | --- | --- | --- | --- |
| `p4 * q2 / (E * E)`, `:903-907` | 速度/运动学核 $p^4q_f^2/E^2$ | $E_{\mathrm{dist}}$ | $E_{\mathrm{kin}}$ | $\sigma$ 各向异性公式 `:119-124` | Thakur et al. Eqs. (37)-(38) | 分母使用普通在壳能量。 |
| `ff = f(1-f)`, `:741-742` | 占据与 Pauli 因子 | 由 $E_{\mathrm{dist}}$ 得到 | 由 $E_{\mathrm{dist}}$ 得到 | 同上 | Thakur et al. Eq. (37) | RS 变形保留在分布中。 |
| $q_f^2$、`τ`、`degeneracy`、角测度 | 电荷、弛豫时间、简并和积分 | 不涉及本次能量拆分 | 不变 | 同上 | 同上 | 本 issue 不改变这些因子。 |

### 5.4 等熵平方核体黏滞 $\zeta$

| 代码位置/表达式 | 数学角色 | 当前能量 | 目标能量 | 仓库公式位置 | 外部文献/方程号 | 结论与理由 |
| --- | --- | --- | --- | --- | --- | --- |
| `E = _energy_for_kernel(...)`, `:989` | `compute_B` 的共用能量 | $E_{\mathrm{dist}}$ | $E_{\mathrm{kin}}$ | $B_a$、$B_{\bar a}$，`:151-173` | Albright-Kapusta Eqs. (111), (138) | $B$ 中的 $E$ 是准粒子色散能量。 |
| `dE_dT = (m/E) * dM_dT`, `:992` | 固定 $p$ 下色散/质量温度响应 | $\frac{M}{E_{\mathrm{dist}}}\frac{dM}{dT}$ | $\frac{M}{E_{\mathrm{kin}}}\frac{dM}{dT}$ | 能量导数与 $B$ 核 | Albright-Kapusta Eqs. (59)-(61), (138) | 对 $E_{\mathrm{kin}}=\sqrt{p^2+M^2}$ 求导；$dM/dT$ 可包含当前 $\xi$ 背景响应。 |
| `dE_dμB = (m/E) * dM_dμB`, `:993` | 固定 $p$ 下色散/质量化学势响应 | $\frac{M}{E_{\mathrm{dist}}}\frac{dM}{d\mu_B}$ | $\frac{M}{E_{\mathrm{kin}}}\frac{dM}{d\mu_B}$ | 能量导数与等熵链式法则 | Albright-Kapusta Eqs. (54), (59)-(61), (138) | 与温度导数同理。 |
| `x = E - μ`, `:1003-1005` | 夸克等熵导数中的 $E-\mu_q$ | $E_{\mathrm{dist}}-\mu_q$ | $E_{\mathrm{kin}}-\mu_q$ | $B_a$，`:157` | Albright-Kapusta Eqs. (54), (61), (138) | 此处属于准粒子热力学组合，不是 RS 分布自变量。 |
| `x = E + μ`, `:999-1001` | 反夸克等熵导数中的 $E+\mu_q$ | $E_{\mathrm{dist}}+\mu_q$ | $E_{\mathrm{kin}}+\mu_q$ | $B_{\bar a}$，`:161` | Albright-Kapusta Eqs. (54), (61), (138) | 同上。 |
| `dx_dT_sigma`, `:1000,1004` | $dE/dT+[dE/d\mu_B\pm b_q](d\mu_B/dT)_\sigma$ | 继承 $E_{\mathrm{dist}}$ 导数 | 继承 $E_{\mathrm{kin}}$ 导数 | 等熵导数展开 | Albright-Kapusta Eqs. (54)-(55), (59)-(61), (138) | $b_q=1/3$ 和 $(d\mu_B/dT)_\sigma$ 不变，只替换色散响应。 |
| `dxt_dT_sigma = dx/T - x/T^2`, `:1008-1009` | $\partial_T[(E\mp\mu)/T]_\sigma$ | 混合使用 $E_{\mathrm{dist}}$ | 全部使用 $E_{\mathrm{kin}}$ | $B_a$、$B_{\bar a}$ | Albright-Kapusta Eq. (138) | 商法则不变，输入的 $x$ 与导数改用普通色散。 |
| `E * dxt_dT_sigma`, `:1011-1012` | $B$ 核第二项 | 外乘和导数均为 $E_{\mathrm{dist}}$ 语义 | 外乘和导数均为 $E_{\mathrm{kin}}$ 语义 | $B_a$、$B_{\bar a}$ | Albright-Kapusta Eqs. (111), (138)；Mykhaylova-Sasaki Eq. (10) | 第二项前保持加号；禁止只替换外乘 $E$ 造成新的混用。 |
| `E = _energy_for_kernel(...)`, `:1022` | 单味积分状态 | $E_{\mathrm{dist}}$ | 同时计算 $E_{\mathrm{kin}}$、$E_{\mathrm{dist}}$ | $\zeta$ 各向异性公式，`:168-180` | Albright-Kapusta Eq. (138) + 仓库完整 RS 约定 | 两种能量在该层显式分流。 |
| `distribution_for_species_from_E(..., E, ...)`, `:1024-1025` | $f_q$、$f_{\bar q}$ | $E_{\mathrm{dist}}$ | $E_{\mathrm{dist}}$ | $\zeta$ 各向异性公式 | RS 分布定义 | 分布继续使用 RS 变形自变量。 |
| `ff_q`, `ff_aq`, `:1027-1028` | $f_\xi(1-f_\xi)$ | 由 $E_{\mathrm{dist}}$ 得到 | 由 $E_{\mathrm{dist}}$ 得到 | 同上 | Albright-Kapusta Eq. (138) 的量子统计因子 + 仓库 RS 替换 | 保留。 |
| `p*p / (E*E)`, `:1036-1038` | 运动学核 $p^2/E^2$ | $E_{\mathrm{dist}}$ | $E_{\mathrm{kin}}$ | $\zeta$ 主公式 `:151-173` | Albright-Kapusta Eq. (138) | 普通在壳能量。 |
| `prefactor` 与各向异性 `/2`, `:980-983,1051-1062` | 径向/角积分归一化 | 不涉及能量选择 | 不变 | $\zeta$ 积分测度 | 仓库完整全角 RS 约定 | 本 issue 不调整归一化。 |

### 5.5 守恒荷扩散矩阵 $\kappa_{XY}$ 与热导率 $\lambda$

| 代码位置/表达式 | 数学角色 | 当前能量 | 目标能量 | 仓库公式位置 | 外部文献/方程号 | 结论与理由 |
| --- | --- | --- | --- | --- | --- | --- |
| `_integrate_species_sum(...)`, `:1128` | 扩散矩阵积分状态 | 单一 `E` 为 $E_{\mathrm{dist}}$（$\xi\ne0$） | 同时提供 $E_{\mathrm{kin}}$、$E_{\mathrm{dist}}$ | $\kappa_{qq'}$ 一般公式 | Das et al. (2022) Eq. (55) | 与 $\eta$、$\sigma$ 共用拆分后的状态，不保留旧混用。 |
| `p4 / (E * E)`, `:1132` | 运动学核 $p^4/E^2$ | $E_{\mathrm{dist}}$ | $E_{\mathrm{kin}}$ | `KappaSeries...` §3-§4 | Das et al. Eq. (55) | 普通在壳准粒子能量。 |
| `_kappa_projection`, `:1087-1103` | $q_X-n_XE/(\epsilon+P)$ | $E_{\mathrm{dist}}$ | $E_{\mathrm{kin}}$ | `KappaSeries...` §3-§4 | Das et al. Eq. (55) | Landau-Lifshitz 投影中的 $E$ 是普通准粒子色散。 |
| `f`, `:1128,1132` | 扩散占据核 | 由 $E_{\mathrm{dist}}$ 得到 | 由 $E_{\mathrm{dist}}$ 得到 | `KappaSeries...` §3.1 + 仓库 RS 约定 | Das et al. Eq. (55) 的 $f_a^{(0)}$；RS 扩展为项目约定 | 各向异性只进入分布。 |
| `kappa_BB/BQ/BS/QQ/QS/SS`, `:1216-1360` | 守恒荷扩散矩阵元素 | 全部继承旧单一 $E$ | 全部继承上述 $E_{\mathrm{kin}}/E_{\mathrm{dist}}$ 分工 | 扩散矩阵一般公式 | Das et al. Eq. (55) | 对角元和非对角元使用同一能量语义。 |
| `lambda_from_kappa_BB`, `:1458-1481` | $\lambda=\kappa_{BB}[(\epsilon+P)/(n_BT)]^2$ | 间接继承旧 $\kappa_{BB}$ | 由修正后的 $\kappa_{BB}$ 重算 | `KappaSeries...` 热导率小节 | legacy Fortran / Landau-Lifshitz 关系 | `energy` 是热力学能量密度 $\epsilon$；不新增单粒子能量选择。 |
| Lorenz / legacy Lorentz / Prandtl | 下游派生诊断 | 间接继承旧 $\sigma$、$\kappa_{BB}$、$\lambda$ | 由修正后的上游量重算 | `KappaSeries...` 派生量小节 | 定义式 | 不单独引入 energy policy。 |

### 5.6 已由作者确认的审计结论

作者确认的完整映射如下：

1. $\eta$、$\sigma$ 的显式能量分母使用 $E_{\mathrm{kin}}$，而 $f_\xi(1-f_\xi)$ 使用 $E_{\mathrm{dist}}$。
2. $\zeta$ 中除分布占据外的全部显式能量角色——$p^2/E^2$、$x=E\mp\mu$、$dE/dT$、$dE/d\mu_B$、$E\partial_T[(E\mp\mu)/T]_\sigma$——都使用 $E_{\mathrm{kin}}$。
3. $M$、$dM/dT$、$dM/d\mu_B$ 可来自当前固定 $\xi$ 的各向异性热力学解；这不等同于采用 $E_{\mathrm{dist}}$ 作为色散。
4. $B$ 核第二项保持加号。
5. `diffusion_coefficient` 和全部 $\kappa_{XY}$ 的 $p^4/E^2$ 与 Landau-Lifshitz 投影使用 $E_{\mathrm{kin}}$，占据分布使用 $E_{\mathrm{dist}}$；$\lambda$ 不直接选择单粒子能量，但由修正后的 $\kappa_{BB}$ 及不变的热力学 $\epsilon,P,n_B,T$ 重新计算。

作者已于 2026-07-15 明确确认以上五项。公式审计 gate 已通过；源码、测试、baseline 和 production 仍按后续里程碑依次推进。

## 6. 确认后的最小实现设计

### 6.1 能量职责拆分

- [x] 引入语义明确的内部 helper，分别得到 $E_{\mathrm{kin}}$ 和 $E_{\mathrm{dist}}$。
- [x] 调整 `_species_transport_state`，让积分调用方显式获得两种能量；避免继续以无角色名称 `E` 传递。
- [x] 保留 provider 的 `prefer_energy_aniso` 作为分布路由兼容字段，并在文档中明确它不是 kernel energy policy。
- [x] 避免在热积分循环中引入不必要分配；使用 `NTuple{5,Float64}` 返回状态，并在 $\xi=0$ 或 provider 无 `energy_from_p_aniso` 时复用已计算的 $E_{\mathrm{kin}}$。

### 6.2 目标调用方

- [x] $\eta$ 的 `p6 / E^2` 改用 $E_{\mathrm{kin}}$，分布保持 $E_{\mathrm{dist}}$。
- [x] $\sigma$ 的 `p4*q2 / E^2` 改用 $E_{\mathrm{kin}}$，分布保持 $E_{\mathrm{dist}}$。
- [x] $\zeta$ 严格按 5.4 表逐项修改，不做全局文本替换。
- [x] `diffusion_coefficient` 的 `p4/E^2` 和 `_kappa_projection` 改用 $E_{\mathrm{kin}}$，分布保持 $E_{\mathrm{dist}}$；全部 $\kappa_{XY}$ 继承该分工。
- [x] $\lambda$ 由修正后的 $\kappa_{BB}$ 重算，并同步重算依赖 $\lambda$ 或 $\sigma$ 的派生诊断量。
- [x] 保持 $\xi=0$ 路径在数值容差内完全退化到原实现。

### 6.3 文档与治理

- [x] 在权威公式文档记录外部来源、方程号、Zotero/BibTeX 映射、“A26”状态和加号勘误。
- [x] 在 transport API 文档链接权威公式，并明确“目标公式”与“当前代码已生效语义”的差异。
- [x] 修正 `docs/api/integrals/OneLoopIntegralsAniso.md` 的 `1/E_{\rm aniso}` 矛盾，使其与 `OneLoopIntegralsAniso.jl` 的普通能量分母一致；若后续证据指向代码应变更，则另开 issue。
- [x] 更新 `TransportCoefficients.jl` docstring，移除把“A26”当成外部公式号的表述。
- [x] 在实现完成后把稳定能量策略同步回 `docs/reference/formula/` 与 `docs/api/`，活动文档只保留过程 gate 和验证记录。

## 7. 测试、回归与验收

### 7.1 Unit 层

- [x] 在 `tests/unit/relaxtime/test_transport_coefficients.jl` 添加能区分 $E_{\mathrm{kin}}$ 与 $E_{\mathrm{dist}}$ 的 deterministic toy provider。
- [x] 验证 $\xi\ne0$ 时 $\eta$、$\sigma$ 的分母使用 $E_{\mathrm{kin}}$，而分布调用实际收到/构造 $E_{\mathrm{dist}}$。
- [x] 对 $\zeta$ 分别覆盖 $p^2/E^2$、$x=E\mp\mu$、$dE/dT$、$dE/d\mu_B$ 和 $B$ 核，避免只测最终数值而遗漏局部混用。
- [x] 验证 $\kappa_{XY}$ 的 $p^4/E^2$ 与两个投影因子使用 $E_{\mathrm{kin}}$，分布使用 $E_{\mathrm{dist}}$。
- [x] 验证 $\lambda$ 严格由修正后的 $\kappa_{BB}$ 和热力学 $\epsilon,P,n_B,T$ 派生。
- [x] 验证 `prefer_energy_aniso=true/false` 在默认 provider 下只改变等价调用路由，不改变目标物理结果。
- [x] 验证 $\xi=0$ 与修改前固定点一致。

### 7.2 Integration 层

- [x] 运行 `tests/integration/relaxtime/test_transport_workflow_smoke.jl` 与 TOML `prefer_energy_aniso` 路由 smoke。
- [x] 验证统一 transport workflow 能同时传递 $E_{\mathrm{kin}}$ 核和 $E_{\mathrm{dist}}$ 分布语义。
- [x] 验证扩散矩阵、$\lambda$、Lorenz number 和 Prandtl number 使用修正后的上游结果。
- [x] 验证 relaxation time、散射率和 `propagator_xi_policy = match_thermo` 在工作流层未改变；源码 diff 未触及这些模块，workflow integration 验证显式传入的 $\tau$ 原样返回，legacy tau guardrail 通过。

### 7.3 Regression 层

- [x] 运行 `tests/regression/relaxtime/test_transport_fixedpoint_regression.jl`。
- [x] 将预期数值漂移限制在 $\xi\ne0$ 的 $\eta$、$\sigma$、$\zeta$、$\kappa_{XY}$、$\lambda$ 及其派生比值；$\tau$、截面、散射率不得因本次修改漂移。
- [x] 不为通过测试而放宽容差；任何 baseline 更新必须记录旧值、新值、相对漂移和物理原因。
- [x] 在作者确认和代码 review 前不得刷新 baseline；作者于 2026-07-15 确认漂移符合预期并批准刷新三个 $\xi\ne0$ 固定点。

### 7.4 Validation 层

- [x] 运行相关 legacy transport guardrail，确认 $\xi=0$ 外部/legacy 对照没有回归。
- [x] 对 $\mu_B=0$ 的 $\zeta$ 核进行 Mykhaylova-Sasaki Eq. (10) 形式交叉核验；该检查验证公式等价，不虚构未提供的文献数值基线。
- [x] 对 $\kappa_{XY}$ 核逐项核验 Das et al. Eq. (55)，并区分外部各向同性公式与仓库完整 RS 分布扩展约定。

### 7.5 验收标准

- [x] 所有目标调用方可从代码命名和测试证据区分 $E_{\mathrm{kin}}$ 与 $E_{\mathrm{dist}}$。
- [x] $\eta$、$\sigma$、$\zeta$、$\kappa_{XY}$ 与表 5.2-5.5 逐项一致，$\lambda$ 使用修正后的 $\kappa_{BB}$。
- [x] `prefer_energy_aniso=false` 不再被描述成 kernel energy policy。
- [x] $s$、$v_{\rm rel}$、$t_\pm$、Pauli blocking、传播子策略、弛豫时间保持已确认语义。
- [x] 角积分有效量和守恒荷矩阵的适用范围在 API/公式文档中明确，不声称完成空间各向异性输运张量分解。
- [x] unit、integration、regression、validation 的目标检查通过，且数值漂移记录可追溯。

### 7.6 PR 1 验证证据（2026-07-15）

- Unit：`tests/unit/relaxtime/test_transport_coefficients.jl` 全部通过；新增能量语义 testset 为 38/38。
- Integration：workflow smoke 80/80；TOML `prefer_energy_aniso` smoke 22/22。
- Validation：公式映射 3/3；legacy transport/tau guardrail 111/111。
- Registry：4 个 `case_slug × mode` 条目通过 required fields、枚举、唯一性、SHA、目录和 manifest 路径校验。
- Regression：刷新前 116/117；没有修改 `rtol=8e-2`，唯一失败是 $\xi=-0.2$ 的 $\zeta$ 相对旧 baseline 漂移 $-10.46\%$。作者于 2026-07-15 确认该漂移是允许的预期物理语义修复，并批准只刷新三个 $\xi\ne0$ 固定点；刷新后 117/117。

| 点位 | $\eta$ 相对漂移 | $\sigma$ 相对漂移 | $\zeta$ 相对漂移 |
| --- | ---: | ---: | ---: |
| $(T,\mu,\xi)=(0.9,0.0,+0.2)$ | $+4.17\%$ | $+3.90\%$ | $+8.01\%$ |
| $(0.9,0.15,+0.2)$ | $+4.19\%$ | $+3.91\%$ | $+8.14\%$ |
| $(0.9,0.0,-0.2)$ | $-5.53\%$ | $-5.13\%$ | $-10.46\%$ |

| 点位 | 字段 | 旧 baseline | 新 baseline | 相对漂移 |
| --- | --- | ---: | ---: | ---: |
| $(0.9,0.0,+0.2)$ | $\eta$ | $7.3798253967065408\times10^{-2}$ | $7.6878503989921562\times10^{-2}$ | $+4.17\%$ |
| $(0.9,0.0,+0.2)$ | $\sigma$ | $1.4110024464777007\times10^{-3}$ | $1.4659685651953233\times10^{-3}$ | $+3.90\%$ |
| $(0.9,0.0,+0.2)$ | $\zeta$ | $3.8817152229423980\times10^{-2}$ | $4.1928048068062808\times10^{-2}$ | $+8.01\%$ |
| $(0.9,0.15,+0.2)$ | $\eta$ | $7.6421177176324587\times10^{-2}$ | $7.9619618296040098\times10^{-2}$ | $+4.19\%$ |
| $(0.9,0.15,+0.2)$ | $\sigma$ | $1.4635989266414187\times10^{-3}$ | $1.5207884350976412\times10^{-3}$ | $+3.91\%$ |
| $(0.9,0.15,+0.2)$ | $\zeta$ | $3.8540178736649888\times10^{-2}$ | $4.1677978356060852\times10^{-2}$ | $+8.14\%$ |
| $(0.9,0.0,-0.2)$ | $\eta$ | $1.3013608965647167\times10^{-1}$ | $1.2293346713344948\times10^{-1}$ | $-5.53\%$ |
| $(0.9,0.0,-0.2)$ | $\sigma$ | $2.4682134925219350\times10^{-3}$ | $2.3416435792456094\times10^{-3}$ | $-5.13\%$ |
| $(0.9,0.0,-0.2)$ | $\zeta$ | $4.9729596683791580\times10^{-2}$ | $4.4530234053340556\times10^{-2}$ | $-10.46\%$ |

全部 $\xi=0$ 固定点只出现浮点重算量级差异，最大相对差约 $4.3\times10^{-11}$。现有 baseline 未覆盖 $\kappa_{XY}$ 和 $\lambda$，因此本 PR 以 term-level unit、workflow 派生链 integration 和 Das Eq. (55) validation 约束其语义，不虚构旧数值基线。

## 8. Production 数据治理

### 8.1 受影响且必须原样保留的正式 case

当前仓库内以下两个 slug 各含 mode A、mode B 的结果和图像：

- `first_canonical_v1_p128_validated_anchored_prod_v1`，source commit `c3f25dde5c9e62e7a3ebcd8ed50aef2b53377c0f`；
- `first_canonical_v1_p128_xi001_validated_anchored_prod_v1`，source commit `700a0845abd09013eb39ff8f8d2993d5574476d3`。

结果根路径为：

- `data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/<case_slug>/`
- `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/<case_slug>/`

图像根路径为：

- `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/<case_slug>/`
- `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/<case_slug>/`

旧 case 内的 README、`PRODUCTION_AUDIT.md`、manifest、CSV、PNG 和 `plot_manifest.json` 受哈希/审计关系约束，必须逐字节保留；不得在旧目录内补标签，也不得以 `--overwrite` 重跑。

### 8.2 外部 registry

作者 gate 通过后，在旧 case 目录外创建：

`data/outputs/results/relaxtime/transport/phase_guided/production_registry.json`

- [x] 定义并校验 registry schema，至少包含 `case_slug`、`mode`、`result_path`、`figure_path`、`source_commit`、`transport_kernel_energy_policy`、`distribution_energy_policy`、`propagator_xi_policy`、`status`、`superseded_by`、`manuscript_eligible`、`audit_manifest_path`、更新时间。
- [x] 把上述旧 case 标记为 `superseded_for_manuscript`、`manuscript_eligible=false`，但不修改旧目录。
- [x] 新 case 尚未生成时允许 `superseded_by=null` 或明确的 pending 状态。
- [ ] 后续新数据 PR 先登记 `current_candidate`；完成收敛、provenance、图像和新旧对比后才晋升为 `approved`、`manuscript_eligible=true`。
- [ ] 回填旧 case 的 `superseded_by`，使 registry 成为论文输入资格的权威状态入口。

### 8.3 新 production 门槛

- [ ] 代码 PR 合并后，从合并提交运行 GitHub Actions；Actions 产物不得直接视为仓库正式数据。
- [ ] 使用新的、语义可辨识的 case slug，不覆盖任何现有正式目录。
- [ ] 先通过参数收敛 gate，再生产全量 CSV 和图像。
- [ ] 比较修改前后 $\eta/s$、$\sigma/T$、$\zeta/s$；若 production workflow 输出 $\kappa_{XY}$、$\lambda$ 或其派生量，也必须纳入差异表。同时证明弛豫时间未因 kernel energy 修复改变。
- [ ] 通过独立数据/图像 PR 导入、审计和晋升。

## 9. 里程碑

### M0：公式来源与审计 gate（当前 Draft PR）

- [x] 确认目标物理约定为“各向异性只进入分布函数，不重新定义准粒子色散关系”。
- [x] 记录 Albright-Kapusta、Mykhaylova-Sasaki、Thakur、Sasaki-Redlich 的精确方程映射。
- [x] 将“A26”降级为 legacy 内部别名并记录加号勘误。
- [x] 创建本逐项审计表，覆盖 $\eta$、$\sigma$、$\zeta$、$\kappa_{XY}$、$\lambda$ 和共享 helper 风险。
- [x] 作者逐行确认 5.6 的五项结论（2026-07-15）。

### M1：代码、测试、稳定文档与 registry（代码 PR）

- [x] 拆分 $E_{\mathrm{kin}}$ 与 $E_{\mathrm{dist}}$。
- [x] 按审计表修改 $\eta$、$\sigma$、$\zeta$、$\kappa_{XY}$，并重算 $\lambda$ 及其派生量。
- [x] 完成 unit、integration、regression、validation 证据。
- [x] 修正 OneLoop API 矛盾和 transport docstring。
- [x] 创建外部 production registry 并标记旧 case 的论文资格。
- [ ] 代码 PR review、CI 和合并完成。

### M2：新数据与图像（独立 production PR）

- [ ] 由代码合并提交运行收敛 gate 和正式 GitHub Actions。
- [ ] 使用新 case slug 导入结果与图像，不覆盖旧产物。
- [ ] 完成新旧差异、provenance、manifest 和论文输入资格审计。
- [ ] 新 case 经 review 后在 registry 晋升为 `approved`。

## 10. Definition of Done

- [x] 作者确认本审计表，确认记录可在 issue/PR 中追溯（2026-07-15）。
- [x] $E_{\mathrm{kin}}$ 与 $E_{\mathrm{dist}}$ 在实现、测试和稳定文档中职责一致。
- [x] $\eta$、$\sigma$、$\zeta$、$\kappa_{XY}$ 的代码逐项映射到已记录的仓库公式和外部方程号，$\lambda$ 的派生链可追溯。
- [x] 非目标的散射、传播子和弛豫时间语义有保护证据。
- [x] 所选测试层通过，数值漂移有物理说明且未通过放宽容差掩盖。
- [x] 旧正式产物未被修改，并由外部 registry 明确标为不再进入当前论文输入包。
- [ ] 新数据和图像使用新 case slug，经独立 production PR 审计和批准。
- [ ] 稳定公式/API 文档已更新；任务完成后本文件按仓库流程归档。

## 11. 风险与回退方案

| 风险 | 触发信号 | 缓解/回退 |
| --- | --- | --- |
| $\zeta$ 某个能量角色仍有公式歧义 | 作者无法确认表 5.4 的某一行 | $\eta$、$\sigma$ 可按独立证据推进；$\zeta$ 保持未修改，production 继续阻塞，不用推测补齐。 |
| 共享 helper 拆分后 $\kappa/\lambda$ 未按目标变化或发生额外漂移 | term-level unit/regression 失败 | 分别检查分母、两个投影因子、分布和 $\lambda$ 派生链；禁止用单一最终数值掩盖局部混用。 |
| $\xi=0$ 发生数值漂移 | fixed-point 或 legacy guardrail 失败 | 停止 baseline 更新，检查 helper 分流和舍入路径；必要时回退实现提交。 |
| provider 兼容性破坏 | toy/custom provider 测试失败 | 保持现有字段兼容，新增内部适配层；不把新策略强加为破坏性 public contract。 |
| 新 production 与旧路径混淆 | 复用了旧 slug 或尝试 `--overwrite` | 立即停止运行；使用全新 slug，从代码合并提交重新生成并走独立数据 PR。 |
| 旧数据状态只能靠聊天记忆 | 论文输入选择未读取 registry | 在新数据前先完成外部 registry；旧目录保持不可变。 |
| 标量结果被过度解释为张量响应 | 论文/API 使用 longitudinal/transverse 完整响应措辞 | 在稳定文档和论文 handoff 中明确本次只覆盖角积分后的标量有效量。 |

回退边界：代码 PR 可整体回退到拆分前实现；旧 production 因始终不被覆盖，无需数据回滚。若新 case 未通过 gate，只保留为 `current_candidate` 或撤销导入，不得晋升为 `approved`。
