# Charged RPA/BU 严格复核与后续实现任务单

更新日期：2026-09-07

当前状态：研究生产入口与配置冻结线验收完成，作者已审核图像；PR收尾中。
当前有效结论以第22.23–22.25节为准，前面各轮记录是历史过程，不覆盖最新验收。
本文承接完整 KMT interaction-kernel 与固定 BQS
quark-only 后处理诊断，目标是重新核对 charged `π^±/K^±` 的 RPA/BU 数值链，
但不在本任务单完成前替换现有 `MesonDensity` 生产语义。`x_min_cut`、低节点数、
未闭合的 retarded/phase 数值 gate 仍只能标记为 diagnostic。

最新研究状态见第22节：固定 GBU partial-yield 方法、同输入比较和稀疏冻结线
诊断已推进；前部“当前链路”保留原默认入口审计上下文，不能读作新因果谱泡
已接管生产。有限 q 截断处方差异仍是论文路线选择前的重要未决项。
第22.18节进一步确认当前真空谱外的负热尾；完整物理粒子谱解释待作者评审，
暂不扩大 Mott/GBU 生产验收。已有冻结线结果保持 diagnostic，不据此改默认。
作者后续明确优先方向B；第22.19节先核对热延伸的静态基础与响应匹配要求，
不将第22.18节的负谱反例扩大成“上游热延伸无效”，也不自动切换方向A。
第22.20节给出现核的固定投影源泛函；第22.21节进一步核验GBU代数与观察量边界。
以上为第22.21节时的历史状态；第22.23–22.24节已完成无限热B路线的
数值验收、显式入口与10点完整配置扫描，不再将这些历史未决项当作当前阻碍。

PR290 从独立的 `origin/main` 基线建立了
[公式路线闭合包](../../reference/formula/relaxtime/ChargedRPA_BU_ProductionRoute.md)
和治理门禁；本任务单的严格数值复核依赖该 candidate 规范，但不会因文档闭合而
自动完成或晋升为 production。

PR289 已以 squash commit `dba1d9b8` 合并，落地独立 charged scalar-RPA
kernel/spec、`K^± -> K45` 映射和第一版 provider 契约。本任务后续按可审核性拆成
三个 PR：ordered retarded 固定点、严格 BU 测度/相位门禁、全链路收敛与冻结线。
second-sheet pole solver 保持为独立非阻塞任务。

## 1. 决策目标与假设

### 1.1 要回答的问题

本审计只回答“完整 charged-RPA/BU 相比当前 scalar BU 诊断，哪些改变是理论上
必需的、哪些只是数值实现差异，以及如何逐项验收”。它不把 `K_{03}/K_{38}`
强行加入 charged `K^±/π^±`：在当前对角平均场背景下，这些项属于中性
`(0,3,8)` 混合块，不能直接进入 `K_{12}` 或 `K_{45}` 的带电标量分母。

### 1.2 固定假设

- 上游仍是 `FixedMuBConservedCharges` 的 BQS quark-only 平衡，约束
  `rho_Q/rho_B=0.4`、`rho_S=0`；不加入 `Omega_M` 或 `Sigma_M` 反馈。
- `phi_f=<bar q_f q_f>`、`G`、`K` 使用 [KMT 平均场到 RPA 核公式](../../reference/formula/relaxtime/couplings/KMT_MFA_to_RPA_QuadraticKernel.md)
  的约定。
- 目标 observables 是 charge-resolved `K^+/pi^+` 与 `K^-/pi^-`；`K^0`、
  `bar K^0` 只作为通道映射的审计参照。
- 本轮不做 strict-support 的凝聚处理，不以 `x_min_cut` 结果授予 production 资格。
- 介子数密度保持四类解耦入口：稳定粒子极限、reduced strict-BW、`q` 依赖复极点
  strict-BW 和 phase-shift BU；最终比较默认使用 `phase_shift_gbu_reference`，
  但四类入口和 phase-shift 的 `current` 对照均保持可调用。
- strict GBU 使用有序 `Pi_us/Pi_su` 的实轴 retarded 泡；`num_s_quark=1` 作为
  Rehberg 有来源的散射正则化处方和课题组旧 Fortran/Cpp oracle 独立保留。

## 2. 当前有效链路（已核对）

当前 charged BU 数值链实际是

```text
PNJL/BQS x_state=(phi_u,phi_d,phi_s,Phi,PhiBar), mu_u/d/s
  -> charge-resolved PolarizationAniso bubble Pi_a
  -> scalar D_a = 2 K_a / (1 - 4 K_a Pi_a)
  -> delta_a = arg(D_a) or -arg(D_a^{-1})
  -> finite (q,omega) BU integration
```

对应源码为 `src/relaxtime/PolarizationAniso.jl`、
`src/relaxtime/MesonPropagator.jl` 和 `src/relaxtime/MesonDensity.jl`。

当前的 charge-resolved 极化输入已经区分味道顺序：

| 通道 | 当前 bubble 顺序 | 完整 KMT 耦合 | 物理组成 |
|---|---|---|---|
| `pi_plus` | `(u,d)` | `K12` | `u bar(d)` |
| `pi_minus` | `(d,u)` | `K12` | `d bar(u)` |
| `K_plus` | `(u,s)` | `K45` | `u bar(s)` |
| `K_minus` | `(s,u)` | `K45` | `s bar(u)` |

但 K 路径当前还会在每个有序通道内部设置 `num_s_quark=1`，把同序 `B0(k0)` 与
`B0(-k0)` 平均。因此“输入顺序已区分”不等于“strict ordered-retarded backend
已实现”。

旧 `EffectiveCouplings` 的 `K123` 代数上等于 `K12`；旧 `K4567` 使用
`H_u=-phi_u`，在 `phi_u!=phi_d` 时等于 `K67`，不是 charged `K45`。
这一映射已经由 `tests/unit/relaxtime/test_meson_interaction_kernel.jl` 的纯代数
测试锁定。`K03/K30/K38/K83` 只在中性 `(0,3,8)` 矩阵中出现；当前没有非对角
夸克平均场或 charged-neutral 混合传播子，因此它们不会直接改变上述四个 charged
标量分母。

## 3. “完整 charged-RPA/BU”必须重新核对的独立层次

### 3.1 相互作用核：`K45` 替换只是必要条件，不是完整 RPA

在固定的 `Pi_a`、相位定义和分母归一化下，`K67 -> K45` 是唯一有意的代数
替换；当前冻结线 A/B 的差异很小，正是这一层的诊断结果。完整 RPA 还必须
证明以下每一层使用同一微观拉氏量的约定：

1. KMT Hartree 收缩得到的 `K12/K45/K67` 与 `phi_f` 的符号；
2. 介子二次作用量的顶角归一化；
3. 夸克泡 `Pi_a` 的 flavor 顺序、正则化和 retarded continuation；
4. 传播子极点、相移和 BU 热权重使用同一归一化。

### 3.2 分母因子：不能把中性矩阵的 `2` 机械搬到 charged 标量道

当前项目并存两种公式：

```math
D_a^{\rm current}=\frac{2K_a}{1-4K_a\Pi_a},
\qquad
D_{0,3,8}^{\rm matrix}=2\mathcal K\,[I-2\mathcal K\Pi]^{-1}.
```

后一个式子是 `MesonRPA` 中实生成元基底的公式。对任一 charged pair，

```math
T_\pm=(\lambda_a\pm i\lambda_b)/\sqrt2=\sqrt2 E_{ij/ji}
```

使矩阵二次作用量中的 charged polarization 成为单个 Rehberg/项目有序泡的两倍：

```math
\Pi^{matrix}_{+-}=2\Pi_{ij},\qquad
\Pi^{matrix}_{-+}=2\Pi_{ji}.
```

所以 `1-2KPi_matrix` 与 `1-4KPi_ordered` 是同一分母，不是两个待数值选择的
候选。该结果不要求 `Pi_ij=Pi_ji`；显式 ladder trace、chiral-limit Goldstone
identity 和错误 `1-2KPi_ordered` 不产生零点的反例已加入
`tests/unit/relaxtime/test_meson_rpa.jl`。低温真空标定区的 pion/kaon pole 继续由
现有 literature/legacy Fortran validation 约束。

### 3.3 极化函数：`Pi_{us}` 与 `Pi_{su}` 不是可自动合并的同一个数组

`PolarizationAniso` 当前通过 `B0(lambda)`、`A1+A2` 和可选的有限宽度修正构造
泡，其中 `lambda=k0+mu1-mu2`，并对 `num_s_quark=1` 做额外的 `k0` 对称平均。
这个开关可追溯到 Rehberg 等 NPA 608 (1996) Eq. (2.22)-(2.23)，并沿用于课题组
旧 Fortran/Cpp；它是有来源的 legacy prescription，而不是任意修补。原文用途是
修复三动量截断下非等质量、有限交换动量散射泡的伪对称破缺，不能据此把同序
`p0` 平均认作反序 `Pi_su`。
严格 charged 路线必须重新确认：

- `K^+` 的 `(u,s)` 与 `K^-` 的 `(s,u)` 是否按 retarded 约定分别取值；
- 在 `mu_u=mu_d=mu_s=0` 时是否满足电荷共轭关系和相同的真空极点；
- 有限 `mu` 时的 Landau 区、阈值和虚部符号是否与相移分支一致；
- `A_f` 与 `B0` 是否使用相同 cutoff、热积分上限和节点。按当前公式的自然单位，
  `A_f` 与 `Pi_a` 的量纲应为 `fm^-2`；`K_a` 与传播子 `D_a` 为 `fm^2`，因此
  `K_a Pi_a` 无量纲。旧 `MesonPropagator` 与 `PolarizationCache` API 页面已同步
  为 `Pi=fm^-2`；这只修正文档契约，不改变当前数值实现。

当前 `phi` 原生适配层只消除了“由 `A_f` 再构造 `H_f`”这一重复步骤；它不等于
已经完成 `Pi_a` 的严格 retarded 解析延拓。`A_f` 仍可能作为夸克泡的 tadpole
输入，这是两个不同的数值角色。

### 3.4 legacy 宽度、retarded 极点与 GBU 实轴对象必须分开

当前旧接口以 `p0=M+iGamma/2` 展开 prefactor，但 `B0` 仍在实轴求值。这个约定
与课题组旧 Fortran/Cpp 一致，并且在 `mu_u=mu_d`、排除 `K45` 错配后可复现旧
质量和正宽度；因此不应全局翻转 `gamma`。它应明确标为
`legacy_real_axis_width`，不能称作完整复平面极点。

若 `q_pole_strict_bw` 声称给出真正的 retarded pole，则每个 charged 通道应由

```math
\Delta_a(z,q)=0,\qquad z=\omega-i\Gamma/2
```

这里采用 `e^{-i\omega t}` 的 Fourier 约定，retarded 共振极点位于跨过 cut 的
second sheet 下半平面。极点求解至少需要：sheet、残差、收敛状态、分支 seed、
阈值及 flavor 顺序。

GBU 主路线只需要物理实轴上 `D^R(omega+i0^+)` 的相位，不需要先求这个复极点。
因此 second-sheet solver 是 `q_pole_strict_bw` 的严格 oracle，不是
`phase_shift_gbu_reference` 的 blocker。

### 3.5 相移与 BU 权重：常数相位、分支和边界都会改变密度

严格路线应从 retarded `D_a^R` 定义

```math
\delta_a(\omega,q)=\operatorname{Im}\log D_a^R(\omega,q)
```

并明确减去高能端常数、处理束缚态的 `pi` 跳变、连续谱相位和 Levinson 边界。
当前 `MesonDensity` 支持 `arg_propagator`/`arg_inverse_propagator`、unwrap 和
`current`/`gbu_reference` 两种权重，但这些是可切换诊断约定，不应在没有公式
固定时混合使用。本任务将两种权重都保留为可调用方案，并把
`gbu_reference` 固定为最终比较默认；这不替代相位边界、Levinson 和收敛验证。
BU 原式及分部积分形式见
[MesonDensity_BU相移公式.md](../../reference/formula/relaxtime/meson_density/MesonDensity_BU相移公式.md)。

特别是分部积分后的密度使用 `delta_a` 本体；若 `delta_a(omega->infty)` 或
低能边界没有正确归一化，会产生与物理相互作用无关的伪密度。因此相位边界和
Levinson 检查必须先于冻结线扫描。

这里的 production gate 不是“已有 atan+unwrap 即通过”，而是逐 `(channel,q)`：
从高能端向低能端连续 unwrap、验证增大 `omega_max` 后端点与积分稳定、独立计数
阈下束缚根，并检查
`delta(omega_thr)-delta(infinity)=pi*n_B`；穿越 Mott 点时还要验证 `n_B` 减一
与阈值相位减少 `pi` 同步。最后对 `omega/q`、`eta`、cutoff 和两种 BU 积分形式
做收敛。

另一个与相位分支无关的常数归一化已由稳定极限固定：本项目对单个 ordered
`K+`/`pi+` 通道使用 `d=1` 和一个 Bose 因子，因此正能量测度必须是
`domega/pi`。束缚态的 `pi` 跳变于是返回一个稳定玻色子；历史
`domega/(2pi)` 会使绝对密度少一半。代码现已提供显式 strict/legacy measure，
且不静默切换旧默认。该共同因子在同口径 `K/pi` 比值中抵消，所以旧 ratio
诊断不因此失效，但四算法绝对密度比较尚未闭合。

### 3.6 Bose 支撑与凝聚：`x_min_cut` 不是严格解

当 `omega<=mu_M` 时，玻色分布在 `omega=mu_M` 有极点。沿化学冻结线先逐点
检查 `mu_M<m_M` 和积分节点的 `omega>mu_M`；若正常相门禁通过，则不引入
介子凝聚零模。当前 `density_policy=:x_min_cut` 把下界移动到
`mu_M+x_min*T`，这是文献复现/异常点诊断的简化截断，不是对零动量凝聚模的处理。
严格路线需要二选一并写入公式合同：

1. 证明所有生产点处于 `mu_M<m_M` 的正常相，并对 `omega>mu_M` 做支撑门禁；或
2. 单独引入凝聚零模、守恒荷和连续激发的分解。

在此之前，节点增加导致的变化不能被解释为物理的 `K/pi` 趋势。

## 4. 预期会造成差异的来源（归因矩阵）

| 变更 | 是否属于 `K45-K67` | 预期影响 | 当前状态 |
|---|---:|---:|---|
| charged 耦合 `K67 -> K45` | 是 | 当前背景下很小；强 `phi_u-phi_d` 时可放大 | 已有 diagnostic A/B |
| charged 分母 `4KPi_ordered` vs `2KPi_matrix` | 否 | 同一二次作用量下严格等价 | 已由 ladder trace + Goldstone 闭合 |
| `Pi_{us}`/`Pi_{su}` 解析延拓与宽度 | 否 | 直接改变 `K^+` 与 `K^-` 的差异和阈值 | 有限 `eta` ordered backend 已实现；全网格收敛待后续门禁 |
| `num_s_quark=1` vs ordered retarded | 否 | 可改变有限 `q/mu` 的 kaon 相位 | 固定 BQS 探针已显示 kaon 传播子可达数十个百分点差异 |
| `A_f/B0` cutoff、热上限、节点 | 否 | 数值漂移或伪峰 | 需收敛门禁 |
| 相位分支、常数边界、Levinson | 否 | 可造成整体密度偏置 | 需单独测试 |
| 单电荷 BU 测度 `1/pi` vs `1/(2pi)` | 否 | 绝对密度差一倍；同口径 ratio 抵消 | 公式已由稳定极限闭合，代码待迁移 |
| Bose 支撑/凝聚处理 | 否 | 接近 `mu_M=m_M` 时可能发散 | 当前仅诊断截断 |
| `Omega_M`/`Sigma_M` 反馈 | 否 | 改变上游 `phi,m,mu`，属于另一条路线 | 暂不实现 |
| 非对角夸克平均场或 charged-neutral 混合 | 否 | 可能让中性 `K03/K38` 间接进入 charged sector | 当前模型未开启 |

因此，不能从现有 `0.0233%` A/B 差异推出完整 charged-RPA/BU 的差异也只有
这个量级；这个小数值只约束“固定当前 `Pi` 和当前 scalar denominator 时的
coupling substitution”。

## 5. 后续实现分层

### Phase A：公式与不变量（当前审计）

- [x] 固定 `K^± -> K45`、旧 `K4567 -> K67` 的味道映射。
- [x] 固定 `K03/K38` 只属于中性 `(0,3,8)` 块的当前模型边界。
- [x] 从 charged ladder 顶角推出 `Pi_matrix=2Pi_ordered`，并以 Goldstone identity
  固定 `1-4KPi_ordered`。
- [x] 记录 `phi` 原生旧耦合适配层与 `A_f` 在泡计算中的不同角色。
- [x] 分离 source-backed `num_s_quark=1`、legacy `+iGamma/2` 与 strict retarded
  GBU/pole 的解析对象。
- [x] 以 `pi` 相移跳变的稳定粒子极限固定单电荷正能量测度为 `domega/pi`；
  `domega/(2pi)` 只保留为文献 ratio adapter。

### Phase B：charged RPA kernel backend（不改旧接口）

- [x] 新增独立的 charged scalar-RPA kernel/spec，显式记录 `pair`、P/S、
  numerator、denominator、retarded convention 和 normalization source。
- [x] 用同位旋对称极限把该 backend 与旧 `MesonPropagator` 做代数/极点 parity
  测试；不以一轮冻结线 A/B 代替该测试。
- [x] 为 `K^+`、`K^-` 分别保留 `Pi_{us}`、`Pi_{su}` 的输入和诊断元数据。

Phase B 由 `src/relaxtime/ChargedRPAKernel.jl` 落地。该模块只提供经过校验的
通道/归一化契约和纯单通道代数；默认 `MesonPropagator`、相移和 BU 密度语义
保持不变。PR290 已用 charged ladder trace 和 Goldstone identity 固定
`D=2K/(1-4KPi_ordered)`，因此本后端不再暴露 `1-2KPi_ordered` 候选。
同位旋对称 parity 和有序 `Pi_{us}/Pi_{su}` 输入由
`tests/unit/relaxtime/test_charged_rpa_kernel.jl` 锁定。

### Phase C：charged bubble provider 与极点

- [x] 新增显式有限 `eta>0` 的复能量 `B0_retarded`，默认
  `:ordered_retarded` 从上半平面求值；历史 `num_s_quark=0` 路线保留为
  `:ordered_legacy_B0`，`num_s_quark=1` 只保留为
  `:legacy_symmetrized_B0` oracle。
- [x] 以 unit gate 验证 ordered 共轭关系，并在 `T=170 MeV, mu_B=240 MeV`
  的 `FixedMuBConservedCharges` quark-only 背景完成四 charged 通道固定点对照。
- [ ] 稀疏冻结线对照推迟到节点、截断、Bose support 和四算法门禁所在的第三个
  后续 PR，不在本 PR 提前消费未认证相位。
- [ ] legacy `Gamma` 固定点与 second-sheet pole solver 独立实施；后者仅在需要
  严格 pole mass/width 时启动，不阻塞实轴 GBU 路线。

固定点脚本为
`scripts/analysis/relaxtime/compare_charged_rpa_ordered_fixedpoints.jl`。它对每个
`pi^±/K^±` 取阈下 `q=0` 与阈上有限 `q` 两个探针，并分别输出 strict eta-coarse、
node-coarse、refined、ordered-legacy 与 symmetrized-legacy 五条处方。原始 CSV
只写入本地 `data/outputs/results/relaxtime/analysis/charged_rpa_ordered_fixedpoint/`，
不纳入版本库，所有行均标记 `production_candidate_status=not_authorized`。

`mu_B=0` charge-symmetric 参考背景残差为 `9.82e-16`；此时
`rho_Q-0.4rho_B=0` 只是零密度退化约束，不将 `Q/B=0.4` 解释为有限比值。
strict refined 在 `q=0` 的正能量 `K^+/K^-` 传播子一致到 `2.64e-10` 相对量级，
有限 `q=0.35 fm^-1` 时仍有约 `1.30e-3` 的 routing/cutoff 差异；symmetrized
legacy 在同一门禁下达到机器精度。该小残差不阻塞 ordered 实现进入后续相位审计，
但必须在 production 评审时与 cutoff/routing 收敛一起解释。

有限 BQS 背景残差为 `1.44e-15`，得到 `mu_Q=-9.22 MeV`、`mu_S=80.27 MeV`。
在该有限探针集合上，ordered legacy 相对 strict refined 的传播子最大差异约
`1.32%`；`num_s_quark=1` symmetrized legacy 的 kaon 传播子最大差异约
`44.8%`，而 pion 的两个 legacy 处方代数上相同。这个结果首次表明旧实现的主要
新增风险不只来自 `K67 -> K45`，还来自 kaon 内部 `k0/-k0` 平均；但它仍不是
相位、密度或冻结线结论，且 strict `eta/N_E` 尚未通过全域收敛门禁。
当前固定点上 `N_E=512 -> 1024` 的传播子最大相对变化约 `1.13e-4`，而
`eta=0.01 -> 0.005 fm^-1` 仍约为 `1.32%`；因此节点固定点已足够支持本 PR 的
实现审核，但 `eta -> 0^+` 外推明确留到后续全域收敛 PR。

### Phase D：strict BU density

- [x] 将 `current` 与 `gbu_reference` 保持为可调用方案，并固定 `gbu_reference`
  为最终比较默认；相移边界和 Levinson 仍待验证后才能成为 production 定义。
- [x] 规定冻结线先执行正常相 Bose 支撑门禁；`x_min_cut` 只能留作 diagnostic，
  若门禁失败再另建凝聚零模分支。
- [x] 将 phase-shift 单电荷实现迁移为显式 strict `domega/pi` 与 legacy
  `domega/(2pi)` measure；旧默认不切换，strict/BW 稳定极限和共同因子由 unit
  regression 锁定。
- [x] 实现从高能端反向 unwrap 的 phase anchor、阈下简单根计数以及可失败的
  Levinson/Mott gate。
- [ ] 将 ordered charged profile 接入上述 gate，并通过 q/omega 节点、`eta`、
  `omega_max` 和 cutoff 收敛门禁。

### Phase E：production candidate review

- [ ] 与至少一个外部文献固定点在相同参数、单位和通道定义下复核。
- [ ] 冻结线只消费通过上述门禁的 quark-only background，输出 `K^+/pi^+`、
  `K^-/pi^-`，并保留失败点和完整 provenance。
- [ ] 单独评审是否需要更改 `Omega_M`/`Sigma_M` 或非对角平均场；不得由本
  charged-RPA backend 自动隐式反馈。

## 6. 验收与非目标

### 必须通过的验收

1. 单位：原始 KMT 六夸克耦合 `[K_KMT]=fm^5`，有效通道耦合
   `[K_a]=fm^2`，极化函数 `[Pi_a]=fm^-2`，传播子 `[D_a]=fm^2`，介子密度
   `[n_M]=fm^-3`，且 `K_a Pi_a` 无量纲。
2. 同位旋对称：`phi_u=phi_d` 时 `K45=K67`、`K03=K38=0`，charged backend
   与旧 scalar 结果在同一归一化下相容。
3. 电荷共轭：`mu=0` 时 `K^+`/`K^-` 的极点和相移满足约定的共轭关系。
4. 密度：单电荷 phase-shift 的窄束缚态极限与 `stable_meson_number_density`
   一致，文献 ratio adapter 单独命名。
5. 宽度：legacy `Gamma` 与 second-sheet pole width 分开命名、分开回归。
6. 数值：所有相移/极点残差、Bose status、节点配置和 fallback 都可追溯。
7. 物理：高能相位边界、Levinson 计数和 Mott 连续性通过后，才可讨论实验趋势。

### 明确非目标

- 不修改 `PNJLCore` 的当前对角平均场求解，不引入非对角夸克自能。
- 不把 `K03/K38` 注入 charged `K/pi`，除非未来先改变平均场/混合基底并完成新
  的二次作用量推导。
- 不实现完整 `Omega_M` 反馈、全 hadronic BQS 或显式 `mu_I` 路线。
- 不提交已有低节点 diagnostic CSV，不更新 production/regression baseline。

### 当前诊断结论与旧实现处置（2026-08-30）

现有冻结线 A/B 的最大 full/legacy 相对差异约为 `0.0233%`。这只约束了
“在同一 quark-only 背景、同一 `Pi`、同一标量分母和同一 BU 数值设置下，把旧
`K4567` 换成 charged `K45`”这一局部替换；它不能外推为严格 charged-RPA/BU
全链路的误差上限。相反，当前旧 scalar BU 的物理比值沿冻结线仍有明显温度/能量
依赖：`K^+/pi^+` 约由 `10.2` 降至 `1.24`，`K^-/pi^-` 约由 `0.0246`
升至 `0.695`。因此目前没有证据表明 `K03/K38` 是 charged ratio 大幅变化的
主因；在当前对角平均场背景中，它们仍只属于中性 `(0,3,8)` 混合块。

代表点 partial-feedback 诊断把 `mu_S` 从 quark-only 的约 `80.3 MeV` 移到约
`1.5--1.7 MeV`，并得到两个 `O(1)` 的 charged ratio。这一结果说明介子守恒荷
外层修正对约束解很敏感，但该路线没有把 `Omega_M`/`Sigma_M` 放入驻点方程，不能
被视为完整热力学反馈或旧实现的 production 修复。

据此，旧实现不需要立即进行破坏性的整体重写，处置分为三层：

1. **保留兼容层。** 旧 `MesonPropagator`、`MesonMass` 和默认 `MesonDensity`
   入口继续保留，避免改变 transport 和已有 baseline 的语义；其 `K4567` 在
   `phi_u != phi_d` 时必须明确标注为旧/legacy 的 `K67` 代理，而不是完整 charged
   `K45`。
2. **定向新增生产候选层。** `ChargedRPAKernel`/`ChargedRPAProvider` 先作为
   显式 `K^\pm -> K45` 的并行后端；严格 retarded `Pi_{us}/Pi_{su}`、极点、相位
   边界和 BU 支撑门禁通过后，再由显式适配器接入生产路线，不静默改写旧接口。
3. **严格复核层。** ordered/legacy 泡的固定点差异、有限宽度/second-sheet
   语义、相移边界、`x_min_cut` 与凝聚处理仍需要单独数值验收。若这些门禁暴露
   系统性偏差，应替换 charged 生产后端；目前没有理由因此修改 `PNJLCore`、
   引入 `K03/K38` 到 charged 标量通道，或重做 transport 主线。

换言之，当前结论是“旧接口存在非对称 charged-kaon 语义缺口，需要被新后端隔离
并最终替代”，而不是“旧 PNJL/BU 全部实现已经被证明错误”。

## 7. 当前证据与引用路径

- 纯代数核：`src/relaxtime/MesonInteractionKernel.jl`、
  `tests/unit/relaxtime/test_meson_interaction_kernel.jl`。
- charged scalar-RPA 契约与 A/B0 provider 适配：
  `src/relaxtime/ChargedRPAKernel.jl`、`src/relaxtime/ChargedRPAProvider.jl`、
  `tests/unit/relaxtime/test_charged_rpa_kernel.jl`、
  `tests/unit/relaxtime/test_charged_rpa_provider.jl`。
- ordered retarded 固定点诊断：
  `scripts/analysis/relaxtime/compare_charged_rpa_ordered_fixedpoints.jl`、
  `tests/integration/relaxtime/test_charged_rpa_ordered_fixedpoint_script_contract.jl`；
  原始 CSV 只在本地保留。
- 中性矩阵 RPA：`src/relaxtime/MesonRPA.jl`、
  `docs/api/relaxtime/propagator/MesonRPA.md`；charged ladder/Goldstone normalization
  gate：`tests/unit/relaxtime/test_meson_rpa.jl`。
- 旧质量固定点：`tests/regression/relaxtime/test_meson_mass_regression.jl`、
  `tests/validation/relaxtime/test_literature_digitized_meson_mass_targets.jl`、
  `tests/validation/relaxtime/test_legacy_fortran_meson_numeric_mu0_targets.jl`。
- 当前 charged BU scalar 诊断：`src/relaxtime/MesonDensity.jl`、
  `docs/api/relaxtime/meson_density/MesonDensity.md`。
- BU 公式与相位边界：[MesonDensity_BU相移公式.md](../../reference/formula/relaxtime/meson_density/MesonDensity_BU相移公式.md)。
- 本轮 full/legacy 冻结线诊断和 `0.0233%` 差异：
  `docs/dev/active/2026-08-29_full-kmt-phase4-phase6.md`，原始 CSV 仅本地保留。

## 8. Phase F：节点/截断、Bose 支撑与冻结线诊断（本 PR）

本阶段新增分析脚本
`scripts/analysis/relaxtime/audit_charged_rpa_bu_convergence.jl`，只复用
`Models` workflow 和 `default` + `baseline_freezeout` 参数化，不修改任何默认
production 入口。脚本按高 `sqrt(s_NN)` 到低能方向续算，并在续算失败时显式记录
`fallback_no_continuation` 或 `failed`，不会静默替换平衡解。

- [x] 复用四类现有密度入口：稳定粒子极限、reduced strict BW、`q` 依赖复极点
  strict BW、`phase_shift_bu`；每条记录同时保留 `K^+/pi^+` 与 `K^-/pi^-`。
- [x] 每个通道执行 `bose_support_gate`，将正常相支撑状态和
  `min(E-μ)` 写入诊断结果；该门禁不等价于凝聚处理。
- [x] 提供代表点 low/refined 节点/截断比较和 `convergence_gate`；四算法失败点
  保留错误文本、数值配置及 `production_candidate_status=not_authorized`。
- [x] 输出路径、`(T,mu_B)`、三味 `mu`、介子质量、夸克/守恒荷密度、算法元数据和
  平衡 residual，供后续复核使用。CSV/README 仅作为本地 diagnostic 产物，不入库。
- [x] 已完成单点 smoke（`sqrt(s_NN)=7.7 GeV`，`p_num=4,t_num=2`，
  `q_nodes=3,omega_nodes=4,stable_q_nodes=8`）。稳定粒子极限的
  `K^+/pi^+` 与 `K^-/pi^-` 均为约 `0.190`，reduced strict-BW 均为约 `0.234`；
  q-pole 与 phase-shift BU 均因 `n_pi<=0` 或非有限而标记
  `status=invalid_density`，不是成功的 ratio 证据。所有通道 Bose 支撑门禁为
  `safe_normal_domain`，平衡 residual 约 `1.3e-15`。带 low/refined 对照的同点运行
  中，稳定和 reduced-BW ratio 分别约由 `0.190 -> 0.219`、`0.234 -> 0.264`，
  在本低节点设置下均未通过默认 `5%` convergence gate；这只说明需要更高节点/截断
  审计，不能作为生产收敛结论。
- [x] 原始 CSV/README 已写入本地
  `data/outputs/results/relaxtime/analysis/charged_rpa_bu_convergence/`，未纳入版本库；
  在本任务单中只保留上述可复现状态，不把低节点结果晋升为 production baseline。

推荐最小运行示例：

```powershell
$env:MESON_CONVERGENCE_SQRTS = "7.7"
julia --project=. scripts/analysis/relaxtime/audit_charged_rpa_bu_convergence.jl
```

若需要三点定性扫描，可将环境变量设为 `3,7.7,200`；这仍属于昂贵的诊断 smoke，
不替代后续全域 `q/omega/eta/omega_max` 收敛和 production candidate review。

本任务单只记录可验证的实现边界和后续 gate，不把任何尚未完成的 charged-RPA/BU
公式或冻结线结果标记为 production。

## 9. Phase G：同一有限-BQS 背景下的四算法与负密度审计（2026-08-31）

新增只读诊断脚本
`scripts/analysis/relaxtime/diagnose_charged_rpa_bu_negative_density.jl`，在同一
`T=170 MeV`、`mu_B=240 MeV`、`rho_Q/rho_B=0.4`、`rho_S=0` 的 quark-only 解上比较
stable、reduced-BW、q-pole 和 phase-shift/GBU。该脚本只写入本地
`data/outputs/results/relaxtime/analysis/charged_rpa_bu_negative_density/`，不改变默认
入口或任何 baseline。

该点的平衡残差为 `1.44e-15`，并得到
`mu_u=0.3742644`、`mu_d=0.4209955`、`mu_s=0.0142260 fm^-1`，
`m_pi=0.6737676`、`m_K=2.1885044 fm^-1`。charged chemical potentials 为
`mu_pi+=-0.0467311`、`mu_K+=0.3600384` 以及相反号的负电荷通道。

### 9.1 四算法结果（同一背景）

| 通道 | stable | reduced strict BW | q-pole strict BW | phase-shift GBU（未锚定） |
| --- | ---: | ---: | ---: | ---: |
| `K+/pi+` | `0.56618` | `0.20221` | `NaN`（Bose 支撑失败） | `2.20354` |
| `K-/pi-` | `0.21112` | `0.07377` | `~6.4e-15` | `0.35023` |

这些算法本来就是不同近似对象，表中数值不要求彼此收敛。q-pole 的 `K+` 失败不是
极点残差问题：其内层 `omega` 网格硬编码从 `0` 开始，而 `mu_K+>0`，所以调用
`bose_distribution` 时遇到 `omega <= mu_K`。同一 helper 的 `K-` 极点四个节点均
接受，残差约 `2.0e-9` 到 `1.7e-12`；`pi` 因零宽度走 stable limit。现有
`bose_support_gate` 能识别该状态，但 q-pole helper 尚未把它作为显式返回状态，
这是后续应单独修复的接口缺口，不能把 `NaN` 当作物理零密度。

### 9.2 BU/GBU 负密度的定位

在本点所有相移积分均使用 `omega_min=0.4 fm^-1`，四个 charged 通道的实际积分
窗口均在正常 Bose 域内；Bose 因子 `g(1+g)` 为正。因此负值不是由
`mu_M=0`、Bose 奇点或 `domega/pi` 与 `domega/(2pi)` 的共同归一化因子造成的。
`arg(propagator)` 与 `arg(inverse_propagator)` 的实现结果逐点相同，故也不是两种
相位对象选择造成的符号差。

真正的证据来自 `phase_shell_breakdown`：

- `K-` 未锚定 GBU 的总积分为 `+0.0239968`；各高-q 壳层仍为正（最低-q 壳层仅有
  `-9.3e-6` 的微小数值项）。
- 启用当前 `high_energy_zero` anchor 后，`K-` 的最高 q 壳层变为
  `-0.0531052`，总和变为 `-0.0375540`，这正对应完整 phase-shift 结果的负密度。
- 同一高能 anchor 在低 q 节点把原始相位约 `-pi` 平移到约 `+pi`，而在高 q 节点
  原始高能相位已是正值（例如 `K-` 的 `1.3423`），于是又施加负的常数平移。不同
  q 的 `tail_before` 从约 `-3.135` 变到 `+1.342`，并非统一的高能极限。

对当前两种权重，
`F_current(delta)=delta`，`F_GBU(delta)=delta-0.5*sin(2delta)`，且
`F_GBU'(delta)=2*sin(delta)^2 >= 0`、`F_GBU(0)=0`。所以 GBU 变换不会自动消除
整体相位符号；当 anchored 相位在某个壳层为负时，正的 Bose 权重会保留该负贡献。
这说明当前 high-energy anchor 尚未满足统一的相位边界/Levinson 条件，不能用来
判定物理负介子密度，也不能通过简单全局翻转相位修复。

需要注意：`phase_shift_meson_number_density` 当前返回的 `status=:ok` 只表示
Bose 能量域通过；它尚未把 `density<0` 作为 helper 内部失败状态。扫描层另行以
`status=invalid_density` 标记非有限或负密度。因此当前“`ok` 但密度为负”是一个
可观测的接口语义缺口，后续应在不破坏旧返回结构的前提下增加显式 density-validity
gate，而不是在数值层静默裁剪。

### 9.3 与旧结果的边界

旧诊断资产 `meson_density_phase_e5_strict_candidates.csv` 已经显示：原始/未锚定
相移密度通常为正，而 `tail_shifted`/高能尾部平移版本可出现大幅负值（例如
`T=208 MeV` 的 pion 为约 `-278`）。因此本轮的有限-BQS `K-` 负值是同一类
相位边界未闭合问题在 ordered charged 背景中的新暴露，并非已经证明旧 PNJL 或
旧 BU 全部错误。当前应保留负值作为 `invalid_density` 诊断证据，先完成相位端点、
Levinson/Mott、节点和截断门禁；在此之前不改 production 默认，不把负值裁剪为零。

## 10. Strict phase backend implementation (2026-09-01)

本阶段新增 `src/relaxtime/ChargedPhaseBackend.jl`，把公式路线中的 strict
real-axis phase/BU 层收束为 solver-independent diagnostic backend。输入是有序
retarded 逆传播子 `Delta^R(omega,q)` 的 callable profile，默认相位约定为
`delta=-arg(Delta^R)`；profile 从高能端反向 unwrap，并显式返回高能端平移量和
tail span。提供 `strict_phase_gate` 组合阈下根计数、Levinson 阈值相位与高能 tail
稳定性；提供 `strict_charged_bu_density` 执行单电荷 `domega/pi` 导数型 BU 积分，
以及 `strict_density_convergence_gate` 比较节点/截断配置。

纯合成路径单测位于
`tests/unit/relaxtime/test_charged_phase_backend.jl`，覆盖相位对象/符号、branch
profile、endpoint tail、单电荷测度、显式 gate 失败和 coarse/refined convergence。
该实现没有接入 `MesonDensity` 默认入口，也没有把有限窗口 anchor 当作已证明的
`delta(infinity)=0`；真实 `ChargedRPAProvider` profile 的 eta、omega_max、节点、
Levinson/Mott 和 Bose-support 收敛仍是未决 production gate。

同时新增 `scripts/analysis/relaxtime/audit_charged_phase_backend.jl`，在固定
`FixedMuBConservedCharges` 背景上通过 `ChargedRPAProvider(:ordered_retarded)` 调用
该后端，输出四个 charged 通道的 gate/测度/节点配置。绑定态数暂以显式 `q -> 0`
诊断输入提供，任何物理 Levinson 失败都保留为 `accepted=false`，不会被当作生产结果。

本阶段的完成状态是“公式实现接口、合成测试和真实 profile 诊断已落地”，不是
production 升格。后续应在真实 profile 通过端点、Levinson/Mott、节点/截断和四算法
对照后，单独进行 production candidate review。

## 11. 真实 ordered profile 诊断（2026-09-01）

新增脚本已在同一有限 BQS quark-only 背景实际运行：
`T=170 MeV`、`mu_B=240 MeV`、`rho_Q/rho_B=0.4`、`rho_S=0`，平衡残差为
`1.44e-15`。脚本对 `pi_plus`、`pi_minus`、`K_plus`、`K_minus` 均使用
`ChargedRPAProvider(:ordered_retarded)`、完整 charged KMT 耦合和 strict
`domega/pi`，并比较了 coarse/refined 的 `q`、`omega` 节点及 `qmax/omega_max`。

低成本复现实验输出为未跟踪文件
`data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/strict_fixed_bqs_t170_mub240_low_refined_v2.csv`。
本次设置为 coarse `(eta=0.01, Pi_nodes=4, qmax=0.5, q_nodes=2,
omega_max=6, omega_nodes=4)`，refined `(eta=0.008, Pi_nodes=6, qmax=0.8,
q_nodes=2, omega_max=7, omega_nodes=6)`，`omega_min=0.5`。

结果边界如下：四个通道的密度均为有限正数，但所有 coarse/refined 行均为
`accepted=false,status=gate_failed`；`convergence_passed=false`。逐 q 诊断显示
阈下 root/Levinson gate 在各通道都失败，部分 profile 的高能 tail 也未稳定。
因此这次运行证明了真实 ordered profile 已接入并能保留失败诊断，但不构成
Levinson/Mott 或节点/截断 production 通过证据；临时 `q -> 0` 束缚态计数不能替代
物理绑定态判定。production 默认与旧 `MesonDensity` 路径均未改变。

## 12. 定向文献综述链接（2026-09-01）

本任务相关的公开文献、公式闭合、来源筛选、代码映射和后续门禁已整理到
`docs/analysis/relaxtime/charged_phase_literature_review_v1/`。该证据包明确区分
文献事实、项目数值约定和 production 授权边界；它不改变本任务的 `in progress`
状态，也不把 strict phase backend 晋升为 production。

## 13. 出版社原文与 PDF provenance（2026-09-01）

在用户明确授权后，使用 Chrome 的 XJTU 网关完成了 APS 访问验证，并对少数需要
出版社版本核对的条目逐篇保留结果。逐项记录见
`docs/analysis/relaxtime/charged_phase_literature_review_v1/tables/publisher_pdf_provenance.csv`：

- Rehberg et al. 1996（PRC 53, 410）的 APS 文章页显示 XJTU 授权，正式 PDF 20 页、
  `pdftotext` 可读，SHA-256 已记录，作为出版社版本证据保留在项目外目录；
- Dashen--Ma--Bernstein 1969 的 APS 文章页授权成功。此前仅按第一页文本首行判断为
  内容错配；重新渲染和逐页检查后确认，页 1 下半页包含目标论文标题、作者、摘要和
  Introduction，页 2--26 为目标论文，只有页 1 顶部残留另一篇液态合金电阻率文章的
  版面片段。因此 S02 可用于公式核对，但引用页 1 时必须避开残留片段并记录该 caveat；
- Hüfner et al. 1994 的 ScienceDirect 文章页显示 XJTU 标识，但明确写明 XJTU 不订阅
  该内容；此前入口触发的机器人 CAPTCHA 未绕过，也未声称取得出版社全文。

版权 PDF 不进入 Git、任务分支或项目分析包；仓库仅保存 DOI、出版社链接、授权路径、
文件 hash/页数、可读性探针和失败原因。上述局部出版社核对不改变 strict phase backend
仍处于探索性、未通过 Levinson/Mott 和收敛门禁的状态。

## 14. S02 重新核对后的公式路线影响（2026-09-01）

对保留的 S02 PDF 逐页渲染后确认：页 1 顶部存在版面残留，但目标论文从页 1 下半页
开始，页 2--26 完整连续。此前基于第一页首行的 `content_mismatch` 记录已修正为
`publisher_pdf_verified_first_page_overlay`，正文可用于 S-matrix/Levinson 公式核对。

本次公式复核得到三项需要纳入路线闭合的限定：

1. DMB 的直接对象是连通 on-shell (S) 矩阵和
   ((4\pi i)^{-1}S^{-1}\overleftrightarrow{\partial_E}S)。单道
   (S=e^{2i\delta}) 时才化为 (\partial_E\delta/\pi)；项目的
   `-arg(Delta^R)` 是否就是这个 δ，仍需独立的 propagator-to-S-matrix 归一化证明。
2. DMB 给出守恒 (B,I,S) 的独立化学势乘子，支持显式记录 μ_B、μ_I/μ_Q、μ_S，
   但不推导 PNJL BQS 闭合或 `mu_s=0.55 mu_u`。
3. DMB 的相对论性 Levinson 推广是基于非相对论结果的合理猜测；因此项目 gate 应标为
   conditional，并继续要求独立束缚态计数、阈值/高能端点和 Mott 补偿检查。S02 不改变
   PV/​i0、ordered charged bubble、KMT 或 Ω_M 反馈的其他文献结论。

## 15. 实轴 PV/cut 复核（2026-09-02）

本阶段在 `codex/charged-profile-gates` 分支增加了显式的实轴诊断适配层：
`OneLoopIntegrals.B0_pv_cut` 与 `ChargedRPAProvider` 的 `:ordered_pv_cut`。
它不修改历史 `B0`、`PolarizationAniso` 或 production 默认，而是复用已有 `B0`
的主值实部，并按四个 `tilde_B0` 复对数项各自的 `i0` 方向计算当前
`e^{-i omega t}` retarded 边界 cut；这不是一个全局虚部翻号：

```math
B_0^{\mathrm{PV+ret}}=\operatorname{PV}B_0+i\operatorname{Im}B_0^R.
```

纯代数/合成测试已锁定：`B0_pv_cut` 的实部与 `B0` 相同，虚部按 cut 区间逐项
确定；不同 `lambda` 区域可能相对历史值翻转，也可能保持同号。其在远离尖锐端点
的有限-`eta` 对照中与 `B0_retarded` 相符；
provider 元数据返回 `analytic_scope=:real_axis_pv_cut`、`eta_inv_fm=0`，并与
`:ordered_retarded` 的 `:upper_half_plane_probe`、两个 legacy oracle 区分。

在同一 `T=170 MeV, mu_B=240 MeV, rho_Q/rho_B=0.4, rho_S=0` 的有限-BQS
quark-only 背景上，低成本网格运行如下：

- `scripts/analysis/relaxtime/audit_charged_phase_backend.jl` 支持通过
  `CHARGED_PHASE_PRESCRIPTION=ordered_pv_cut` 与
  `CHARGED_PHASE_OMEGA_MEASURE` 选择诊断处方；CSV 额外保留
  `polarization_prescription`、Levinson residual、阈值相位范围和高能 tail span。
- PV-cut 结果保存在未跟踪的
  `data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/`，不纳入版本库。
- 在 `qmax=0.4/0.5 fm^-1`、`omega_max=7/8 fm^-1` 的 low/refined 对照中，四个
  `pi^±/K^±` 密度均为负且 `status=invalid_density`；对应有限-`eta` 运行给出相同
  量级的正密度，但仍为 `status=gate_failed`。两者均有 `failed_q_count=2`、
  `tail_failed_q_count=2`，Levinson/root gate 未通过。
- 这不是 Bose 支撑失败：所有运行都在正常相 `omega>mu_M` 窗口；也不是
  `domega/pi` 与 `domega/(2pi)` 的共同倍数。PV/retarded 切支差异会改变相位的
  连续谱方向，但在高能端点和束缚态计数未闭合前，不能把任一符号的密度当作物理
  结果，更不能做全局相位翻转或把负值裁剪为零。

先前同背景的 `charged_rpa_bu_negative_density` shell breakdown 进一步显示，
历史 GBU 路径在高能尾端施加常数平移时，K- 的最高 q 壳层可由未锚定的
`+0.0107` 变为 `-0.0531`；K+ 也出现由 `+0.0414` 变为 `-0.5447` 的壳层。
这只说明旧的非线性 GBU 权重对相位分支非常敏感，不能直接等同于 strict
`d(delta)/domega` 积分的误差。strict PV 运行的负值还必须结合 retarded/PV
切支、端点和独立 Levinson 计数共同判定。当前真实 profile 没有独立的
Mott 前后态对和物理束缚态计数，因此 Mott 只由合成单元测试覆盖，未宣称真实
Mott gate 已通过。

当前结论是：旧 `B0` 的实轴主值结构可以作为复核基础，但其 cut 符号不能在没有
Fourier 边界说明时直接称为本项目 retarded 值。严格 production 前仍需对
`eta->0+`、PV/cut、相位端点、独立束缚态计数、Levinson/Mott、节点、截断和
`omega/q` 测度做同一 profile 的联合门禁；本阶段不修改 PNJLCore、Omega_M、
旧密度入口或 production baseline。

## 16. 严格 phase 路线实施更新（2026-09-02）

### 16.1 显式 phase--S-matrix 适配

`PhaseNormalization` 新增 `propagator_phase` 和 `propagator_to_s_matrix`。它们把
当前诊断约定 `delta = -arg(Delta_R_inverse)`（或显式选择 propagator/符号）映射为
`S_diagnostic = exp(2im*delta)`，并返回所选对象、符号和映射标签。该适配只闭合
代数因子二，不声称任意 off-shell 传播子已经等于 on-shell 物理 S 矩阵；连续分支、
端点和 DMB 测度仍由独立 gate 负责。`tests/unit/relaxtime/test_phase_normalization.jl`
新增反向映射和对象选择测试。

### 16.2 独立束缚态计数与 q continuation

`BUPhaseGates.count_bound_states` 对给定 `q` 独立采样阈下逆传播子并复用简单实根
bracket 计数；`continue_bound_state_counts` 对显式 `q_values` 每点重新计数，返回
前一点差值。该 API 不从 phase unwrap 推断束缚态，也不把有限虚部丢弃；超过阈值时
返回 `status=:complex_subthreshold`、`passed=false`。strict density backend 现在
接受整数或此类计数结果，并在每个 `q_profile` 保留 `bound_state_diagnostic`。

真实 fixed-BQS ordered profile 已改为调用该独立计数器，而非旧的 `q -> 0` 占位。
在低成本设置（`eta=0.01/0.005`、`Pi_nodes=16/32`、`q_nodes=2`、
`BOUND_STATE_NODES=64/96`）下，四个通道均出现 `bound_state_status=complex_subthreshold`，
且所有 profile 仍为 `gate_failed`；这说明有限 eta 下无法认证物理束缚态数，不能把
符号变号计数写成 Levinson 的 `n_B`。

### 16.3 真实 Mott 前后配对

新增分析脚本 `scripts/analysis/relaxtime/audit_charged_mott_profiles.jl`，在两个
显式温度上配对真实 ordered profile，输出质量--阈值差、两侧独立计数状态、q 点计数
范围、Levinson 和 Mott gate。高温 pion 阈值低于普通 `omega_min` 时，脚本只在该诊断
profile 内自动降低窗口下界以覆盖阈值，并记录实际设置。

低成本 `T_before=170 MeV`、`T_after=230 MeV` 运行结果保存在未跟踪本地文件
`data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/strict_mott_profile_pair_low.csv`。
该运行显示 pion 的 `m_M-(m_q1+m_q2)` 从负值变为正值，具有阈值跨越迹象；但两侧
阈下计数均为 `complex_subthreshold`，K 通道的符号变号计数也随 q 改变。因此两侧
Levinson/Mott gate 均不能通过，不能将这组温度称为已定位的物理 Mott 转变。

### 16.4 联合收敛合同

`BUPhaseGates.joint_convergence_gate` 新增对多个样本的 pairwise 数值比较，并同时
要求 finite、显式 `accepted` 和（若提供）`tail_stable`。它保留 `eta`、q/omega 节点
和 cutoff 元数据，适用于下一轮真实 profile 的 PV--finite eta--节点--截断--端点联合
门禁。`audit_charged_phase_backend.jl` 现在同时输出
`joint_convergence_passed` 与 `joint_convergence_endpoint_stable`。当前真实 profile
因复阈下、Levinson 和 tail 失败，联合 gate 仍为 false；没有修改 production 默认或
旧 `MesonDensity`、PNJLCore、Omega_M 反馈。

### 16.5 当前阶段判定

本次代码与诊断只完成了严格路线的可审核基础设施和失败可见性：

- [x] 显式 phase--S-matrix 代数适配与纯测试；
- [x] 独立阈下计数、q continuation API 与脚本接入；
- [x] 真实 Mott profile 配对脚本与低成本诊断；
- [x] 多轴联合收敛 gate；
- [ ] 真实 profile 的 PV/有限 eta 极限一致性；
- [ ] 物理束缚态计数、Levinson/Mott 补偿；
- [ ] 节点、截断、端点联合通过；
- [ ] production candidate 评审。

所有新增数值 CSV 均保持为未跟踪本地诊断产物；当前任务状态继续为 `in progress`。

## 17. 阈值坐标与逐项 PV cut 修正（2026-09-02）

历史记录限定：本节的外频率平移正确，但 `E_1(q)+E_2(q)` 有限总动量阈值错误，
已由第 19 节修正；不能继续把本节当作有限 q 运动学闭合证据。

本阶段在独立分支 `codex/charged-phase-coordinate-cut-fix` 中完成两项最小闭合：

1. `ChargedRPAProvider.charged_pair_continuum_thresholds` 明确返回内部
   `lambda_thr=E_1(q)+E_2(q)` 和外部 `k0_thr=lambda_thr-(mu_1-mu_2)`。
   strict phase backend 采样的 `omega` 是外部 `k0`，审计脚本与 Mott profile 脚本
   已改用 `k0_thr`；CSV 同时保留内部阈值和 `threshold_coordinate=external_k0`。
2. `OneLoopIntegrals.B0_pv_cut` 不再采用未经证明的
   `Complex(Re(B0),-Im(B0))` 全局翻号，而是对四个复对数项逐项计算解析 cut，
   `k=0` 使用主值极点的残数边界。历史 `B0`、有限 `eta` 的 `B0_retarded`、
   legacy provider 和 production 默认均保持不变。

纯单元测试覆盖阈值平移、正负 `lambda` 的不同 cut 符号，并在有限 `eta` 探针远离
尖锐端点的范围内做数值对照。该修正尚未证明真实 charged profile 通过
`eta`、节点、截断、端点或 Levinson/Mott gate；若真实 profile 仍出现负密度，必须
继续检查相位对象/端点和独立束缚态计数，不能把负值裁剪为零。

## 18. Fig.2 同构负密度相位可视化（2026-09-02）

历史记录限定：下述原始数值保留用于追溯，不再作为已验证的壳层量纲/GBU 权重。
第 19 节独立复核发现额外 `1/T` 和跳变链式法则错误，原 CSV/PNG 未被覆盖。

为直接定位负密度来自相位对象、分支还是 q 壳层，新增分析脚本
`scripts/analysis/relaxtime/plot_charged_phase_negative_density_fig2_like.jl`。
脚本在同一个有限 BQS quark-only 背景
`T=170 MeV`、`mu_B=240 MeV`、`rho_Q/rho_B=0.4`、`rho_S=0` 上，使用
`ChargedRPAProvider(:ordered_pv_cut)`、完整 charged KMT 耦合和
`delta=-arg(Delta_inverse^R)`，对 `pi^±/K^±` 在
`q={0,0.25,0.5,0.75,1.0} fm^-1` 上输出密集频率剖面。`fold_0_pi` 只用于
Fig.2 风格的显示，不进入 BU 导数；真正积分分支仍是
`raw -> unwrapped -> high_energy_zero anchored`。

本次正式诊断的本地产物（均未跟踪、未纳入 production）为：

- `data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/negative_density_phase_fig2_like/charged_phase_profile_detail.csv`
- `data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/negative_density_phase_fig2_like/charged_phase_q_shell_summary.csv`
- `data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/negative_density_phase_fig2_like/charged_phase_fig2_like.png`
- `data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/negative_density_phase_fig2_like/charged_phase_negative_density_attribution.png`
- `data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/negative_density_phase_fig2_like/plot_manifest.json`

关键可复核事实：

1. q=0 的壳层因相空间前因子 `q^2` 恒为零，不能单独解释总密度；负值来自有限 q。
   在 q=1 fm^-1，current anchored 壳层为
   `pi+=-1.758e-2`、`pi-=-2.105e-2`、`K+=-1.163e-1`，而 `K-` 为
   `+2.335e-3 fm^-3 per fm^-1`。GBU 权重将其压低到
   `-1.592e-3`、`-2.231e-3`、`-6.197e-3` 和 `+2.598e-3`，但没有统一变号。
2. 在 K+、q=0 的 PV 曲线中，`Re Delta_inverse` 在阈值前跨过零，
   `Im Delta_inverse=0`；principal phase 在负实轴段取 `-pi`，阈值后取 `+pi`。
   反向 unwrap 后低能端约为 `2pi`，连续谱再降至高能端 0；折叠显示才呈现
   Fig.2 中熟悉的 `0 -> pi -> 0` 形状。
3. `unwrapped` 与 `anchored` 只相差常数，因此它们的相位导数逐点相同；高能
   anchor 本身不能消除 current 的负导数。raw 与连续分支在阈值跳变附近不同，
   但这属于相位边界/束缚态补偿问题，而不是可以直接全局翻转的符号错误。
4. 低节点有限-`eta` 对照另存于
   `.../negative_density_phase_fig2_like_ret_smoke/`；其 K+ 相位出现窄而密的
   branch oscillation，说明有限 `eta` 仍只能作为边界探针，不能替代 PV 的端点、
   Levinson/Mott 和节点收敛门禁。

因此，本图把当前负密度定位为“真实 PV 相位在有限 q 上的连续谱负导数与尚未闭合的
阈值/束缚态相位补偿共同造成的 diagnostic failure”。它排除了 `mu_M=0` 或 Bose
支撑奇点作为当前主因，但没有证明哪一个分支是最终物理 production 分支；下一步
仍应先完成物理束缚态计数、Levinson/Mott 补偿以及 PV/有限 `eta`、节点和截断联合
收敛，禁止裁零或切换 production 默认。

## 19. 独立接手复核与 gap/split 诊断（2026-09-05）

### 19.1 范围与审查结论

工作树固定 `D:\w\jrt-ord`，分支 `codex/charged-phase-coordinate-cut-fix`，
起点 `fa34a3975d858ca5c517f3490476178720f0d23e`。task-ledger preflight 通过，
主 track 仍为 `formula-route-closure`；本轮属于该 track 的 `required_follow_up`
及 physical-closure `research`，不改变其他任务、PR base 或授权状态。
PR310 是 OPEN 的探索性分支，原 CI 不覆盖本轮未提交改动。
最终在线核对：base 仍为 `codex/charged-profile-gates`，远端 head 仍为上述 SHA，
14/14 个 check run 为 SUCCESS，mergeStateStatus=CLEAN。暂存区为空，本轮没有
提交、推送或合并；CLEAN 不代表物理闭合或 production 授权。

独立发现并修正的代码问题：

1. `ChargedPhaseBackend` 和 Fig2 脚本把导数型 `g_B delta'` 积分再除以 T。
   这导致 density 具有 `fm^-2` 而非 `fm^-3`，壳层也少一个能量量纲。
   已去掉多余因子，加入全部能量同时放大 2 时密度放大 8 的合成测试。
2. 有限总动量阈值应为 `hypot(q,m1+m2)`，而不是 `hypot(q,m1)+hypot(q,m2)`。
   外频率平移仍为 `k0_thr=lambda_thr-(mu1-mu2)`。新增极小化及“真实阈值与旧阈值
   之间已存在非零 cut”的测试；q=0 与旧接口默认不变。
3. Fig2 脚本将 `2 sin(delta)^2 delta'` 用到未解析的 pi 跳变，两端正好落在
   整数 pi 时会漏掉整个 GBU 离散贡献。改为先对 `delta-sin(2delta)/2` 求差；
   新 split API 则直接加认证根的 Bose 权重，完全不让 unwrap 决定根数。
4. strict gate 原先没有核对采样根数与显式根数，也未消费独立计数结果的 `passed`。
   现已补拒绝路径；Mott helper 拒绝不同相位对象/符号/anchor target 的配对。
   联合 convergence helper 遇到 NaN 现在返回失败诊断，不在已检测非有限后继续抛错。

原 Fig2 23 项全为存在性/字符串契约。接手时独立执行四组 tests 得到 92/92；
这说明旧记录中的通过数不能代替数值语义验证。本轮新增了可执行 shell 代数测试，
不加载绘图库、不调用平衡 solver。绘图脚本现在拒绝覆盖非空输出目录，修正自定义
q 列表时仍硬写 q=0 的标签，并记录量纲、阈值、源码 hash 和 PV 实际固定 16/32 节点。

### 19.2 文献事实、历史实现与数学推导

- 文献定位：Blaschke 2020 的广义替换仍为 `delta-sin(2delta)/2`，见公式包所列
  DOI `10.3390/particles3010014`。本轮没有重新取得或重新判读论文 PDF，不新增
  文献关于有限 BQS、有限 q 或 positive density 的主张。
- temp7 当前 CSV 有 8 组、每组 241 点；全部 principal 与 forward-unwrapped
  完全相同。`fold_0_pi` 是最后一步展示变换，不是连续物理积分分支。
  其配置为 T=90 MeV、q=0、mu_u=mu_d=mu_q、mu_s=0.2 mu_q、cutoff=602 MeV；
  图示 mu_q=350 MeV。`paper_2020` 标签下 K+=(s,u)、K-=(u,s)，与当前项目标签相反。
  temp7 使用 `atan(-Im inverse_from_legacy_B0,Re)`，也不等同于当前逐项 PV-ret cut。
- 数学推导：解析简单根附近 `F(z)=a(z-omega_b)`，
  `-Im d log F/domega` 的离散部分为 `+pi delta_D(omega-omega_b)`，与非零实斜率 a
  的符号无关。零 eta 实轴的 signed-zero/atan 不能自行选择这个定向。
- 项目诊断约定：`certify_gap_roots` 只认证外部指定的解析 gap 内简单根；
  `count_scope=provided_analytic_gaps`、`completeness_certified=false`。
  `continue_gap_roots` 给出位置、ID 与未恢复事件，不把 count change 自动升级为 Mott。
- anchor 限定：普通 BU 的 `d delta` 对任意常数平移不变；GBU 的 `dF(delta)`
  只对整数倍 pi 平移不变。合成连续谱测试明确验证一般常数 anchor 会改变 GBU，
  因此不能把普通 BU 的 anchor 不变性推广至 GBU。
- 公式缺口：provider 的 `lambda=k0+mu_M` 与热权重必须同步转换。
  `g_B(k0;0,T)=g_B(lambda;mu_M,T)`；历史 `g_B(k0;mu_M,T)` 再减一次化学势。
  新审计同时输出两者，旧 Fig2 入口只明确标记历史口径，未静默改默认。
  完整 charged KMS/正负频率组织、physical S 映射和多 cut Levinson 仍待评审。

### 19.3 新保留证据与数值边界

新入口：`scripts/analysis/relaxtime/audit_charged_gap_spectrum.jl`，读取原 Fig2
CSV/manifest，显式重建 A（64 节点、pmax=16 fm^-1），用每通道 24 个分散探针
检查逆传播子，最大绝对差为 **0.0**。除生成 T=230 的新 BQS snapshot 外，gap
搜索、split spectrum、温度配对和 temp7 overlay 均不调用平衡 solver。
T=230 snapshot 的 quark-only 残差为 `4.743e-14`，不是 production 背景收敛证明。

所有新产物保存在未跟踪的
`data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/`，不暂存：

- `gap_spectrum_T170_review_v3/`、`gap_spectrum_T230_review_v3/`：根、分段相位、
  有符号壳层、输入 hash、执行时源码快照和 manifest。
- `gap_mott_pair_review_v3/`：按相同 channel/q 唯一配对的 conditional Mott 表。
- `temp7_controlled_overlay_20260905/`：未拟合平移的 principal/unwrapped/folded
  受控 overlay，逐曲线标出 flavor 与 q，manifest 列出两套温度/化学势与相位规则。
- 原 `negative_density_phase_fig2_like/` 及所有历史 CSV/PNG/manifest 保持不变。
  `gap_spectrum_audit_20260905/` 和 `gap_spectrum_endpoint12_20260905/` 是早期
  窗口诊断，不取代带源码快照的 review_v3；v1/v2 同样保留而不覆盖。

v3 校验：T170/T230 分别有 36/0 行根、72/72 行壳层和 23508/25668 行连续谱；
无根 CSV 仍有完整表头。两组各 5 份执行时源码快照均与 manifest 及当前源码的
SHA-256 匹配。Mott 表为 36 个唯一 channel/q 配对，全部保留
`physical_mott_certified=false` 与 `production_authorized=false`。

T=170 时 `q=0:0.125:1 fm^-1` 的四通道共 36 个 gap 搜索均找到一个简单实根；
初步 128 节点与 256 节点计数一致。没有在该 q 范围观测到根消失/进入 cut。
根的外频率位置为：

| channel | q=0 的 k0 [fm^-1] | q=1 的 k0 [fm^-1] |
|---|---:|---:|
| pi+ | 0.7258110844 | 1.3603251997 |
| pi- | 0.6383697121 | 1.2706461080 |
| K+ | 2.1682351339 | 2.4679001872 |
| K- | 2.8520675806 | 3.0894346860 |

这些是 gap 内根，不是“所有阈下 root 已完备认证”。unitary 阈值相位检查在
omega_max=8 时 K- 的 q>=0.25 失败；延伸到 12 后全部 36 个 conditional gate
通过，最大残差约 `3.294e-3 rad`。因此前者至少包含有限高能窗口错误，不能只
用固定末四点 tail span 判断。这里的 endpoint 检查服务于 Levinson 定义，
不等于已启动或通过 eta/节点/cutoff 联合收敛。

T=230 的相同 36 个 gap 搜索均未找到简单根，unitary 阈值相位接近 0。
T=170/230 的配对给出 gap count `1->0` 与阈值相位下降约 pi；这只通过条件
计数/端点合同，**未定位真实根进入 unitary cut 的温度，也未认证物理 Mott 转变**。

在原低端 `k0_min=0.3620384126599023 fm^-1`、高端 12、q=1 下，新分离式
`g(k0;0)` 的 current 部分和依次为 `0.0096052, 0.0134741, 0.00173166,
0.00641242 fm^-2`（pi+,pi-,K+,K-）。这是截断谱窗的正部分和，不是总密度。
沿用历史 `g(k0-mu_M)` 时 K+ 部分和仍为 `-0.0283861 fm^-2`。
将低端另取 0.05 的独立窗口诊断时，Landau 区域可重新压过正的离散项。
因此不能把一组正部分和宣称为“负密度已修复”。正常 Bose 支撑通过不等于
低能端点积分已闭合；不可通过抬高下界、fold、全局翻号或裁零隐藏这个差别。

### 19.4 验证、分类与停止条件

- focused unit：522/522，含新增 gap/simple-pole 排除、q continuation、离散 pi 的
  current/GBU 稳定权重、能量坐标平移、量纲缩放与旧 oneloop/meson density 测试。
- 旧 meson mass 与 density regimes regression：81/81；没有改 baseline 或旧默认。
- docs、formula-route、script governance、data-output guard、active-doc governance
  均通过；`git diff --check` 与 task-ledger preflight 也通过。
- 图生成成功，但本轮图像工具未向模型提供可判读图像，不能声称完成主观版式或
  与论文像素细节的视觉验收；CSV 数值与背景/标签合同已独立检查。
  overlay 为 1800x590，T230 两图为 1640x480 和 1480x500；抽样像素分别有
  358/329/289 种颜色，均非空白。这不代替文字重叠与曲线细节的视觉检查。

最终分类：额外 1/T、有限 q 阈值、GBU 离散跳变和 fail-open gate 属于代码缺陷；
KMS 坐标及多 cut/physical-S 归一化属于公式缺口；omega_max 与低能 Landau 支撑
属于未收敛/未闭合数值问题；quark-only、有限 cutoff/routing 及非正相关谱权重
是否构成最终模型限制，当前证据不足以定论。

下一步责任属于本 strict-audit track：先按受控温度 continuation 定位根进入 cut
的事件，并补正负频率/KMS 与 Landau 多 cut 计数/边界公式。只有这一层闭合后，
才启动 PV 与 eta->0+ 的分布/积分极限及 omega/q、cutoff、端点联合门禁；
不要求 PV 与有限 eta 在尖锐端点逐点相等。PR310 保持探索性 OPEN；本轮不自动
提交/推送或合并，留待审阅这些超出原 Fig2 契约测试的新缺陷及条件证据。

## 20. 两线谱与低能边界推进（2026-09-05）

### 20.1 本轮范围与来源

用户授权继续查文献、运行和修改现有实现；本轮仍只修改 `D:\w\jrt-ord`，
HEAD 保持 `fa34a397`，主 track 保持 `formula-route-closure`。
新工作属于 `required_follow_up`，并识别出先于 phase/Mott 的 analytic-cut blocker。
未运行新的 BQS 平衡、未加入 Omega_M/反馈、未改变旧 MesonDensity 或 baseline，
未暂存 data/outputs、未提交、推送或合并。台账只调整当前任务的 next_action。

联网状态已更正：内置搜索工具的 gateway 不兼容，不等于不能联网。直接访问
`https://arxiv.org/abs/1912.13162` 返回 HTTP 200，INSPIRE 公开 JSON API 也成功。
本轮直接下载并读全文的来源保存在 `work/charged_phase_low_energy_20260905/`：

| arXiv | 使用依据 | PDF SHA-256 |
|---|---|---|
| [1305.3907v3](https://arxiv.org/abs/1305.3907v3) | Eq. (82)-(90)、(105)、Appendix A.3；频率、log-propagator BU、计数和两线谱 | `196fbd1c9b07705f6f5642f725e0d09a89b18d8617d8241fc8328fb8621bf09d` |
| [1912.13162](https://arxiv.org/abs/1912.13162) | q=0 外推假设、full/no-anomalous 比较及广义替换 | `dc5c4ef922053e19c849759ac1cdc63a016b673577189b1cf77550d43f0737a5` |
| [2301.09882](https://arxiv.org/abs/2301.09882) | 有限 q LD/QP 压强分离与 Appendix A 谱性质 | `b54777da443d241b2c208cd8a1be0587dab6fc530b15658743b6e2bc60a2df9a` |

2014 论文页 14 已渲染，但本轮模型仍不能判读返回图像；方程核对使用正文提取及
独立代数/极限测试，不声称做完图像视觉验收。

重要文献修正：2014 Gaussian BU 直接来自 `Tr log D^-1`，Eq. (90) 后明确说
该相位不是可观测 on-shell 散射相移。因此第 19 节把 physical-S 映射列为必需
阻塞项过强；本轮撤下该前置要求，但保留真正的复泡解析性、热权重及态数检验。
2014 的 Gaussian BU 与 2020 的 `delta-sin(2delta)/2` 比较权重不能仅因同称
generalized BU 就视为同一推导；后者不自动修复前者的非物理 cut。

### 20.2 先于相移的失败：真空 spacelike 与零外频率

新入口 `scripts/analysis/relaxtime/audit_charged_phase_low_energy.jl` 完全使用冻结
Fig2 CSV/manifest；A64、pmax16 重建后，每通道 24 个逆传播子探针最大差仍为 0.0。
新输出位于同一 analysis 根下 `low_energy_cut_audit_review_v2/`，包含 4 行真空
检查、96 行低能探针、96 行边界积分，以及六份执行时源码快照和输入哈希。
早期 `low_energy_cut_audit_20260905/` 保留，未被覆盖。

等质量真空探针取 `m1=m2=q=1`、`mu1=mu2=0`、`T=0.001 fm^-1`：

| lambda | 区域 | 旧 PV Im B0 | 旧有限 eta 探针 | 独立两线 Im B0 |
|---:|---|---:|---:|---:|
| 0.1 | spacelike | 0.3141592654 | 0.3141734766 | 0 |
| 0.5 | spacelike | 1.5707963268 | 1.5722077678 | 0 |
| 0.9 | spacelike | 2.8274333882 | 2.8252273462 | 0 |
| 3.0 | unitary | 2.2214414691 | 2.2098013786 | 2.2214414691 |

有限 eta 这里只是 `eta=1e-4`、512 节点的反例旁证，不是 eta 极限收敛认证。
无热占据的真空不应有 Landau occupation difference；旧结果却包含 `pi*lambda/q`
量级的虚部。因此 PV 与有限 eta 相互接近只证明共享表达一致，不能证明物理正确。
不能通过对此项直接减一个拟合函数来修复，必须重新闭合同一 regulator 的复泡。

冻结 T170 BQS 的 `k0=0,q=1`，旧 Im B0 为 pi+/- 的 `-/+0.1421148963`，
K+/- 的 `-/+0.6755933468`；新两线占据差检验均为 0。
on-shell 条件 `sE1-tE2=k0+mu1-mu2` 在 k0=0 要求两线的带化学势能量相同，
所以相同 PNJL 背景分布的差应抵消。此检验不涉及相位折叠或根数选择。

第 18/19 节的 gap 根和 conditional Mott 数值仍是保留表达式的代数诊断，但不能
在 provider 物理 cut 检验失败后继续把它们称为已认证的物理谱。
当前“负密度是否是模型限制”尚无定论；先定位为 finite-q analytic continuation /
regulator closure 缺陷。q=0 的 temp7 图形复现不检验这个有限 q 真空反例。

### 20.3 新 oracle 与代码修复

- `OneLoopIntegrals.B0_spectral_cut`：从两线 Matsubara residue sum 出发，独立
  解析角向 delta 支撑并积分 occupation difference；返回 pair/Landau/total 虚部。
  两条线都限制在 pmax 内，明确不同于旧 shift 后的一线 cutoff。
- 该 oracle 只有虚部，未实现对应实部，不能直接替换旧 B0，不能把它与旧 PV
  实部、旧 tadpole A 拼接后声称新 RPA 已闭合。下一步须从同一 regulated loop
  同时推导完整 Pi 的实虚部与 tadpole/contact 边界，验证 q=0 匹配及必要校准。
- `ChargedRPAProvider` 返回 `physical_cut_certified=false`，明确未认证状态；
  该字段是上游诊断元数据，不能代替调用方的实际 production 授权检查。
- `bu_phase_integral_parts` 显式报告体项、边界项和离散乘积恒等式残差。
- `propagator_phase` 拒绝精确复零点；不再让 atan(0,0) 给出伪相位。
- 修正 `ChargedPhaseBackend` 顶部旧公式注释残留的额外 1/T。

### 20.4 省略边界项为何会得到伪正值

对固定 K+、q=1、原 provider 的低能段，采用 `g(k0;0)`：

| 窗口 [fm^-1] | 权重 | 导数积分 | 分部积分体项 | 边界项 |
|---|---|---:|---:|---:|
| [0.01,0.05] | current | -0.329906 | 3.070675 | -3.400581 |
| [0.001,0.05] | current | -0.815905 | 40.933153 | -41.749059 |
| [0.01,0.05] | GBU | -0.012268 | 0.040278 | -0.052546 |

表值是初轮 257 节点的无量纲能量积分，不是总密度；v2 用 513 个对数节点复核。
旧相位低端 `delta(0)=0.1553871`，非零斜率使 current/GBU 都有负的对数下界依赖。
这说明“保持正相位、使用正的分部积分体项”若没有证明边界项消失，就可以得到
与原导数积分符号相反的结果。不能将这种正值用作 strict 密度已修复的证据。

条件性数学推导：若正确 provider 在正常相满足 `delta=c*k0+...`，则
`F_GBU=(2/3)*c^3*k0^3+...`，其 `g*F_GBU'` 低能项为 O(k0)，不像一般
`g*delta'` 那样有 O(1/k0) 项。这提示 GBU 的低能极限可能改善，但前提是先
验证正确的零频率相位与对应近似，不能用本轮已失败的旧 provider 代替验证。

### 20.5 验证与下一步

最终九组 focused unit 为 **1107/1107**，包括等质量有限 T 的 Fermi-Dirac
Landau 解析积分及双线 cutoff 端点项；旧 meson mass/density regimes regression
为 **81/81**。docs、formula-route、script governance、data-output guard、
active-doc governance、task-ledger preflight 和 `git diff --check` 均通过。

v2 的 6 份源码快照与 manifest、当前源码 SHA-256 全部一致，行数为 4/96/96；
257 到 513 节点的导数积分最大绝对差为 `4.732e-4`，离散分部积分残差最大
`1.596e-11`（大体项与大边界项的浮点抵消）。这只验证反例稳定性，不是生产收敛。
PR310 在线核对仍 OPEN、base=`codex/charged-profile-gates`、head=`fa34a397`，
14/14 远端 checks 成功但不覆盖本轮 dirty 修改。暂存区为空。
未做新的 production baseline、完整 q 积分或联合 eta/cutoff 数值生产。

下一步优先级由“继续调 unwrap / 继续 Mott 扫描”调整为：

1. 同一 cutoff/loop routing 下建立完整 causal Pi，而不拼接不同 regulator 的实虚部；
2. 验证真空无 spacelike cut、零外频率 occupation cancellation、正负频率关系；
3. 在该 provider 上重做 q=0 文献点、有限 q 根/连续谱、低能边界及 Mott；
4. 最后比较论文 q=0 外推近似和直接 finite-q 路线，不将 folded/no-anomalous
   结果升级为 strict 物理闭合；完整收敛与 production 授权仍独立。

## 21. 完整因果谱泡与研究密度路线（2026-09-05）

### 21.1 选择及物理解释

用户明确要求继续推进到“文献级研究可行方案”。选择为：固定 BQS quark-only
背景，使用同一双线 regulator 的完整 ordered P/S 谱泡，以配套 Cauchy/PV
重建实轴，独立计数解析 gap 根，将稳定粒子项、unitary 和 Landau 连续谱全部
计入 GBU partial yield。公式推导见路线包 §5.2.1；不修改旧入口或生产 baseline。

“可行”指该固定背景后处理近似有明确公式、可运行数值及逐层验证，不指
已获得唯一介子守恒数、完整热力学反馈或全冻结线 production 授权。
普通单正频率 BU 在有限 q Landau 区可有红外对数依赖，不将它与 GBU 一起
宣称为已修复。2020 的广义替换为本候选的明确近似；其 q=0 Lorentz 外推可以
作为另一个比较模型，但不能覆盖本轮直接 finite-q 的 regulator 选择和误差。

文献本轮核对：arXiv 直连 HTTP 200；内置 web gateway 仍返回不兼容错误，
没有可调用的 INSPIRE MCP。使用已有 hash-bound 2014/2020 正文，2014 的
A.18/A.36 支撑顶角投影和 occupation difference，2020 保留 PDF 的 Eq. (14)/(16)
支撑 partial density/GBU 替换；原路线表中的开放版本 Eq. (18)/(20) 编号已改正。
没有新作出版社访问或声称完成图片视觉判读。

### 21.2 实现和受控失败

- `CausalSpectralBubble.jl`：完整两线 Pi、B0、导数及同域 C(q)，并提供独立
  解析 cut 与该谱的分段线性精确 Cauchy/PV 变换。q=0 的 C 与旧 A1+A2 匹配；
  有限 q 不把旧常数 A 强塞回新泡。
- 两个 cutoff 显式保留：真空两线 Lambda，热两线 L_th。默认同截断；本轮
  BQS 研究诊断使用 L_th=16 fm^-1，以匹配上游的真空截断/热项延伸。
  双线硬截断导致有限 q 的模型依赖，不能靠增大数值节点消除。
- `OneLoopIntegrals.B0_spectral_cut` 加入 full/vacuum/thermal 分量及精确运动学
  排除；旧 B0/PV/log 数值语义不变。
- `audit_causal_charged_bu.jl`：冻结输入，不调用平衡 solver；独立 gap 根、
  分离壳层、上下边界、q 求积、输入/源码哈希及运行快照。
- 原 v1 在落盘时暴露局部变量遮蔽 `cp`，已修正；仅留下部分空目录，不作证据。
  v2 保留了两处 `complex_gap` 失败，不作为闭合积分。
- v2 的失败来自阈值坐标往返舍入和截断端点的支撑歧义：微小 cut 值被插值
  进整个解析 gap。v3 改为内部 lambda 直接生成谱节点，并先执行精确运动学
  包络排除，没有放宽 gap 容差。近等质量的 gap 回归测试固定了此问题。
- 首轮测试中任意“有限 eta 虚部足够小”断言不成立，改测无热 cut 的实轴零值
  和有限 eta 尾部的线性衰减；不是要求有限 eta 与 PV 逐点相等。

### 21.3 初步数值证据

所有产物位于未跟踪的 `data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/`。
`causal_bu_review_v3/` 为修正后的 mesh64 诊断；`causal_bu_review_v4/` 使用
mesh128、直接泡 p160/x80、omega1200，omega_max=40 fm^-1。旧文件均保留。

在 T170、q=0/0.5/1 的 12 个通道点，独立 gap 根均为一个；阈值相位接近 pi，
高端超出完整截断谱支撑后为零。v4 的直接二维泡与独立谱重建在 eta=0.4 探针
上的最大绝对差为 `1.926e-4 fm^-2`，主要在 q=0 的硬截断跳变，仍需 mesh 收敛。
这是真正独立表示的比较，不是旧 PV/log 的相互 parity。

q=1、低端 1e-4、上端40 时的 GBU 壳层（`fm^-2`）为：

| pi+ | pi- | K+ | K- |
|---:|---:|---:|---:|
| 0.002372896 | 0.002643042 | 0.001103055 | 0.0004933954 |

全部包含精确离散 Bose 权重、有符号 unitary 以及完整正外频率 Landau 区，
没有 fold、全局翻号、裁零或删除 anomalous。低端从 1e-2 降到 1e-4，GBU 壳层
变化最高约 `1.54e-10 fm^-2`；普通 BU 仍有可见对数增长，例如 pi+ 壳层
`0.0042474 -> 0.0073100`。因此“新 provider 已修复虚假真空 cut”与“普通
BU 的单正频率数目积分仍不闭合”是同时成立的两个结论。

第一轮 q8、qmax4 积分得到四通道正 GBU 密度，且逐 q 条件态数 gate 通过；
这还不是 q 节点/尾部收敛。有限 q 根在约 q=3.05 到3.59之间由一个变为零，
阈值相位同时由 pi 变为零；该事件依赖当前三动量 regulator，不能直接叫作
无 cutoff 依赖的热 Mott 温度。更细 q 积分与 T230 冻结对照的验收见第 21.5 节。

### 21.4 未决边界

核心需保持：真实数值只作 research diagnostic；相位部分数密度不是新上游
热力学闭合；双线硬 cutoff/routing 是模型近似，不是可通过误差容差掩盖的项。
完整冻结线、无 Bose/不稳定域、全谱根完备性、热 Mott 定位及跨 regulator
比较仍需独立证据。PR310 不合并、不晋升 production，本轮源码尚未提交。

### 21.5 本轮完成的代表点验收

`causal_bu_q8_qmax4/` 与 `causal_bu_q12_qmax4/` 隔离比较 q 节点；GBU 密度
最大相对变化约 0.10%。`causal_bu_q12_qmax6/` 扩大 q 上界，最高变化约
0.88%，说明 qmax4 尚有可测尾部影响，不能只看 q 节点稳定。

最终 `causal_bu_q16_refined/` 使用 q16/qmax8、谱 mesh256、omega2400/
omega_max48、热 cutoff20（真空 Lambda=3.0522942112582996 fm^-1 不变）。
相对于 q12/qmax6、mesh128、omega1200/max40、热 cutoff16 的联合设置：

| channel | GBU density [fm^-3] | 联合变化，相对最终值 |
|---|---:|---:|
| pi+ | 0.005380468682 | 0.2055% |
| pi- | 0.005985903182 | 0.1994% |
| K+ | 0.002961524084 | 0.2119% |
| K- | 0.001337759592 | 0.1807% |

这些是低端1e-4的 research diagnostic，未注册为 production baseline。
联合变化不是完整误差上限：它同时改变多轴，尚未逐项量化模型 regulator
系统误差、所有 q 接近阈值事件的网格误差或全参数域误差。最终组包含256行
因果对照、28个 gap 根、384行壳层、305120行 profile；64个 channel/q 的
条件态数 gate 均通过。该组直接二维泡与独立谱重建的最大复值差为
`1.238e-5 fm^-2`（eta=0.4 探针），不是有限 eta 的实轴逐点收敛声明。

`causal_eta_integral_review/` 在 q=1、同一 mesh256 谱上比较完整有限 eta
GBU 积分与独立“PV连续谱+离散根”。四通道各取6个 eta，从0.02降到
0.000625 fm^-1，24行结果的每通道误差均单调下降；末点相对误差分别为
`4.734e-5,4.793e-5,5.628e-5,5.293e-5`。必须先在固定正低端做 eta 极限，
再在 PV GBU 中移低端；有限 eta 可使 delta(0) 非零，不能交换极限后强行
要求其红外形式与 eta=0 相同。该检查不要求尖锐跳变或阈值逐点相等。

`causal_bu_T230_review/` 复用已有 T230 quark-only snapshot，没有新求解。
与 T170 的 q=0/0.5/1 对应，12个 gap count 均 `1 -> 0`，阈值相位
`pi -> 0`，T230 阈值相位/pi 最大约 `8.18e-8`；q=1 的四个 GBU 壳层
仍为正。它支持两端点的条件 Mott 补偿，不定位中间温度、不认证根数全谱完备。

最终测试：10组 focused unit **1331/1331**，其中新核 **163/163**；旧
meson mass/density regimes regression **81/81**。首轮 regression 命令
误写文件名，在执行测试前失败；改用实际 `_regression.jl` 文件后通过，未改
baseline 或容差。docs、formula-route、relaxtime-script、data-output、active-doc、
task-ledger preflight 和 `git diff --check` 均通过。七组有效主诊断各9份源码
快照、eta诊断10份快照，全部与各自 manifest 的 SHA-256 匹配。

最终判断：旧 finite-q analytic-cut 缺陷与此前量纲/阈值/离散计数缺陷可定性为
代码/解析实现问题；完整因果谱泡与 GBU 分离式已形成可运行且有代表点验证的
研究候选。普通单正频率 BU 的红外、partial yield 与守恒荷导数的区别、有限 q
regulator 依赖属于必须公开的公式/近似边界，不因数值变正而消失。下一步以
该候选做精确 Mott continuation、regulator/论文外推的受控比较及冻结线覆盖；
无需继续用 unwrap/fold/anchor 试图修补已知坏 cut。PR310 仍 OPEN、base不变、
远端fa34a397、14/14原checks成功但不覆盖本轮修改；未暂存、提交、推送或合并。

## 22. 方法 v1、受控比较与冻结线前置验证（2026-09-05）

### 22.1 本轮定义与授权

作者批准按“固定方法、同输入比较、Mott/误差、稀疏冻结线”顺序推进。
新增研究合同 `config/models/pnjl/charged_gbu_research_v1.toml` 和
[方法文档](../../reference/formula/relaxtime/ChargedGBU_ResearchMethod_v1.md)。
observable 为固定 quark-only BQS 背景上的 GBU partial yield，不是含介子反馈的
完整热力学导数；真空双线 cutoff、热项延伸、ordered flavor、external k0/Bose
坐标、独立根及有符号连续谱均显式固定。没有改变 PNJLCore、旧默认或 baseline。
本工作仍是同一 strict track 的 required_follow_up，不晋升 production。

计算入口为 `audit_causal_gbu_comparison.jl`、`audit_causal_gbu_validation.jl`、
`audit_causal_gbu_routing.jl`，公共契约/围道计数在 analysis-only
`causal_gbu_research_utils.jl`。所有输出都在既有 charged_rpa_phase_backend
analysis 根下的新 `method_v1_*` 目录，拒绝覆盖并留存执行时源文件快照。

### 22.2 同背景 q=0 与外推比较

`method_v1_comparison/` 固定 T=170 MeV、muB=240 MeV 的保留 BQS 背景，
mu_u/d/s=(0.3742643949,0.4209954887,0.0142259822) fm^-1；四通道使用
(u,d)/(d,u)/(u,s)/(s,u)，耦合采用当前 full-KMT K12/K45，不使用 temp7 的 GS。
这不是论文 T90、mu_q=350、mu_s=0.2mu_q 的原背景图形重现。

只读加载 temp7 原 `audit_fig2_phase_shift.jl` 的 i1/i2 函数；i2 调用当前保留
`OneLoopIntegrals.B0`，不是声称重建了 temp7 历史二进制依赖。归一化固定为
Pi_project=pol_temp7/2，再使用同一 D=1-4K*Pi，未拟合符号或 fold。
在同真空/热 cutoff 下，164 个 q=0 探针的最大实部差 `4.66592082e-4 fm^-2`、
虚部差 `1.50725404e-4 fm^-2`。这是有限网格的一致性证据，不是精确相等证明。
单独延伸热项的 Pi 差最高 `0.01618525 fm^-2`，不能把它算作代码误差。

同一新核进一步比较 direct finite-q 和 lambda-invariant timelike q0 reference：

| channel | direct density [fm^-3] | q0 reference density [fm^-3] |
|---|---:|---:|
| pi+ | 0.005380468682 | 0.05718618896 |
| pi- | 0.005985903181 | 0.06447132038 |
| K+ | 0.002961524084 | 0.02719845842 |
| K- | 0.001337759592 | 0.01203024323 |

两路线使用相同 q16/qmax8、mesh256、omega2400/max48、热 cutoff20。
128 个 channel/q/route 条件检查均通过。K+/pi+ 为 0.5504212 对 0.4756124，
K-/pi- 为 0.2234850 对 0.1865984。密度相差约一个数量级，不是节点误差。
reference 是文献 Lorentz 外推思想在当前内部 lambda 坐标的项目实现，不是
temp7 外频率公式的逐字移植；其 spacelike 区没有建模，不能用零值删除 direct LD。
参考边界相位不为零时会报告 gate 失败。

### 22.3 独立根计数、热 Mott 与动量消失

在 normal kinematic gap 内，根扫描 128/256 节点并与复围道零点数比较。
合成 profile 测试包含偶重根、近邻根、域外根、复零点、边界零点和近边界根；
围道绕数只计数，不为密度选择相移分支。真实结果只覆盖所围解析域，不能证明
额外 cutoff gap、cut 内低能结构及远离实轴的不稳定零点不存在。

`method_v1_mott/` 使用 76 个 quark-only 背景，muB=240、p24、q=0/0.5/1、
mesh128/256。24 个热 bracket 均通过分离端点的根数 1->0 与阈值相位减少 pi。
mesh256 的 q=0 bracket 为：

| channel | T bracket [MeV] |
|---|---|
| pi+ | [203.6328125, 203.75] |
| pi- | [203.28125, 203.3984375] |
| K+ | [201.7578125, 201.875] |
| K- | [204.1015625, 204.21875] |

区间宽 0.1171875 MeV；阈下偏移固定 1e-7 fm^-1，另记录 1e-9 对照。
计数/相位端点位于 bracket 外约 0.5 MeV，不等于在无限接近阈值处完成计数。
K-、q=0.5 的两种 mesh bracket 相差一个网格区间，其他 bracket 重合。

原 q=0:0.25:4 的 68 个检查中 67 个通过，K-、q=3.5 围道分辨率不足。
`method_v1_q_refine/` 用自适应围道节点加倍修复该诊断分辨率问题，没有改容差；
该点根数/围道数均为 1，根与阈值距离约 1.28246e-4 fm^-1。
八个 mesh/channel 动量 bracket 的根数和相位端点检查均通过。mesh256 时，
pi+/pi-/K+/K- 的 q 消失区间分别是 [3.130371,3.131104]、[3.121582,3.122314]、
[3.461426,3.462158]、[3.535400,3.536133] fm^-1；mesh128->256 位置变化
0.00366--0.00513 fm^-1，大于二分区间宽。不得把区间宽当作总误差。

### 22.4 数值误差与 regulator 依赖分账

`method_v1_backgrounds/` 比较 T170/muB240 的 p8/24/48/96：p8 有可见背景
积分偏差，p24->48->96 质量/化学势在约 1e-12 fm^-1 内稳定，残差约 1e-15。
既有 FixedMuBConservedCharges 执行器在收敛候选间按残差挑选；没有据此证明
多解区域的全局自由能分支选择。该独立问题不通过修改 PNJLCore 来规避。

`method_v1_axes/` 使用保留 T170 背景；基准 mesh128、omega1200、q12/qmax8。
32 个数值轴 channel 结果和四个 all-hard 模型探针均通过条件谱检查。
相对同基准的四通道最大密度变化为：

| 改变项 | 最大相对变化 |
|---|---:|
| 谱 mesh128->256 | 0.004325% |
| omega1200->2400 | 0.044814% |
| q12->16 | 0.229468% |
| qmax8->10，q 节点数不变 | 0.588307% |
| omega 上端48->56，omega 节点数不变 | 0.018775% |
| 下端1e-4->1e-5 | 0.000002935% |
| 热 cutoff20->24 | 0.000543% |
| 热项也硬截断在真空 Lambda（模型改变） | 4.512--9.017% |

改变上端且节点数不变同时改变了节点分布，不能称为纯尾部误差；逐轴小变化也
不是全域误差上限。谱虚部的一维 energy_nodes 默认64，与直接二维泡 p/x 节点
不同：仅提高后者不会加密用于密度的 Cauchy 谱。其独立检查另行记录。
本轮旧 axes 的 manifest 绑定了源码和输出，但没有直接记录冻结背景输入哈希；
来源可由执行快照追溯，不把这一点夸称为完备输入 provenance。

`method_v1_routing/` 的 80 个 gap 探针比较双线球交集和相对动量中心球 cutoff，
q=0/0.01/0.02/0.05/0.1/0.25/0.5/1/2/3，p96/192；只算根，不算替代密度。
centered 根扫描与围道数均通过、两种 p 节点稳定，但两处方根的位置差很大。
例如 pi+ q=1 的外频率根为 2.615743346 对 1.343208900 fm^-1；q=0.01
也已为 0.766128736 对 0.726863103。q=0 的约 6.47e-4 fm^-1 差来自谱插值
与直接 gap 求积的表示差；有限 q 的大差异不能这样解释。

小 q 双线球交集有线性 |q| 接触边界项，中心球以 q^2 起始，推导和纯数值代数
极限测试见方法文件。这是影响色散的显著模型系统项，不是 phase bug，也不是
提高节点或调容差能消除的项。此前口头将 routing 差异称为小量的表述撤回。
不得自动把 centered 探针晋升新密度方法，亦不得宣称 direct finite-q 优于文献外推。

### 22.5 稀疏冻结线的通过项与保留失败

`method_v1_freezeout/` 在解析环境参数时失败，`method_v1_freezeout_v2/`
的六个背景虽已求出，但 q_nodes=4 不满足最低8节点，24个密度求值全部失败。
这两份目录原样保留；stage completed 字段不表示物理/数值接受。参数现在在求解
背景之前做有类型的解析和合法性检查，测试覆盖这些错误。

`method_v1_freezeout_v3/` 使用合同默认 mesh256、energy64、omega2400、
q16/qmax8、热 cutoff20，仓库现有 default freezeout 参数化，三个点各算 p24/48
背景。六背景残差均小于 1.2e-14；24个密度都是正的，但仅22个通过全部条件 gate。
384个 channel/q 检查中382个通过，两个失败都是 sqrt(sNN)=3 GeV 的 pi-，
q=3.61994996065 fm^-1，分别对应 p24/48 背景。不能因总密度为正而忽略它们。

| sqrt(sNN) [GeV] | T [MeV] | muB [MeV] | K+/pi+，p48 | K-/pi-，p48 | v3 条件 gate |
|---:|---:|---:|---:|---:|---|
| 200 | 165.9230564 | 23.5251799 | 0.3614574 | 0.3320646 | 四通道均通过 |
| 7.7 | 139.6121195 | 421.6498501 | 0.5441895 | 0.1410776 | 四通道均通过 |
| 3 | 79.9569438 | 719.0764156 | 0.5061299 | 不接受 | pi- 围道未通过 |

这些只是当前处方的 diagnostic 比值，未作实验比较；尤其不能凭三个点判断
horn、定量拟合或旧结果的物理优劣。尚无逐点全部谱轴收敛及全谱计数完备性证明。

失败点根扫描128/256一致，p24根 k0=5.17118355187、距阈值9.25624e-6 fm^-1；
静态检查、分支检查及当前相位容差均通过。旧均匀围道加到每边4096节点后，
两种绕数仍约为1，但 max_step=1.57852482 略高于固定 pi/2 门槛。
缩小围道高度的只读定位探针可通过，表明是近边界采样问题的证据；正式修正不
缩小矩形，而在原矩形的实轴交点聚集节点，保持同一解析域与相角门槛。
新增合成根距围道边界1e-8的测试，另作 `sparse_refine` 针对性数值复核。
v3 manifest 仍保持 conditional_gates_passed=false，不追溯改写失败记录。

`method_v1_sparse_refine/` 验证聚集围道修正：两背景下 contour-only、energy128、
omega4800 均通过；mesh512 的独立根数/围道数也均为1，但固定阈上偏移1e-7的
相位/pi 残差为 -0.00645311，超过原0.005门槛，记录 failed_gates=levinson。
进一步只读端点序列为：

| 阈上偏移 [fm^-1] | mesh256 残差 | mesh512 残差 |
|---:|---:|---:|
| 1e-7 | -0.00426636 | -0.00645311 |
| 1e-8 | -0.000437598 | -0.000665164 |
| 1e-9 | -0.0000439123 | -0.0000668137 |
| 1e-10 | -0.00000439316 | -0.00000668528 |

这表明固定偏移不等于 Levinson 的阈值极限；当前离散谱的极限仍趋于 pi。
正式壳层检查因此保留整个四点序列，要求末两点稳定、阈值非零，再比较独立根数。
频率积分也补入这三个较近端点，未折叠相位或改原容差；这不证明谱网格的极限
可与端点极限任意交换，仍保留不同 mesh 的独立验证。新增13个合成端点测试
覆盖0/pi极限、阈值零点、非有限输入、振荡不收敛及不可分辨的偏移。

`method_v1_axes_refine/` 已补齐 energy64->128、q12->24 和独立 q=[8,10]
尾部带。12个 channel/variant 条件检查均通过；energy64->128 密度差在浮点
舍入量级，q12->24 最大变化0.207096%，尾部带占原总密度最高0.013155%。
因此原 qmax8->10 的0.588307%主要不能解释为纯尾部贡献，还包含节点重分布。
这里仍只验证T170代表背景；不是所有冻结线点的绝对误差上界。

`method_v1_sparse_refine_v2/` 在原两个失败壳层上重新比较 contour-only、energy128、
mesh512、omega4800，8/8通过独立根/围道/端点相位检查。最大末点相位/pi 残差
约6.69e-6，末两点变化约6.02e-5；原phase门槛0.005及新增极限稳定要求0.0005
均未放宽。原始失败目录仍保留。此结果解决针对性诊断失败，不是全谱计数授权。

### 22.6 可复核运行入口

所有命令在 `D:\w\jrt-ord` 执行。每次 `GBU_RESEARCH_OUTPUT` 必须是尚不存在的
新目录；不能再次使用表中已留存目录。`audit_causal_gbu_validation.jl` 接受：

| GBU_RESEARCH_STAGE | 内容 | 是否新求背景 |
|---|---|---|
| backgrounds | 背景p8/24/48/96对照 | 是 |
| mott | 热Mott和粗q continuation | 是 |
| q_refine | 动量消失区间及旧失败点 | 是 |
| axes / axes_refine | 原代表点逐轴/补充求积与尾部检查 | 否 |
| freezeout | 三个能量点，各p24/48背景与密度 | 是 |
| sparse_refine | 读取旧稀疏输出，仅复核失败壳层 | 否 |
| saved_sparse | 读取p48背景，联合加密三个能量点 | 否 |

示例（output后缀由执行者选择唯一值）：

```powershell
$env:GBU_RESEARCH_STAGE = 'saved_sparse'
$env:GBU_RESEARCH_INPUT = 'D:\w\jrt-ord\data\outputs\results\relaxtime\analysis\charged_rpa_phase_backend\method_v1_freezeout_v3'
$env:GBU_RESEARCH_OUTPUT = 'D:\w\jrt-ord\data\outputs\results\relaxtime\analysis\charged_rpa_phase_backend\method_v1_saved_sparse_new'
julia --project=. scripts/analysis/relaxtime/audit_causal_gbu_validation.jl
```

同背景temp7比较需显式设置 `GBU_TEMP7_SOURCE` 指向外部原
`audit_fig2_phase_shift.jl`，执行 `audit_causal_gbu_comparison.jl`；routing极点探针
执行 `audit_causal_gbu_routing.jl`。三个入口均只生成 diagnostic，且禁止覆盖。
运行中不修改src/config/研究脚本，否则结束时源哈希门禁失败。

### 22.7 后续物理取舍的依据

等待数值运行期间复核已留存的 `2301.09882.txt` Appendix A.4 和正文：其两味
mu=0 有限 q 谱路线采用 Lambda_4^2(q)=4(Lambda^2+m^2)+q^2 的色散上界，
不是当前两个夸克动量球交集。外层 Landau 动量上界的敏感性又是另一项。
因此文献不能用来证明本项目双线交集是唯一严格选项；也不能把该文的压力计算
直接当作 finite-BQS 的 GBU charged yield。源PDF哈希与第20.1节相同。

建议后续比较真空不变质量谱截断、中心球及当前双线交集，先在同输入下验证
q=0匹配、真空小q色散、因果实虚部及热谱，再比较密度，不按实验拟合好坏选择
切支/截断。方法文件已记录这一候选；本轮没有擅自将任何处方替换为论文主路线。

### 22.8 联合加密收尾、验证与保留边界

`method_v1_saved_sparse/` 已正常结束，复用 v3 的三个 p48 背景，没有调用新的
平衡求解。联合设置为 mesh512、energy128、omega4800、q24/qmax8、
k0 下端1e-5/上端56、热 cutoff24，所有频率与动量单位均为 fm^-1。
12/12 密度通过当前条件门禁及相对变化目标，288/288 channel/q 检查通过，
失败求值数为0。与 v3 同背景原网格相比，最大密度变化为 **0.087847896%**；
这不是绝对误差上界，也不能替代逐点正则化、尾部和全谱检查。

| sqrt(sNN) [GeV] | n(pi+) [fm^-3] | n(pi-) [fm^-3] | n(K+) [fm^-3] | n(K-) [fm^-3] |
|---:|---:|---:|---:|---:|
| 200 | 0.00483144823140 | 0.00488548375657 | 0.00174358812006 | 0.00161973459659 |
| 7.7 | 0.00216377090305 | 0.00247590072581 | 0.00117682988366 | 0.000349097496445 |
| 3 | 0.000117875921441 | 0.000131432873664 | 0.0000596315986622 | 0.00000147969435578 |

对应 K+/pi+ 为 0.360883121696、0.543879151899、0.505884475245；K-/pi-
为 0.331540268537、0.140998180099、0.0112581754818，能量顺序同表。
这些是当前候选处方的 diagnostic partial-yield 比值，没有实验符合性结论。
旧 v3 的3 GeV pi- 失败仍保留；新的通过记录不回写旧 manifest。

阈值相位/pi 与独立根数的最大残差为1.95925e-5，末两偏移的最大相位/pi变化为
1.76025e-4，分别满足原0.005及极限稳定0.0005门槛。总密度正值并不意味着
每个壳层都为正：3 GeV K+ 的8个高q壳层仍保留负贡献，q约6.18--7.98，
shell 约 -1.75e-18 至 -3.49e-20 fm^-2。未裁零，也没有把这种极小尾部的符号
单独认证为可靠物理预测。离散根、unitary、Landau 原始分量均在 CSV 中保留。

最终 manifest 绑定4份输入、295份源码/配置快照、2份输出；独立 SHA-256
核验全部相符。连同本节前述10个已完成结果目录，共11份 manifest 的源码快照、
输出和方法合同哈希均通过核验。旧 comparison 的输入目录未直接写入 manifest，
本次按执行快照解析路径后另核对2份输入和外部 temp7 原脚本；旧 axes/routing
未完整记录输入哈希的限制不因这次检查而消失。核验以各次执行快照为准，
不要求早期运行与之后改过的当前源码相同。保留原错误与失败目录，没有清理数据。

本轮所选 focused unit 已通过1312/1312（含新 research contract 的118项），
旧 meson mass/density regimes regression 为81/81。中断恢复后的收尾另重跑
research contract 118/118 和这两组旧 regression 81/81，不把重复运行数累加。此前一次误写
regression selector 的命令在选文件阶段失败，修正后才得到81/81；没有改 baseline。
docs、formula-route、script governance、data-output guard、active-doc governance、
task-ledger preflight 和 diff-check 已通过；收尾改文档和 ledger 后的相关复查也通过。
没有运行完整 CI suite，远端 CI 只覆盖已提交的 HEAD，不覆盖这些 dirty changes。

结论分账如下：

- **已修正的代码/诊断问题**：环境参数类型与合法性前置检查、近边界围道分辨率、
  将固定阈上偏移误当阈值极限。修正均有合成测试，没有放宽原判据。
- **数值证据已增加但未全闭合**：三点密度联合加密稳定；热 Mott/动量事件仅有
  明确范围内的 bracket 与计数证据，网格误差不等于二分区间宽。
- **显著模型/公式选择仍未解决**：双线硬截断的线性小q边界项及大幅色散/密度
  依赖；真空不变质量谱截断候选尚未完成不等质量 contact/减法及因果配对推导。
- **独立验收缺口**：额外解析 gap、cut 内结构及不稳定模的全谱计数，背景多解
  区的全局稳定性，以及实验产额口径。不能以当前稀疏正值替代这些验证。

因此本轮完成研究合同、受控比较、限定范围的根/Mott和稀疏数值验证，但不冻结
完整论文曲线。下一步先解决22.7的正则化比较，再决定是否扩展冻结线。
两处 ledger next_action 已据此更新；历史 accepted 字段不是本次数值或
production 授权，未改历史审核状态。

2026-09-05 收尾核对：分支 `codex/charged-phase-coordinate-cut-fix`，HEAD
`fa34a3975d858ca5c517f3490476178720f0d23e`；PR310 OPEN，base仍为
`codex/charged-profile-gates`，14/14已提交 checks成功。本轮未暂存、提交、
推送、合并或改变 production 默认，所有 data/outputs 仍是未跟踪 diagnostic。

### 22.9 冻结线 direct finite-q 与 q=0 reference 受控比较（2026-09-05 恢复后）

新增 scripts/analysis/relaxtime/audit_causal_gbu_freezeout_comparison.jl，并以纯代数
契约测试 tests/unit/relaxtime/test_causal_gbu_freezeout_comparison_contract.jl
固定汇总规则。输入为同一组 FixedMuBConservedCharges quark-only BQS p48 背景；
direct 读取 method_v1_saved_sparse/，reference 读取
method_v1_freezeout_q0_reference_v3/。每个 bound/unitary/Landau 分量先乘

$$
q^2/(2\pi^2),
$$

再使用与原运行相同的 24 点 Gauss--Legendre q 权重积分；没有把三个分量重复使用
shell_inv_fm2。v1 汇总目录保留为脚本错误的审计痕迹，v2 才是修正后的结果。

修正后的输出位于
data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/method_v1_freezeout_comparison_v2/，
包含 component_totals.csv、ratio_comparison.csv、historical_context.csv、
freezeout_route_audit.png 与输入/输出 SHA manifest。12 个 direct 与 12 个 reference
总密度均逐行闭合到原始 density CSV；两条路线都标为 production_authorized=false。

| sqrt(sNN) [GeV] | direct K+/pi+ | q=0 reference K+/pi+ | direct K-/pi- | q=0 reference K-/pi- | reference/direct pion density |
|---:|---:|---:|---:|---:|---:|
| 200 | 0.3608831 | 0.2989259 | 0.3315403 | 0.2740314 | 11.62--11.63 |
| 7.7 | 0.5438792 | 0.4154218 | 0.1409982 | 0.1056556 | 13.77--13.93 |
| 3 | 0.5058845 | 0.2731274 | 0.0112582 | 0.0060339 | 33.93--34.09 |

因此此前“约 9--11 倍”只适用于 T=170、muB=240 的代表点；沿这三个冻结线点，
低能 pion 密度因子明显更大。ratio 变小不是由单一 K 通道异常造成，而是 reference
路线同时抬高 pion 与 K 密度，且 pion 抬升更强。

分量账显示两条路线的总密度主要由阈下离散根贡献，unitary 连续谱给出负补偿，
Landau 项在这些加密点通常远小于前两者。例如 direct K+ 的 bound/unitary/landau
分量约为：200 GeV (0.00278673, -0.00104321, 6.93e-8) fm^-3；7.7 GeV
(0.00161549, -0.00043866, 9.10e-11) fm^-3；3 GeV
(6.29316e-5, -3.30002e-6, 5.81e-19) fm^-3。负的 unitary 分量是当前
d\delta/d\omega diagnostic 的连续谱贡献，不是脚本裁零或缺少 Bose support；
它与离散 bound 项相消后仍得到正总密度。故“总密度为正”不能被解释为每个谱分量
都正，也不能把负 unitary 直接称为代码错误。

历史 trho_asymmetric_kplus_piplus_scan_v1 同步生成了按温度的上下包络，但明确标记
为 historical_context_unmatched_background：它使用 FixedAsymmetricRho、
rho_u/rho_d=0.876、x_min_cut 与 bose_x_min=0.01，而当前冻结线使用
FixedMuBConservedCharges、rho_Q/rho_B=0.4、rho_S=0，且不允许裁零或
扩展 Bose 域。历史 GBU K+/pi+ 在 T=200 MeV 的 rho 扫描范围约为 0.35--1.63，
当前三点 ratio 不能与其作同点定量匹配；图中只作分面上下文，未叠加成同一路线。

现有 method_v1_routing/routing_gap_comparison.csv 仍只证明正则化处方对 gap 极点
色散的敏感性，manifest 明确 density_computed=false。代表点 q=1 fm^-1 的双线球
交集与相对动量中心球根相差约 1.1--1.35 fm^-1，故 finite-q/q=0 密度差异不能
仅归因于相位 unwrap；但在完成同一 regulator 下的实部、虚部、contact/subtraction
和因果性闭合前，也不能把任一处方授权为文献级 production 路线。

本节新增契约测试通过 7/7。本次修正只改变 diagnostic 汇总脚本和文档，未修改
PNJLCore、旧 MesonDensity 默认、production baseline、PR 合并状态或任何
data/outputs 正式产物。

### 22.10 FIG.4-like ratio 图（2026-09-05）

根据 Blaschke 等 2020 年论文 FIG.4 的展示结构，新增
scripts/analysis/relaxtime/plot_causal_gbu_fig4_like.jl。输出目录为
data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/fig4_like_freezeout_ratio_20260905/，
包含当前 BQS 面板、历史路线上下文面板、合并图、README 和 manifest。

当前面板的横轴是 sqrt(s_NN)_GeV，但目前只保留 200、7.7、3 GeV 三个
加密冻结线点；实线/圆点为 direct finite-q，虚线/菱形为 q=0 lambda-reference，
黑色表示 K+/pi+，红色表示 K-/pi-。这是 FIG.4-like 诊断展示，不含论文实验点，
也不声称已经得到完整冻结线或 horn 拟合。

历史面板不能使用相同横轴：保留的
trho_asymmetric_kplus_piplus_scan_v1 是 FixedAsymmetricRho 的 T-rho 扫描。
因此该面板以 T 为横轴，显示 rho 扫描的 K+/pi+ 中位数和包络，并显式标注
rho_u/rho_d=0.876、x_min_cut、bose_x_min=0.01。它是历史处方的温度上下文，
不是沿当前化学冻结线的 sqrt(s_NN) 曲线；没有把其点插值或伪匹配到当前 BQS 结果。

论文 FIG.4 本身比较的是不同扫描线（包括临界线/常温线）与实验比值；因此当前
图只借用坐标和图例语义，不把论文的扫描线、实验点或 anomalous-mode 阴影
移植到当前 strict GBU 结果。两张图均 production_authorized=false。

### 22.11 十点冻结线 diagnostic 与执行复核（2026-09-05）

作者明确要求先补充当前 FIG.4-like 点、再进行 regulator/routing 闭合检查。
本轮据此扩展的是可观察趋势的 diagnostic，不是绕过第22.8节的论文/production
冻结门禁。背景、通道、GBU 权重与模型默认均未改变，没有修改 PNJLCore 或
旧 MesonDensity；真空 cutoff 仍取模型 Lambda，BQS 仍为 quark-only。

有效结果目录为
`data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/fig4_like_freezeout_ratio_dense_20260905_v3/`。
其中 `plots/current_bqs_fig4_like.png` 为1200×700加密图，`ratio_comparison.csv`
为比值，`audit/` 保存分量账及独立哈希/行数审核。原三点图未覆盖，其 SHA-256
仍为 `fcf3b1aba660389f3a4d4b0a4766473c72c680e8701ac5370e928fe531a93937`。
早期无后缀空目录、串行 `_v2` 的中断检查点也保留，不作为十点完成结果。

计算入口 `audit_causal_gbu_dense_parallel.jl` 只读复用已验收三个点的哈希绑定
背景/密度，另求七个 p48 背景。源代码和设置匹配先于复用；新背景沿能量下降
续算，残差最大3.207e-14。设置保持 mesh512、energy128、omega4800、q24/qmax8、
k0=1e-5..56、thermal cutoff24（内部频率/动量单位fm^-1）。新加点不是插值。

串行先实测约15分钟完成6个密度求值，随后改用四通道并行、单写入者检查点。
`causal_gbu_dense_execution.jl` 的 reference 适配器只复用 q=0 谱，省去原算法
随后丢弃的有限q谱和围道；13项合成逐字段比较与原慢算法完全一致。
130 GeV 的新并行 direct pi+ 也与串行检查点逐值相同；单个新增能量点实测约
200--292秒。未减少积分节点、改容差、使用fold或裁零。

| sqrt(sNN) [GeV] | direct K+/pi+ | q0 reference K+/pi+ | direct K-/pi- | q0 reference K-/pi- |
|---:|---:|---:|---:|---:|
| 3 | 0.5058845 | 0.2731274 | 0.0112582 | 0.0060339 |
| 5 | 0.5364496 | 0.3726560 | 0.0821954 | 0.0561975 |
| 7.7 | 0.5438792 | 0.4154218 | 0.1409982 | 0.1056556 |
| 11.5 | 0.5224628 | 0.4165374 | 0.1846516 | 0.1444813 |
| 19.6 | 0.4737229 | 0.3874937 | 0.2328538 | 0.1875820 |
| 27 | 0.4450840 | 0.3665644 | 0.2570979 | 0.2090965 |
| 39 | 0.4177728 | 0.3454075 | 0.2801826 | 0.2294399 |
| 62.4 | 0.3922877 | 0.3249242 | 0.3024475 | 0.2489020 |
| 130 | 0.3687292 | 0.3054899 | 0.3240666 | 0.2676136 |
| 200 | 0.3608831 | 0.2989259 | 0.3315403 | 0.2740314 |

在本次采样点中，direct K+/pi+ 以7.7 GeV最高，reference以11.5 GeV最高；
后者与7.7 GeV差很小，不能从十点折线认证精确峰位或峰位偏移。K-/pi-在这十点
随能量上升而增加。没有实验点、实验拟合、完整物理谱认证或论文曲线冻结。

`summarize_causal_gbu_dense_outputs.jl` 独立核对302份运行源码/配置快照、7份输入
哈希、80个密度、1920个壳层；全部当前条件gate通过。bound/unitary/Landau分量
与密度的最大闭合差2.082e-17 fm^-3；8个微小负尾壳层全部保留。
最大Levinson相位/pi残差1.960e-5，阈值末两偏移变化1.761e-4；这不是逐点全轴
误差上界，也不是regulator独立性。图中折线仅引导视线；失败/缺失点在CSV保留，
绘图必须断线，有限但失败的值另画橙色叉号。

PNG生成、尺寸、输入/输出哈希及数据合同已核对；当前图像接口向模型返回
“不支持图像输入”，因此没有宣称已完成主观版式视觉验收。

### 22.12 同背景 regulator/routing 闭合与数值缺陷定位（2026-09-05）

本节属于 research / required_follow_up，不变更 primary track、上游BQS、
v1物理处方或production授权。新增 `causal_gbu_regulator_checks.jl` 和
`audit_causal_gbu_regulator_closure.jl`，在保留的 T=170 MeV、muB=240 MeV
quark-only BQS背景上比较双线球交集（two_line）与相对动量中心球（centered）。
两者使用相同质量、化学势、Polyakov、耦合、真空Lambda和热cutoff24 fm^-1；
不是把有限q强制拟合到q=0外推，也不使用实验数据挑选处方。

**数学推导与实现**：设两线动量为 p1=P+q/2、p2=P-q/2，则
P^2=(p1^2+p2^2)/2-q^2/4。角向delta满足 E2=t(sE1-lambda)，中心球限制
P^2<L^2因此化为E1上的二次不等式。辅助代码以它截取运动学cut区间；另在
中心坐标直接二维求积得到B0、Pi、两个多项式矩及同域contact。没有删除Landau
结构。两种有限域的能量极值给出几何谱支撑；解析gap不是通过“小虚部”猜测。

首轮输出 `regulator_routing_closure_20260905/` 含64条比较、320个复数探针、
146个gap记录和48个弱eta记录。独立审核核验302份源码快照、4份输出、2份背景
输入哈希，无求值异常；但**联合门禁仅4/64通过**。不能把无异常或gap计数大多
通过误报为整体通过。预设目标始终为谱重建差<1e-4 fm^-2、独立二维节点变化
<1e-6 fm^-2、contact恒等式残差<1e-10 fm^-2；未调整容差。

多数首轮失败来自二维节点128/64到256/128未达到1e-6目标，而不是发现因果关系
被破坏。q=0的谱mesh128/256/512重建差另呈约四倍下降；这两个离散轴分开记录。
此外 `K_plus centered q=3 mesh256` 的完整gap稳定性失败，不能静默忽略。

**已定位的代码缺陷**：`audit_causal_gbu_closure_failure.jl` 对同一分段线性谱做
256-bit独立解析积分，发现相距约8.88e-16 fm^-1的节点有约4.5e7的大斜率。
原实轴表达式的斜率乘距离与对数相消，造成根附近逆传播子约1.1--1.8e-8的噪声。
64/128/256/512根扫描中128节点未获得认证根，其余根的位置也有不稳定偏差。
这是数值求值缺陷，不是额外物理束缚态或改变相位规范的证据。

修正 `CausalSpectralBubble.cauchy_transform` 的单元求值形式。令h=b-a、
t=h/(a-z)、dy=y(b)-y(a)，解析积分写为

```math
y(a)\log(1+t)+dy\,H(t),\qquad H(t)=1-\frac{\log(1+t)}t.
```

小|t|使用H=t/2-t^2/3+...的机器精度展开；实轴仍使用原PV减法与rho虚部。
这只是代数等价的稳定计算，不合并谱节点、不改变cut或丢掉边界项。
新合成窄单元/远场24项测试在修正前12项失败，修正后全部通过；另有三角谱、
正负谱、实轴/上半平面和独立loop测试。真实mesh256根探针修正后与256-bit
结果最大差3.331e-15，四种扫描所得根约3.78107820499 fm^-1，最大跨度
4.605e-11 fm^-1。旧失败目录保留，新证据为
`regulator_closure_failure_probe_fixed_20260905/`。
另在 `regulator_closure_failure_legacy_stable_20260905/` 保留旧/新表达式的
全部8条扫描状态：旧128节点明确为 `uncertified_zero`，未省略失败行；新四种
分辨率均认证1根。该对照特意保持旧扫描域，避免把端点margin修正混入算术比较。

同次修正还统一了support-gap实根扫描与围道的端点margin：原辅助代码先缩gap，
再被扫描器默认margin缩一次，实际排除2e-6而记录1e-6。现只施加一次1e-6，
并用距端点1.5e-6的合成根防止重犯。没有缩小物理失败门槛；排除域与元数据一致。

**加密诊断证据**：`regulator_routing_closure_refined_20260905/` 对四通道、
q=0/1/3/6 fm^-1、两处方，使用谱mesh512、cut能量128，二维节点比较
256/128到512/256。32/32联合条件检查通过；160个复数探针、70个解析gap、
48个弱eta记录均完整。独立 `summarize_causal_gbu_regulator_outputs.jl` 重新
计算gate和计数，核验305份源码快照、4份输出、2份背景哈希。

- 最大谱重建差2.181e-5 fm^-2；最大二维节点变化6.086e-7 fm^-2。
- 最大contact恒等式残差5.879e-12 fm^-2。
- 同一mesh的q=0四通道根在两处方间逐值相同。
- 若错误复用C(0)，q=1的逆传播子偏移：双线pi约-0.217、K约-0.190；
  中心球pi约-0.00849、K约-0.00658。实际计算均使用各自C(q)，此列是反事实
  一致性探针，不是说当前实现正在错误复用C(0)。

| 通道 | q=1 two_line根 [fm^-1] | q=1 centered根 [fm^-1] |
|---|---:|---:|
| pi+ | 2.615734887 | 1.343213048 |
| pi- | 2.527476128 | 1.252935573 |
| K+ | 3.234165646 | 2.400628067 |
| K- | 3.918437428 | 3.083012101 |

q=6时two_line四通道正常gap无根，centered仍各有一个。不是据此定位精确Mott，
也不能只凭更接近自由色散就选择centered。加密后根差仍约0.83--1.27 fm^-1
（q=1），远大于当前离散误差；这是有限q处方依赖的证据，不能全部归于unwrap、
遗漏contact、二维节点不足或q=0外推近似。后两路线的密度差仍需分账研究。

**弱极限与范围**：在q=1，用Poisson权重的半群恒等式检查整个实轴bubble谱的
弱eta极限；eta=0.2到0.00625 fm^-1，8组误差均单调下降，末次减半误差比为
1.997--2.001。末点绝对弱误差约0.00167--0.00375 fm^-2，不能说已经为零。
这不是GBU相位的Bose加权密度eta收敛；没有要求PV与有限eta在阈值逐点相等。

所有计数只覆盖正外频率[0,64] fm^-1内、margin=1e-6的几何解析gap；
不认证cut内嵌零点、端点margin内零点、UHP不稳定零点或无限窗口外全谱。
因此 `full_spectrum_certified=false`、`density_computed=false`：本节没有计算
centered替代总密度，没有以32个条件通过点授权全q积分或完整生产闭合。

下一步是将同样的support/根/端点门禁扩展到实际密度q壳层，并继续比较文献支持
的真空不变质量谱截断候选；完成不等质量的同谱实虚部/contact配对后，才比较
centered/谱截断与v1总密度。截断选择属于待解决的模型系统项，不用实验拟合
或正密度判据替代这个选择。精确Mott、全谱/稳定性以及完整冻结线收敛仍未认证。

### 22.13 Cauchy修正的冻结线漂移回归与交付边界（2026-09-05）

`audit_causal_gbu_cauchy_density_drift.jl` 使用十点v3中保留的3、7.7、200 GeV
p48背景，完全不调用BQS求解器，保持mesh512、ne128、nw4800、nq24/qmax8、
k0=1e-5..56、thermal24。对两条路线、四通道重算24个完整密度和576个壳层，
输出至 `cauchy_stability_density_drift_20260905/`。

24/24密度和576/576壳层通过原条件gate；8个微小负尾壳层继续保留。相对原v3的
最大密度变化6.658e-10（3 GeV pi- reference），direct最大约2.052e-13。
这些代表点不支持将大幅finite-q/reference密度差归因于本次消减误差。
未在修正后重新跑全部十点；十点图保留原v3源码快照及数值，三个代表点的全积分
漂移回归支持其图示趋势不因本次修正改变，不把这等同于其他点的逐点回归证明。

漂移回归独立核验305份源码快照、4份输入与输出哈希，全部相符。manifest SHA-256：
`e488fb4e1a0f836adb14e5f21c271fdcdfe42b30b7ba313f1965e7f6e3805d60`。
加密regulator运行manifest SHA-256：
`9d16eb983cc4f3545567e949076edffc2f56e5819169b4e853c5997310cbee10`。
十点v3主manifest SHA-256：
`da2caecd08342ae2c1ca1b2f5182d5b298a3052d08d97f1b8dee0bafe2944c6e`；
十点当前图SHA-256：
`1b7e62cac22c80dc8a688716aedd7ff59a3eec86aab40a2b0613b4bc5fc76794`。
原三点图哈希仍不变，所有原失败/中断产物保留，没有覆盖正式数据。
另由独立分量汇总重新q积分，24个bound/unitary/Landau分量账与总密度最大差
1.388e-17 fm^-3。

本次最终所选9组focused unit为6716/6716，旧meson mass/density regimes
regression为81/81；此前更大选择集合的7285/7285与这次不重复累加。
窄单元测试首次存在Julia数字小数点语法错误，修正语法后才观测到12/24数值失败；
稳定表达式修正后24/24通过。refinement轻量入口首次docstring绑定include失败，
没有开始数值运行；实际成功运行通过显式include closure模块后调用
`main(configurations=[(q,512) for q in (0.,1.,3.,6.)],np=512,nx=256)`。
后续轻量入口已改为有main和PROGRAM_FILE保护的模块，单独验证可include且不启动扫描。

本轮交付分账：

- 代码缺陷：实轴Cauchy消减和辅助gap重复margin已修正，并有纯合成测试。
- 数值未收敛：首轮独立loop节点不足，经加密在本节有限探针集合过原门槛；
  不宣称全积分域均已通过。
- 模型系统项：两个自配套regulator在有限q仍有明显色散/计数差；尚未选择论文
  优选处方，也没有替代centered密度曲线。
- 公式与验收缺口：全谱、稳定性、端点根、精确Mott、全部实际q壳层及文献支持的
  谱截断推广仍未闭合；bubble弱eta证据不能冒充GBU密度弱极限。

只在D:\w\jrt-ord工作。PR310在线核对仍OPEN，base为codex/charged-profile-gates，
HEAD为fa34a3975d858ca5c517f3490476178720f0d23e，远端14/14 checks成功只覆盖
已提交HEAD，不覆盖dirty changes。本轮没有暂存、提交、推送、合并或production晋升。

### 22.14 冻结线路径叠加历史T-rho热图（2026-09-05）

作者要求先快速叠加背景轨迹，暂不继续新的regulator密度计算。新增
`compute_freezeout_trho_overlay.jl` 与 `plot_freezeout_on_historical_trho.py`，
复用十点v3保存的p48状态与三味化学势，经Models.model_rho补算净夸克密度，
映射rho_B/rho0=(rho_u+rho_d+rho_s)/(3rho0)，rho0=0.16 fm^-3。
不调用平衡求解器，不调用介子质量/相移/密度积分，不改变保存的背景分支。

坐标采用原p48求积，与p96/p192对照，并以p192的quark-antiquark分布积分独立
核验压力导数。十点全部通过：p48/p192坐标最大相对变化2.750e-15，分布对照
最大绝对差8.154e-17 fm^-3，原BQS约束残差最大1.641e-15 fm^-3。
14项纯代数/域覆盖测试通过；数值漂移仅检查坐标，不声称旧介子结果已同点匹配。

| sqrt(sNN) [GeV] | T [MeV] | rho_B/rho0 | 历史网格覆盖 |
|---:|---:|---:|---|
| 3 | 79.9569 | 0.00290349 | 温度、密度均在范围外 |
| 5 | 118.5236 | 0.01451084 | 温度、密度均在范围外 |
| 7.7 | 139.6121 | 0.03189367 | 密度在范围外 |
| 11.5 | 151.5934 | 0.04692547 | 密度在范围外 |
| 19.6 | 160.0084 | 0.05298796 | 范围内 |
| 27 | 162.5747 | 0.04877581 | 密度在范围外 |
| 39 | 164.2385 | 0.04017734 | 密度在范围外 |
| 62.4 | 165.2674 | 0.02839003 | 密度在范围外 |
| 130 | 165.8213 | 0.01476026 | 密度在范围外 |
| 200 | 165.9231 | 0.00977933 | 密度在范围外 |

历史网格为T=120..220 MeV、rho_B/rho0=0.05..1；九个点落在低密度边界之外。
“理论上可以映射轨迹”不代表历史网格实际覆盖整条轨迹，更不代表可以用热图
插值提取完整旧方法冻结线ratio。原u/d=0.876，当前严格BQS对应0.875，图中均
明确标记；不存在把两者当作完全相同背景的定量认证。

输出目录：
`data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend/freezeout_on_historical_trho_20260905/`。
`freezeout_trho_coordinates.csv`保存坐标、净密度、检查结果与域覆盖标记。
`plots/historical_gbu_with_freezeout.png`为GBU全图与低密度放大图；
`plots/historical_all_regimes_with_freezeout.png`为历史四面板叠加图，另保存GBU SVG。
底色仅为旧ratio，线仅为当前背景路径；没有叠加新ratio或插值旧ratio。橙色菱形
表示历史范围外点，外部浅色表示无历史数据；只裁切最外半格的绘图足迹到实际
采样矩形，防止暗示范围外仍有数据。旧log色标0.310601..260.526881保持不变。

291份源码/配置快照与脚本、输入、输出哈希保留；原历史CSV、PNG、SVG和manifest
在绘图前后哈希一致。主GBU图2610×1368，四面板2430×1980；图像接口仍未提供
模型可见像素，因此只报告生成/尺寸/数据与哈希验收，不宣称主观视觉验收完成。
首次动态include遇到Julia world-age错误，在数值计算前修正为invokelatest入口；
最终实际执行成功。所有新产物为diagnostic，不覆盖历史正式产物，不提交data/outputs。

### 22.15 同域独立坐标、完整有符号实gap与regulator分账（2026-09-06）

**作者范围与执行决定**：先完成同背景同regulator的q=0匹配、实虚部/contact/
因果表示、全部解析gap根计数，再判断finite-q差异。停止手动冻结线加密；现有
ratio足作诊断。待模型处方评审后，将检查与失败保留接入研究生产脚本，再统一
生产完整曲线。本节为当前primary track的required_follow_up，不改变旧默认、
PNJLCore、BQS背景、baseline或PR310探索性状态。

**新增可重复入口**：`causal_gbu_coordinate_oracles.jl` 在固定参数下改用
(p1,p2)积分，q>0测度为p1*p2/(q*E1*E2) dp1 dp2；同一two_line或centered
域在解析边界上划分。原loop坐标不复用该划分。真空/热谱、B0、contact和两个
多项式矩分别累计；关闭热占据保留当前质量，不调用新的真空平衡求解器。

`causal_gbu_gap_completeness.jl` 不读unwrap，检查每个紧支撑线性谱的全部正、负、
内部及外侧实gap。符号变化单元精确写为两个端点三角形，不裁掉微小权重，也不
插入浮点数无法表示的零交叉节点。符号定号单元的Cauchy值在实gap内单调，用
端点范围之和排除无根区；不能排除时同时用两种分辨率实根扫描和辩值原理计数。
旧margin遗漏的端点条带单独检查，不通过则缩小条带以解析根，不能忽略该区。
设支撑半径U、M>=integral|rho|/pi，取R=U+max(1,8|K|M)，则外侧
|4K Pi|<=4|K|M/(R-U)<=1/2，排除无穷尾的逆传播子零点。

这些范围是显式带舍入余量的Float64数值界，不是对原连续核的区间算术证明。
计数范围是**每个插值谱的全部实gap**；cut内嵌零点、UHP不稳定模、精确Mott端点、
连续核完备性与物理极点留数仍需分别验收。两种谱网格计数相同和根漂移需同时
报告，不能把计数稳定等同于根位置精确。正负频率两根不是同一电荷密度加两次。

`audit_causal_gbu_full_closure.jl` 自动执行上述检查，保留全部失败并由单写入者
输出检查点；最多4个Julia线程只并行独立case。串行pilot的每case成本约5--17秒，
现有Gauss节点缓存使用锁。默认q列表复用24个[0,8]密度GL节点，加0、0.01、0.02、
1、3、6，共30个q；不求新背景或密度。独立入口
`summarize_causal_gbu_full_closure.jl` 验证快照/输入/输出哈希、case覆盖及各原始
探针到门禁的归约，任何数据缺失或伪通过均报错。

**保留失败与实现修正**：

- `full_signed_regulator_pilot_20260906/` 首轮q=0/1共16个case，15个完成；
  pi- two_line q1在新增审计器中因“无法表示的谱零交叉节点”失败。改为端点
  三角形代数分解后解决，未修改物理谱；原失败保留。
- `full_signed_regulator_probe_v2_20260906/` 完成32个case，其中5个支撑门禁失败：
  pi± centered q3/q6及pi- two_line q1。根扫描/围道本身通过，但几何支撑不匹配。
  centered热pair外端点残值约2.8e-17，two_line Landau端点约1e-28；插值使其
  延伸到gap，mesh512最大越界宽度9.412e-6 fm^-1。该问题不能用幅值裁零或
  放宽几何门禁掩盖。
- 修正 `B0_spectral_cut`：先用同一两球交集的pair/Landau各符号精确极值排除
  支撑外及零测集端点，再求二次式允许域。`build_bubble_dispersion` 补齐这些
  极值节点；centered诊断核用其自身域的极值做同样排除。支撑内部的极小热权重
  保留；不改B0/B0_retarded/B0_pv_cut旧入口，不翻相位、不删Landau。
- 四q修正前后448个同坐标复频探针的Cauchy Pi最大变化8.6465e-11 fm^-2，
  说明端点伪支撑虽破坏严格gap语义，但该修正不解释原数量级密度差。本轮未
  重算完整密度，不能把这一bubble探针对照冒充ratio全积分漂移回归。
- 新的几何检查同时限制谱超出和遗漏支撑，不能把数值缺失伪装成物理gap。
  纯合成测试覆盖1e-200谱、近端点根、偶重根、近重合节点、域外q和假通过报告。

**四q回归证据**：`full_signed_regulator_probe_v3_20260906/`，同一冻结
T=170 MeV、muB=240 MeV quark-only BQS，四通道×两处方×四q=32/32通过。
直接loop np512/nx256，独立径向np256/512对照；谱mesh256/512、cut ne128。

| 检查量 | 观测最大误差/差异 | 原门槛 |
|---|---:|---:|
| q0两regulator的Pi | 0 | 1e-12 |
| q0 contact与独立A1+A2 [fm^-2] | 3.197e-14 | 1e-10 |
| q0旧B0与同一全硬截断B0 | 3.169e-9 | 1e-8 |
| 同regulator两loop坐标Pi [fm^-2] | 6.229e-12 | 1e-6 |
| 独立径向节点变化 [fm^-2] | 1.877e-11 | 1e-6 |
| contact代数残差 [fm^-2] | 5.879e-12 | 1e-10 |
| mesh512 Cauchy与直接loop [fm^-2] | 1.045e-5 | 1e-4 |
| mesh512实轴虚部插值误差 [fm^-2] | 3.574e-6 | 1e-4 |
| flavor反射虚部误差 [fm^-2] | 6.661e-16 | 1e-10 |

q0旧B0对照限定同一全硬截断域；不能把它与thermal24相比较后要求相等。
有限q不再借用q0的A。独立同域loop与Cauchy一致，连同同域矩恒等式，支持
当前Pi实/虚部及contact配套的因果表示；不等同于已证明完整RPA传播子无UHP极点。
mesh512支撑越界宽度为0，支撑缺失最大8.882e-16 fm^-1，属于已声明坐标舍入。
64份插值谱共272个gap逐项通过；最大根网格漂移4.896e-4 fm^-1，不能写成零。
308份源码/配置快照、2个背景输入、6个输出哈希核验相符；主manifest：
`9c2fd5759ea01d4b181bb4480689c79f88acbf94761d4d19ff829573d016a293`。

**模型分账证据**：取内部lambda的z0=2+0.6i fm^-1，比较关闭热占据后
Pi_vac(sqrt(z0^2+q^2),q)-Pi_vac(z0,0)。q=1时：

| 质量对 | two_line差 [fm^-2] | centered差 [fm^-2] |
|---|---:|---:|
| pion | -0.238583 - 0.014222i | -0.009267 - 0.000502i |
| kaon | -0.212057 - 0.010467i | -0.007335 - 0.000338i |

真空contact的处方差分别约5.55335、5.12348 fm^-2；热contact差约1e-11 fm^-2。
因此大幅处方差在没有热occupation时也存在，不能全部称为物理介质效应。
这与双线球交集的|q|边界项相符；centered也不精确保持真空boost不变，不能
凭根更接近自由色散就选定它。该探针不是对总密度差的百分比归因，也不否认
热介质的真实有限q依赖。q1正k0根的处方差仍为0.8335--1.2745 fm^-1，远大于
当前网格漂移；q6 two_line无实gap根而centered各电荷方向仍有根，不能据此
声称已精确定位Mott。

**实际密度q节点扩展**：`full_signed_regulator_shells_20260906/` 完整运行及
独立报告均通过240/240，无evaluation失败。仍是同一T170/muB240背景，不是
240个冻结线点。480份插值谱、1960个有符号实gap、3360个复频探针和8280个
实轴cut探针均纳入覆盖检查；308份源码/配置快照、2个背景输入、6个输出哈希
一致。主manifest SHA-256：
`ec3c17d7960abd32893f25d617d4eccf675e16eee9848251e3c8d919462d890d`。

全部q上的最大同域坐标差6.694e-12、径向节点变化7.398e-11、contact恒等式残差
6.380e-12，单位均fm^-2；fine Cauchy/loop差仍为1.045e-5 fm^-2。实轴cut
插值误差最大3.701e-6、flavor反射误差最大2.798e-14 fm^-2。几何支撑越界为0，
缺失宽度最大2.665e-15 fm^-1，小于显式舍入余量；根网格漂移最大4.896e-4 fm^-1。
全部使用既定门槛，未通过放宽容差或删谱修复。

fine谱中centered的120个case均有一个正k0根和一个负k0根；two_line的62个case
各有正负一根，另58个case无实gap根。没有在其他实gap发现额外根。两种谱网格
计数一致；这仍不认证cut内、UHP或连续核全谱。pion两电荷的two_line根消失
被现有q节点夹在[3,3.235525] fm^-1，kaon夹在[3.235525,3.743772] fm^-1；
只是固定T下的节点区间，不是精确Mott位置或热Mott温度。centered在检查上界
q=7.980749 fm^-1仍有根，并不据此判定其物理优越。

**验证状态**：最终所选四组谱/完整gap/loop单元测试7148/7148，另四组相位、
provider及研究合同测试344/344；台账单元测试185/185；旧meson mass/density
regimes回归81/81。台账测试初次因硬编码旧SHA而184/185，依据实际HEAD同步
预期值后通过，未改变台账状态接受规则。
旧生产默认及模型接口未变；没有运行完整全仓单元测试或新的冻结线密度回归。
初始新测试误把T=0.8615的热尾预估成<1e-15，实际为5.01e-12；将合成测试温度
明确取0.4 fm^-1后验证极小热谱仍保留，没有改数值接受阈值或真实背景。
docs、active-doc、script governance、data-output guard、formula-route、task-ledger
preflight及git diff --check均通过。PR310在线核对仍OPEN，base仍为
codex/charged-profile-gates，已提交HEAD fa34a397的14/14 checks成功不覆盖本轮
dirty changes；没有暂存、提交、推送、合并或production晋升。

**生产衔接边界**：先由作者评审正则化定义/候选同谱减法及真空性质，再把这些
可调用检查接入统一研究生产入口。该入口须保留背景→谱→全gap/稳定性→
Levinson/Mott→积分收敛→ratio的逐点状态，失败不能转成正常曲线点。此处不
创建替代生产默认，不重跑冻结线，不授权完整生产、baseline更新或PR合并。

### 22.16 Pereira 同域公式核对与 PNJL/热项延伸分账（2026-09-06）

**作者授权与范围**：以 Pereira 2024 为直接公式依据，核对双线积分域、配套
单线项和 retarded，再单独评估 PNJL 占据与热项延伸。不以删掉第二条夸克线
cutoff 条件为修复目标。本轮为 required_follow_up；保留原有 dirty changes，
不修改 src、PNJLCore、旧默认或 baseline。受控诊断不重求平衡态，不计算
新密度/冻结线；例行旧 mass/density 回归仍按既有测试流程调用其原求解器。

**文献与公式判定**：本地57页 PRC 109,025206 的 Eq.(26)--(28)/(60) 明确
区分平衡态单球 A 与双线域继承的 mathcal A(q)。Appendix C18/C19 明确
外频率 k0+i0 相对于 mass-i0 的虚部符号修改。双球交集、q>=2Lambda 时
全硬截断泡消失均属于该论文处方，不应单凭这些现象判为实现错误。
A/B 分解仍然适用，但两种项必须同域；有限q不能借用平衡态常数。
同域高能矩保证 B0=C(q)/z^2+O(z^-3)，Pi 的常数项消去，没有额外任意减法。
完整信源页码、数学推导、项目推广及边界分别写入
[Pereira公式审查](../../reference/formula/relaxtime/Pereira2024_RegularizationAudit.md)。

**独立实现**：

- `causal_gbu_pereira_reference.jl`：柱坐标透镜求积 pair/scattering 和两个
  单线项，再按 Eq.(26) 重构 P/S。透镜平移只是同域换坐标，不是 centered
  单球 regulator；不复用当前角向域、四留数、cut 区间或 Cauchy 实现。
- 真空单线项另以径向原函数闭式和解析透镜体积核验；实轴用 (E,epsilon)
  允许域及 C18/C19 分项核验 pair/Landau。q=0 使用独立 Jacobian，
  等质量零频零动量退化点保留未认证标记。
- `audit_causal_gbu_pereira.jl`：非覆盖、哈希绑定的受控入口。普通 Fermi
  全硬截断只作为原始统计/域公式对照，质量与mu仍为冻结BQS输入，不称为
  重新求出的 NJL 平衡。PNJL 占据替换与热尾8/12/16/24分别记录。

**保留失败**：

1. 新 oracle 初次345项测试338通过、7失败；原因是 `-4f2(L)` 被 Julia
   解析成 Float32 字面量与隐式乘法，而非 -4*f2(L)。修正新增参考核语法后
   345/345通过，未改容差，也未改原物理核。runner 初次加载的 `10.+1im`
   语法错误同样在新代码中修复；后续加入 runner 加载与合成适配器测试。
   新增高能矩测试又捕获两项失败：测试表达式 `-2E2` 被作为科学计数法，
   不是 -2*E2；补显式乘号后复核。两次均是新增oracle/测试代码的问题，
   不是当前泡的contact残差，也没有据此改动受控数值输出。
2. `pereira_formula_pilot_20260906/`，q=0/1.4、48案例，46通过。两项失败为
   pi±、q0、thermal24 的192/96节点变化1.165294e-6，超过原1e-6门槛；
   两实现已经接近，而粗网格未收敛。保留失败manifest并以非零状态退出。
3. 新目录 `pereira_formula_checks_20260906/` 使用384/192节点，144/144通过。
   同一门槛下节点差降到约1e-10量级。此次显式调用 nodes=384；验收后把
   runner 默认从192改为384，仅改变默认求积分辨率。运行时源码以目录内
   source_snapshot 为准，后续台账/默认值编辑不回写已绑定快照。

**受控数值证据**：均复用 T=170 MeV、muB=240 MeV quark-only BQS；
Phi=0.278902886、PhiBar=0.289038273、Lambda=3.052294211 fm^-1。
四通道×q={0,0.8,1.4,2.7,3.2,6.2}×六变体=144个case；720个复频探针、
1152个实轴cut探针、600项推广差异记录和7项弱极限记录。
没有evaluation失败或伪装成通过的缺失行。

| 验证对象 | 观测最大绝对差 | 原门槛 |
|---|---:|---:|
| 独立参考与当前 Pi [fm^-2] | 4.323824e-12 | 1e-6 |
| 独立参考与当前 B0 [无量纲] | 2.496785e-12 | 1e-6 |
| 配套 C(q) [fm^-2] | 8.846257e-13 | 1e-6 |
| 真空 C(q) 与独立解析式 [fm^-2] | 4.511946e-13 | 1e-6 |
| pair/Landau Im B0 分项 [无量纲] | 2.060574e-13 | 1e-8 |

节点门禁分别比较 Pi/C 的 fm^-2 数值和无量纲 B0 数值，并非将不同单位相加。
独立核与当前核的逐量最大节点变化数值分别为1.366522e-10、2.418909e-11，
均低于对应1e-6门槛。结果不是对所有参数/近实轴频率的统一误差上界。

**PNJL/热延伸的模型分账**：

| 受控改变 | 五个复频探针、四通道、六q中 max abs(dPi) [fm^-2] | max abs(dC) [fm^-2] |
|---|---:|---:|
| Fermi hard -> PNJL hard，仅占据改变 | 0.116254596 | 1.74482738 |
| PNJL hard -> PNJL thermal24 | 0.040206670 | 1.20574185 |
| thermal8 -> thermal12 | 0.010120635 | 0.18979336 |
| thermal12 -> thermal16 | 1.338135e-4 | 3.590620e-3 |
| thermal16 -> thermal24 | 1.791037e-6 | 5.013112e-5 |

上述改变明显大于节点/实现误差，是定义改变或截断尾部误差，不能靠收紧
unwrap 或强迫两公式点值相等消除。thermal16->24 的最大差出现在检查中的
高q区域；q=1.4时 max abs(dPi)=4.055952e-8、max abs(dC)=9.769601e-7。
这不是介子密度/ratio 的相对误差估计，也没有证明 thermal24 已是精确无限尾。
q=6.2>2Lambda 时两种 hard 变体均为零，因此其占据替换差为零并无特殊
物理含义；thermal延伸仍有贡献，不能误称整个泡必须在2Lambda消失。

**retarded 弱极限**：K+、q1.4，内部lambda中心2、Poisson宽度0.8 fm^-1，
独立实轴 Im B0 加权积分0.05096516530；柱坐标复泡给出0.05096522425，
表示差5.895656e-8<1e-4。gamma从0.2逐次减半到0.00625 fm^-1时，
相应Poisson半群加权值与PV的差从0.0125460单调降到0.000409554。
这是 B0 的积分意义证据，不宣称最后一个有限gamma已与PV相等，更不把它
当作非线性 GBU 密度的 eta 收敛认证。未要求阈值点值逐点相等。
此外144项 k0=0 探针的 pair 均为零、Landau 最大绝对值2.115803e-17，
支持所列背景下的占据差静态消去；不是零动量退化极限或完整 Ward/Goldstone
恒等式的替代认证。

**provenance**：独立读取核对310份源码/配置快照、2份背景输入、5份CSV
输出和原论文PDF哈希，均相符。主manifest SHA-256：
`b80f76f98c5c59cbf2f44b3cb7d50aa8fe1fb5fca5c04a0c45ae361df44d9781`。
两轮目录均保留在 charged_rpa_phase_backend 下，全部 production_authorized=false。

**当前判定与下一步责任**：

- 代码层：当前同域完整泡未发现 Pereira 域/单线项/retarded 符号不一致。
  本轮修复的是新增 oracle/runner 语法，不是通过改变原物理公式修复密度。
- 数值层：384/192节点与所列cut探针已通过；首轮失败属于热尾大区间上的
  求积分辨率不足。分布检查范围有限，不外推到完整密度收敛。
- 公式/模型层：双线交集有直接文献依据，其有限q边界效应仍属于处方依赖；
  文献支持不等于真空 Lorentz 性质、PNJL 热力学闭合和实验可用性已认证。
  全硬截断与当前 thermal-extended 上游也不能混称同一热力学模型。
- 后续由本 track 结合22.15证据形成作者的 regulator/热处方评审；在选定
  版本上分别补 PNJL 静态/Goldstone 一致性、物理极点留数与UHP稳定性、
  精确Mott/Levinson和密度积分收敛，再接入统一研究生产脚本。另一种
  regulator可作明确模型对照，不要求有限q结果相同，也不因Pereira通过
  自动授权删第三条件或择模。现有冻结线图继续保留为诊断，不再手动加密。

**最终验证与交付**：Pereira新增测试364/364；与完整谱泡、regulator、
全gap、研究合同及台账一起运行7220/7220。旧介子质量/密度regimes回归81/81，
未更新baseline。独立重新归约原门槛并逐项核验144 case覆盖通过。
docs、formula-route、active-doc、relaxtime script governance、script-entrypoints、
data-output guard、task-ledger preflight、git diff --check均通过；没有运行全仓
完整测试或新的冻结线积分。本轮没有修改 src 中任何已有实现。

在线核对PR310仍OPEN，base为codex/charged-profile-gates，head为fa34a397；
已提交部分14/14 checks成功，不覆盖本轮dirty源码/测试/文档。没有暂存、提交、
推送或合并，data/outputs未进入暂存区，生产默认与授权状态保持不变。

### 22.17 PNJL gap、极点权重与单球 routing 探针（2026-09-06）

**作者授权与范围**：按已对齐计划继续 PNJL 一致性与极点稳定性审查，然后
再接续 Mott/Levinson 与密度积分。主候选仍为双线交集；另用第一条线单球
探针澄清“仅两个条件”不等同于历史代码。不改 PNJLCore、src 已有实现、
旧默认、baseline、背景或耦合，不重求平衡、不积分新密度、不手动加密冻结线。
本轮为 required_follow_up，仍只在 D:/w/jrt-ord 工作并保留既有 dirty changes。

**方法与实现**：

- `causal_gbu_pnjl_stability.jl` 新增同谱解析导数、gap 内有限差分对照、
  D=2K/F 的实极点权重 W=-2K/F'、逐线性单元的 k0*ImPi 符号检查与 UHP
  矩形计数。实极点权重不乘入每模一个 pi 的 GBU 粒子计数。
- `audit_causal_gbu_pnjl_stability.jl` 复用冻结 T170/muB240 BQS 背景，
  四通道、q=0/1.4/3.2/6.2，分别检查 all_hard 与 thermal24。
  background_gap 重算三个味 A/凝聚/mass 与 K12/K45，不调用 solver 或回写。
- 单球探针同时用角向与两径向坐标；两 A 与 B 始终继承该单球域。探针不
  自动平均 flavor，不重写旧 shifted/log B0，也不宣称大截断近似已经可用。
- 新 [公式文件](../../reference/formula/relaxtime/ChargedGBU_PNJLStability.md)
  区分 Pereira 文献事实、项目谱稳定性充分条件推导、数值诊断与生产授权；
  方法/API 文档同步。`summarize_causal_gbu_pnjl_stability.jl` 独立核验保留
  哈希、案例覆盖及输出中的门槛归约，不冒称重新证明连续核或底层 gap 排除界。

**保留失败及审计器修正**：

1. 首次加载新 runner 的 `-1.+0.6im` 语法在 main 前失败，修为明确浮点数；
   后续测试加载 runner。初版139项通过但有测试局部函数覆盖警告，改用闭包后
   141/141，无需改物理核。
2. `pnjl_routing_stability_pilot_20260906/` 完成16 case，但 pi± q0 的最低
   eta 围道未满足相角步长门槛。固定 anchor 到根的区间不随均匀网格加密，
   最大步长略高于 pi/4。改为细分合并后的完整围道；保留原门槛和失败目录。
   这不是检出不稳定极点。原谱端点的斜率-截距重建还可能制造约1e-14的
   舍入负值，符号检查改用原端点及重心插值，不裁零或修改任何谱权重。
3. 新增围道回归测试初版选择的合成函数不需要加密，只有“必须加密”的
   断言失败；改用确有同向背景相角变化的实根函数。新增 gap 代数测试曾
   对 Tuple 调用不支持的 isapprox，改为向量比较。两者均为测试构造问题。
4. `pnjl_routing_stability_v2_20260906/` 完成32/32案例；30个通过当前
   联合诊断，pi± q3.2 all_hard 的有限差分导数对照失败。极点权重符号正常，
   最大相对导数误差3.052008e-4，超过原1e-5。近阈值有限差分原步长占
   gap 裕量的1/20，尚有可见截断误差。后续加入 gap 内逐次减半，要求
   两次差分变化及解析导数误差都<1e-5，不改变根、谱、物理公式或容差。
   合成 log 端点函数验证减半收敛与用尽预算后的失败状态。

**四动量数值证据（v2，np384/192，谱128/256，cut ne128）**：

| 检查量 | 结果与边界 |
|---|---|
| 同热处方 q0 contact 与独立 A1+A2 | 8/8通过，最大差1.421085e-14 fm^-2 |
| 单球两坐标 Pi/contact | 最大差7.887025e-13 fm^-2；节点差最大1.168083e-13 fm^-2 |
| 实 gap 极点权重 | 两种谱网格共88条记录全部符号通过；每个有根 case 各正负一根 |
| 全实 gap 计数 | 64份插值谱检查通过，未发现同一电荷方向的额外实 gap 根 |
| UHP 矩形 | 192/192解析并给出零绕数；最低 eta=0.003 fm^-1，近轴条带仍未认证 |
| 实根位置 | 最大谱网格漂移0.003342487 fm^-1，发生在 q0；不能称为根位置已精确收敛 |
| 全硬谱被动性 | 16个 fine case 数值条件均通过，少量约1e-37负值不升级为严格数学证明 |
| thermal24 谱被动性 | 16个 fine case 均不满足，min(k0*ImPi)最低约-0.559388 fm^-3 |

fine all_hard：q0/1.4/3.2的12个case各有正负一根，q6.2的4个case无实根；
fine thermal24：q0/1.4的8个case和q3.2的两个kaon各有正负一根，其余无实根。
同一背景下，pi在q3.2的根是否仍在gap已经对热处方敏感；这不是精确Mott定位。
无实gap根也不等于没有cut内共振。

热延伸的有符号高能谱是模型处方风险，不应裁去；它与 GBU 连续密度可负
是不同命题。谱被动性条件不满足不证明 UHP 不稳定；当前有限围道也不证明
完整传播子稳定。该充分条件检查不被静默提升为不可豁免的生产门槛。

**单球 routing 的反例证据**：q0四通道与双线核一致。q1.4时 Pi 的
flavor交换/反射差，pion约3.804829e-5、kaon约0.01525943 fm^-2；
四q最大差0.05123030 fm^-2，明显大于求积误差。A1仍等于平衡A1，但A2
相对平衡值已有q依赖：q1.4 pi+约0.863681、K+约0.539142 fm^-2。
因此“删第三条件即可沿用常数A1+A2/旧路线”不成立。这不等于证明所有
单球处方不可用；交换对称平均或 centered 是其他处方，须单独定义和审查。

**冻结背景的单线一致性**：np192/384重算的三个味 gap 质量残差，
全硬最大0.02166023 fm^-1（约4.27 MeV）；上游默认热上限10 fm^-1最大
2.593114e-5 fm^-1（约0.00512 MeV）；thermal24最大2.669308e-5 fm^-1
（约0.00527 MeV）。对应K12变化约+3.452987e-4、-2.151322e-7、
-8.253685e-7 fm^2。同一热上限的两套求积节点变化远小于这些残差。
原冻结manifest记录了背景求解残差，但未绑定完整solver网格；此处不能把
残差差异全部归为某个假定网格。没有重查Polyakov/守恒荷驻点，也没有
真实手征极限平衡；合成Goldstone代数通过不能代替这些物理验收。

**provenance与生产状态**：v2独立报告核对312份源码/配置快照、2个背景输入、
7份CSV哈希与行覆盖通过。主manifest SHA-256：
`9b26a6587e2cd2f49428eaaca4653b1303d0c67d9035bfa9d5d52b8b175f5c2d`。
所有旧失败目录保留；本轮没有修改src已有实现，没有暂存/提交/推送或合并。
PR310在线核对仍OPEN，base仍codex/charged-profile-gates；已提交fa34a397的
14/14 checks成功不覆盖dirty源码。完整生产、baseline更新与论文预测仍未授权。

**定向复核与误差归因**：
`pnjl_routing_stability_refinement_20260906/` 只取pi±、q0/3.2、all_hard，
np384/192、谱256/512。4/4诊断通过；16条极点、24条围道全部通过对应检查。
同一mesh256根的初始导数差仍能复现v2失败，步长减半3--4次后误差降到
不超过2.925491e-6，两次差分变化不超过8.775599e-6，均低于原1e-5。
证明该项失败来自有限差分验证步长，而非解析导数或极点权重符号缺陷；
没有重新定义物理谱或通过调容差隐藏差异。

q0根的谱256->512漂移为0.000840243 fm^-1，约为128->256漂移的四分之一；
q3.2漂移为2.484924e-6 fm^-1，同样约缩小四倍。这支持所检查根的网格收敛，
不是对所有热处方/冻结线或Mott端点的统一误差界，也不作为GBU密度收敛。
独立报告再次核对312源码快照、2输入、7输出哈希及原门槛归约。
主manifest SHA-256：
`237377e7592943a97ad5adc14ed51b43117c411f9994e6d470be7a43b8518804`。
子集复核不回写v2的30/32，也不冒充完整矩阵已在最新网格重新运行。

**最终验证与后续边界**：新增测试157/157；与谱泡、Pereira、regulator、
全gap、研究合同和台账共同运行7377/7377。225份src文件与本轮pilot快照
哈希无变化；本轮未重跑旧mass/density回归或全仓测试，因为没有改其有效
实现、默认与数据。生产数值回归仍须在后续接入研究入口时重新执行，
不沿用旧通过记录作为本轮的新证据。
docs、active-doc、formula-route、script governance、script-entrypoints、
data-output guard均通过；台账next_action同步本阶段边界。

下一阶段优先处理近轴条带及原连续核稳定性、热延伸的物理性质评审，然后
再做精确Mott/Levinson和非线性GBU的eta/PV及外层积分联合验收。本轮并未
完成这些项目，也没有产出新的density/ratio。单球探针的routing反射失败
阻止将其直接称为已验收的第二生产路线；不能自行以flavor平均或更换regulator
消除失败。未来可以保留多条具名研究路线，但各自须明确公式、状态和适用域。

### 22.18 热谱适用性、原连续核与近轴探针（2026-09-06）

**授权及范围**：按作者同意的“热处方物理适用性、连续核、近轴，然后
Mott/Levinson 与 GBU 积分”顺序执行，分类 required_follow_up。仅在
D:/w/jrt-ord；启动台账 preflight 通过，HEAD仍fa34a397，tracked dirty20、
untracked58的既有改动保留。本轮不修改src已有实现、PNJLCore、旧默认或baseline，
不求新背景、不计算密度、不加密冻结线、不自动变更regulator或flavor平均。

**原文与推导**：

- 重核 Pereira 第9页 Sec.II C：原文选择截断作用于全部发散和收敛项，
  热延伸不是其原样处方。INSPIRE DOI元数据核对 arXiv:2310.05749；
  内置网页网关404，改用公开INSPIRE API/直接HTTP，未绕过访问限制。
- Laine/Vuorinen arXiv:1701.01554 的保留PDF Eq.(8.2)--(8.21) 给出
  Wightman/KMS/retarded谱表示。正热密度矩阵、同一算符及其共轭、
  H-sum(mu*N)演化下，k0*ImG>=0 是完整物理关联函数解释的必要条件。
  带电频率应用、PNJL推广边界及当前D的符号均在新公式文件明确。
- q0正lambda pair可独立写成 ImB0=(2pi*p/lambda)*
  [1(p<Lambda)-1(p<Lth)*(n1+nbar2)]。Lambda<p<Lth是开放区间，
  权重为负，与角度分支或截断端点赋值无关。P通道投影在此为正。
  同时 ImD=8K^2*ImPi/abs(1-4KPi)^2，实contact不能消除这个符号。
- 不把必要谱条件与上一节排除UHP根的充分条件混同；负谱不证明已发现
  UHP不稳定，也不同于GBU连续产额可负。不宣称所有热项不截断方案均不可用。
  目前停止在完整物理粒子谱解释待审，而不是擅自把观察量改成减法相关量。

完整公式、文献/推导/约定/诊断/授权分账见
[热谱适用性](../../reference/formula/relaxtime/ChargedGBU_ThermalAdmissibility.md)。
PDF转文本与渲染中间件保留在tmp/thermal_spectral_review_20260906；
图像读取工具返回不支持image input，未声称已完成视觉验收，未交付新PDF图件。

**新增实现**：

- `causal_gbu_thermal_admissibility.jl`：PNJL正多项式/方差导数、精确q0
  Jacobian、独立Pereira原cut、解析常数减法的PV/近轴Cauchy求积、
  RPA谱符号恒等式及拒绝将采样误差冒充Rouche统一界的标记。
- `audit_causal_gbu_thermal_admissibility.jl`：冻结T170/muB240 BQS；
  四通道、Fermi/PNJL、q0/1.4/3.2/6.2，各3个真空pair支撑外探针，共96条。
  另对pi+/K+、q0/3.2、all_hard/thermal24的gap/unitary/tail，
  eta=0/.003/.001/.0003 fm^-1共96条原连续核探针。
- 原cut取32/64节点，Cauchy Gauss阶数8/12、全局误差目标1e-7/1e-8；
  双配置变化原门槛1e-6 fm^-2。与谱128/256/512分别比较。
  分区使用现有几何支撑提示，谱值不用插值；不是独立全支撑证明。
- `summarize_causal_gbu_thermal_admissibility.jl` 独立检查快照、输入、
  原文、输出哈希及行覆盖/门槛归约，不重新认证连续核物理性质。

**失败保留和原因**：

1. 新占据导数测试初版h=1e-4产生二阶差分截断误差，4个断言失败；
   h=1e-5后通过原atol=1e-9，未改公式/容差。另一个`2.+0.1im`
   Julia字面量歧义修正为`2.0+0.1im`；与物理核无关。
2. `thermal_admissibility_20260906/` 的96个cut全部通过，但96个连续
   探针的收敛布尔全部false。首版给每个区间分配正比其宽度的误差预算，
   对极窄运动学分区反复细分至深度限制，虽总误差估计已小于目标仍拒绝通过。
   保留该目录；这是求积验收器过严的局部预算，不是已发现积分不收敛。
3. 求积改为对全部叶区间误差估计总和验收、优先细分最大误差区间；
   总atol、最大深度和物理谱不变。合成窄区间、内部硬跳变、解析近轴谱、
   预算不足返回失败均测试覆盖。不把估计差提升为严格误差界。

**第二轮结果**：`thermal_admissibility_v2_20260906/`。

| 检查 | 数值证据 |
|---|---|
| 独立cut与现核 | 96/96通过，最大差2.810252e-16 fm^-2 |
| q0解析式与独立cut | 24条，最大差1.942890e-16 fm^-2 |
| 真空外负热尾 | 96/96确认，包括普通Fermi；ImPi范围约[-0.213571,-0.00345472] fm^-2 |
| 原连续核双配置 | 96/96通过；最大变化1.385022e-8 fm^-2，低于1e-6 |
| 细配置总误差估计 | 最大9.951723e-9 fm^-2，低于1e-8；不是严格上界 |
| 插值谱相对原连续核 | mesh128/256/512最大差1.707785e-3 / 4.273650e-4 / 1.068252e-4 fm^-2 |
| RPA虚部恒等式 | 最大差5.551115e-17 fm^2 |
| 原网格contact代数 | 最大差8.886713e-14 fm^-2；此项在eta>=.03计算，不冒充近轴直接二维积分收敛 |

最大插值误差在q0真空cutoff跳变附近的PV实部作用中出现，网格加倍后约缩小四倍；
因此旧插值不是原连续核的逐点数值真值。所测finite-q子集mesh512最大误差
2.662103e-6 fm^-2，不能外推为密度误差或全域uniform bound。

四条PNJL PV tail探针（lambda是真空pair上界加0.25 fm^-1）：

| 通道、q [fm^-1] | ImPi [fm^-2] | ImD [fm^2] |
|---|---:|---:|
| pi+，0 | -0.0525858142 | -0.00180115090 |
| pi+，3.2 | -0.0627948418 | -0.0158892654 |
| K+，0 | -0.0528288115 | -0.00165053043 |
| K+，3.2 | -0.0622369563 | -0.0136394581 |

这说明完整RPA函数保留负谱，不是单纯B0符号或相移后处理现象。
本轮没有算GBU积分，不能据此声称历史ratio偏大/低能趋势全由负尾造成。

**证据与状态**：独立报告核对314份源码/配置快照、2个背景输入、
2份原文PDF、2份CSV及两组96条覆盖，全部通过。v2主manifest SHA-256：
`c192e72824915bc5c5f196efc688b1382569759616786b5b18a3c96e2d642059`。
报告在该目录audit/manifest.json。所有输出仍untracked diagnostic。
近轴只做有限探针，`continuous_UHP_count_certified=false`；
没有复用采样误差作Rouche统一界，没有新的全UHP、Mott/Levinson或GBU积分认证。

**结论及作者决策点**：

- 代码：当前cut的负号由独立实现和解析式复现，未定位为相位/retarded实现错误；
  本轮只修正新增测试语法、差分验证及求积预算算法。
- 数值：原连续核可在所测PV/近轴探针稳定求值；旧谱插值仍有明确网格误差。
  全连续核根计数及非线性GBU的积分收敛未完成。
- 公式/模型：上游热延伸兼容性不足以认证完整物理谱。当前组合在真空支撑外
  只留负Pauli修正，需要明确其低能/减法解释或重新推导处方与观察量。
- 下一步由作者评审：Pereira同硬截断可作为有文献依据的模型候选，但不能
  冒充与现有thermal-extended背景同一驻点；带匹配/减法的热延伸需要新推导。
  本轮不自动选择、重求背景、改PNJLCore或将结果重新标成粒子产额。
  暂停扩大Mott/GBU生产验收，保留现有冻结线诊断及多路线规划。

**最终验证与Git交付**：新增合成/契约测试343/343；连同PNJL稳定性、
Pereira、谱泡、regulator、全gap、研究契约及台账共7720/7720。
台账测试打印的missing-track失败是预期异常测试；独立台账检查和preflight
均通过。docs、formula-route、active-doc、script governance、script-entrypoints
（258个Julia文件）、data-output guard、git diff --check全部通过。
225份src哈希与v2快照一致；本轮未修改其实现，因此未重跑旧mass/density数值
回归或全仓测试，亦未更新baseline。后续接入生产入口时仍需重新执行相应回归。

在线刷新PR310仍OPEN，base codex/charged-profile-gates，head fa34a397；
14/14已提交CI成功，不覆盖本地dirty。没有暂存、提交、推送或合并，
暂存区为空，data/outputs未提交；最终tracked dirty20、untracked63。
本阶段交付的是处方解释的可复查反例及数值误差分账，不是论文级完整链路验收。

### 22.19 方向B前置：静态势一致性与减法反例（2026-09-06）

**作者请求与边界**：先尝试解释/构建热延伸的自洽响应，再视结果考虑方向A。
本轮分类required_follow_up，只在D:/w/jrt-ord；开始preflight通过，
HEAD fa34a397，tracked dirty20、untracked63；全部既有改动保留。
不修改PNJLCore/src、上游热上限、旧默认或baseline，不求平衡、不算新介子密度，
不裁谱、不改相移、不切换A或进行production授权。

**具体执行**：

1. 从有效Models入口核对实际使用的是pnjl_physics/PNJLModel/Core，
   而非旧同名core。上游动能势分真空Lambda与热Lth；
   OneLoopIntegrals.A的质量/分布/上限可与该势的一阶导数相配。
2. INSPIRE公开API进行3组针对性检索并核对精确记录；无专用INSPIRE MCP，
   内置web仍404，arXiv直连成功。保留并阅读2105.14323v2、
   2102.02844、hep-ph/9509363的PDF，另重核已有2301.09882。
   首版元数据打印假定每条有DOI/arXiv字段，遇缺项报错；改为可空数组后
   核对选定记录，不据首次不完整列表作“无相关文献”的结论。
3. 新 `causal_gbu_direction_b.jl` 提供独立每味PNJL动能势、
   凝聚/净密度及合成两次减法谱；新增纯代数/合成测试136/136通过。
   `audit_causal_gbu_direction_b.jl` 调用实际Models动能势做AD一/二阶偏导，
   固定质量和mu，只作Fermi/PNJL、热Lambda/10/24、节点64/128的对照。
4. 新[方向B公式文件](../../reference/formula/relaxtime/ChargedGBU_DirectionB.md)
   划分文献事实、项目推导、静态证据和未实现的动态方案。方法/API/台账同步，
   不更改方法TOML的物理定义；新核仍未选定。

**文献评审结论**：

- Maslov/Blaschke 2301.09882 Eq.(4)--(7)明确从真空截断/热全动量的PNJL势
  求均场，直接支持“热延伸上游有文献先例”，不是charged-BQS完整谱认证。
- Lopes等2102.02844 Sec.II A讨论TRS热不截与MSS介质分离；不能将MSS
  等同于现有泡函数。Xue等2105.14323报告其两味mu0模型的高温凝聚异常，
  这是需检查适用温区的反面证据，不是本项目背景已失败。
- Broniowski等hep-ph/9509363 Secs.9--10区分regulator高能负cut、复极点
  及额外围道项。正则化辅助响应不自动等于全频正物理谱；若用于GBU，须从
  同一泛函核对实轴公式及所有项，不能仅通过重新命名绕过审查。
- 已渲染Maslov/Lopes关键页，图像工具不支持image input，未声称视觉验收；
  本轮以原文文本和方程定位为依据，不交付论文图形复现。

**冻结背景静态结果**：
`direction_b_static_20260906/` 共36条（2占据*3热上限*2节点*3味），全部通过。

| 同一静态处方下的检查 | 最大绝对差 |
|---|---:|
| Models动能势与独立积分 | 5.627499e-12 fm^-4 |
| dOmega/dM 与 Nc*M*A/(4pi²) | 6.321610e-13 fm^-3 |
| -dOmega/dmu 与独立夸克净密度 | 3.608225e-16 fm^-3 |
| dOmega/dM 与独立分布积分 | 6.330492e-13 fm^-3 |
| d²Omega/(dM dmu) 的交换对称 | 2.775558e-17 fm^-2 |

各量使用其声明单位下1e-8绝对门槛。这些是同一给定势的偏导恒等式，
不是完整驻点/最低自由能验证、动态charged Hessian或不同热处方的等价性。
没有重新计算Polyakov场/守恒荷，也没有建立有限q响应；固定热上限数值误差
与无限尾模型极限仍分开。结果反对“仅由负热尾即可判上游不能热延伸”。

**合成减法反例**：
对rho(t)=c*t*theta(t-S)，DeltaPi(s)=-(c*s/pi)log(1-s/S)，
常数和一次项在s0都为零，正cut虚部为c*s。但s=-Q²下补项正增长，
最简单contact-RPA F=1-4KDeltaPi必有UHP零点。
示例c=3/(8pi)、S=49 fm^-2、K=.24 fm²得Q=6.5800650547 fm^-1，
残差4.440892e-16。**这不是实际项目极点，也不是被采用的新核。**
它只说明恢复正谱和静态低阶匹配不充分；不能任意把补项加回真实泡。
对真实charged核，k0静态点与内部lambda0亦不相同，匹配坐标必须先推导。

**provenance**：独立检查316份源码/配置快照、2个输入、3份新PDF、
2个输出哈希，以及36个唯一key和原门槛归约，全部一致。主manifest：
`07aae6bfe9f5b58aa0743455766cd9c2706d23e557e6e20b8dc717feac0787a4`。
另重核的历史Maslov原PDF SHA：
`b54777da443d241b2c208cd8a1be0587dab6fc530b15658743b6e2bc60a2df9a`；
它未被冒充为本轮runner绑定的3份新PDF之一。
数据均为untracked diagnostic；静态通过和合成极点标记分别保留，不改变生产资格。

**下一步**：优先B1——明确带charged外源的真空regulator和介质差泛函，
从同一泛函推导二阶响应与contact/减法，而非仅指定热上限。之后B2检查解析
结构和有效域，B3核对GBU观察量及额外围道/边界项，B4才接计数与密度验收。
普通Gaussian BU不自动推导GBU的sin(2delta)权重，须单独说明。
若B在目标域内无法闭合、需要无依据匹配系数或必须改变观察量/上游，提交
作者评审A或范围调整；无需先普遍证明“热延伸平衡不可用”。当前B未排除，
也未宣称已有可生产的B路线。

**最终验证与交付**：新增测试136/136，连同前序谱/热解释/PNJL/Pereira/
regulator/全gap/研究契约/台账focused合计7856/7856。台账测试中的missing-track
输出是预期异常断言，独立台账检查通过。docs、formula-route、active-doc、
script governance、script-entrypoints（260个Julia文件）、data-output guard与
git diff --check通过。225份src当前哈希与本轮快照一致；未改有效生产实现，
本轮未重跑旧mass/density数值回归或全仓测试，不沿用旧通过冒充本轮结果。
PR310在线刷新仍OPEN，base codex/charged-profile-gates、head fa34a397；
未暂存、提交、推送或合并，所有diagnostic未进入暂存区。

### 22.20 方向B1：固定投影的同源减法泛函与charged Hessian（2026-09-06）

**授权与边界**：作者同意按B优先顺序执行；本轮required_follow_up仍在
D:/w/jrt-ord，preflight通过，HEAD fa34a397，启动tracked dirty20、
untracked67。保留全部旧改动；不改src、PNJLCore、模型TOML或旧默认，
不求新平衡、不积分介子密度、不裁谱、不切换A、不作生产授权。

**推导交付**：见[方向B第6节](../../reference/formula/relaxtime/ChargedGBU_DirectionB.md)。
将regulator明确为与外源J无关的空间动量投影P_L，
先在ran(P_L)限制fermion Dirac算符再取logdet。具名候选为
Gamma_vac,Lambda[J]+Gamma_full,Lth[J]-Gamma_vac,Lth[J]；
真空参考在同一J、质量、mu协变频率取值，去掉Dirac海的源无关mu线性项，
不是重求真空背景。此处只在|mu_i|<M_i和固定Polyakov背景作二阶展开。

1. 固定P_L的二阶变分在两条线各生成一个投影，正好给出双球交集；
   这里delta P=0，所以没有被漏掉的regulator源变分顶角。
   这不认证源依赖截断、局域规范Ward恒等式或对q求导的边界项。
2. 独立4x4 Dirac矩阵投影的P/S迹与现核壳上分子代数相同，且在同域
   分解出C(q)。因而当前完整核可以来自该减法泛函；不需任意补contact。
3. 独立8x8 charged源行列式检验Gamma_q''=-Pi(k0=0,q)，源规范为
   tr_flavor(T²)=1。另有同味退化q0手征曲率恒等式的合成测试。
   这是源Hessian检查，不是相移导数/BU分部积分检查。
4. 复化color holonomy的多项式恒等式可导出项目PNJL占据；
   Phi/PhiBar独立不等实值一般不是unitary SU(3)，不据此证明正热迹。
5. **负热尾仍在**：完整核是vacLambda+fullLth-vacLth的有符号组合，
   泛函存在不等于Lehmann正性或GBU部分粒子产额已有推导。
   本轮没有发现上游热延伸不可用的普遍证明，也没有构造新正谱修正。

**文献复核**：INSPIRE专用工具不存在，内置web网关仍404；公开INSPIRE API
精确查arXiv2310.05749成功，核对本地Pereira文本Sec.II B的determinant
展开和Sec.II C逐线截断/全硬处方。固定投影减法热延伸是本项目明确提出的
解释候选，不能冒充Pereira原文或GBU推导。本轮不交付PDF图形复现。

**新增实现**：analysis-only
causal_gbu_source_functional.jl、audit_causal_gbu_source_functional.jl，
及test_causal_gbu_source_functional.jl；同步unit runner、API、方法和台账。
独立oracle不用现核的B0分子、四留数或谱插值求Pi；使用Dirac矩阵迹和
柱坐标双球求积。静态另直接对源矩阵特征值logdet作二阶差分。

**保留失败与修复**：

- 首轮合成测试140/196通过；测试中的2e1被Julia解析为科学计数常数，
  而非2*e1，修正显式乘号后196/196通过。没有修改物理核/容差。
- 第一轮背景输出direction_b_source_20260906保留：动态32/48，
  静态0/16。动态失败来自32/64内部节点变化，最大2.880414e-5 fm^-2，
  与fine独立核差仅约4e-9不是同一误差。静态汇总新代码中的2f0被解析为
  Float32字面量，改成2*f0；补上积分级源曲率回归，不能仅靠局部block测试。
- 改正新审计算式后，动态使用64/128、静态使用32/64内部节点；
  源步长仍0.004/0.002/0.001 fm^-1。绝对门槛保持1e-6 fm^-2，
  没有放宽容差或更改src。最终新测试208/208通过。

**第二轮真实背景证据**：
direction_b_source_v2_20260906（T170 MeV、muB240 MeV的保留quark-only BQS）。
动态覆盖pi+/K+及其负通道反射、P/S、Fermi/PNJL、q=0/1.4/3.2 fm^-1、
热上限Lambda/10 fm^-1，共48组，每组3个external k0 UHP探针：
0+0.8i、2+0.8i、6+1.1i fm^-1。48/48通过。
静态是pi+/K+、P、Fermi/PNJL、q=0/1.4、两热上限，共16/16通过。
Fermi控制只改变占据，不重求mass、mu或Polyakov平衡。

| 同一具名处方的检查 | 最大绝对差（fm^-2） |
|---|---:|
| 独立Dirac源核与现核 | 7.212045e-13 |
| 动态两套内部节点变化 | 2.102986e-9 |
| flavor交换+频率反射 | 3.997027e-14 |
| 现核同域contact分解恒等式 | 9.991968e-14 |
| 静态特征值源曲率与-Pi | 2.780913e-8 |
| 最后两源步长的曲率变化 | 6.610253e-8 |
| 静态内部节点变化 | 6.681622e-11 |

**provenance**：第二轮结束后、更新台账前，独立验证318份源码/配置与快照
和当前文件一致；2个冻结输入、2个CSV输出哈希一致，48/16行唯一key及
原门槛归约一致。manifest SHA256：
25bb30ca513d356df76a85d5fec53c8bc032e2c4f2e32196b405ab58855ffd2a。
其后的台账更新是治理状态记录，不改本轮绑定的数值实现；所有diagnostic
仍untracked，首轮失败产物未覆盖。

**结论及下一步**：B1将问题从“现双球/热延伸核是否必然漏源项”
缩小到“有符号减法响应能否推导所需GBU观察量”。现核的fermion二阶
来源得到代数和有限探针验证，不认证完整PNJL+KMT驻点、全连续UHP结构、
Mott/Levinson、密度或正物理热迹。下一步B2/B3先检查减法辅助核的
解析结构、额外围道/端点和GBU关联泛函来源；普通Gaussian BU不自动
推出GBU权重。若须新UV定义/匹配参数或换观察量，先作者评审，不自行实施。

**测试**：相关focused共8064/8064通过（新208及前序7856，本轮均重跑）；
台账missing-track FAILED仍是预期异常断言，非仓库门禁失败。
本轮不改有效数值src，未重跑全仓/历史mass-density生产回归；
诊断脚本通过不冒充完整生产回归通过。治理和最终Git核对见本节收尾记录。

**收尾（2026-09-07中断恢复）**：docs consistency、formula-route、
active-doc governance、relaxtime script governance、script entrypoints
（262个Julia文件）、data-output guard、task-ledger与git diff --check均通过。
225份src相对B0快照未变。B1运行后只有台账文件与绑定的源码/配置版本不同，
原因是记录本轮结论；数值实现没有运行中修改。
恢复后preflight仍通过，tracked dirty20、untracked70，暂存区为空；
PR310仍OPEN，base codex/charged-profile-gates，head fa34a397。
未提交、推送或合并，未暂存data/outputs；所有既有改动与首轮失败证据保留。

### 22.21 方向B2/B3：辅助响应、GBU权重与部分产额语义（2026-09-07）

**授权与范围**：作者要求继续优先验证B，分类仍required_follow_up。
完整阅读仓库/技能约束并通过task-ledger preflight，启动时tracked dirty20、
untracked70、stage为空。保留全部既有改动；只新增analysis测试/runner与文档，
不改src、PNJLCore、方法TOML、旧默认或baseline，不重求上游、不算新密度。

**文献与推导**：见[GBU观察量闭合条件](../../reference/formula/relaxtime/ChargedGBU_ObservableClosure.md)。
使用research-engineer区分文献事实、独立数学恒等式与项目定义；julia-pro要求
纯合成代数测试；doc-coauthoring保留未决物理状态。INSPIRE专用MCP缺失，
依mcp-inspire技能改用公开INSPIRE API及arXiv。新读2512.03876v1与1612.09556，
复核1305.3907和1912.13162。pdf技能要求渲染核对，但当前image-input工具不支持
查看渲染页，未声称视觉验收；引用使用已核对文本与页码。

1. 1305.3907的Gaussian logdet虽称generalized BU，其推导仍是delta形式；
   不能仅据题名推导正弦项。2020 Eq.(16)是目标GBU部分密度的明确文献先例。
2. 两loop Phi路线将相移与自能项组合成W=delta-sin(2delta)/2；完整热力学
   解释还需同一Phi的驻点与抵消，不能由B1单独认证。
3. 项目D=2K/(1-4KPi)对应介子D的自能Sigma_M=2Pi，不是介子诱导的夸克
   自能Sigma_q^(M)。独立代数证明Im(Sigma_M)*Re D=sin(2delta)/2，
   不要求ImPi非负；Im(Sigma_M D)及Born相位都不是目标修正。
4. 减法先作用于Pi，再统一RPA；各片分别重求和再相减一般不同。
5. 有符号Cauchy核可以开UHP解析，但F=0仍可产生RPA不稳定极点。
   新增带条件的外圆无零点界、近零分支尾项界；真实核没有统一误差/总变差
   认证，状态仍false。有限probe通过不冒充全UHP证明。
6. 固定profile的部分记账偏导与包含profile变化的总导数不等价。新合成
   测试分账端点与额外导数，不给上游添加新化学势或反馈。B保持已约定的
   固定quark-only BQS GBU部分产额候选，不静默改成另一观察量。

**新增实现**：causal_gbu_observable_closure.jl、
audit_causal_gbu_observable_closure.jl、test_causal_gbu_observable_closure.jl及unit入口；
同步API、方法说明和台账。实际runner不取相移导数、不数新根、不积分密度。
新模块纯代数/合成测试最终841/841通过；首次702/703，失败断言把合成正响应
误要求为>0.01，而实际0.0030205715610241。按证明改为>0，并额外验证其等于
total-tag；这是反例断言错误，不是物理核/物理门槛调整。随后新增外圆与分支
条件测试；测试负谱、负部分贡献和不稳定辅助极点均保留。

**真实背景诊断**：direction_b_observable_20260907；只读T170 MeV、muB240 MeV
quark-only BQS，PNJL占据、Lth=10 fm^-1。pi+/K+、q=0/1.4/3.2，分别在Landau、
gap、pair、负热尾取局部probe。内部cut节点32/64，Cauchy求积8/12节点与
绝对误差估计目标1e-7/1e-8 fm^-2；对照门槛仍1e-6，cut一致性1e-9。
24/24通过数值和代数检查。直接PV实部不依赖谱插值，虚部来自独立Pereira cut。

| 检查 | 本轮证据 |
|---|---:|
| 最大两节点Pi变化 | 1.4126884e-8 fm^-2 |
| 最大GBU权重恒等式差 | 5.5511151e-17 |
| 最大光学代数恒等式差 | 2.7755576e-17 fm^2 |
| 三个external UHP k0点、两通道三q的log非交换对照 | 18条 |
| 上述combined与piecewise复log差的模 | 0.0164766--3.2296105 |

6个thermal-tail探针的ImPi均约-0.053至-0.063 fm^-2；局部主值delta约
-0.013至-0.045 rad，W约-1.5e-6至-5.9e-5。**负值保留，GBU只压低小角
幅度，不改符号。** 这不是整段尾积分上界，也不是物理密度分支认证；不得
据此宣称负尾已可忽略、总密度变正或已解释旧冻结线全部差异。

**provenance**：320份src/config/analysis源码与source_snapshot及运行后当前
文件一致；2个冻结输入、2个CSV、4份PDF哈希匹配。manifest SHA256：
6af5a0d2784ca1b2d074321951e6d6c8af4067cffb70ffade069db18cc5b716c。
台账在运行结束后更新，不改其绑定的数值实现。输出均untracked diagnostic。

**阶段结论与下一步**：B2/B3的局部代数和观察量边界更清楚，B未被排除，
但还没有完整驻点推导、全连续核UHP零点认证、Mott/Levinson或密度验收。
下一步优先原连续核的目标域近轴/外圆稳定性，随后将独立根计数、Bose端点、
有符号连续项及尾误差接入既定partial-yield积分；成功后才接统一研究生产。
不按“更像实验”择模，不自动切A；若需新UV参数、改观察量或加入反馈先请作者评审。

**收尾验证**：本轮新测试841项及相关focused共2743/2743通过，覆盖源泛函、
方向B静态、热谱、PNJL稳定性、gap/full-closure、Pereira、因果谱泡、研究合同
与台账；missing-track FAILED输出是预期异常测试。未重跑全仓生产回归，
原因是本轮没有更改src、默认或新增生产计算；不能将focused通过替代该回归。
docs consistency、formula-route、active-doc governance、relaxtime script governance、
script entrypoints（264个Julia文件）、data-output guard、task-ledger preflight
及git diff --check通过。CSV唯一key为24，最大PV虚部/原cut差4.4408921e-16 fm^-2。
运行后只有台账更新造成一份配置与manifest不同，225份src及分析实现未变。
最终分支/HEAD不变，tracked dirty20、untracked74，stage为空；PR310 OPEN，
base codex/charged-profile-gates。未提交、推送、合并或暂存data/outputs。

### 22.22 前三步联合推进：连续核端点障碍、源延续与正常Mott（2026-09-07）

**范围**：作者要求尽可能一次完成剩余四步或先前三步。本轮按解析结构→计数→
积分的依赖推进；并未把“执行检查”当作“验收通过”。启动preflight通过，
tracked dirty20、untracked74、stage为空。全部工作在D:/w/jrt-ord；
不改src、PNJLCore、方法TOML、旧默认或baseline，不裁谱/密度、不fold积分。
使用research-engineer区分数学证明与数值证据，julia-pro补纯合成回归，
doc-coauthoring记录失败与授权边界；没有启动子代理。

**本轮没有完成前三步整体验收**。完成了原连续核的一个明确端点反例、
解析外圆界、跨单味onset的源参考延续，以及正常低能Mott的局部检查。
全连续核圆内/近轴的全部零点、包含辅助UV结构的Levinson和真实GBU密度
仍未闭合；第四步研究生产接入未执行。B未因此被否定，也没有自动切A。

#### 22.22.1 硬热端点不能用谱插值计数替代

推导见[端点与源延续](../../reference/formula/relaxtime/ChargedGBU_EndpointClosure.md)。
投影Dirac迹给出谱总变差的解析上界，从而排除一个显式外圆之外的F零点；
未给出圆内/近轴无零点证明，也未作有向舍入认证。

q0、Lth>Lambda、T>0时，正频热支撑端点rho(S-)<0。分离
c/pi*log[(S-z)/(S-h-z)]后余项有限，F(S+0)=-infinity、F(infinity)=1，
故**至少一个cut上方正频辅助实零点**存在。它不同于正常正谱束缚态，
不是新发现的UHP不稳定。线性谱插值把端点平滑到零，会移除对数奇性，
故前序插值核的“全gap通过”不能晋升为原连续核全谱证明。

direction_b_endpoint_closure_20260907：只读T170/muB240保留BQS，四通道乘
Lth=10/20/24，共12/12端点证据通过；正常低能gap仍每例1个正根。
log10[(r-S)/(fm^-1)]的渐近估计如下，不是Float64直接求到的根或残差：

| 通道 | Lth=10 | Lth=20 | Lth=24 |
|---|---:|---:|---:|
| pi+ | -6525 | -1.580e8 | -1.119e10 |
| pi- | -6174 | -1.496e8 | -1.059e10 |
| K+ | -6343 | -1.477e8 | -1.038e10 |
| K- | -10185 | -2.307e8 | -1.615e10 |

q0无q²相空间贡献，有限q端点通常连续趋零，不能直接外推此定理至全部q。
单根Bose上界不是根数上界或总密度界，也不能据本证据归因旧的负密度。
合成有符号band模型独立验证：连续谱W下降-pi与负谱留数辅助根的+pi跳变
相互补偿；只数正常正谱根会漏一单位。新增带权积分测试保留离散根和连续项，
检验分部积分与eta→0加权极限，不要求端点或极点处逐点相等。

#### 22.22.2 原cut Mott检查及保留的诊断程序失败

causal_gbu_continuation_gate.jl使用原Pereira cut，无谱插值，按物理断点作
cosine映射；Float64无可表示内部节点的窄panel改用128-bit坐标重算同一cut。
不把节点夹到端点、不丢panel，Gauss权仍为Float64，故不是全核任意精度证明。

第一轮direction_b_continuation_20260907保留：5条normal行、12条finite-q行，
其他案例因未解析窄panel/nonfinite或局部变量error遮蔽Base.error而失败。
这是新诊断代码问题，不是物理模型失败。合成测试第一次409通过/3错误是
source_domain误用Symbol索引普通Vector；改为NamedTuple后通过。
另一个normal=false原因是阈上1e-7探针32/64节点差2.819e-6 fm^-2，
超过原1e-6门槛；不是改变相位符号或裁负尾可以解决的问题。

第二轮direction_b_continuation_v2_20260907：修正上述诊断实现，Cauchy用64/128，
阈上偏移用1e-9/1e-11；门槛未改。8例均有正常根1→0、相位减少约pi，
完整normal门禁初次6/8通过：K+ q1保留背景区间宽0.46875 MeV，
K- q1正常根128/256内部节点漂移5.914e-7 fm^-1，均如实保留失败。

finite-q检查为四通道、q=0.01/1.4/3.2、Lth=10/24，共24/24通过；
两阶Pi最大差1.478e-14 fm^-2，端点F实部最小1.0013764。
只证明这些端点探针，不证明所有外部实根或UHP根不存在。

direction_b_mott_refinement_20260907只补一个T=196.015625 MeV、muB240、
p24 quark-only BQS背景（不是冻结线加密），残差1.646e-14；质量/Phi位于
保留邻点之间。K- q1改用256/512内部节点，漂移降至2.843e-10 fm^-1。
结合具名复用的前序成功案例，正常Mott共8/8通过；所有区间宽<=0.25 MeV，
相位减少量除pi在0.9999723--0.9999936之间：

| 通道 | q0温度区间（MeV） | q1温度区间（MeV） |
|---|---|---|
| pi+ | [203.6328125,203.75] | [194.140625,194.375] |
| pi- | [203.28125,203.3984375] | [193.7890625,193.90625] |
| K+ | [201.7578125,201.875] | [195.78125,196.015625] |
| K- | [204.1015625,204.21875] | [197.890625,198.0078125] |

这是Lth10、muB240、p24背景上的正常局部Mott；不是Lth收敛、全冻结线或
全谱Levinson验收。q1个别区间与前序插值阈值结果不同，不能覆盖旧数据假装一致。

#### 22.22.3 源参考跨夸克起始面：不改现有核

旧B1 oracle用|mu_i|<M_i保证grand-energy<0等于负Dirac能带，
但Mott附近部分味已越过该限制。新增causal_gbu_source_continuation.jl
将参考明确为固定最低四个Dirac海能带，full仍用PNJL，thermal=full-sea；
全动量参考隙min(2M1,2M2,sqrt(q²+(M1+M2)²)-|mu1-mu2|)>0时可延续。
源步长还满足sqrt(2)|J|<gap。异能带隙关闭仍拒绝，旧oracle限制保留。
这是局部二阶源Hessian延续，不是整个PNJL/KMT驻点或全GBU泛函证明。

direction_b_source_continuation_20260907使用T205.625/muB240保留BQS，
四通道乘q0/1，共8例，均已超出旧pair证明范围。64/128节点初次6/8通过；
pi± q0内部节点差最大2.363e-6 fm^-2，但fine现核parity约3e-11。
未改变容差或热处方，direction_b_source_continuation_refined_20260907
针对两例用128/256，二者均通过：node差<7.75e-11、动态parity<2.30e-12、
源曲率差<2.16e-8、源步长差<3.32e-8 fm^-2。最终覆盖8/8。
8个Mott事件的参考能带隙也均开启。此结果补上原证明的覆盖缺口，
不把“mu>M”误判为PNJL背景不合法，也不消掉负热尾或端点辅助结构。

#### 22.22.4 provenance、验收与下一步

全部CSV、manifest、source_snapshot为untracked diagnostic，没有提交data。
独立重算六个产物目录的源码快照、输入、输出哈希，无不一致：

| 目录后缀（均在charged_rpa_phase_backend下） | 源快照数 | manifest SHA256 |
|---|---:|---|
| direction_b_endpoint_closure_20260907 | 322 | 5bab9fbb0966eed36730affd721fefa3643c93043e44a4eb4af07ab7d19f7e78 |
| direction_b_continuation_20260907 | 324 | 500d3de493eb7182a5f5e40808cd9c8c8389ae9bab040055e99c765306386d21 |
| direction_b_continuation_v2_20260907 | 326 | 4d5c36ae3909c707fb84aa8643700984789b9ddb04948019592082f281c40d8c |
| direction_b_source_continuation_20260907 | 326 | 0bb9b83dd4e8362c30c955345b5eeda3056a407e546329a8cde4976ceaa4e0e0 |
| direction_b_source_continuation_refined_20260907 | 327 | 60cb0a63e461790e74fe7a169e4170045a2f75e2025f3884f2de1441cd45b6df |
| direction_b_mott_refinement_20260907 | 327 | c13ebbe6785758e1bec4b2a4a2c0b00a9fb95fa3fa65e419a50aec1a3823218c |

新增三个模块、四个audit入口与三个unit文件；所有真实密度积分仍未执行。
原225份src与本轮端点快照完全相同。唯一新平衡解为上述Mott定位点，
没有修改上游模型；其余真实探针复用保留背景，旧失败输出不覆盖。

**当前验收状态**：步骤1有外圆界，但全圆内/近轴计数未完成；步骤2正常Mott
局部通过，但完整有符号Levinson未通过；步骤3只完成合成带权/eta回归，
真实密度的联合端点/截断/节点收敛未执行；步骤4不接入生产。
下一步优先显式端点对数分离与辅助零点/连续谱补偿，或对原本意图的无限
热积分建立可控尾极限；先不换观察量或UV定义，不删辅助根、不自动切A。
研究生产必须在该障碍及目标域原核计数通过后才接入统一入口。

**测试**：本轮新增纯合成测试514项（端点409、原cut求积53、源延续52），
相关focused含前序源泛函/GBU代数/热谱/PNJL/全gap/Pereira/因果泡/研究合同/
台账，共3257/3257通过。missing-track FAILED是预期异常断言。
未跑全仓生产回归：本轮没有改src或默认，也没有运行生产密度；focused结果
不能替代该回归。

**收尾**：docs consistency、formula-route、active-doc governance、relaxtime
script governance、script entrypoints（271个Julia文件）、data-output guard、
task-ledger preflight、git diff --check均通过。最后preflight为tracked dirty20、
untracked85；暂存区为空。分支/HEAD保持codex/charged-phase-coordinate-cut-fix /
fa34a397，PR310仍OPEN、base codex/charged-profile-gates。
未提交、推送、合并或暂存data/outputs；既有修改和首轮失败产物均保留。

### 22.23 无限热目标、加权误差与研究生产阻碍（2026-09-07）

**作者定义与范围**：有限Lth是数值上限，理论热积分到无穷；真空Lambda仍是
模型参数。固定quark-only BQS上的GBU部分产额，不要求先完成介子反馈驻点。
本轮属于当前charged任务的blocker/required_follow_up，不切A、不改PNJLCore、
默认或baseline。只在D:/w/jrt-ord工作，启动preflight通过，既有dirty保持。
使用research-engineer、julia-pro和doc-coauthoring；未启动子代理。

#### 22.23.1 不再把数值硬热边当作模型边界

新增原cut的无限内部动量积分、独立径向Dirac投影对照，以及保留真空硬边的
快速Cauchy表示。无限尾采用半无限坐标映射，数值分区之外的尾仍实际积分，
不是换一个更大的硬Lth。公式见
[无限热积分](../../reference/formula/relaxtime/ChargedGBU_InfiniteThermal.md)。

快速表示解析分离unitary平方根项，余项插值。q0真空端点保留左右极限；
有限q的连续节点两侧共用同一数值，不制造两条不同的外推极限。
已知运动学gap及其连续端点严格为零，避免平方不变量舍入造出伪谱支撑。
这是解析域定义，不删除负热尾或Landau。

本轮保留的新诊断实现失败及修正：

- +0/-0被unique视为不同值，产生零宽panel；规范化有符号零。
- 新径向对照把`2E1`写成Julia科学计数常量20；改为显式`2*E1`。
- 很长但有限的Landau能量区间用线性Gauss时漏掉热边界层；改用保留有限终点的有理映射。
- 对已是Float64的插值器传播宽坐标，造成不必要的BigFloat全cell计算；固定表示精度。
- 在连续finite-q节点分别外推左右值，触发伪jump；改为共享连续值。

`direction_b_infinite_20260907`保留4条shell和部分probe，因性能诊断中断并补
partial manifest；不是完整验收。后续完整证据与该partial输出分开。
单元测试先后出现132/134、172/179，通过修正上述实现后179/179；加入计数、
分账及全gap测试后325/325。没有改变物理核符号、相位fold或门槛。

#### 22.23.2 原核、实轴近邻与带权eta极限

`direction_b_infinite_acceptance_20260907`初次2例成功、2例因上述连续节点伪jump
失败，原数据保留。`direction_b_infinite_acceptance_v2_20260907`使用保留
T170/muB240 BQS、pi+/K+、q=1.4/3.2，4/4通过。

原cut与快速核对照门槛1e-6 fm^-2不变，最大差1.48e-7 fm^-2。独立正负gap
根、根邻域围道及实轴Nyquist均通过；外半径28/30、根缩进1e-5/3e-6、
两阶节点一致给出UHP计数0。根缩进中的区域由独立根围道检查，不只检查
Im(z)>=某个有限eta的矩形。仍是数值计数，不是有向舍入的全频证明。

| 通道、q（fm^-1） | Lth10与无限目标的shell相对差 | eta=0.0003与PV的带权相对差 |
|---|---:|---:|
| pi+、1.4 | 7.49e-5 | 5.00e-5 |
| pi+、3.2 | 2.51e-4 | 3.14e-3 |
| K+、1.4 | 5.64e-5 | 5.07e-5 |
| K+、3.2 | 1.54e-4 | 1.50e-3 |

eta使用0.003/0.001/0.0003；Bose下端为0.001/0.0003，均为fm^-1。
比较保留`gW`两端边界项，先固定IR端点取eta极限，再检查IR变化；不要求
极点、阈值或Bose原点逐点一致。束缚态在有限eta整体积分中不再另加一次。

例如pi+、q1.4，bound=2.81625e-3、Landau=2.33888e-7、pair=-3.81524e-4，
总shell=2.43496e-3 fm^-2。负连续贡献保留。该结果不等于积分后的总密度，
也不能据此声明所有冻结线点已通过。

#### 22.23.3 无限目标的正常温度Mott检查

`direction_b_infinite_mott_20260907`将原8个q0/q1事件转移到无限目标。
初次7/8：K+ q1在旧下端T195.78125的F(U)=+1.64113e-4，旧区间不再夹根；
前后正常根1→0及相位补偿仍通过。这是数值热上限改变后的Mott定位微移，
不是新增物理相或上游失效。

`direction_b_infinite_mott_refined_20260907`只读两个保留背景，改用
[195.6640625,195.78125] MeV；F(U)分别为-4.42361e-3、+1.64113e-4。
没有重求或插值平衡态。具名复用其余7个事件及原前后计数，最终8/8通过。
全部区间宽<=0.25 MeV，根数减少一，相位下降/pi在
0.99997224--0.99999356之间。最大原cut节点差8.52e-12 fm^-2。

#### 22.23.4 自动外层验收与保留的失败

`direction_b_infinite_density_20260907`在同一背景上完成pi+/K+的首轮外层探针。
未分段q网格8/16阶下，总密度较接近，但bound与pair的分账收敛慢；原因是
Mott动量处离散态边界项变化，两项误差会相互抵消。
q>约14时固定omega_max=24触发原UV相位门槛，7个尾区失败保留。
不把这些失败补零，也不输出看似完整的积分值。

后续`audit_causal_gbu_infinite_production_gate.jl`从原核F(U,q)=0自动定位Mott
动量并划分q积分；同时划分2Lambda和尾区，omega上端随q提高。
对bulk每个节点独立检查正负解析gap及cut实轴交点，不由unwrap数根。
较近Mott的固定频率偏移触发相位极限未收敛，失败点继续保留。
其解析补救是利用rho(U+u)=A sqrt(u)+o(sqrt(u))及非零实F(U)求单侧相位极限，
不把独立根数作为相位极限的输入。阈值/窗口纯代数及数值类型合同最终137/137通过。

首轮四通道外层输出406行、10个固定偏移失败，另有4行围道验收失败。
后者包括低q围道分辨率不足，以及K+ q=3.53857009的负频根被固定1e-7边距
排除。增加节点只能解决前者，不能解决后者；按根到gap端点距离调整围道
窗口才恢复正确计数。保留原非零边界、pi/2相邻辐角门槛，不放宽容差。
合成near-endpoint根回归明确演示旧窗口计数0、新窗口计数1。

`direction_b_infinite_gate_recovery_20260907`验证旧输入、输出哈希及源码变化
白名单后，复用402个原先已通过的行，自动重算14个缺失/计数失败行。
允许的旧源码变化仅为Y的解析阈值guard和Q的围道窗口/加阶；物理核、
密度公式未改。最终416/416行齐全、0失败；其中288个bulk节点明确做了
双gap及cut拓扑检查，128个尾区节点只参加积分/尾部验收，不能合称416个
全谱认证。近Mott原始失败CSV不覆盖、不补零。

同一保留T170 MeV、muB240 MeV、quark-only BQS背景的积分结果为：

| 通道 | 无限热密度（fm^-3） | q 8→16阶总密度相对变化 | L10正常分支加权绝对差/无限总密度 |
|---|---:|---:|---:|
| pi+ | 0.005377208397 | 8.47e-5 | 1.30e-4 |
| pi- | 0.005981917841 | 9.34e-5 | 1.27e-4 |
| K+ | 0.002964868724 | 1.20e-7 | 8.82e-5 |
| K- | 0.001339314252 | 7.32e-8 | 1.13e-4 |

以上不是新冻结线扫描，也不应与实验最终产额直接等同。最大bound/pair/Landau
分账变化/总密度为9.34e-5，故总数稳定不是仅靠两项误差抵消。
q>8的尾区贡献占总数3.65e-5至1.72e-4；四个尾带递减，最后24--32带
占比不超过5.77e-11。omega条件尾界约1.13e-12 fm^-3。q>32采用连续尾带
数值收敛判据，不是全区间解析上界；omega界仍以该分支控制成立为条件。

**L10解释限制**：该列比较同一正常分支/有限频率窗，不包括对有限硬热边
所有辅助根的重新计数。因此不能把约0.013%的差称作旧硬端点全谱算法的
严格误差界，更不能据此保留漏根实现并宣布无影响。可行方案是直接计算
无限热目标；这保留了热尾，而不是删除辅助根后沿用原有限Lth算法。

#### 22.23.5 独立原核与最终自动验收

`direction_b_infinite_oracle_20260907`补成完整的4个q0/finite-q案例：
pi+/K+、q=0/1.4，20/20原cut/快速表示探针通过；其中8个复频率点的
独立径向Dirac积分与原cut最大差1.27e-14 fm^-2，快速表示最大差
3.16e-7 fm^-2。非零q的mesh/omega节点加阶密度相对差不超过6.42e-7。
q0密度因相空间因子为零，不能用其零差替代相位核对照。
保留CSV中阈值pi曾以Julia常量文字序列化；随后明确返回Float64，新增16个
类型断言。该修正只影响输出数值类型，不改极限、计数或产额。

新增`audit_causal_gbu_infinite_readiness.jl`统一验证证据目录的源码快照、
输出和背景哈希，再检查完整行数、eta/IR序列及各阶段结果，输出
`direction_b_infinite_readiness_20260907/checks.csv`：**5/5通过**。
五项为无限原核对照、代表点近轴/根计数、PV/eta加权极限、正常温度Mott、
四通道外层积分。空记录、缺eta/IR组合及非收敛序列不会被all([])误判通过。

| 完整证据目录 | manifest SHA256 |
|---|---|
| direction_b_infinite_oracle_20260907 | 417e979f8cfa4a7f6db042021a1996c1ac047d7a1e66311471b51492a50ebe1f |
| direction_b_infinite_acceptance_v2_20260907 | fb2badde3acde908f7314372c5073a88f658f52aaf802c746d434dfe0b18eae9 |
| direction_b_infinite_mott_refined_20260907 | 7df40d0c1363387185a0f24bff8902ddd538b42caaac829e99023975b93096cf |
| direction_b_infinite_gate_recovery_20260907 | 80bbfafde5b1f51096384f23da7cf5465218f149c35f82381cabc8ea244214c4 |
| direction_b_infinite_readiness_20260907 | 2e8bbde4095990fa175ab11042473d5b7e533619068d5f00ce2589e8e1ea7a9e |

#### 22.23.6 明确结论与生产边界

1. **无限热B路线：在固定quark-only BQS的GBU部分产额定义和已验收背景内，
   数值上可行。** 正常根独立计数、Mott补偿、近轴积分极限与外层分账收敛
   已有真实证据，不再以“有限Lth辅助根存在”作为该目标不可行的依据。
2. **旧有限硬端点/平滑插值路线：不可直接作为全谱计数正确的生产方案。**
   端点辅助根与插值原核不等价的诊断继续有效；未证明其全谱误差可忽略。
3. **本轮不授权替换默认或冻结完整论文曲线。** 这是有符号辅助响应上的固定
   背景GBU部分产额，不是全频正Lehmann谱、完整介子反馈驻点或实验最终产额。
   数值围道/采样不是全频区间算术证明，完整冻结线也尚未逐点验收。
4. 下一项工程工作应是将无限热表示和自动失败门禁接到**显式选择的研究生产
   provider/脚本**，由该入口自动验收新背景并保留失败原因；旧默认不变。
   出现额外拓扑、Bose不安全或误差超标时拒绝给出产额，不能静默沿用当前
   单正常根公式。此覆盖问题不要求重开完整介子反馈推导，也不靠手动加密
   冻结线解决。正式论文方法选定和production晋升仍由作者评审。

相关focused 3708/3708；新增最终归并器8/8及阈值Float64新增16项另行验证。
按不重复计数共3732项。全仓生产回归未运行，因为本轮225份src相对端点
快照全部未改、默认和baseline未动；focused不能冒充生产全仓回归。
此轮技能用于方法/数值/证据边界，未改变用户授权范围。

**收尾核对**：docs consistency、formula-route、active-doc、relaxtime script、
script entrypoints（285个Julia文件）、data-output guard、task-ledger preflight
和git diff --check通过；台账更新后另跑185/185合同测试。
tracked dirty20、untracked107，暂存区为空。PR310仍OPEN，base为
codex/charged-profile-gates，分支/HEAD仍为codex/charged-phase-coordinate-cut-fix /
fa34a397。未提交、推送、合并或暂存data/outputs；既有src改动保留，本轮未改。

### 22.24 显式研究生产入口与完整冻结线扫描

作者授权将已打通的无限热GBU链路接入生产入口并产出完整冻结线ratio。
本轮为required_follow_up；仍仅在D:/w/jrt-ord工作，preflight与dirty核对后实施。
使用research-engineer、julia-pro、doc-coauthoring，未启动子代理。

新增`Models.run_charged_gbu_freezeout_scan`及
`scripts/relaxtime/run_charged_gbu_freezeout_scan.jl`。统一workflow负责冻结线
映射、BQS求解、四通道、Mott分段、加阶、门禁、checkpoint和图像；CLI不手拼
物理链路。旧入口/旧regime默认不变。已验收核暂保留analysis单一源码，由
显式研究适配器调用，非悄悄晋升全部analysis API。入口合同见
[ChargedGBUResearchWorkflow](../../api/relaxtime/meson_density/ChargedGBUResearchWorkflow.md)。

方法配置`charged_gbu_infinite_v1.toml`与旧research_v1分开，明确热目标无穷。
默认完整能量网格为3、5、7.7、11.5、19.6、27、39、62.4、130、200 GeV；
沿项目default冻结线由高能到低能延续quark-only BQS，48节点、残差<=1e-7。
每个通道q1.4重新验收原cut/径向/快速核、两围道及eta/IR序列，再进行q积分。
bulk8/16阶不满足预设门槛则32阶，不修改容差。由于wide-coordinate fallback
会使用setprecision，通道并行使用独立进程而非共享线程。

输出目录为`data/outputs/results/relaxtime/meson_density/charged_gbu_infinite/`
下新建的`freezeout_20260907`；419份源码/config/script绑定到run.json快照。
checkpoint带SHA256，失败JSON用null、ratio CSV用NaN，不用零替代失败。
初次合同测试17通过、1个测试语法错误（Dict误用.get），修正测试后18/18。
相关unit485/485；新入口集成11/11，旧统一入口/external-path集成11/11，
旧介子密度regime和旧冻结线GBU路径回归389/389，未改baseline。
入口/文档/公式/active-doc/script/data-output治理均通过。

完整扫描已由sysimage wrapper以fallback启动；无可用sysimage时使用普通Julia。
初次启动出现Julia1.12惰性模块binding的world-age警告，当前版本可运行；
不在源码绑定的长计算过程中改写加载器。10个背景全部成功，残差<=3.21e-14。
最终ratio、失败点及manifest验收在本节后续记录。

首轮200 GeV四通道及130 GeV两个pion通过后，62.4 GeV两个pion触发PV节点
不一致。独立复核在U=3.832563530306745发现余谱节点为0与6.02e-9：
`a+(b-a)*sinpi(1/2)^2`将右端舍入到U+1ULP，再追加U产生反向小区间。
这不是物理cut不连续。原运行停止，保留run.json/source_snapshot/checkpoint
及interruption.json；不把partial目录声称为完整验收。

修复`panel_nodes`：先精确指定左右端，再去重排序并断言严格正宽。
未调整PV容差、占据、截断、相位或GBU公式。原失败探针复查ordered=true、
threshold_cut=0、F(U)=-0.7821659422660217、shell=0.0020695302513674 fm^-2。
合成端点/相关profile/qgate/入口合同289/289通过，入口集成11/11重跑通过。
惰性入口改为在invokelatest闭包内访问新binding，消除Julia1.12警告。
新目录`freezeout_20260907_v2`通过同一生产CLI完整重算，四进程隔离；
没有跨版本复用首轮通过的密度。10个背景重新求解与首轮一致。

#### 22.24.1 完整生产入口运行结果

`freezeout_20260907_v2`正常退出，manifest状态为
`complete_research_curve_accepted`：配置网格3–200 GeV的10/10个能量点、
40/40个通道通过，0个失败点。此处“完整”指完整遍历上述配置网格，
不是连续区间内所有背景的认证或精确峰位定位。

| sqrt(sNN) / GeV | K+/pi+ | K-/pi- |
|---:|---:|---:|
| 3 | 0.505861 | 0.011258 |
| 5 | 0.536345 | 0.082180 |
| 7.7 | 0.543734 | 0.140961 |
| 11.5 | 0.522328 | 0.184603 |
| 19.6 | 0.473653 | 0.232814 |
| 27 | 0.445055 | 0.257074 |
| 39 | 0.417779 | 0.280178 |
| 62.4 | 0.392319 | 0.302463 |
| 130 | 0.368773 | 0.324100 |
| 200 | 0.360928 | 0.331578 |

正电比值在采样点中以7.7 GeV最大，向3 GeV仅缓慢回落；不能据此声称复现
实验低能快速下降或已经定量匹配实验。负电比值在当前网格随能量单调增加。
量的定义仍是固定quark-only BQS背景上的GBU部分产额，非最终强子产额。

共4160条壳层记录（含低/高阶比较），所有通道最终16阶通过。独立汇总的
最大q加阶相对差为0.00256011（0.256011%），分账差/总数为0.00256011；
原核绝对差2.95560e-7 fm^-2，独立径向差1.01134e-14 fm^-2，
eta序列末端最大相对差1.24494e-5，内层加阶差3.45757e-7。
最后q尾带占比最大1.72471e-12，omega条件界占比4.84542e-10。
这些分别是数值收敛与有相位条件的尾界，不是全频解析证明。

独立读取CSV核对10行ratio、40行density，两个比值均与对应通道相除完全一致。
运行结束、任何收尾修改之前，419份现场源码/config与source_snapshot及登记
哈希全部一致；146个输出哈希和52个校验侧文件全部匹配。

| 文件（均位于上述v2目录） | SHA256 |
|---|---|
| manifest.json | 2c180374b66ce5a0563a44be1f697341cdbb90910bc700e1be0ce6b926be2078 |
| ratios.csv | f3855cbf02c2bb341ee041fac37ef93463af878ee39e02368840a213450d416f |
| freezeout_ratios.png | 86068ffb82dc1b7c797e84c857f6c0ad320c33199ad765093b14d57fc3977b0e |

入口同时导出PNG/PDF；绘图脚本和输入已核对，失败不以零代替。当前会话图像
工具返回不支持图像输入，不能声称已完成像素级视觉审核，图面仍待作者确认。
不覆盖原图，不把数据验收等同于论文图表冻结。

#### 22.24.2 收尾代码与授权边界

完整运行后修复一项仅涉及续算manifest的缺陷：旧`manifest.json.sha256`
不能进入新manifest的output_hashes，否则重写manifest会使其自引用校验过期。
新helper同时排除manifest与该侧文件，并加4项纯合成测试。此修复不改数值核、
配置或上述产物；v2的首次输出不存在该自引用，其原快照及哈希保持不变。
随后更新任务台账。当前源码身份因此与v2快照有意不同，不能将其伪装成同一
源码运行或绕过resume身份检查；数值复现以v2保存的完整source_snapshot为准。

运行后复核profile/yield/qgate/入口unit301/301，入口integration11/11。
此前旧regime与旧冻结线GBU路径回归389/389；未跑全仓回归，不改baseline。
旧冻结线CLI与两份PNJLCore无diff，暂存区为空。PR310仍OPEN，未提交、推送或
合并；新的研究生产入口可显式使用，但旧默认、全反馈驻点与论文production
晋升未获本次自动授权。后续应审核图面和方法适用边界，不再手动拼接诊断链路。

收尾manifest排除测试与入口/台账合同207/207通过（含入口22项、台账185项，
不与上述301项重复累加）。models-entry、287脚本入口、relaxtime分类、
data-output、docs consistency、formula-route、active-doc、task-ledger preflight
及git diff --check均通过。最终tracked dirty24、untracked116，既有改动保留。
v2快照与当前现场差异仅为研究workflow的manifest修复和任务台账更新；
419份原快照与所有输出仍保持原哈希。PNG尺寸1100x680，PDF文件21345字节。

### 22.25 作者图像审核与依赖PR收尾

作者已审核v2图像，接受其与之前有限q直算趋势基本吻合，并授权按依赖顺序
收尾OPEN PR，以最新研究生产结果更新过时内容。因此第22.24节“待作者确认”
属于已完成的历史待办。该审核不将有限Lth正常分支的一致性扩大成全谱伪根
误差定理，也不授权替换旧默认或加入介子反馈。

本次required_follow_up仅在D:/w/jrt-ord执行。PR顺序为307、308、309、310；
最终物理方法和数值证据以第22.23–22.24节及ChargedGBU_InfiniteThermal为准。
不删除既有失败/诊断数据，不暂存data/outputs、tmp或work，不重写已审核产物。
前置PR的旧CI与审查意见须先处理，不能借后续PR的通过状态跳过验收。
