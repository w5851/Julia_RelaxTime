# Charged RPA/BU 严格复核与后续实现任务单

更新日期：2026-08-31

当前状态：in progress。本文承接完整 KMT interaction-kernel 与固定 BQS
quark-only 后处理诊断，目标是重新核对 charged `π^±/K^±` 的 RPA/BU 数值链，
但不在本任务单完成前替换现有 `MesonDensity` 生产语义。`x_min_cut`、低节点数、
未闭合的 retarded/phase 数值 gate 仍只能标记为 diagnostic。

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
