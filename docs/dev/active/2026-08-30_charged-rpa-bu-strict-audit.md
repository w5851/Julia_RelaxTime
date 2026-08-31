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
`domega/pi`。束缚态的 `pi` 跳变于是返回一个稳定玻色子；当前代码的
`domega/(2pi)` 会使绝对密度少一半。该共同因子在同口径 `K/pi` 比值中抵消，
所以旧 ratio 诊断不因此失效，但四算法绝对密度比较尚未闭合。

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
- [ ] 将 phase-shift 单电荷实现迁移到 `domega/pi`，并新增窄束缚态与 stable/BW
  绝对密度回归；确认旧 `K/pi` ratio 只发生共同因子抵消。
- [ ] 实现高能 phase anchor、束缚态计数和 Levinson/Mott gate，并通过 q/omega
  节点、`eta`、`omega_max` 和 cutoff 收敛门禁。

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

本任务单只记录可验证的实现边界和后续 gate，不把任何尚未完成的 charged-RPA/BU
公式或冻结线结果标记为 production。
