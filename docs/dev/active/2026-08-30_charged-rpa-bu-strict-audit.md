# Charged RPA/BU 严格复核与后续实现任务单

更新日期：2026-08-30

当前状态：in progress。本文承接完整 KMT interaction-kernel 与固定 BQS
quark-only 后处理诊断，目标是重新核对 charged `π^±/K^±` 的 RPA/BU 数值链，
但不在本任务单完成前替换现有 `MesonDensity` 生产语义。`x_min_cut`、低节点数和
未验证的归一化结果仍只能标记为 diagnostic。

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

后一个式子是 `MesonRPA` 中 `(0,3,8)` 基底、矩阵归一化下的公式；它不能在没有
重新定义 `Pi_a` 和场归一化的情况下直接替换前一个式子。charged 严格实现至少
要比较两个候选标量形式

```math
D_a^{(A)}=\frac{2K_a}{1-4K_a\Pi_a},\qquad
D_a^{(B)}=\frac{2K_a}{1-2K_a\Pi_a},
```

并用同位旋对称极限、真空介子质量/Goldstone 条件和已固定的旧传播子定义确定
唯一可接受的因子。仅改变分母因子就可能移动极点、Mott 阈值和相移，不能把这种
差异归因于 `K45-K67`。

### 3.3 极化函数：`Pi_{us}` 与 `Pi_{su}` 不是可自动合并的同一个数组

`PolarizationAniso` 当前通过 `B0(lambda)`、`A1+A2` 和可选的有限宽度修正构造
泡，其中 `lambda=k0+mu1-mu2`，并对 `num_s_quark=1` 做额外的 `k0` 对称平均。
严格 charged 路线必须重新确认：

- `K^+` 的 `(u,s)` 与 `K^-` 的 `(s,u)` 是否按 retarded 约定分别取值；
- 在 `mu_u=mu_d=mu_s=0` 时是否满足电荷共轭关系和相同的真空极点；
- 有限 `mu` 时的 Landau 区、阈值和虚部符号是否与相移分支一致；
- `A_f` 与 `B0` 是否使用相同 cutoff、热积分上限和节点。按当前公式的自然单位，
  `A_f` 与 `Pi_a` 的量纲应为 `fm^-2`；`K_a` 与传播子 `D_a` 为 `fm^2`，因此
  `K_a Pi_a` 无量纲。现有部分旧 API 页面仍把 `Pi` 写成 `fm^2`，这与公式页和
  `PolarizationAniso` 的实际量纲不一致，列为后续文档修复项，不在本轮静默改写。

当前 `phi` 原生适配层只消除了“由 `A_f` 再构造 `H_f`”这一重复步骤；它不等于
已经完成 `Pi_a` 的严格 retarded 解析延拓。`A_f` 仍可能作为夸克泡的 tadpole
输入，这是两个不同的数值角色。

### 3.4 极点与宽度：先固定定义，再谈 BU 相移

每个 charged 通道应由

```math
\Delta_a(z,q)=0,\qquad z=\omega+i\Gamma/2
```

定义复极点（具体半平面符号必须与 retarded 传播子核对）。束缚态、连续谱和
Mott 解离要分别记录；不能用一个负质量或仅靠有限 `eta` 的相位峰替代极点
定义。极点求解至少需要：残差、收敛状态、分支 seed、阈值
`m_1+m_2`、以及 `K^+`/`K^-` 的 flavor 顺序。

### 3.5 相移与 BU 权重：常数相位、分支和边界都会改变密度

严格路线应从 retarded `D_a^R` 定义

```math
\delta_a(\omega,q)=\operatorname{Im}\log D_a^R(\omega,q)
```

并明确减去高能端常数、处理束缚态的 `pi` 跳变、连续谱相位和 Levinson 边界。
当前 `MesonDensity` 支持 `arg_propagator`/`arg_inverse_propagator`、unwrap 和
`current`/`gbu_reference` 两种权重，但这些是可切换诊断约定，不应在没有公式
固定时混合使用。BU 原式及分部积分形式见
[MesonDensity_BU相移公式.md](../../reference/formula/relaxtime/meson_density/MesonDensity_BU相移公式.md)。

特别是分部积分后的密度使用 `delta_a` 本体；若 `delta_a(omega->infty)` 或
低能边界没有正确归一化，会产生与物理相互作用无关的伪密度。因此相位边界和
Levinson 检查必须先于冻结线扫描。

### 3.6 Bose 支撑与凝聚：`x_min_cut` 不是严格解

当 `omega<=mu_M` 时，玻色分布在 `omega=mu_M` 有极点。当前
`density_policy=:x_min_cut` 把下界移动到 `mu_M+x_min*T`，这是明确的诊断截断，
不是对零动量凝聚模的处理。严格路线需要二选一并写入公式合同：

1. 证明所有生产点处于 `mu_M<m_M` 的正常相，并对 `omega>mu_M` 做支撑门禁；或
2. 单独引入凝聚零模、守恒荷和连续激发的分解。

在此之前，节点增加导致的变化不能被解释为物理的 `K/pi` 趋势。

## 4. 预期会造成差异的来源（归因矩阵）

| 变更 | 是否属于 `K45-K67` | 预期影响 | 当前状态 |
|---|---:|---:|---|
| charged 耦合 `K67 -> K45` | 是 | 当前背景下很小；强 `phi_u-phi_d` 时可放大 | 已有 diagnostic A/B |
| charged 分母 `4KPi` vs `2KPi` | 否 | 可移动极点、相移和密度，可能显著 | 未决，必须推导 |
| `Pi_{us}`/`Pi_{su}` 解析延拓与宽度 | 否 | 直接改变 `K^+` 与 `K^-` 的差异和阈值 | 未决，必须审计 |
| `A_f/B0` cutoff、热上限、节点 | 否 | 数值漂移或伪峰 | 需收敛门禁 |
| 相位分支、常数边界、Levinson | 否 | 可造成整体密度偏置 | 需单独测试 |
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
- [x] 记录 `2/4` 归一化不能机械互换的原因。
- [x] 记录 `phi` 原生旧耦合适配层与 `A_f` 在泡计算中的不同角色。

### Phase B：charged RPA kernel backend（不改旧接口）

- [x] 新增独立的 charged scalar-RPA kernel/spec，显式记录 `pair`、P/S、
  numerator、denominator、retarded convention 和 normalization source。
- [x] 用同位旋对称极限把该 backend 与旧 `MesonPropagator` 做代数/极点 parity
  测试；不以一轮冻结线 A/B 代替该测试。
- [x] 为 `K^+`、`K^-` 分别保留 `Pi_{us}`、`Pi_{su}` 的输入和诊断元数据。

Phase B 的第一项已由 `src/relaxtime/ChargedRPAKernel.jl` 落地。该模块只提供
经过校验的通道/归一化契约和纯标量代数；默认 `MesonPropagator`、极化 provider、
相移和 BU 密度语义保持不变。其两个分母候选分别标记为
`:legacy_scalar_diagnostic` 与 `:neutral_matrix_candidate`，尚未选择生产定义。
同位旋对称 parity 和有序 `Pi_{us}/Pi_{su}` 输入已由
`tests/unit/relaxtime/test_charged_rpa_kernel.jl` 锁定；这仍不代表 Phase C 的
retarded bubble provider 已实现。

### Phase C：charged bubble provider 与极点

- [x] 从现有 `PolarizationAniso` 提取可替换的 charged bubble provider 契约，复用
  相同 `A/B0` 正则化；实现位于 `src/relaxtime/ChargedRPAProvider.jl`。
- [ ] 完成有限宽度严格 retarded 解析延拓验证；当前 `:finite_width` 仍只是实轴
  兼容 proxy，不能作为该门禁的通过证据。
- [ ] 实现极点残差/阈值/Mott 分支记录，覆盖 `mu=0` 共轭关系和有限同位旋点。

### Phase D：strict BU density

- [ ] 固定相移边界、Levinson 检查和 generalized-BU 权重的唯一生产定义。
- [ ] 把正常相 Bose 支撑与凝聚零模明确分支；`x_min_cut` 只能留作 diagnostic。
- [ ] 通过 q/omega 节点、`eta`、`omega_max` 和 cutoff 收敛门禁。

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
4. 数值：所有相移/极点残差、Bose status、节点配置和 fallback 都可追溯。
5. 物理：高能相位边界、Levinson 计数和 Mott 连续性通过后，才可讨论实验趋势。

### 明确非目标

- 不修改 `PNJLCore` 的当前对角平均场求解，不引入非对角夸克自能。
- 不把 `K03/K38` 注入 charged `K/pi`，除非未来先改变平均场/混合基底并完成新
  的二次作用量推导。
- 不实现完整 `Omega_M` 反馈、全 hadronic BQS 或显式 `mu_I` 路线。
- 不提交已有低节点 diagnostic CSV，不更新 production/regression baseline。

## 7. 当前证据与引用路径

- 纯代数核：`src/relaxtime/MesonInteractionKernel.jl`、
  `tests/unit/relaxtime/test_meson_interaction_kernel.jl`。
- charged scalar-RPA 契约与 A/B0 provider 适配：
  `src/relaxtime/ChargedRPAKernel.jl`、`src/relaxtime/ChargedRPAProvider.jl`、
  `tests/unit/relaxtime/test_charged_rpa_kernel.jl`、
  `tests/unit/relaxtime/test_charged_rpa_provider.jl`。
- 中性矩阵 RPA：`src/relaxtime/MesonRPA.jl`、
  `docs/api/relaxtime/propagator/MesonRPA.md`。
- 当前 charged BU scalar 诊断：`src/relaxtime/MesonDensity.jl`、
  `docs/api/relaxtime/meson_density/MesonDensity.md`。
- BU 公式与相位边界：[MesonDensity_BU相移公式.md](../../reference/formula/relaxtime/meson_density/MesonDensity_BU相移公式.md)。
- 本轮 full/legacy 冻结线诊断和 `0.0233%` 差异：
  `docs/dev/active/2026-08-29_full-kmt-phase4-phase6.md`，原始 CSV 仅本地保留。

本任务单只记录可验证的实现边界和后续 gate，不把任何尚未完成的 charged-RPA/BU
公式或冻结线结果标记为 production。
