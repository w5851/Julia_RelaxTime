# Charged-RPA/BU 公式路线闭合包（试点）

route_id: charged_rpa_bu_quark_only
状态：candidate（未授权 production）
初始基线：`origin/main` @ `bc9b2990bcfe3b8c32d2ec0f00066b52b4cf800b`
更新日期：2026-08-30

本文件不是“已经完成的完整 charged-RPA/BU 数值实现”，而是为 PR290 建立的、
可交给独立审阅者复核的公式规范。它把微观 PNJL/KMT 模型、固定 BQS
quark-only 平衡、带电介子 RPA/BU 以及输出比值连接起来，并明确哪些层仍未
完成。任何未决项都必须在升格 production 前单独关闭。

## 1. 范围与计算目标

目标观测量是 charge-resolved

```text
K^+/pi^+  和  K^-/pi^-
```

计算背景固定为：给定 `(T, mu_B)`，先求 `rho_Q/rho_B=0.4`、`rho_S=0` 的
quark-only PNJL 平衡；然后在该背景上构造介子传播子和 BU 数密度。

本路线不包含：

- 介子压强或介子守恒荷反馈到 `Omega_PNJL` 驻点方程；
- 非对角夸克平均场、自能或带电-中性介子混合；
- `K^0`/`bar K^0` 的生产输出；
- 显式 `mu_I` 路线或 pion condensate；
- strict-support 凝聚零模、复平面极点求解和已认证的数值收敛。

因此它是“完整公式路线的 candidate 规范 + 当前 diagnostic 后端的审查边界”，
不是新的稳定 `Models` 入口。

### 1.1 四类可调用的介子数密度算法

本项目把介子数密度算法保持为解耦的四类入口；本路线不把它们合并成单一
“BU 数密度”实现：

| 算法族 | 当前入口/语义 | 本路线边界 |
|---|---|---|
| `stable_particle_limit` | 稳定粒子极限 | 可调用；用于稳定极限基线 |
| `reduced_strict_bw` | reduced strict-BW 单通道积分 | 可调用；用于有限宽度过渡诊断 |
| `q_pole_strict_bw` | `q` 依赖复极点 strict-BW | 可调用；极点状态必须保留 |
| `phase_shift_bu` | 相移双积分 BU | 可调用；`current` 与 `gbu_reference` 均保留 |

四类算法共享上游参数和输出契约，但不共享未经证明的归一化。最终比较的默认
方案固定为 `phase_shift_gbu_reference`（广义 BU）；`phase_shift_current` 仍可
被显式调用作为对照，不能因为比较默认而删除或改写其它三类算法。这个选择只
规定比较口径，不授予 strict charged-RPA/BU production 资格。

## 2. 公式闭合链

本路线必须按下列顺序解释，不能跳过中间层：

```text
PNJL/KMT 微观拉氏量
  -> Hartree/mean-field: phi_f, M_f, Phi, PhiBar
  -> 固定 BQS 守恒荷的 quark-only 平衡
  -> KMT Hartree 收缩后的 K_ab
  -> 有序夸克泡 Pi_ff'（同一正则化和解析延拓）
  -> charged RPA 分母、极点和 retarded 相移
  -> BU 统计权重和介子密度
  -> K^+/pi^+、K^-/pi^- 及完整 provenance
```

这条链中每一层的变量、单位、近似和边界都必须能追到下一节的公式或代码/测试
映射；“公式看起来相同”不能代替归一化和符号转换。

## 3. 微观模型与约定

### 3.1 PNJL/KMT 起点

在项目的自然单位制中，三味 PNJL 的作用量可写成

```math
\begin{aligned}
\mathcal L={}&\bar q(i\gamma^\mu D_\mu-\hat m_0)q
 +G\sum_{a=0}^{8}\left[(\bar q\lambda^a q)^2
 +(\bar q i\gamma_5\lambda^a q)^2\right]\\
&-K\left[\det\bar q(1+\gamma_5)q
 +\det\bar q(1-\gamma_5)q\right]
 -\mathcal U(\Phi,\bar\Phi;T).
\end{aligned}
```

`D_mu` 只表示 PNJL 背景中的规范协变导数；若退回 NJL，则取普通导数。`lambda_0`
采用 `sqrt(2/3) * I_3`。`G` 是四夸克耦合，`K` 是 KMT 六夸克耦合；内部单位
分别为 `[G]=fm^2`、`[K]=fm^5`。

KMT 行列式的整体正负号必须和配置及质量方程一起读取。若外部论文把 KMT 项的
符号吸收到 `K` 中，路线包必须记录转换，不能仅比较 `K` 的数值。

### 3.2 平均场凝聚与质量

本路线把

```math
\phi_f\equiv\langle\bar q_f q_f\rangle,\qquad f=u,d,s
```

作为基本凝聚变量；在当前热积分约定下 `phi_f` 通常为负。平均场质量方程是

```math
\begin{aligned}
M_u&=m_{0u}-4G\phi_u+2K\phi_d\phi_s,\\
M_d&=m_{0d}-4G\phi_d+2K\phi_u\phi_s,\\
M_s&=m_{0s}-4G\phi_s+2K\phi_u\phi_d.
\end{aligned}
```

与之相容的巨热力学势（省略已由项目 `PNJLCore` 展开的夸克对数项）为

```math
\Omega_{\mathrm{PNJL}}
=2G(\phi_u^2+\phi_d^2+\phi_s^2)
 -4K\phi_u\phi_d\phi_s
 +\mathcal U(\Phi,\bar\Phi;T)+\Omega_q.
```

上游平衡只对 `Omega_PNJL` 求驻点，不含 `Omega_M`。因此这里的
“RPA/BU”是固定平均场背景上的后处理，不是介子热力学反馈。

### 3.3 BQS 化学势与守恒荷

项目采用正的奇异度化学势 `mu_S`；奇异夸克的奇异度为 `S_s=-1`，故

```math
\begin{aligned}
\mu_u&=\frac{\mu_B}{3}+\frac{2\mu_Q}{3},\\
\mu_d&=\frac{\mu_B}{3}-\frac{\mu_Q}{3},\\
\mu_s&=\frac{\mu_B}{3}-\frac{\mu_Q}{3}-\mu_S.
\end{aligned}
```

三味夸克密度 `rho_u,d,s` 到 B/Q/S 的映射是

```math
\rho_B=\frac{\rho_u+\rho_d+\rho_s}{3},\qquad
\rho_Q=\frac{2\rho_u-\rho_d-\rho_s}{3},\qquad
\rho_S=-\rho_s.
```

本路线的平衡约束为 `rho_Q/rho_B=0.4` 和 `rho_S=0`。最后一个负号来自夸克的
奇异度量子数，而不是数密度定义的正负任意选择。若将介子密度加入总守恒荷，
必须另建含 `Omega_M`/hadronic charge 的路线，不能在本 candidate 中隐式加入。

在当前 quark-only 约定下，`rho_S=0` 等价于 `rho_s=0`。只要
`rho_B != 0`，无除法的 BQS 残差可以写成

```math
r_Q = \rho_Q-0.4\rho_B=0,
\qquad
r_S = \rho_s=0.
```

在 `rho_s=0` 的子域内，`rho_Q/rho_B=0.4` 进一步给出
`rho_u/rho_d=7/8=0.875`。当 `rho_B=0` 时，比值语义本身未定义；项目 solver
实际求解的是 affine residual `r_Q`，因而在零密度点具有连续代数延拓，但该点
不得被解释为已经测得或验证了 `rho_Q/rho_B=0.4`。

## 4. KMT Hartree 收缩与完整相互作用核

令 `epsilon_P=+1`、`epsilon_S=-1`，并把
`phi=(phi_u,phi_d,phi_s)` 直接代入 KMT 的一次收缩结果，则

```math
\begin{aligned}
K_0^{(X)} &= G+\epsilon_X\frac{K}{3}(\phi_u+\phi_d+\phi_s),\\
K_{12}^{(X)}=K_3^{(X)} &=G-\epsilon_X\frac{K}{2}\phi_s,\\
K_{45}^{(X)} &=G-\epsilon_X\frac{K}{2}\phi_d,\\
K_{67}^{(X)} &=G-\epsilon_X\frac{K}{2}\phi_u,\\
K_8^{(X)} &=G-\epsilon_X\frac{K}{6}(2\phi_u+2\phi_d-\phi_s),\\
K_{03}^{(X)}=K_{30}^{(X)}
 &=-\epsilon_X\frac{K}{2\sqrt6}(\phi_u-\phi_d),\\
K_{08}^{(X)}=K_{80}^{(X)}
 &=-\epsilon_X\frac{\sqrt2K}{12}(\phi_u+\phi_d-2\phi_s),\\
K_{38}^{(X)}=K_{83}^{(X)}
 &=\epsilon_X\frac{K}{2\sqrt3}(\phi_u-\phi_d).
\end{aligned}
```

其中 `X=P` 是赝标量、`X=S` 是标量；所有 `K_ab` 的单位为 `fm^2`。

### 4.1 味道通道映射

自旋零介子双线性按下表进入 charged 或 neutral 通道：

| 物理通道 | 夸克组成（示意） | 相互作用核 |
|---|---|---|
| `pi+` | `u bar(d)` | `K12` |
| `pi-` | `d bar(u)` | `K12` |
| `K+` | `u bar(s)` | `K45`（d spectator） |
| `K-` | `s bar(u)` | `K45` 的共轭有序通道 |
| `K0` | `d bar(s)` | `K67`（u spectator） |
| `bar K0` | `s bar(d)` | `K67` 的共轭有序通道 |

因此 `K^\pm` 不能在一般 `phi_u != phi_d` 时使用旧合并字段 `K4567`。
旧 helper 的纯代数关系为

```math
K_{123}^{\mathrm{legacy},\pm}=K_{12}^{P/S},\qquad
K_{4567}^{\mathrm{legacy},\pm}=K_{67}^{P/S}.
```

只有在 `phi_u=phi_d` 时才有 `K45=K67=K4567`。`K03/K38` 是中性 `(0,3,8)` 基底
的混合元；在当前对角平均场、无 charged-neutral 混合的路线里，它们不直接进入
`K^\pm` 或 `pi^\pm` 的标量分母。

## 5. 有序夸克泡与极化函数

### 5.1 Matsubara 起点

对组成味道 `f,f'` 和通道顶角 `Gamma_X`，先定义有序泡

```math
\Pi_{ff'}^X(i\nu_m,\mathbf q)
 =T\sum_n\int_{|\mathbf p|<\Lambda}\frac{d^3p}{(2\pi)^3}
 \operatorname{Tr}_{D,c}
 \left[\Gamma_X S_f(i\omega_n,\mathbf p)
 \Gamma_X S_{f'}(i\omega_n-i\nu_m,\mathbf p-\mathbf q)\right].
```

整体符号由下式的项目约定固定。`K^+` 使用 `(u,s)`，`K^-` 使用 `(s,u)`；在
有限化学势下两者不应自动合并为同一个缓存数组。

将每条夸克线写成正/负能量谱投影，`a,b in {+1,-1}`，Matsubara 求和后的
有序 retarded 泡具有通用结构

```math
\Pi_{ff'}^{R,X}(\omega,\mathbf q)
=\int_{|\mathbf p|<\Lambda}\frac{d^3p}{(2\pi)^3}
\sum_{a,b=\pm1}\mathcal C_{ab}^{X}(\mathbf p,\mathbf q)
\frac{n_F(bE_{f'}-\mu_{f'})-n_F(aE_f-\mu_f)}
{\omega+\mu_f-\mu_{f'}-aE_f+bE_{f'}+i0^+}.
```

`C_ab^X` 包含 Dirac/color trace、P/S 顶角和本项目整体 loop 符号；PNJL 时
`n_F` 替换为同一 `Phi,Phibar` 背景下的广义分布。这个式子直接固定了外能量
只通过 `omega+mu_f-mu_f'` 进入。因此

```text
Pi_us:  omega + mu_u - mu_s,  masses/order=(u,s)
Pi_su:  omega + mu_s - mu_u,  masses/order=(s,u)
```

二者在有限化学势下不是同一个函数。严格 provider 必须用相同 cutoff 和变量
变换验证交换不变量

```math
\Pi_{ff'}^R(\omega,q;\mu_f,\mu_{f'})
\stackrel{?}{=}
\left[\Pi_{f'f}^R(-\omega,q;\mu_{f'},\mu_f)\right]^*,
```

并在 `mu_f=mu_f'=0` 恢复正能量 charged-conjugate 极点/相移一致性；问号表示
这是后续实现必须通过的验证关系，而不是当前 legacy B0 已经证明的性质。

### 5.2 `A_f`/`B0` 组合

本路线采用 `e^{-i omega t}` Fourier 约定；retarded 能量从上半平面
`z=omega+i 0^+` 延拓，有限 `eta` 只是实轴数值探针。项目目前采用的代数形式为

```math
\Pi_{ff'}^{P,S}(z,q)
=-\frac{N_c}{8\pi^2}\left\{
A_f+A_{f'}+\left[(m_f\mp m_{f'})^2
-(z+\mu_f-\mu_{f'})^2+q^2\right]B_0^{ff'}(z,q)\right\}.
```

上式中的 `upper P` 取 `m_f-m_f'`，`lower S` 取 `m_f+m_f'`。量纲为

```text
[A_f] = fm^-2,  [B0] = 1,  [Pi_ff'] = fm^-2,
[K_ab] = fm^2,  [K_ab Pi_ff'] = 1.
```

一个与当前凝聚符号相容的 tadpole 关系是

```math
\phi_f=\frac{N_c}{4\pi^2}m_f A_f,
\qquad H_f\equiv-\phi_f
=-\frac{N_c}{4\pi^2}m_f A_f.
```

所以旧 `calculate_G_from_A` 返回的是历史 helper `H_f`，不是另一个物理凝聚。
新路线应优先复用平衡解中已经得到的 `phi_f`；`A_f` 仍可作为夸克泡本身的积分
输入，二者不能因为命名相似而混为同一计算角色。

### 5.3 `num_s_quark=1` 的边界

Rehberg 等的 Eq. (2.22)–(2.23) 讨论了三动量截断下非等质量 kaon 泡在变量
平移后出现的 `p_0` 对称性伪破缺，并给出同一有序表达在 `p_0` 与 `-p_0` 的平均
形式。这个平均是针对其散射计算正则化伪影的处方，不自动等价于有限 `mu` 下的
严格 retarded 解析延拓。

因此本路线把两种模式分开：

- `ordered_retarded`：生产候选目标，分别保留 `Pi_us` 与 `Pi_su`，须完成实轴/复
  平面延拓和共轭检查；
- `legacy_symmetrized_B0`：兼容诊断，对应当前 `num_s_quark=1` 的显式平均，不能
  在没有额外证据时作为严格生产定义。

当前代码的 legacy 平均具体是

```math
B_0^{legacy}(k_0;f,f')=\frac12\left[
B_0(k_0+\mu_f-\mu_{f'};m_f,m_{f'})+
B_0(-k_0+\mu_f-\mu_{f'};m_f,m_{f'})\right].
```

它既没有交换 `m_f,m_f'`，也没有把化学势差改成
`mu_f'-mu_f`，所以不能等同于 `Pi_f'f`；这正是 `num_s_quark=1` 只能保留为
legacy diagnostic 的代数边界。

## 6. Charged RPA 传播子、极点与相移

### 6.1 标量 charged 分母的项目规范

在上述同一 `Pi_ff'` 定义下，本路线选择 Rehberg 的标量形式

```math
D_{ff'}^X(z,q)=\frac{2K_{ff'}^X}
 {1-4K_{ff'}^X\Pi_{ff'}^X(z,q)},
\qquad
F_{ff'}^X(z,q)=1-4K_{ff'}^X\Pi_{ff'}^X(z,q).
```

它与项目旧 `MesonPropagator`/`MesonDensity` 的 charged 标量接口连续；
`K^\pm` 用 `K45^P`，`pi^\pm` 用 `K12^P`。这里的 `2`/`4` 不是可脱离泡归一化
单独比较的物理常数。

另一个常见矩阵写法是

```math
\mathcal M=2\mathcal K\,[I-2\mathcal K\Pi]^{-1}.
```

它用于 `(0,3,8)` 基底的矩阵核；只有在同时给出场归一化和 `Pi` 的因子转换时，
才可与上面的 scalar 形式比较。当前 candidate 明确禁止把矩阵的 `2` 机械替换
charged 分母的 `4`。

单通道的代数桥接条件可以明确写成

```math
\Pi_{\mathrm{matrix}}=2\Pi_{\mathrm{charged}}
\quad\Longrightarrow\quad
2K\,[1-2K\Pi_{\mathrm{matrix}}]^{-1}
 =\frac{2K}{1-4K\Pi_{\mathrm{charged}}}.
```

因此 `2`/`4` 的差异等价于泡或介子场归一化相差一个明确的因子 2；只有从同一
二次作用量证明 `Pi_matrix=2 Pi_charged`，该桥接才可用于严格实现。当前项目
先把这个等式作为 normalization gate 和单通道回归测试目标，不把它当作已经
完成的 charged-RPA 证明。

### 6.2 极点定义

retarded 传播子的共振极点定义为

```math
F_{ff'}^X(z_p,q)=0,
\qquad z_p=M_{ff'}(q)-i\,\Gamma_{ff'}(q)/2,
```

其中 `Im(z_p)<0` 是 retarded 约定。每个极点记录 `z_p`、残差、求根状态、
flavor 顺序、连续谱阈值 `m_f+m_f'` 和 Mott 分支。当前 `PolarizationAniso` 的
公共输入仍主要是实数 `k0`，所以这一步尚未实现为 production pole solver；有限
`eta` 的实轴相位峰不能冒充复极点证据。

### 6.3 相移规范

在实轴上先构造 `D^R(omega,q)=D(omega+i eta,q)`，定义

```math
\delta_{ff'}(\omega,q)
 =\operatorname{unwrapArg}D_{ff'}^R(\omega,q)
 -\lim_{\omega\to\omega_{\max}}\operatorname{unwrapArg}D_{ff'}^R(\omega,q).
```

等价地，在远离分支常数时可以使用 `-Arg(F^R)`。必须显式记录 unwrap、低能
阈值、高能归一化和束缚态 `pi` 跳变；不能把 `arg(propagator)`、
`arg(inverse propagator)` 与任意常数相位混用。

## 7. BU 数密度与带电化学势

### 7.1 介子化学势

按“夸克减去反夸克”的有序组成，

```math
\mu_{\pi^+}=\mu_u-\mu_d,\quad
\mu_{\pi^-}=\mu_d-\mu_u,
\quad
\mu_{K^+}=\mu_u-\mu_s,\quad
\mu_{K^-}=\mu_s-\mu_u.
```

单个 charged spin-0 物种的简并度取 `d=1`；若改用 `pi` 或 `K` 聚合通道，必须
显式使用其聚合简并度，不能与单个 charged ratio 混用。

### 7.2 BU 原式

在相移边界已经固定且正常相支撑成立时，

```math
n_M=d_M\int\frac{d^3q}{(2\pi)^3}
\int_{\omega_{\min}}^{\omega_{\max}}
\frac{d\omega}{2\pi}
g_B(\omega;\mu_M,T)\,\frac{\partial\delta_M(\omega,q)}{\partial\omega},
```

其中

```math
g_B(\omega;\mu_M,T)=\frac{1}{\exp[(\omega-\mu_M)/T]-1}.
```

若低能/高能边界项按零处理，分部积分形式为

```math
n_M=\frac{d_M}{T}\int\frac{dq\,q^2}{2\pi^2}
\int\frac{d\omega}{2\pi}
g_B(\omega;\mu_M,T)[1+g_B(\omega;\mu_M,T)]\,\delta_M(\omega,q).
```

生产实现必须保存采用的是导数形式还是分部积分形式，并以同一 phase branch、
阈值和边界验证两者；当前项目优先使用相移本体的分部积分诊断形式。

### 7.2.1 四类算法与 BU 比较口径

四类顶层数密度算法由各自入口解耦实现：稳定粒子极限、reduced strict-BW、
`q` 依赖复极点 strict-BW，以及 phase-shift BU。它们可以在同一固定 BQS
背景上分别调用和记录，不得用一个算法的数值替代另一个算法的验证。

phase-shift BU 内部保留两种相移权重：

```math
F_{\mathrm{current}}(\delta)=\delta,
\qquad
F_{\mathrm{gBU}}(\delta)=\delta-\frac12\sin(2\delta).
```

两种权重都属于可调用的分析/诊断方案；本路线把 `F_gBU` 注册为最终比较的
默认方案，因为它是当前文献比较中采用的广义 BU 形式。这个默认不改变稳定、
两类 strict-BW 或 `F_current` 的调用能力，也不绕过相位边界、Levinson、
支撑和数值收敛门禁。

### 7.3 Bose 支撑

正常相要求积分支撑满足 `omega>mu_M`。沿当前化学冻结线，本路线首先检查
`mu_M < m_M` 和所有积分节点的 `omega-mu_M > 0`；在该正常相门禁通过时，
不引入介子凝聚零模。经验上的“通常不超过有效质量”只能作为待验证假设，
不能替代逐点检查。

`density_policy=:x_min_cut` 把下界改为 `mu_M+x_min*T`，保留为文献复现和
异常点诊断的简化处理。它没有引入零动量凝聚模，也不能被解释为严格 Bose
凝聚处理。若任何冻结线点满足 `mu_M >= m_M`、节点触及 `omega <= mu_M`，
或正常相门禁失败，则该点标记为 `unsafe_bose_domain`，停止当前 route 的
production 升格，并另建“凝聚零模 + 连续谱”路线；在此之前不得把
`x_min_cut` 的有限结果解释为凝聚后的物理密度。

## 8. 公式 → 代码 → 测试映射

| 公式层 | 当前代码/入口 | 当前证据与边界 |
|---|---|---|
| PNJL 平均场、质量、BQS | `src/models/pnjl_physics/PNJLCore.jl`、`src/models/solver/spec/ConstraintModes.jl`、`src/models/solver/spec/Conditions.jl` | `FixedMuBConservedCharges` 代数/求解测试（含 `tests/unit/models/test_fixed_mub_conserved_charges.jl`）；不含 `Omega_M` |
| KMT `K_ab` | `src/relaxtime/MesonInteractionKernel.jl` | `tests/unit/relaxtime/test_meson_interaction_kernel.jl`；纯代数 backend |
| 中性 `(0,3,8)` 矩阵 | `src/relaxtime/MesonRPA.jl`、`MesonRPAAdapter.jl` | `tests/unit/relaxtime/test_meson_rpa.jl`、`test_meson_rpa_adapter.jl`；不代表 charged 混合 |
| `A_f/B0` 泡 | `src/relaxtime/OneLoopIntegrals.jl`、`OneLoopIntegralsAniso.jl`、`PolarizationAniso.jl` | `tests/unit/relaxtime/test_oneloopintegrals*.jl`、`test_polarization_aniso.jl`；`num_s_quark=1` 仍属兼容语义 |
| 四类介子数密度 | `stable_meson_number_density`、`strict_bw_meson_number_density`、`strict_bw_qpole_meson_number_density`、`phase_shift_meson_number_density` | `tests/unit/relaxtime/test_meson_density.jl` 及对应 workflow tests；四类均可调用，当前仍是有限 `eta`/数值相位诊断 |
| route registry | `config/governance/formula_route_closure.toml`、`scripts/dev/check_formula_route_closure.jl` | `tests/unit/config/test_formula_route_closure.jl`；只检查闭合包完整性 |

当前实现中 `FullKMTInteraction` 可以把 `K12/K45` 注入既有 charged BU 入口，但这
并不等于已经完成本文件第 5、6、7 节的严格 retarded、极点、相位和支撑门禁。

## 9. 外部来源与项目转换

| 来源 | 使用内容 | 与项目规范的转换/限制 |
|---|---|---|
| Rehberg, Klevansky, Hüfner, *Nucl. Phys. A* 608 (1996) 356–388, DOI [10.1016/0375-9474(96)00247-3](https://doi.org/10.1016/0375-9474(96)00247-3) | Eq. (2.1) NJL/KMT 拉氏量；Eq. (2.2)–(2.5) 平均场耦合与 `A`；Eq. (2.7) 质量方程；Eq. (2.8)–(2.14) 传播子/极化；Eq. (2.22)–(2.23) kaon 泡对称化 | 采用项目内部 `fm` 自然单位；把其 `G_f` 记为 `H_f=-phi_f`；Eq. (2.23) 仅作 legacy 诊断，不能直接宣称是有限 `mu` 的 retarded 定义 |
| Rehberg & Klevansky, *Ann. Phys.* 252 (1996) 422–457, DOI [10.1006/aphy.1996.0140](https://doi.org/10.1006/aphy.1996.0140)；开放版本 [hep-ph/9510221](https://arxiv.org/abs/hep-ph/9510221) | arXiv v2 Eq. (19) 的 flavor-order 结构、Eq. (84)–(89) 的一圈积分/解析结构；出版物方程编号需按具体版本复核，不再写成不存在的 Eq. (4.9)/(6.5)–(6.10) | 保留有序 `(f,f')`；项目 `Pi` 的 `fm^-2` 约定与 `K Pi` 无量纲检查必须同时使用 |
| Tian et al., *Phys. Rev. D* 114 (2026) 034012, DOI [10.1103/d7nm-y2vp](https://doi.org/10.1103/d7nm-y2vp) | Eq. (2)–(3) 完整 KMT 有效耦合（含 `K03/K30/K38/K83`）；Eq. (20)–(22)、(26) 中性 `(0,3,8)` RPA 矩阵组织 | 该文使用 NJL、外磁场和 Pauli–Villars；本项目只采用耦合/矩阵结构，把 PNJL/零磁场/BU 数值另行闭合 |
| Blaschke et al., *Phys. Rev. D* 96 (2017) 094008, DOI [10.1103/PhysRevD.96.094008](https://doi.org/10.1103/PhysRevD.96.094008) | PNJL + BU 的 bound/resonant/scattering 统一、Mott 跳变、Levinson/高能相位边界 | 只采用方法结构；具体相位分支、正则化和节点仍按本项目公式包冻结 |
| Blaschke et al., *Particles* 3 (2020) 169–177, DOI [10.3390/particles3010014](https://doi.org/10.3390/particles3010014)；arXiv [1912.13162](https://arxiv.org/abs/1912.13162) | 开放版本 Eq. (18) 标准 BU 部分密度、Eq. (20) `delta -> delta - sin(2delta)/2` 的广义替换；最终比较采用广义形式 | 四类数密度算法均保持可调用；最终路线比较默认 `phase_shift_gbu_reference`，单 charged 物种取 `d=1`；实验趋势只作参照 |

这些来源分别支撑“微观模型/传播子”“有序泡与归一化”“BU 方法结构”和
“`pi/K` 带电比值工作流”。来源之间若采用不同的 `Pi`、相移或 KMT 符号，必须
先执行本文件的转换表和不变量检查，不能拼接公式后再靠数值调参消除差异。

## 10. 选择表与不变量

| 项目 | 本 candidate 选择 | 必须保持的不变量 |
|---|---|---|
| KMT 凝聚输入 | 直接使用平衡 `phi_f`；`H_f` 仅为 legacy helper | `H_f=-phi_f`，`[phi]=[H]=fm^-3` |
| charged 核 | `pi^\pm -> K12`，`K^\pm -> K45` | `phi_u=phi_d` 时 `K45=K67`；任意背景旧 `K4567=K67` |
| 中性混合 | `K03/K38` 只留在 `(0,3,8)` 矩阵 | `phi_u=phi_d` 时 `K03=K38=0` |
| charged 分母 | `D=2K/(1-4KPi)` | `K Pi` 无量纲；不能与矩阵 `1-2KPi` 混用 |
| bubble 顺序 | 目标为 `ordered_retarded` | `mu=0` 时 `K^+/-` 真空极点共轭/相同；有限 `mu` 保留方向 |
| 相移 | retarded `D` 的 unwrap 相位并归一化高能常数 | 束缚态跳变、Levinson 计数和高能 `delta -> 0` |
| BU 支撑 | 正常相 `omega>mu_M` | `x_min_cut` 只能标 diagnostic；接近凝聚必须分支 |
| 数密度算法 | 四类入口均可调用 | 最终比较默认 `phase_shift_gbu_reference`；不得删除其它算法 |
| 上游背景 | quark-only BQS | 不把后处理介子荷写成完整热力学平衡 |

## 11. 生产边界与升格条件

在以下条件全部满足前，route status 保持 `candidate`：

1. `K12/K45/K67` 的符号、单位和味道 spectator 由纯代数测试覆盖；
2. charged `Pi_us/Pi_su` 的正则化、解析延拓和 `num_s_quark` 语义有可复核推导；
3. `2/4` 归一化与至少一个真空或同位旋对称固定点相容；
4. 四类数密度入口均可复现；phase-shift 的 `current`/`gbu_reference`、极点、
   Mott 阈值、相位分支、Levinson 边界和 BU 两种积分形式通过验证；
5. `q/omega` 节点、`eta`、截断和 Bose 支撑的变化低于预先登记的误差预算；
6. 至少一个外部文献固定点按相同参数、单位、简并度和通道定义复核；
7. 独立审阅者确认公式 → 代码 → 测试映射，并明确接受所有近似；
8. 另行决定是否建立 `Omega_M` 反馈或非对角平均场路线，不由本 backend 自动开启。

即使以上公式/数值条件通过，`production_authorized=true` 仍须由人工在 route
registry 和专题 SOP 中显式批准。当前 `production_authorized=false`。

## 12. 未决项与审查问题

供后续独立审查（包括计划中的 5.6Sol 审核）的最小问题清单：

1. 项目 `A_f/B0/Pi` 符号是否能从同一 Matsubara 定义逐行推到第 5 节；
2. `Pi_us`、`Pi_su` 与 Rehberg 开放版本 Eq. (19) 的 flavor-order 交换在有限 `mu` 下如何
   映射到 retarded 共轭，而不是误用同序 `p_0` 平均；
3. charged 标量 `1-4KPi` 与中性矩阵 `I-2KPi` 的场归一化转换是否完整；
4. `z_p=M-iGamma/2` 的半平面、虚部符号和有限 `eta` 相位是否一致；
5. 相移常数、束缚态跳变和 BU 分部积分边界是否在有限窗口内闭合；
6. `x_min_cut` 何时会删去物理支撑、何时必须引入凝聚零模；
7. 当前 BQS quark-only 背景是否足以支撑实验趋势比较，还是需要独立 hadronic/BQS
   路线；
8. 哪些结论可以晋升为 production，哪些必须继续标为 diagnostic。

在这些问题没有书面回答前，不得把本文件的 `candidate` 改为
`production_authorized`，也不得更新正式 baseline。

## 13. 当前诊断结果的解释边界

PR287/既有诊断显示，在当前 BQS 背景上 `K45-K67` 很小，full/legacy 冻结线 A/B
差异约为 `0.0233%` 以内；这只约束“固定当前泡、分母和积分口径时的耦合替换”。
它不能推出严格 charged-RPA/BU 的差异也只有这个量级，因为有序泡、分母归一化、
极点/相位边界、Bose 支撑和节点收敛仍可能改变结果。

同理，`K03/K38` 在当前对角、无中性-带电混合的路线中不直接进入 `K^\pm/\pi^\pm`
分母；未来若改变平均场基底或二次作用量，必须新建 route id 并重新闭合。

## 14. 审阅快照

- 公式包基线：`bc9b2990bcfe3b8c32d2ec0f00066b52b4cf800b`
- route registry：`config/governance/formula_route_closure.toml`
- 治理检查器：`scripts/dev/check_formula_route_closure.jl`
- 纯契约测试：`tests/unit/config/test_formula_route_closure.jl`
- 当前状态：`candidate` / diagnostic-only
- 未提交数值产物：任何低节点 CSV、冻结线 CSV 或临时 benchmark 均不属于本公式包

本文件的目的，是让后续实现可以在一个冻结、可追踪、可逐项质询的规范上继续，
而不是把尚未验证的选择包装成唯一的物理真理。
