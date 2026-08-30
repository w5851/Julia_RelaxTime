# KMT 平均场背景到完整介子相互作用核

状态：Phase 0/1.5/2 的公式合同与诊断实现。本文不是完整 RPA 数值接入，也不授予 production 资格。

## 1. 层次边界

三味 PNJL/NJL 的微观拉氏量可按下列顺序使用：

```text
L_4q + L_KMT
  -> Hartree/mean-field background (phi_u, phi_d, phi_s, M_f, Phi, PhiBar)
  -> contract one quark pair in L_KMT
  -> effective four-fermion interaction K_ab^± on that background
  -> quadratic meson action / RPA polarization matrix
  -> propagator, pole or phase shift
```

因此本文件中的 `K_ab^±` 是“在给定平均场背景上生成 RPA 二次核”的输入，不是把介子热力学势反馈回平均场驻点方程。`MesonInteractionKernel` 只实现中间的纯代数一步；RPA 的夸克泡图和极化函数仍由后续阶段接入。

## 2. 约定

本阶段锁定以下实现约定：

| 项目 | 约定 |
| --- | --- |
| 中性基底 | `(lambda_0, lambda_3, lambda_8)` |
| `lambda_0` | `sqrt(2/3) * I_3` |
| 凝聚 | `phi_f = <bar q_f q_f>`，直接作为公式中的 `sigma_f` 使用 |
| 四夸克耦合 | `G`，单位 `fm^2` |
| KMT 六夸克耦合 | `K`，单位 `fm^5`，保留调用方传入的符号 |
| P 通道 | 记为上标 `+`，对应赝标量；取 `s_channel=+1` |
| S 通道 | 记为上标 `-`，对应标量；取 `s_channel=-1` |

不同论文可能把 `i tr S^f`、`-<bar q_f q_f>` 或 KMT 耦合本身定义成相反符号。因此下面的正负号只在这张表的约定下成立，不能脱离约定逐项比较。

## 3. 来源对齐后的完整耦合公式

令 `s=+1`（P）或 `s=-1`（S），并写 `phi_u=u`、`phi_d=d`、`phi_s=sigma_s`。Phase 1.5 按 Tian 等，Phys. Rev. D 114, 034012 (2026)，Eq. (3) 对齐下列公式。这里采用 `phi_f = sigma_f`，并把上标 `+` 解释为 P（赝标量）、上标 `-` 解释为 S（标量）：

```math
\begin{aligned}
K_0^{(s)} &= G + s\frac{K}{3}(u+d+\sigma_s),\\
K_3^{(s)} &= G - s\frac{K}{2}\sigma_s,\\
K_{45}^{(s)} &= G - s\frac{K}{2}d,\\
K_{67}^{(s)} &= G - s\frac{K}{2}u,\\
K_8^{(s)} &= G - s\frac{K}{6}(2u+2d-\sigma_s),\\
K_{03}^{(s)}=K_{30}^{(s)} &= -s\frac{K}{2\sqrt{6}}(u-d),\\
K_{08}^{(s)}=K_{80}^{(s)} &= -s\frac{\sqrt{2}K}{12}(u+d-2\sigma_s),\\
K_{38}^{(s)}=K_{83}^{(s)} &= s\frac{K}{2\sqrt{3}}(u-d).
\end{aligned}
```

中性 `(0,3,8)` 矩阵按

```math
K^{(s)}_{0,3,8}=
\begin{pmatrix}
K_0 & K_{03} & K_{08}\\
K_{30} & K_3 & K_{38}\\
K_{80} & K_{83} & K_8
\end{pmatrix}
```

组装；带电/味道通道中，`K_12` 对应 `pi^±`，`K_45` 对应 `u\bar s/s\bar u` 的 charged kaon，`K_67` 对应 `d\bar s/s\bar d` 的 `K^0/\bar K^0`。

在 `u=d` 时，`K_03=K_30=K_38=K_83=0`，中性矩阵退化为旧的 `0/8` 块；同时 `K_45=K_67`，恢复旧接口把所有 kaon charged/neutral 味道通道合并成一个 `K_4567` 的极限。

### 3.1 从 KMT 行列式到 charged 二次核的可复核代数步骤

为避免把非对称背景下的通道映射当作文献记号的直接替换，定义

```math
J_{ij}^{R/L}=\bar q_i(1\pm\gamma_5)q_j,
\qquad
\det J^{R/L}=\frac{1}{6}\epsilon_{ijk}\epsilon_{lmn}
J_{il}^{R/L}J_{jm}^{R/L}J_{kn}^{R/L}.
```

在对角平均场中只收缩一个对角双线性，`J_{ff}^{R/L} -> phi_f`，并将剩余两个
双线性投影到 `J_a^P=bar(q)i gamma_5 lambda_a q` 或
`J_a^S=bar(q)lambda_a q`。把一次收缩后的 KMT 项按 `J_a^X J_b^X` 的系数
定义为 `Delta K_ab^(X)`，则两种 determinant orientation 的和以及 P/S 投影给出

```math
\begin{aligned}
S_1^2+S_2^2 &= 4(\bar u d)(\bar d u),\\
P_1^2+P_2^2 &= 4(\bar u i\gamma_5 d)(\bar d i\gamma_5 u),\\
\mathcal L_{KMT}^{(ud;\phi_s)}
 &=2K\phi_s(\bar u d)(\bar d u)
   -2K\phi_s(\bar u i\gamma_5 d)(\bar d i\gamma_5 u).
\end{aligned}
```

这里前两式使用项目的 Gell-Mann 生成元归一化；将最后一式重新写回
`S_1^2+S_2^2`、`P_1^2+P_2^2`，便得到 S/P 通道的 `1/2` 系数。循环交换
`(u,d,s)` 后，

```math
\begin{aligned}
\Delta K_{12}^{(X)}=\Delta K_3^{(X)}&=-\epsilon_X\frac{K}{2}\phi_s,\\
\Delta K_{45}^{(X)}&=-\epsilon_X\frac{K}{2}\phi_d,\\
\Delta K_{67}^{(X)}&=-\epsilon_X\frac{K}{2}\phi_u.
\end{aligned}
```

因此每个 charged flavor pair 的 spectator 是未出现在该 pair 中的味道：
`(12,s)`、`(45,d)`、`(67,u)`。对角投影则给出 `K_0`、`K_8`，而 `u-d` 的
反对称部分只出现在 `K_03/K_38`；这两个量都正比于 `phi_u-phi_d`，故在
`phi_u=phi_d` 时同时消失。将 `G delta_ab` 加回上述 KMT 二次系数，得到第 3 节
列出的完整 `K_ab^(X)`。这个步骤明确说明：`K4567` 的 legacy 修正使用
`phi_u` 时只能对应 `K67`，不能在非对称背景下重解释为 `K45`。

该推导只使用 KMT 行列式、对角平均场和项目的 P/S 投影约定；它不证明
`Pi_ff'` 的 retarded 解析延拓、极点 sheet 或 BU 热力学权重。后者仍需独立的
泡函数和数值验证。

## 4. 与当前代码的兼容关系

`src/models/pnjl_physics/PNJLCore.jl` 当前的标量平均场质量方程已经允许 `phi_u`、`phi_d`、`phi_s` 分开出现，但在实际平衡状态和旧传播子接口中，u/d 介子耦合仍被合并。本路线不改变这些上游方程，只把三味凝聚作为输入生成一个可审计的完整核。

当 `phi_u=phi_d=phi_l` 时，旧 `EffectiveCouplings.calculate_effective_couplings` 使用 `G_u=-phi_l`、`G_s=-phi_s`。在本文件约定下，它的 `K0/K123/K4567/K8/K08` 与新后端的 P/S 对应分量相等；单元测试锁定这一兼容性。这个等式是 API/符号映射的测试，不是对所有有限同位旋物理的证明。

### 4.1 非对称背景下的旧接口映射

旧接口只有一个 `K4567` 字段，且其修正项使用 `H_u=-phi_u`。因此在
`phi_u != phi_d` 时，纯代数对应关系是

```math
\begin{aligned}
K_{123}^{\mathrm{legacy},\pm} &= K_{12}^{P/S},\\
K_{4567}^{\mathrm{legacy},\pm} &= K_{67}^{P/S},\\
K^\pm &= K_{45}^{P/S}.
\end{aligned}
```

这里 `K_{45}` 是 `u\bar s/s\bar u` 的物理 charged-kaon 通道（`K^+`、`K^-`），
而 `K_{67}` 是 `d\bar s/s\bar d` 的中性通道（`K^0`、`\bar K^0`）。只有在
`phi_u=phi_d` 时才有 `K_{45}=K_{67}=K_{4567}^{\mathrm{legacy}}`；旧函数签名
保持不变，但非对称 charged-kaon 计算不得把 `K4567` 直接当成 `K45`。

所有凝聚在这里均为 `phi_f=<bar q_f q_f>`，单位 `fm^-3`；旧 helper `H_f`
仅是 `-phi_f` 的兼容命名，`K H_f` 与完整核中的 `K phi_f` 修正项均为 `fm^2`。

## 5. Phase 3 数值桥接合同

`MesonRPAAdapter` 将当前 `PolarizationAniso` 接到上述矩阵公式，但只采用
三次同味夸克泡作为输入：

```text
Pi_u = Pi_P/S(k0, k; m_u, m_u, mu_u, mu_u, A_u, A_u)
Pi_d = Pi_P/S(k0, k; m_d, m_d, mu_d, mu_d, A_d, A_d)
Pi_s = Pi_P/S(k0, k; m_s, m_s, mu_s, mu_s, A_s, A_s)
```

随后调用方可以把 `(Pi_u, Pi_d, Pi_s)` 传入 `neutral_polarization_matrix`，
或直接使用 `neutral_rpa_from_quark_params` 完成矩阵组合。这里的
`Pi_03`、`Pi_08`、`Pi_38` 是味道对角泡在 `(0,3,8)` 基底下的投影，不代表
adapter 已经实现了非对角夸克传播子或非对角平均场自能。

所有新入口的 `k0_inv_fm`、`k_norm_inv_fm`、`m_f`、`mu_f`、`T` 和 `A_f`
遵循项目内部自然单位；`num_s_quark` 默认逐味为零，只有显式设置为 `1`
时才使用当前 `PolarizationAniso` 中的 `k0` 对称平均。缺少 `A` 时的自动
补值来自 `AFieldBuilder`，不会触发 gap solver。

这一步只验证中性数值接线和输入合同，不证明当前极化的正则化、retarded
延拓或外部数值固定点，也不产生极点、相移或介子热力学反馈。charged 顶角的
纯代数归一化在第 6.1 节独立闭合。相应的诊断 API 边界见
`docs/api/relaxtime/propagator/MesonRPAAdapter.md`。

## 6. 文献证据与未决项

完整 `0/3/8` 结构可追溯到 KMT 六费米项的 Hartree 收缩和三味 RPA 组织。项目调研中使用的主要线索是 Rehberg、Klevansky、Hüfner（Phys. Rev. C 53, 410, DOI: [10.1103/PhysRevC.53.410](https://doi.org/10.1103/PhysRevC.53.410)）、Mei 等（Phys. Rev. D 107, 074018, [arXiv:2212.04778](https://arxiv.org/abs/2212.04778)）以及本轮提供的 Tian 等（Phys. Rev. D 114, 034012 (2026)）。Tian 等 Eq. (3) 明确给出了 `K_{03}`、`K_{38}` 的相对符号；当前实现已按该来源和本项目的 `phi_f=sigma_f` 约定对齐。

在调用方提供三味 flavor-diagonal 泡 `Pi_u, Pi_d, Pi_s` 后，Tian 等 Eq. (26) 的中性极化矩阵元为

```math
\begin{aligned}
\Pi_0 &= \frac{2}{3}(\Pi_u+\Pi_d+\Pi_s), &
\Pi_3 &= \Pi_u+\Pi_d, &
\Pi_8 &= \frac{1}{3}(\Pi_u+\Pi_d+4\Pi_s),\\
\Pi_{03}=\Pi_{30} &= \frac{\sqrt{6}}{3}(\Pi_u-\Pi_d), &
\Pi_{08}=\Pi_{80} &= \frac{\sqrt{2}}{3}(\Pi_u+\Pi_d-2\Pi_s), &
\Pi_{38}=\Pi_{83} &= \frac{\sqrt{3}}{3}(\Pi_u-\Pi_d).
\end{aligned}
```

中性 RPA 传播子按 Eq. (20) 使用

```math
\mathcal M^P=2\mathcal K^P\left[1-2\mathcal K^P\Pi^P\right]^{-1},
```

其中矩阵乘法顺序必须保留；它不能在一般同位旋不对称背景下被替换成逐元素标量公式。
这里的 RPA 不是“只取作用量二阶项”的一次微扰修正，而是在固定平均场背景
上保留二次介子涨落并对夸克泡链作无限重求和。相应的 `Pi` 还依赖
`T,mu_f,Phi,\bar Phi`、正则化和外部运动学；不能简化成只依赖 `m_f` 与 `q^2`。

### 6.1 charged scalar 与 matrix 归一化

对 charged pair `(i,j)`，令 `lambda_a,lambda_b` 为对应的两个实 Gell-Mann
生成元，并定义

```math
T_+=(\lambda_a+i\lambda_b)/\sqrt2=\sqrt2 E_{ij},\qquad
T_-=(\lambda_a-i\lambda_b)/\sqrt2=\sqrt2 E_{ji}.
```

若 `Pi_ij` 是项目 `PolarizationAniso`/Rehberg 定义的单个 ordered flavor 泡，
则二次作用量中的 flavor trace 直接给出

```math
\Pi^{matrix}_{+-}=2\Pi_{ij},\qquad
\Pi^{matrix}_{-+}=2\Pi_{ji}.
```

所以在 charge basis 中

```math
2K[1-2K\Pi^{matrix}_{+-}]^{-1}
=\frac{2K}{1-4K\Pi_{ij}}.
```

这个因子 2 来自两个 charged ladder 顶角，不依赖 isospin 对称，也不是事后
匹配。`tests/unit/relaxtime/test_meson_rpa.jl` 以显式 3x3 生成元、任意且不相等的
`Pi_ij/Pi_ji` 和 chiral-limit Goldstone identity 锁定该结论。它闭合的是场/泡
归一化；ordered retarded 泡的解析和数值实现仍由 charged-RPA/BU 路线验收。

接入极化矩阵前仍需逐项核对：

1. `phi_f` 与文献 `sigma_f` 或 `i tr S^f` 的符号；
2. P/S 上标与 `K^±` 的对应关系；
3. `K` 的配置符号和单位；
4. `lambda_0` 归一化；charged 通道必须使用第 6.1 节的 ladder 转换，不能把单个
   ordered `Pi_ij` 直接代入矩阵 `1-2KPi`。

只有完成上述核对、并在同位旋对称极限和一个外部文献固定点通过验证后，才允许把后续完整 RPA 数值结果纳入生产候选。当前符号对齐只证明代数约定，不证明有限温密度 PNJL 极化或介子数密度的物理正确性。
