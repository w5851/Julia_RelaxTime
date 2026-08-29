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

## 4. 与当前代码的兼容关系

`src/models/pnjl_physics/PNJLCore.jl` 当前的标量平均场质量方程已经允许 `phi_u`、`phi_d`、`phi_s` 分开出现，但在实际平衡状态和旧传播子接口中，u/d 介子耦合仍被合并。本路线不改变这些上游方程，只把三味凝聚作为输入生成一个可审计的完整核。

当 `phi_u=phi_d=phi_l` 时，旧 `EffectiveCouplings.calculate_effective_couplings` 使用 `G_u=-phi_l`、`G_s=-phi_s`。在本文件约定下，它的 `K0/K123/K4567/K8/K08` 与新后端的 P/S 对应分量相等；单元测试锁定这一兼容性。这个等式是 API/符号映射的测试，不是对所有有限同位旋物理的证明。

## 5. 文献证据与未决项

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

接入极化矩阵前仍需逐项核对：

1. `phi_f` 与文献 `sigma_f` 或 `i tr S^f` 的符号；
2. P/S 上标与 `K^±` 的对应关系；
3. `K` 的配置符号和单位；
4. `lambda_0` 归一化以及 RPA 矩阵中的额外 2 因子。

只有完成上述核对、并在同位旋对称极限和一个外部文献固定点通过验证后，才允许把后续完整 RPA 数值结果纳入生产候选。当前符号对齐只证明代数约定，不证明有限温密度 PNJL 极化或介子数密度的物理正确性。
