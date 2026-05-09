# 介子热力学：BU / off-shell / Landau damping EOS 公式主线

本文档用于把 `Maslov & Blaschke, Phys. Rev. D 107, 094010 (2023)` 的 mesonic off-shell correlation / BU equation-of-state 主线，整理成适合当前仓库的系统公式文档。

它的目标不是逐字复刻论文，而是回答三件事：

1. 文献的 EOS 主公式到底是什么；
2. 它应如何映射到当前项目的 `propagator / polarization / phase shift` 链；
3. 按当前仓库标准，哪些量属于 `relaxtime` 公式层，哪些量应交给 `Models` 的统一 `Omega -> AD thermodynamics` 流程。

---

## 1. 适用边界与当前不完全同构处

先明确当前项目与该文献的同构边界：

1. 文献模型是 `N_f = 2` PNJL，主通道是 `pi / sigma`。
2. 结合当前项目的 meson 命名，这里的 `sigma` 更接近轻味标量伙伴 `sigma_pi`，而不是 strange 赝标量 `K`。
3. 当前项目已落地的 meson thermo 第一版主通道是 `pi / K`，还没有把 `sigma_pi` 纳入正式 EOS workflow。
4. 因此，本文档只把这篇文献作为：
   - **方法结构参考**
   - **BU / off-shell / Landau damping 物理口径参考**
   - **未来 `pi/sigma_pi` validation 路径的上游公式依据**
5. 当前项目不能据此声称“已经严格复现该文献数值结果”。

换句话说，本文档固定的是**公式主线与实现语义**，不是当前版本的完全 validation 声明。

---

## 2. 总 EOS 结构

文献主线的总压强结构是：

```math
P_{\mathrm{total}}(T,\mu)
=
P_{\mathrm{MF}}(T,\mu)
+ \sum_M P_M(T,\mu).
```

对当前项目，应保持同样的分层：

- `P_MF`：均匀平均场 / quark + mean-field + Polyakov 背景压强；
- `P_M`：介子相关压强；
- `P_total`：二者相加后的总压强。

与巨热力学势的关系是：

```math
P = -\Omega,
```

因此总巨热力学势必须写成：

```math
\Omega_{\mathrm{total}}
=
\Omega_{\mathrm{MF}}
- \sum_M P_M.
```

这条符号关系在当前项目里非常关键，因为后续所有 EOS 派生量应从 `\Omega_total` 统一导出，而不是在 workflow 层各写一套手工差分逻辑。

---

## 3. 介子压强的 BU 相移公式

文献使用的介子压强主线是广义 Beth-Uhlenbeck 表达：

```math
P_M
=
d_M \int \frac{d^3 q}{(2\pi)^3}
\int_0^\infty \frac{d\omega}{\pi}
\frac{\delta_M(\omega,q;T)}{e^{\omega/T}-1}.
```

对当前项目的实现语义，应理解为：

1. `\delta_M(\omega,q;T)` 不是外加输入，而是来自当前项目的传播子相位；
2. `P_M` 是**介子相关压强**，不是简单理想玻色气结果；
3. `stable` 与 `strict BW reduced` 只能作为：
   - 稳定粒子极限参考；
   - reduced resonance 对照口径；
   - 数值 sanity check；
   不能取代 BU / off-shell 主 EOS 口径。

当前项目若要贴近该文献，正式 EOS 主口径应优先使用：

- `phase_shift_current`
- `phase_shift_gbu_reference`

而不是 `stable_meson_pressure` 或 `strict_bw_stage1_reduced_pressure`。

---

## 4. 相移与当前项目传播子的映射

文献中的介子相关压强建立在 RPA 重求和传播子上。对当前项目，保持已有 `MesonMass / MesonPropagator / PolarizationAniso` 口径即可。

当前 simple meson channel 的项目传播子写成：

```math
\mathcal D_M(\omega,q;\xi)
=
\frac{2K_M}{1 - 4K_M \Pi_M(\omega,q;\xi)}.
```

当前 BU / phase-shift 主线应定义：

```math
\delta_M(\omega,q;\xi)
=
\arg \mathcal D_M(\omega + i\eta, q; \xi).
```

因此本项目的 `meson thermo` 与 `meson density` 在公式层必须共享同一条链：

```text
polarization -> propagator -> phase shift -> BU integral
```

而不能在 pressure 侧另造一套与已有传播子不同的记号系统。

---

## 5. Landau damping 必须作为一等公民

这是该文献相对当前第一版实现最重要的物理强调点。

文献把介子压强按频率-动量区域拆为：

```math
P_M = P_M^{\mathrm{QP}} + P_M^{\mathrm{LD}}.
```

其中：

- `QP`（quasipole / timelike）区：通常对应 `\omega > q`
- `LD`（Landau damping / spacelike）区：对应
  ```math
  \omega^2 - q^2 < 0,
  \qquad 0 < \omega < q
  ```

因此更接近文献的分拆写法应为：

```math
P_M^{\mathrm{QP}}
=
d_M \int \frac{d^3 q}{(2\pi)^3}
\int_{\omega \in \mathrm{QP}} \frac{d\omega}{\pi}
\frac{\delta_M(\omega,q;T)}{e^{\omega/T}-1},
```

```math
P_M^{\mathrm{LD}}
=
d_M \int \frac{d^3 q}{(2\pi)^3}
\int_{\omega \in \mathrm{LD}} \frac{d\omega}{\pi}
\frac{\delta_M(\omega,q;T)}{e^{\omega/T}-1}.
```

对当前项目，这意味着文档上必须先固定：

1. `LD` 区不是“误差项”，而是正式物理贡献；
2. `LD` 与 `QP` 应在合同层可分开导出；
3. 若对 LD 引入动量截断，应显式写成治理参数，而不是隐藏在脚本里。

建议的目标字段：

- `P_pi_qp`, `P_pi_ld`
- `P_K_qp`, `P_K_ld`
- `P_meson_qp`, `P_meson_ld`
- `ld_cutoff`
- `ld_threshold_mode`

当前仓库第一版 `phase_shift_*_meson_thermo` 还没有做到这一步；本文档把它定义成**后续公式与合同收口的硬约束**。

---

## 6. cutoff / threshold 的文献方向

文献对 LD 区的一个核心结论是：

- `P_M^{LD}` 对 cutoff / threshold 很敏感；
- 这种敏感性会直接反映到 mesonic pressure 与 trace anomaly 上。

因此对当前项目，凡是后续进入 canonical scan 的 BU / LD 热力学结果，都应记录：

- `qmax`
- `omega_max`
- `ld_cutoff`
- `scheme`
- `phase_shift_variant`

并把 `cutoff sensitivity` 作为正式 plot-review 资产的一部分，而不是只在交互式分析里临时观察。

---

## 7. 与现有 meson density 文档的关系

当前仓库已有：

- `MesonDensity_BU相移公式.md`

它固定的是**数密度**主线。本文档固定的是**压强 / EOS** 主线。

二者共用同一个相移对象 `\delta_M(\omega,q)`，但权重不同：

### 7.1 数密度主线

当前优先实现口径可写成：

```math
n_M(T)
=
\frac{d_M}{T}
\int \frac{dq\,q^2}{2\pi^2}
\int_0^\infty \frac{d\omega}{2\pi}
g_M(\omega)\bigl[1+g_M(\omega)\bigr]\delta_M(\omega,q).
```

### 7.2 压强主线

本文固定的压强口径是：

```math
P_M(T)
=
d_M \int \frac{d^3 q}{(2\pi)^3}
\int_0^\infty \frac{d\omega}{\pi}
\frac{\delta_M(\omega,q;T)}{e^{\omega/T}-1}.
```

因此在工程分层上：

- `MesonDensity` 模块负责 `n_M`
- `MesonThermodynamics` 模块负责 `P_M`
- 二者共享 `polarization / propagator / phase-shift` 底座

---

## 8. 派生 EOS 量不在本层手工定义

文献关心的常见 EOS 派生量包括：

- entropy density
- energy density
- trace anomaly

但对当前项目，公式层应明确遵守以下约束：

1. `relaxtime/meson_thermo` 层只负责给出 `P_M` 或等价地给出 `\Omega_M = -P_M`
2. 总压强 / 总巨热力学势写成：
   ```math
   \Omega_{\mathrm{total}}
   =
   \Omega_{\mathrm{MF}}
   - \sum_M P_M
   ```
3. `s, n, \epsilon, I=(\epsilon-3P)/T^4` 等派生量，统一交给 `Models` 主域热力学流程，通过自动微分从 `\Omega_total` 导出

也就是说，**系统公式标准**应是：

```math
P_{\mathrm{total}} = -\Omega_{\mathrm{total}},
```

```math
\rho_i = - \frac{\partial \Omega_{\mathrm{total}}}{\partial \mu_i},
```

```math
s = - \frac{\partial \Omega_{\mathrm{total}}}{\partial T},
```

```math
\epsilon = -P_{\mathrm{total}} + T s + \sum_i \mu_i \rho_i,
```

```math
I(T) = \frac{\epsilon - 3P_{\mathrm{total}}}{T^4}.
```

这正是你指出的项目标准流程：  
只要把介子压强并入总巨热力学势，就不应再把 `entropy_meson` / `trace_anomaly` 写成 workflow 层的私有差分定义。

---

## 9. 当前项目的推荐实现顺序

按本文档对应的系统主线，推荐顺序是：

1. 先固定 `phase-shift` 为正式 meson EOS 主口径；
2. 在 `MesonThermodynamics` 中显式拆出 `QP / LD`；
3. 把介子相关项提升为 `\Omega_total` 的组成部分；
4. 复用 `src/models/thermo_kernel.jl` 的统一 AD 派生流程；
5. 再做 canonical `mu_B = 0` 温度扫描、cutoff sensitivity、LD/QP 分拆图。

---

## 10. 当前版本与文献的关系总结

截至当前仓库状态，更准确的描述应是：

1. 已经具备“把介子压强并入 EOS”的实现骨架；
2. `phase-shift` 路径已经在方法结构上对齐该文献；
3. 但尚未完成该文献最关键的 `LD` 一等公民化：
   - 还没有 `QP / LD` 显式分拆输出
   - 还没有 `LD cutoff sensitivity` 正式资产
4. 当前正式 workflow 仍是 `pi/K`，而文献方向更接近 `pi/sigma_pi`
5. 派生量的系统标准应回到：
   - `Omega_total`
   - `Models.model_thermo`
   - `ForwardDiff`

这也是本文档最重要的收口结论。

补充的 `QP / LD` 分区与 cutoff 治理口径见：

- [MesonThermo_QP_LD_Cutoff_Governance.md](MesonThermo_QP_LD_Cutoff_Governance.md)
