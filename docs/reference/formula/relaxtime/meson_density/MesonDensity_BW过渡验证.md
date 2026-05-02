# 介子数密度：BW（Breit-Wigner）过渡验证

本文固定当前介子数密度主线中的 BW（Breit-Wigner）过渡验证层，用于在稳定粒子极限和 BU 全相移公式之间建立可检查的中间台阶。

相关文档：

- [MesonDensity_稳定粒子与KPi比值.md](./MesonDensity_%E7%A8%B3%E5%AE%9A%E7%B2%92%E5%AD%90%E4%B8%8EKPi%E6%AF%94%E5%80%BC.md)
- [MesonDensity_BU相移公式.md](./MesonDensity_BU%E7%9B%B8%E7%A7%BB%E5%85%AC%E5%BC%8F.md)
- [MesonMass_RPA_Pole.md](../propagator/MesonMass_RPA_Pole.md)

## 1. 引用口径与文献分工

本文中的 BW 需要同时处理“极点从哪里来”和“为什么它可以作为 BU 的中间近似层”两个问题，因此文献分工必须拆开：

1. `Hatsuda:1994pi`
   - 角色：NJL/RPA 介子传播子与极化函数结构的综述性背景；
   - 作用：解释为什么可以围绕传播子极点做共振近似。
2. `Rehberg:1996da`
   - 角色：当前项目简单介子道传播子/极点方程的直接 SU(3) NJL 背书；
   - 作用：支撑
     ```math
     \mathcal{D}_M=\frac{2K_M}{1-4K_M\Pi_M},
     \qquad
     \mathcal{G}_M(z,q)=1-4K_M\Pi_M(z,q)=0
     ```
     这一 strict BW 的项目底座。
3. `Maslov:2023wul`
   - 角色：phase-shift / off-shell 热力学主线背书；
   - 作用：说明 BW 在本项目中不是独立主理论，而是位于 stable 与 full BU 之间的受控极点近似层。

因此，本文的 strict BW 应理解为：

- **极点定义**对齐 `Hatsuda:1994pi` + `Rehberg:1996da`；
- **热力学位置**对齐 `Maslov:2023wul` 的 BU / off-shell 主线；
- **当前 D3 最小代理写法**仍可参考 `Blaschke:2020bzh` 周边的 BW 表达，但不再把它当作当前传播子底座的直接背书。

## 2. BW 的含义与当前角色

本文中的 `BW` 指 `Breit-Wigner` 近似。

在当前语境下，它的具体作用是：

1. 用 `Breit-Wigner` 形式近似共振极点附近的相移；
2. 等价地，用一个 Lorentzian 形式近似相移导数；
3. 以此构造稳定粒子极限与完整 BU 相移实现之间的最小过渡层。

因此，这里的 `BW` 不是一个泛称，而是对相移/相移导数所采用的**具体近似公式**。

## 3. 理论上的 BW 极点近似

文献中对窄宽度共振采用：

```math
z_M = E_M - i\Gamma_M/2.
```

对应的相移写为：

```math
\delta_M(\omega)
= -\arctan\frac{\Gamma_M/2}{\omega-E_M}.
```

### 来源

- `Blaschke:2020bzh` Eq. (26) 周边

## 4. 相移导数

由上式可得：

```math
\frac{d\delta_M}{d\omega}
= \frac{1}{2}\frac{\Gamma_M}{(\omega-E_M)^2+\Gamma_M^2/4}.
```

### 来源

- `Blaschke:2020bzh` Eq. (26)

## 5. BW 数密度公式

把 BW 相移导数代入 BU 原式，可得：

```math
n_M(T)
= d_M \int \frac{d^3q}{(2\pi)^3}
\int \frac{d\omega}{2\pi}
g_M(\omega)\frac{\Gamma_M}{(\omega-E_M)^2+\Gamma_M^2/4}.
```

### 来源

- `Blaschke:2020bzh` Eq. (27)

## 6. 与当前项目传播子口径的关系

`z_M = E_M - i\Gamma_M/2` 不是一个与传播子完全无关的独立定义。

它本质上来自对某个给定传播子在共振极点附近的局部近似，因此：

1. `E_M` 与 `\Gamma_M` 的定义依赖于你采用什么传播子；
2. 若项目中的传播子分母结构与文献不同，那么“如何从传播子提取 `E_M,\Gamma_M`”这一步会受影响；
3. 但一旦你已经从当前项目自己的传播子口径中提取出了 `(M_M,\Gamma_M)`，再把它们代入一个 Lorentzian / Breit-Wigner 型近似，本身仍然是允许的近似步骤。

因此，真正与传播子强绑定的是：

- `M_M` / `\Gamma_M` 的提取方式；

而不是“写出一个带宽度的 Lorentzian 核”这一步本身。

## 7. 与当前项目的对接方式

当前项目不直接把文献中的 BW 参数当作独立输入。更合适的做法是：

1. 先从当前项目主极点条件
   ```math
   1-4K_M\Pi_M(p_0,q;\xi)=0
   ```
   求得 `M` 与 `\Gamma`
2. 再把它们代回 BW 近似表达式
3. 用作稳定极限与 BU 之间的过渡验证层

这样做可以保证 BW 近似仍服从当前仓库的传播子和极点口径。

### 7.1 复极点定义不等于 BW 一阶展开

这里需要把三个层级明确分开：

1. **复极点定义**
   ```math
   \mathcal{G}_M(z,q)=0,
   \qquad
   z_p(q)=E_M(q)-i\Gamma_M(q)/2
   ```
   这是对传播子极点位置的定义；
2. **把 `z_p(q)` 参数化成 `E_M(q)-i\Gamma_M(q)/2`**
   只是把复根拆成实部与虚部，仍然属于极点定义层；
3. **BW 一阶展开**
   ```math
   \mathcal{G}_M(z,q)
   \approx
   \mathcal{G}_M'(z_p(q),q)\,[z-z_p(q)]
   ```
   这一步才是真正把传播子压缩成 Breit-Wigner 极点近似。

因此：

- `z_M = E_M - i\Gamma_M/2` 本身不自动等于 BW；
- `1-4K_M\Pi_M(M_M+i\Gamma_M/2, q)=0` 若按复变量完整求值，也仍属于“复极点定义”；
- strict BW 则是在已经得到该复极点后，再对分母做极点邻域一阶线性化。

## 8. 基于当前项目传播子的 strict BW 推导

对当前项目简单介子道，传播子写为：

```math
\mathcal{D}_M(\omega,q)
= \frac{2K_M}{1-4K_M\Pi_M(\omega,q)}.
```

定义分母函数：

```math
\mathcal{G}_M(\omega,q)
\equiv 1-4K_M\Pi_M(\omega,q).
```

若沿实轴写

```math
\Pi_M = \Pi_M^R + i\Pi_M^I,
```

则有

```math
\mathcal{G}_M(\omega,q)
= A_M(\omega,q)+i B_M(\omega,q),
```

其中

```math
A_M(\omega,q)=1-4K_M\Pi_M^R(\omega,q),
\qquad
B_M(\omega,q)=-4K_M\Pi_M^I(\omega,q).
```

### 8.1 复极点定义

strict BW 不是额外假设一套与传播子无关的 `z_M`，而是从当前传播子的复极点出发。

在固定 `q` 下，设共振极点为

```math
z_p(q)=E_M(q)-i\Gamma_M(q)/2,
```

其定义满足

```math
\mathcal{G}_M(z_p(q),q)=0.
```

因此，`E_M(q)` 与 `\Gamma_M(q)` 的定义由当前项目自己的传播子分母决定。

### 8.2 极点邻域展开

在极点附近对分母做一阶展开：

```math
\mathcal{G}_M(z,q)
\approx
\mathcal{G}_M'(z_p(q),q)\,[z-z_p(q)].
```

于是传播子可写为

```math
\mathcal{D}_M(z,q)
\approx
\frac{2K_M}{\mathcal{G}_M'(z_p(q),q)}
\frac{1}{z-z_p(q)}
\equiv
\frac{Z_M(q)}{z-E_M(q)+i\Gamma_M(q)/2},
```

其中 residue 为

```math
Z_M(q)=\frac{2K_M}{\mathcal{G}_M'(z_p(q),q)}.
```

这就是在当前项目传播子上的 strict BW 极点近似。

### 8.3 strict BW 相移

沿实轴取 `z=\omega+i0^+`，有

```math
\mathcal{D}_M(\omega+i0^+,q)
\approx
\frac{Z_M(q)}{\omega-E_M(q)+i\Gamma_M(q)/2}.
```

因此相移写为

```math
\delta_M^{BW}(\omega,q)
=
\arg Z_M(q)
-\arctan\frac{\Gamma_M(q)/2}{\omega-E_M(q)}.
```

若 `Z_M(q)` 近似为实正量，或只关心导数形式，则常数相位 `\arg Z_M(q)` 可忽略，于是得到熟悉的 BW 相移：

```math
\delta_M^{BW}(\omega,q)
\approx
-\arctan\frac{\Gamma_M(q)/2}{\omega-E_M(q)}.
```

### 8.4 strict BW 相移导数

由上式可得

```math
\frac{d\delta_M^{BW}}{d\omega}
=
\frac{\Gamma_M(q)/2}
{\left[\omega-E_M(q)\right]^2+\Gamma_M^2(q)/4}.
```

因此，严格地说，BW 导数的 Lorentzian 形状是极点邻域展开的结果，而不是与传播子无关的先验模板。

### 8.5 基于当前传播子的 strict BW 数密度

将上式代回 BU 原式：

```math
n_M^{BW}(T)
= d_M \int \frac{d^3q}{(2\pi)^3}
\int \frac{d\omega}{2\pi}
g_M(\omega)
\frac{\Gamma_M(q)/2}
{\left[\omega-E_M(q)\right]^2+\Gamma_M^2(q)/4}.
```

若进一步写成一维动量积分形式：

```math
n_M^{BW}(T)
= d_M \int_0^\infty \frac{dq\,q^2}{2\pi^2}
\int_0^\infty \frac{d\omega}{2\pi}
g_M(\omega)
\frac{\Gamma_M(q)/2}
{\left[\omega-E_M(q)\right]^2+\Gamma_M^2(q)/4}.
```

### 8.6 对应到当前项目实现层级

把上述 strict BW 落到当前仓库时，还要再区分两个实现层级：

1. **Stage 1 reduced strict BW**
   - 只使用 `q=0` 处 meson workflow 已给出的 `(M_M,\Gamma_M)`；
   - 在积分中采用
     ```math
     E_M(q)\approx \sqrt{q^2+M_M^2},
     \qquad
     \Gamma_M(q)\approx \Gamma_M(0).
     ```
2. **Stage 2 q-pole strict BW**
   - 在 `q` 网格上逐点求解
     ```math
     \mathcal{G}_M(z_p(q),q)=0,
     \qquad
     z_p(q)=E_M(q)-i\Gamma_M(q)/2,
     ```
   - 再把该 `E_M(q),\Gamma_M(q)` 回填到 BW 双积分。

因此，本项目里：

- `z_p(q)=E_M(q)-i\Gamma_M(q)/2` 属于**复极点定义层**；
- 在极点邻域把传播子分母线性化，并得到 Lorentzian 型相移/导数，才属于**BW 近似层**；
- Stage 1 与 Stage 2 的差异，不在于“有没有复极点”，而在于是否真正保留了 `q` 依赖极点信息。

这就是当前项目传播子口径下应实现的 strict BW 目标公式。

### 8.6 与当前最小冻结近似的关系

若暂时只掌握 `q=0` 的 workflow 输出 `(M_M,\Gamma_M)`，可再引入额外近似

```math
E_M(q)\approx \sqrt{q^2+M_M^2},
\qquad
\Gamma_M(q)\approx \Gamma_M.
```

则得到一个“半严格”的 reduced BW 公式：

```math
n_M^{BW,red}(T)
= d_M \int_0^\infty \frac{dq\,q^2}{2\pi^2}
\int_0^\infty \frac{d\omega}{2\pi}
g_M(\omega)
\frac{\Gamma_M/2}
{\left[\omega-\sqrt{q^2+M_M^2}\right]^2+\Gamma_M^2/4}.
```

这比 `BW-proxy` 更接近 strict BW，但仍不是完全由 `q` 依赖极点决定的最终口径。

## 9. 当前仓库里的实际 BW 实现不是严格极点 BW

当前仓库脚本级 `BW` 实现位于：

- [validate_meson_density_bw_minimal.jl](/C:/Users/Wmzx/.codex/worktrees/346c/Julia_RelaxTime/scripts/analysis/relaxtime/validate_meson_density_bw_minimal.jl)

它目前并没有直接实现：

```math
n_M(T)
= d_M \int \frac{d^3q}{(2\pi)^3}
\int \frac{d\omega}{2\pi}
g_M(\omega)\frac{\Gamma_M}{(\omega-E_M)^2+\Gamma_M^2/4}.
```

当前实际采用的是一个更弱的 `BW-proxy`：

```math
n_M^{BW-proxy}(T)
= d_M \int_0^\infty dm\,\rho_M^{BW}(m; M_M,\Gamma_M)\,n_M^{stable}(T;m),
```

其中 `\rho_M^{BW}` 是按 `(M_M,\Gamma_M)` 构造的 Lorentzian 质量展宽核。

这意味着：

1. 当前 `BW` 结果主要反映“有限宽度展宽”这个效应；
2. 它不是严格从当前传播子相位或传播子极点残数直接推出的完整 BW 数密度；
3. 它的职责是 `Phase D3` 过渡验证，而不是最终主线物理结果。

## 10. strict BW 的实现前提

若未来要在项目中真正实现 strict BW，对当前简单介子道至少还需要：

1. 在固定 `q` 下，从 `\mathcal{G}_M(z,q)=0` 提取复极点
   `z_p(q)=E_M(q)-i\Gamma_M(q)/2`
2. 判断 residue `Z_M(q)` 的相位是否可近似视为常数
3. 在 `q` 网格上稳定计算 `E_M(q)` 与 `\Gamma_M(q)`
4. 再将它们代入本页的 strict BW 双积分公式

因此，strict BW 是当前项目中合理的后继对照层，但它不同于现有 `BW-proxy`。

## 11. 稳定极限退化关系

当

```math
\Gamma_M \to 0
```

时，BW 数密度应退化回稳定粒子极限。

这条退化关系是当前 `Phase D` 的关键一致性检查之一：

1. 若退化不成立，则实现、单位或归一化存在问题；
2. 若退化成立，则 BW 层可作为过渡验证台阶继续保留。

## 12. 当前用途

本文只服务于：

1. `Phase D3`：是否需要 BW 过渡验证的判断
2. `Phase D4`：把 BW 明确标记为“验证台阶而非最终目标”
3. 后续与 BU 结果的最小差异对照
4. 后续 strict BW 正式实现时的公式基线

本文不把 BW 视为当前主线最终物理结果。
