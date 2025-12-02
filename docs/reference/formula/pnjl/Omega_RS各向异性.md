下面我把 PNJL 模型在**动量各向异性（momentum-space anisotropy）**情形下的修改整理成完整的 Markdown 文档。结论先说在前面，然后把所有修改（带数学表达）写清楚、并给出数值实现的要点与注意事项。

**结论（简短）**

* 如果你的“各向异性”是指**非平衡的动量空间变形（例如 Romatschke–Strickland 型变形）**，且你只希望在热（有限温度/化学势）部分反映该各向异性 —— **只需把出现在分布函数（热项、含 (T) 的那部分）里的能量 (E=\sqrt{p^2+m^2})** 替换为带各向异性项的 (E_{\rm aniso}(p,\hat t,\xi))（或等价地在分布函数中替换 (p\to p_{\rm eff})）。
* **真空项（零温的 Dirac sea、即带截断的项）通常保留原来的各向同性能量 (E=\sqrt{p^2+m^2})**，不作修改，因为真空项来源于基态并非热介质的动量分布变形；这是多数文献和数值实现中常用且物理上合理的做法。
* 但若你认为各向异性会改变准粒子色散（例如自能引入方向依赖），那必须同时修改真空项与能隙方程（见下“例外/警告”）。

下面是详细的 `.md` 文档，已把所有必要公式、替换说明、以及数值实现要点写好（并保留自动微分友好的原始偏导定义）。

> 实现说明：`src/pnjl/solvers/AnisoGapSolver.jl` 已按照上述约定拆分能量贡献，`calculate_energy_sum` 始终使用各向同性能量，而热项 `calculate_log_sum` 才读取 `ξ`。若未来需要让真空项也带各向异性，请同步更新该模块并在此文档注明新的假设。

---

# PNJL（各向异性）模型：巨热力学势与序参量求解（Markdown）

> 基于论文中 PNJL 部分（式 2.23–2.39，2.26 等）。原始无各向异性表达见论文。

---

## 1. 约定与符号

* ( \mathbf{p} ) — 三维动量矢量，模长 (p=|\mathbf{p}|)。
* ( \hat{t} )（或 (\mathbf{n})）— 规定的单位矢量，表示各向异性方向（你用的符号是 (t)，文中我写成 (\hat t)）。
* ( \xi ) — 各向异性参数（实数）。常见情形：(\xi>0) 表示沿 (\hat t) 压缩/拉伸的不同类型变形，根据你采用的变形约定符号可能不同；在文中按你给出的形式使用。
* ( m_i )、(M_i)、( \phi_i ) 等保持与无各向异性时相同的定义（夸克流质量、准粒子质量与凝聚）。

---

## 2. 各向异性下的“有效能量”用于分布函数

将分布函数中出现的能量 (E=\sqrt{p^2 + M_i^2}) 替换为：

```math
E_{i}^{\rm aniso}(\mathbf p;\xi,\hat t)
= \sqrt{\, p^2 + M_i^2 + \xi\,(\mathbf p\cdot\hat t)^2 \, }.
```

（你给出的形式 (\sqrt{p^2 + m_i^2 + \xi (p\cdot t)^2}) 与上式等价；注意这里 (M_i) 表示准粒子质量（由能隙方程给出））

物理含义：在分布函数的指数中用 (E_i^{\rm aniso})（或等价形式）反映动量分布的方向性。

---

## 3. 巨热力学势 Ω（含各向异性，仅在热项中替换）

以论文式（2.26）为基础（参见原文），写出带各向异性的巨热力学势。**仅对“有限温度的热项”使用 (E_i^{\rm aniso})，真空（零温、带截断）项保持原式（各向同性）**。

```math
\Omega(\phi_i,\Phi,\bar\Phi;T,\mu,\xi)
=
2G \sum_{i=u,d,s} \phi_i^2
- 4K\,\phi_u\phi_d\phi_s
- 2N_c \int_{|\mathbf p|<\Lambda}\frac{d^3p}{(2\pi)^3}\sum_{i=u,d,s} E_i(p)           % vacuum term: E_i = sqrt(p^2 + M_i^2)
- 2T \sum_{i=u,d,s} \int\!\frac{d^3p}{(2\pi)^3}\, \big( \mathcal{Q}_1^{\rm aniso} + \mathcal{Q}_2^{\rm aniso} \big)
+ U(\Phi,\bar\Phi;T).
```

说明：

* 真空积分里的 (E_i(p)=\sqrt{p^2+M_i^2}) **不替换**（保留各向同性）。
* 热项（含温度因子 (T)）中的 (\mathcal{Q}_{1,2}^{\rm aniso}) 使用 (E_i^{\rm aniso}(\mathbf p)) 代替 (E_i)。

> 真空部分的三维积分仍可套用 `Omega_各向同性.md` 中给出的解析结果 (球坐标下一维积分)，实现时只需将该封闭表达与热项的各向异性处理进行拼合即可。

---

## 4. 各向异性下的 (\mathcal{Q}_1,\mathcal{Q}_2) 与有效分布

把式 (2.27),(2.28) 的 (\mathcal{Q}) 改为：

```math
\mathcal{Q}_1^{\rm aniso}(\mathbf p)
= \ln\!\Big(1 + 3\Phi e^{-(E_i^{\rm aniso}-\mu_i)/T}
  + 3\bar\Phi e^{-2(E_i^{\rm aniso}-\mu_i)/T}
  + e^{-3(E_i^{\rm aniso}-\mu_i)/T}\Big)
```

```math
\mathcal{Q}_2^{\rm aniso}(\mathbf p)
= \ln\!\Big(1 + 3\bar\Phi e^{-(E_i^{\rm aniso}+\mu_i)/T}
  + 3\Phi e^{-2(E_i^{\rm aniso}+\mu_i)/T}
  + e^{-3(E_i^{\rm aniso}+\mu_i)/T}\Big)
```

同时有效费米分布函数（自动微分/数值时可能会用到）也替换能量参数：

```math
n_+^{\rm aniso}(E_i^{\rm aniso}-\mu_i)
=
\frac{\Phi e^{-(E_i^{\rm aniso}-\mu_i)/T} + 2\bar\Phi e^{-2(E_i^{\rm aniso}-\mu_i)/T} + e^{-3(E_i^{\rm aniso}-\mu_i)/T}}
{1 + 3\Phi e^{-(E_i^{\rm aniso}-\mu_i)/T} + 3\bar\Phi e^{-2(E_i^{\rm aniso}-\mu_i)/T} + e^{-3(E_i^{\rm aniso}-\mu_i)/T}}
```

以及对应的 (n_-^{\rm aniso}(E_i^{\rm aniso}+\mu_i))。
（原式见论文式 2.27–2.39，已在此替换为各向异性能量。）

---

## 5. 能隙方程（Gap equations）与最小化条件

巨热力学势极小化条件仍然成立（原始导数形式，适合自动微分）：

```math
\frac{\partial\Omega}{\partial \phi_i} = 0,\quad i=u,d,s,
```

```math
\frac{\partial\Omega}{\partial \Phi} = 0,\quad
\frac{\partial\Omega}{\partial \bar\Phi} = 0.
```

**实现要点**：在计算这些偏导时，请确保 (\Omega) 的热项中使用的是 (E_i^{\rm aniso}(\mathbf p))（对 (\phi_i) 的隐含依赖通道仍通过 (M_i(\phi)) 反映），而真空项中仍用 (E_i(p))。用自动微分求导会自动把热项中 (E_i^{\rm aniso}) 对 (\phi_i) 的依赖（通过 (M_i)）带入。

---

## 6. 热力学量（仍以 Ω 的偏导定义）

**（保持原来的偏导定义，适合自动微分）**

* 压强：

  ```math
  P = -\Omega
  ```
* 夸克数密度：

  ```math
  \rho_i = -\frac{\partial \Omega}{\partial \mu_i}
  ```
* 熵密度：

  ```math
  s = -\frac{\partial \Omega}{\partial T}
  ```
* 能量密度：

  ```math
  \epsilon = Ts + \sum_i \mu_i\rho_i - P
  ```

注意：由于热项使用了 (E_i^{\rm aniso})（含方向依赖），在做偏导（尤其是对 (\mu) 和 (T)）时，自动微分会正确考虑方向角积分下的贡献，无需把展开形式写出。

---

## 7. 数值实现要点（积分与坐标）

因为 (E_i^{\rm aniso}) 含 ((\mathbf p\cdot\hat t)^2)，原先可直接用“仅对 (p) 做一维 radial 积分”的做法**不再充足**。推荐做法：

1. 使用球坐标，令 (\hat t) 为极轴（即把坐标系旋转使 (\hat t) 指向 (\theta=0)），则
   [
   (\mathbf p\cdot\hat t)^2 = p^2 \cos^2\theta.
   ]
   因此
   [
   E_i^{\rm aniso} = \sqrt{ p^2 + M_i^2 + \xi p^2 \cos^2\theta }.
   ]

2. 热项为双重或三重积分（角向不可退化）：
   [
   \int \frac{d^3p}{(2\pi)^3}(\cdots)
   = \int_0^\infty \frac{p^2 dp}{2\pi^2}\int_{-1}^{1}\frac{d\cos\theta}{2}(\cdots)\times 2\pi
   ]
   其中积分核中含 (E_i^{\rm aniso}(p,\cos\theta))。（即数值上至少需要对 (p) 和 (\cos\theta) 做二维积分；(\phi) 角可直接积分给 (2\pi)）

3. 真空项仍可沿用原来的一维径向积分（上限为 (\Lambda)）：
   [
   \int_{|\mathbf p|<\Lambda}\frac{d^3p}{(2\pi)^3} E_i(p)
   = \int_0^\Lambda \frac{p^2 dp}{2\pi^2} E_i(p).
   ]

4. 积分方案建议：

   * 对 (p) 使用合适的高斯节点或自适应积分（注意 (p) 范围通常从 0 到较大值；热项在高温时尾部贡献更大）；
   * 对 (\cos\theta) 使用 Gauss-Legendre；
   * 将二维积分矢量化/并行化以提高速度。

---

## 8. 代码/自动微分实现提示（简要）

* 把 (\Omega) 实现为函数 (\Omega(\phi,\Phi,\bar\Phi;T,\mu,\xi))；热项为带二维积分的数值函数。
* 使用自动微分（例如 Zygote、ForwardDiff、ReverseDiff 等——取决于你的语言）对 (\Omega) 求偏导得到 gap 方程和热力学量（(\rho,s)）。
* 注意：自动微分对包含数值积分的函数有两种实现路径：

  1. **对被微分的数值积分使用可微的积分操作**（例如将积分实现为权重和节点的加权和，节点与权重固定），这样 AD 能直接传播导数；或
  2. **使用手动推导的雅可比/微分（例如对参数做有限差分）**（不推荐，如果你想保留 AD 的精度与可靠性）。
* 建议把积分用固定节点（如 Gauss-Legendre/Gauss-Laguerre 根据变换）离散化，然后用 AD 在离散和上求导（这样 AD 能工作且稳定）。

---

## 9. 例外与物理/数值警告

* 上述做法是假定**各向异性仅是动量空间分布/相空间变形**（非对准粒子自能/色散产生方向依赖）的情况下。这是常见且保守的做法。
* 若你认为介质的自能（self-energy）或准粒子谱本身会变成方向依赖（例如介质极化使质量 (M_i) 或动能项方向依赖），则：

  * 必须把真空项、能隙方程、以及(E_i(p)) 的定义一并修改；
  * 那样就不再是“仅修改分布函数”的情形，且理论上更复杂（需要重新考察重整化/截断的一致性）。
* 数值上，若 (\xi) 很大，积分变得更加尖锐（角向依赖更强），需要提高角向积分分辨率。

---

## 10. 与原论文公式的对应与引用

* 无各向异性时的巨热力学势与 (\mathcal{Q}_1,\mathcal{Q}_2) 的原式见论文式 (2.26)–(2.28)。我在上面将其中的热项（(\mathcal{Q})）里的 (E_i) 全部替换为 (E_i^{\rm aniso})，真空项保留原式。
* 原论文中用于定义 (n_\pm) 的表达（式 2.38–2.39）在本处相应替换为 (n_\pm^{\rm aniso})。

---

## 11. 更新后的关键公式汇总（供复制粘贴）

（只列出关键替换，不展开热项的完整展开，以便与 AD 兼容）

* 各向异性能量（用于热项）：

  ```math
  E_i^{\rm aniso}(\mathbf p) = \sqrt{p^2 + M_i^2 + \xi (\mathbf p \cdot \hat t)^2 }.
  ```
* 热项中 (\mathcal{Q}_{1,2})：

  ```math
  \mathcal{Q}_1^{\rm aniso}(\mathbf p)
  = \ln\!\Big(1 + 3\Phi e^{-(E_i^{\rm aniso}-\mu_i)/T} + 3\bar\Phi e^{-2(E_i^{\rm aniso}-\mu_i)/T} + e^{-3(E_i^{\rm aniso}-\mu_i)/T}\Big)
  ```

  ```math
  \mathcal{Q}_2^{\rm aniso}(\mathbf p)
  = \ln\!\Big(1 + 3\bar\Phi e^{-(E_i^{\rm aniso}+\mu_i)/T} + 3\Phi e^{-2(E_i^{\rm aniso}+\mu_i)/T} + e^{-3(E_i^{\rm aniso}+\mu_i)/T}\Big)
  ```
* Ω（总式，热项用 aniso，真空项原式）：

  ```math
  \Omega = 2G \sum_i \phi_i^2 - 4K\phi_u\phi_d\phi_s
    - 2N_c \int_{|\mathbf p|<\Lambda}\frac{d^3p}{(2\pi)^3}\sum_i E_i(p)
    - 2T \sum_i \int\frac{d^3p}{(2\pi)^3}(\mathcal{Q}_1^{\rm aniso} + \mathcal{Q}_2^{\rm aniso})
    + U(\Phi,\bar\Phi;T).
  ```

---

## 12. 参数（与无各向异性时一致）

（见之前整理的 PNJL 参数；若需要我可以把它们复制到此处；原参数来自论文表 2.1/2.2）

---


