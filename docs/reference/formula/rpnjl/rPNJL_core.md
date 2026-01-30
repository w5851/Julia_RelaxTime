## 【提取报告：rPNJL模型拉氏量与热力学势】

### 📍 源信息
- **文献**：强相互作用物质相变与输运性质的研究.pdf (HWB博士学位论文)
- **位置**：第3.3.1节 "rPNJL模型简介"
- **提取时间**：2024-01-20 15:45:00

---

### 🔢 目标公式（按原文顺序）

#### **【公式 3.26】rPNJL模型拉氏量**
**【类型】** 量子场论拉氏密度  
**【内容】**  
\[
\begin{aligned}
\mathcal{L}_{rPNJL} &= \sum_{f=u,d,s} \bar{\psi}_f \gamma_\mu i D^\mu \psi_f - \sum_f m_f \bar{\psi}_f \psi_f + \sum_f \mu \gamma_0 \bar{\psi}_f \psi_f \\
&\quad + \frac{g_S}{2} \sum_{a=0,\dots,8} \left[ (\bar{\psi}_f \lambda^a \psi_f)^2 + (\bar{\psi}_f \gamma_5 \lambda^a \psi_f)^2 \right] \\
&\quad - \frac{g_D}{2} \left[ \det \bar{\psi}_f (1+\gamma_5) \psi_{f'} + \det \bar{\psi}_f (1-\gamma_5) \psi_{f'} \right] \\
&\quad + \frac{g_1}{2} \left[ (\bar{\psi}_f (1-\gamma_5) \psi_m)(\bar{\psi}_m (1+\gamma_5) \psi_f) \right] \\
&\quad + g_2 \left[ (\bar{\psi}_f (1-\gamma_5) \psi_m)(\bar{\psi}_m (1+\gamma_5) \psi_f)(\bar{\psi}_f (1-\gamma_5) \psi_k)(\bar{\psi}_k (1+\gamma_5) \psi_f) \right] \\
&\quad + U(T,\Phi,\bar{\Phi})
\end{aligned}
\]

**📖 上下文关联**  
- **引入文本**："我们的计算中，采用了一个实际的PNJL（rPNJL）模型。该模型在NJL拉氏量中加入八夸克相互作用以解决由六夸克相互作用引起的真空不稳定问题，并且在Landau-Ginzburg型胶子场有效势中考虑了Vandermonde项[193-194]。"  
- **解释文本**：该拉氏量包含四费米子、六费米子、八夸克相互作用项以及Polyakov圈有效势 \(U(T,\Phi,\bar{\Phi})\)，用于描述有限温度密度下的手征对称性破缺和退禁闭相变。

---

#### **【公式 3.27】胶子场有效势**
**【类型】** 有效势表达式  
**【内容】**  
\[
\frac{U}{T^4} = -\frac{b_2(T)}{2} \Phi\bar{\Phi} - \frac{b_3}{6} (\Phi^3 + \bar{\Phi}^3) + \frac{b_4}{4} (\Phi\bar{\Phi})^2 - \kappa \ln[J(\Phi,\bar{\Phi})]
\]

**📖 上下文关联**  
- **引入文本**："胶子场有效势取如下形式："  
- **解释文本**："\(U/T^4\) 的最后一项是 Vandermonde 项，其中的 \(J(\Phi,\bar{\Phi})\) 是从 Wilson 线到 Polyakov 圈的 Jacobi 行列式[195]："

---

#### **【公式 3.28】系数 \(b_2(T)\) 的温度依赖**
**【类型】** 参数温度依赖关系  
**【内容】**  
\[
b_2(T) = a_0 + a_1 \frac{T_0}{T} \exp\left( -a_2 \frac{T}{T_0} \right)
\]

**📖 上下文关联**  
- **引入文本**："系数 \(b_2(T)\) 是温度依赖的："

---

#### **【公式 3.29】Vandermonde 项 Jacobi 行列式**
**【类型】** Jacobi 行列式表达式  
**【内容】**  
\[
J(\Phi,\bar{\Phi}) = \frac{27}{24\pi^2} \left[ 1 - 6\Phi\bar{\Phi} + 4(\Phi^3 + \bar{\Phi}^3) - 3(\Phi\bar{\Phi})^2 \right]
\]

**📖 上下文关联**  
- **引入文本**："\(J(\Phi,\bar{\Phi})\) 是从 Wilson 线到 Polyakov 圈的 Jacobi 行列式[195]："

---

#### **【公式 3.30】巨热力学势（平均场近似）**
**【类型】** 热力学势积分表达式  
**【内容】**  
\[
\begin{aligned}
\Omega &= g_S \sum_{f=u,d,s} \sigma_f^2 - \frac{g_D}{2} \sigma_u \sigma_d \sigma_s + 3\frac{g_1}{2} \left( \sum_f \sigma_f^2 \right)^2 + 3g_2 \sum_f \sigma_f^4 \\
&\quad - 2N_c \int_0^\Lambda \frac{d^3p}{(2\pi)^3} \left( \sum_f E_f \right) - 2T \sum_{f} \int \frac{d^3p}{(2\pi)^3} (\bar{Q}_1 + \bar{Q}_2) \\
&\quad + U(T,\Phi,\bar{\Phi})
\end{aligned}
\]
其中 \(E_f = \sqrt{p^2 + M_f^2}\)。

**📖 上下文关联**  
- **引入文本**："在平均场近似下可以得到巨热力学势："
- **补充说明**：公式中的 \(\bar{Q}_1\) 和 \(\bar{Q}_2\) 是耦合到Polyakov圈的夸克统计分布函数，其标准形式为：
  \[
  \bar{Q}_1 = \ln\left(1 + 3\Phi e^{-(E_f-\mu_f)/T} + 3\bar{\Phi} e^{-2(E_f-\mu_f)/T} + e^{-3(E_f-\mu_f)/T}\right)
  \]
  \[
  \bar{Q}_2 = \ln\left(1 + 3\bar{\Phi} e^{-(E_f+\mu_f)/T} + 3\Phi e^{-2(E_f+\mu_f)/T} + e^{-3(E_f+\mu_f)/T}\right)
  \]
  它们分别对应夸克（正化学势）和反夸克（负化学势）的贡献，体现了夸克禁闭效应对统计分布的修正。

---

#### **【公式 3.31】夸克能隙方程**
**【类型】** 自洽质量方程  
**【内容】**  
\[
M_f = m_f - 2g_S \sigma_f + \frac{g_D}{4} \sigma_{f+1} \sigma_{f+2} - 2g_1 \sigma_f \left( \sum_{f'} \sigma_{f'}^2 \right) - 4g_2 \sigma_f^3
\]

**📖 上下文关联**  
- **引入文本**："以及夸克的能隙方程："

---

### 🔧 公式解析

#### **变量说明表**
| 符号 | 物理意义 | 单位/备注 |
|------|----------|-----------|
| \(\psi_f\) | 夸克场（味 \(f=u,d,s\)） | 旋量场 |
| \(m_f\) | 流夸克质量 | MeV |
| \(\mu\) | 化学势 | MeV |
| \(g_S, g_D, g_1, g_2\) | 四费米子、六费米子、八夸克耦合常数 | 无量纲（与截断Λ组合） |
| \(\sigma_f\) | 夸克凝聚序参量 | MeV |
| \(\Phi, \bar{\Phi}\) | Polyakov圈及其复共轭 | 无量纲 |
| \(T\) | 温度 | MeV |
| \(N_c\) | 色数 | 通常=3 |
| \(\Lambda\) | 动量截断 | MeV |
| \(E_f\) | 单夸克能量 | MeV |
| \(\bar{Q}_1\) | 耦合到Polyakov圈的夸克统计权重对数 | 对应夸克分布 \(f^+\) 的积分形式 |
| \(\bar{Q}_2\) | 耦合到Polyakov圈的反夸克统计权重对数 | 对应反夸克分布 \(f^-\) 的积分形式 |

#### **\(\bar{Q}_1\) 和 \(\bar{Q}_2\) 的详细说明**
在PNJL/rPNJL模型中，\(\bar{Q}_1\) 和 \(\bar{Q}_2\) 是描述有限温度密度下夸克统计分布的关键函数：
- **物理意义**：它们来源于夸克传播子与Polyakov背景场的耦合，反映了在有限温度下夸克受到色禁闭效应的修正。
- **与分布函数的关系**：对 \(\bar{Q}_1\) 和 \(\bar{Q}_2\) 关于化学势求导，可得到对应的夸克和反夸克分布函数：
  \[
  f^+(E_f,\mu_f) \propto \frac{\partial \bar{Q}_1}{\partial \mu_f}, \quad f^-(E_f,\mu_f) \propto \frac{\partial \bar{Q}_2}{\partial \mu_f}
  \]
- **渐近行为**：
  - 在高温退禁闭相（\(\Phi \rightarrow 1\)），\(\bar{Q}_1\) 和 \(\bar{Q}_2\) 简化为通常的费米-狄拉克分布对数形式。
  - 在低温禁闭相（\(\Phi \rightarrow 0\)），夸克分布被强烈压制，体现了禁闭效应。

---

### 📊 参数表格（从原文提取）

#### **表3.1 rPNJL模型参数**
| 参数 | 数值 | 单位 | 说明 |
|------|------|------|------|
| \(m_{u,d}\) | 5.5 | MeV | u,d夸克流质量 |
| \(m_s\) | 183.468 | MeV | s夸克流质量 |
| \(\Lambda\) | 637.72 | MeV | 动量截断 |
| \(g_S \Lambda^2\) | 2.9147 | 无量纲 | 四费米子耦合 |
| \(g_D \Lambda^5\) | 5.968 | 无量纲 | 六费米子耦合 |
| \(g_1\) | \(2.193 \times 10^{-21}\) | MeV\(^{-8}\) | 八夸克耦合1 |
| \(g_2\) | \(-5.890 \times 10^{-22}\) | MeV\(^{-8}\) | 八夸克耦合2 |

#### **表3.2 胶子场有效势参数**
| 参数 | 数值 | 单位 | 说明 |
|------|------|------|------|
| \(T_0\) | 175 | MeV | 纯规范临界温度 |
| \(a_0\) | 6.75 | 无量纲 | |
| \(a_1\) | -9.8 | 无量纲 | |
| \(a_2\) | 0.26 | 无量纲 | |
| \(b_3\) | 0.805 | 无量纲 | |
| \(b_4\) | 7.55 | 无量纲 | |
| \(\kappa\) | 0.1 | 无量纲 | Vandermonde项系数 |

---

### 📚 相关补充材料（推荐参考文献）

以下文献与本模型构建和物理背景直接相关，建议进一步阅读：

1. **Osipov, A. A., Hiller, B., da Providência, J. (2006)** – 引入八夸克相互作用以稳定真空，为rPNJL模型提供理论依据。  
   *Phys. Lett. B 634, 48* [193]

2. **Bhattacharyya, A., et al. (2010)** – 研究八夸克相互作用对PNJL模型相图和热力学性质的影响。  
   *Phys. Rev. D 82, 014021* [194]

3. **Ilgenfritz, E. M., Kripfganz, J. (1985)** – 给出Polyakov圈有效势中Vandermonde项的Jacobian行列式表达式。  
   *Z. Phys. C 29, 79* [195]

4. **Bhattacharyya, A., et al. (2017)** – rPNJL模型参数化研究，使模型在零化学势下与格点QCD的手征和退禁闭临界温度一致。  
   *Phys. Rev. D 95, 054005* [187]

5. **Fukushima, K. (2004)** – PNJL模型的原始论文，详细推导了耦合到Polyakov圈的夸克分布函数。  
   *Phys. Lett. B 591, 277-284*

6. **Ratti, C., et al. (2006)** – PNJL模型在QCD相图研究中的应用，包含了\(\bar{Q}_1\)和\(\bar{Q}_2\)函数的详细讨论。  
   *Phys. Rev. D 73, 014019*

---

### 💎 总结
本报告提取了rPNJL模型在平均场近似下的核心公式组，包括：
- 拉氏量（式3.26）
- 胶子场有效势及其温度依赖（式3.27–3.29）
- 巨热力学势（式3.30），特别补充说明了\(\bar{Q}_1\)和\(\bar{Q}_2\)的显式表达式及其物理意义
- 夸克能隙方程（式3.31）

模型通过引入八夸克相互作用和Vandermonde项改进了传统PNJL模型，能更准确地描述QCD相图。公式中的\(\bar{Q}_1\)和\(\bar{Q}_2\)函数是耦合Polyakov圈效应的关键，体现了夸克禁闭对统计分布的修正。