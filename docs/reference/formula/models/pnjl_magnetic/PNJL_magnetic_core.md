## 【提取报告：外磁场下的PNJL模型相关公式】

### 📍 源信息
- **文献**：高雪艳博士论文《强相互作用物质相变与重子数涨落的研究》
- **位置**：第五章“磁场对QCD一阶相变的影响”，第5.1节“外磁场下的PNJL模型”
- **提取时间**：2026-2-20

---

## 🔢 核心公式组

### **公式1：外磁场下的PNJL模型拉格朗日密度**
**【标签】** 式(5-1)
**【类型】** 拉格朗日密度
**【内容】**
\[
\begin{aligned}
\mathcal{L}=&\bar{q}\left(i\gamma^{\mu}D_{\mu}+\gamma_{0}\hat{\mu}-\hat{m}_{0}\right)q \\
&+G\sum_{k=0}^{8}\left[\left(\bar{q}\lambda_{k}q\right)^{2}+\left(\bar{q}i\gamma_{5}\lambda_{k}q\right)^{2}\right] \\
&-K\left[\operatorname{det}_{f}\left(\bar{q}\left(1+\gamma_{5}\right)q\right)+\operatorname{det}_{f}\left(\bar{q}\left(1-\gamma_{5}\right)q\right)\right] \\
&-U(\Phi[A],\bar{\Phi}[A],T)-\frac{1}{4}F^{\mu\nu}F_{\mu\nu}
\end{aligned}
\]

**【上下文关联】**
- **引入文本**：“在外部磁场存在的情况下，夸克物质的拉格朗日密度具有如下形式”
- **关键修改**：协变导数增加电矢势项 \(D^{\mu} = \partial^{\mu} - iA^{\mu} - iq_{i}A_{EM}^{\mu}\)
- **磁场设定**：沿z方向静态恒定磁场 \(A_{EM}^{\mu} = \delta^{\mu 2}x_{1}B\)

---

### **公式2：磁场下PNJL模型的巨热力学势**
**【标签】** 式(5-2)
**【类型】** 巨热力学势表达式
**【内容】**
\[
\Omega = \sum_{f = u,d,s}\left(\Omega_{f}^{0} + \Omega_{f}^{T}\right) + 2G\left(\phi_{u}^{2} + \phi_{d}^{2} + \phi_{s}^{2}\right) -4K\phi_{u}\phi_{d}\phi_{s} + U(\Phi ,\bar{\Phi},T)
\]

---

### **公式3：真空项 \(\Omega_{f}^{0}\)**
**【标签】** 式(5-3)、(5-7)
**【类型】** 真空贡献的积分表达式
**【内容】**
\[
\Omega_{f}^{0} = -N_{c}\frac{|q_{f}|eB}{2\pi}\sum_{n = 0}^{\infty}\alpha_{n}\int_{-\infty}^{\infty}\frac{dp_{z}}{2\pi} E_{f,n}
\]

**平滑截断正则化形式**（式5-7）：
\[
\Omega_{f}^{0} = -N_{c}\frac{|q_{f}|eB}{2\pi}\sum_{n = 0}^{\infty}\alpha_{n}\int_{-\infty}^{\infty}\frac{dp_{z}}{2\pi} f_{\lambda}^{2}(p_{f})E_{f,n}
\]
\[
f_{\Lambda}(p) = \sqrt{\frac{\Lambda^{2N}}{\Lambda^{2N} + p^{2N}}}
\]
其中 \(N=10\) 在数值计算中采用。

---

### **公式4：温度相关项 \(\Omega_{f}^{T}\)**
**【标签】** 式(5-4)
**【类型】** 温度贡献的积分表达式
**【内容】**
\[
\Omega_{f}^{T} = -T\frac{|q_{f}|eB}{2\pi}\sum_{n = 0}^{\infty}\alpha_{n}\int_{-\infty}^{+\infty}\frac{dp_{z}}{2\pi}\left(Z_{f}^{+} + Z_{f}^{-}\right)
\]

**其中：**
\[
\begin{aligned}
Z_{f}^{+} &= \ln \left(1 + 3\Phi e^{-\frac{E_{f,n} - \mu_{f}}{T}} + 3\bar{\Phi} e^{-2\frac{E_{f,n} - \mu_{f}}{T}} + e^{-3\frac{E_{f,n} - \mu_{f}}{T}}\right) \\
Z_{f}^{-} &= \ln \left(1 + 3\bar{\Phi} e^{-\frac{E_{f,n} + \mu_{f}}{T}} + 3\Phi e^{-2\frac{E_{f,n} + \mu_{f}}{T}} + e^{-3\frac{E_{f,n} + \mu_{f}}{T}}\right)
\end{aligned}
\]

---

### **公式5：准粒子能量（朗道能级）**
**【标签】** 式(5-6)
**【类型】** 能量本征值
**【内容】**
\[
E_{f,n} = \left(2n|q_{f}|B + p_{z}^{2} + M_{f}^{2}\right)^{1 / 2}
\]
其中 \(n = 0,1,2,\ldots\) 表示朗道能级指标。

---

### **公式6：数密度公式**
**【标签】** 式(5-10)、(5-11)
**【类型】** 夸克数密度表达式
**【内容】**
\[
\rho_{f} = \frac{\partial\Omega}{\partial\mu_{f}} = \sum_{n = 0}^{\infty}\rho_{f,n}
\]
\[
\rho_{f,n} = 3\frac{q_{f}B}{2\pi}\alpha_{n}\int_{-\infty}^{\infty}\frac{\mathrm{d}p_{z}}{2\pi}\left[f_{f,n}^{+} - f_{f,n}^{-}\right]
\]

**其中分布函数：**
\[
\begin{aligned}
f_{f,n}^{+} &= \frac{\Phi e^{-\frac{E_{f,n} - \mu_{f}}{T}} + 2\bar{\Phi} e^{-2\frac{E_{f,n} - \mu_{f}}{T}} + e^{-3\frac{E_{f,n} - \mu_{f}}{T}}}{\left(1 + 3\Phi e^{-\frac{E_{f,n} - \mu_{f}}{T}} + 3\bar{\Phi} e^{-2\frac{E_{f,n} - \mu_{f}}{T}} + e^{-3\frac{E_{f,n} - \mu_{f}}{T}}\right)} \\
f_{f,n}^{-} &= \frac{\bar{\Phi}e^{-\frac{E_{f,n} - \mu_{f}}{T}} + 2\bar{\Phi} e^{-2\frac{E_{f,n} + \mu_{f}}{T}} + e^{-3\frac{E_{f,n} + \mu_{f}}{T}}}{\left(1 + 3\Phi e^{-\frac{E_{f,n} + \mu_{f}}{T}} + 3\bar{\Phi} e^{-2\frac{E_{f,n} + \mu_{f}}{T}} + e^{-3\frac{E_{f,n} + \mu_{f}}{T}}\right)}
\end{aligned}
\]

---

### **公式7：低温近似下的数密度**
**【标签】** 式(5-13)
**【类型】** 近似表达式
**【内容】**
\[
\rho_{f,n}\approx 3\frac{q_{f}B}{2\pi^{2}}\alpha_{n}\int_{0}^{\infty}\frac{\mathrm{d}p_{z}}{(2\pi)}\left[\frac{1}{1 + e^{3\left(\frac{E_{f,n} - \mu_{f}}{T}\right)}}\right]
\]

---

### **公式8：磁场依赖的耦合常数**
**【标签】** 式(5-9)
**【类型】** 经验公式（用于描述IMC效应）
**【内容】**
\[
G(B) = G_{0}\frac{1 + a\zeta^{2} + b\zeta^{3}}{1 + c\zeta^{2} + d\zeta^{4}}
\]
其中 \(\zeta = eB / \Lambda_{QCD}^{2}\)，\(\Lambda_{QCD} = 300\mathrm{MeV}\)

---

## 🔧 变量说明表

| 符号 | 物理意义 | 单位/量纲 | 备注 |
|------|----------|-----------|------|
| \(\Omega\) | 巨热力学势 | MeV⁴ |  |
| \(\Omega_{f}^{0}\) | 真空贡献项 | MeV⁴ | 每味夸克 |
| \(\Omega_{f}^{T}\) | 温度贡献项 | MeV⁴ | 每味夸克 |
| \(q_{f}\) | 夸克电荷 | e | u:+2/3, d:-1/3, s:-1/3 |
| \(eB\) | 磁场强度 | MeV² |  |
| \(B\) | 磁感应强度 | G | 1 G ≈ 1.95×10⁻²⁰ GeV² |
| \(n\) | 朗道能级指标 | 无量纲 | n=0,1,2,... |
| \(\alpha_{n}\) | 简并度因子 | 无量纲 | n=0时为1，n≠0时为2 |
| \(E_{f,n}\) | 准粒子能量 | MeV | 依赖于朗道能级 |
| \(M_{f}\) | 组分夸克质量 | MeV | 由于征对称性破缺产生 |
| \(\Phi,\bar{\Phi}\) | Polyakov圈及其共轭 | 无量纲 | 禁闭序参量 |
| \(U(\Phi,\bar{\Phi},T)\) | Polyakov圈有效势 | MeV⁴ | 对数形式 |
| \(G\) | 四费米子耦合常数 | MeV⁻² | 磁场依赖 |
| \(K\) | 六费米子耦合常数 | MeV⁻⁵ | t'Hooft相互作用 |
| \(\phi_{i}\) | 夸克凝聚 | MeV | i=u,d,s |
| \(f_{\Lambda}(p)\) | 截断函数 | 无量纲 | 平滑截断用 |
| \(\Lambda\) | 动量截断 | MeV | 602.3 MeV |
| \(\Lambda_{QCD}\) | QCD能标 | MeV | 300 MeV |
| \(Z_{f}^{\pm}\) | 统计分布的对数形式 | 无量纲 | 包含Polyakov圈耦合 |
| \(f_{f,n}^{\pm}\) | 夸克/反夸克分布函数 | 无量纲 | 含Polyakov圈修正 |

---

## 📊 参数表格

**表5-1：磁场依赖耦合常数G的相关参数**
| 参数 | 数值 | 单位 | 说明 |
|------|------|------|------|
| a | 0.108805 | 无量纲 | 式(5-9)参数 |
| b | -1.0133 × 10⁻⁴ | 无量纲 | 式(5-9)参数 |
| c | 0.02228 | 无量纲 | 式(5-9)参数 |
| d | 1.84558 × 10⁻⁴ | 无量纲 | 式(5-9)参数 |
| \(G_{0}\) | - | MeV⁻² | 零磁场时的耦合常数 |
| \(\Lambda_{QCD}\) | 300 | MeV | QCD能标 |

**PNJL模型固定参数（第二章）**
| 参数 | 数值 | 单位 | 说明 |
|------|------|------|------|
| \(\Lambda\) | 602.3 | MeV | 三动量截断 |
| \(G\Lambda^{2}\) | 1.835 | 无量纲 | 四费米子耦合 |
| \(K\Lambda^{5}\) | 12.56 | 无量纲 | 六费米子耦合 |
| \(m_{u}\) | 5.5 | MeV | u夸克流质量 |
| \(m_{d}\) | 5.5 | MeV | d夸克流质量 |
| \(m_{s}\) | 140.7 | MeV | s夸克流质量 |

---

## 💡 数值实现建议

1. **朗道能级求和**：
   - 对于给定的化学势 \(\mu_{f}\)，填充的最高朗道能级为 \(n_{\mathrm{max}} = \mathrm{Floor}\left(\frac{\mu_{f}^{2} - M_{f}^{2}}{2|q_{f}|B}\right)\)
   - 求和需包含足够多的能级以确保收敛

2. **真空项处理**：
   - 使用平滑截断函数 \(f_{\Lambda}(p)\) 避免截断不连续性
   - 建议采用N=10以获得良好的平滑效果

3. **低温数值处理**：
   - 低温时（\(T<20\) MeV），\(\Phi \approx \bar{\Phi}\approx 0\)，可简化分布函数
   - 需注意数值积分在费米面附近的奇异性

4. **多重一阶相变**：
   - 低温高磁场下可能出现多个一阶相变结构
   - 求解能隙方程时需小心处理多解情况

5. **IMC效应实现**：
   - 需采用磁场依赖的耦合常数 \(G(B)\) 式(5-9)
   - 参数a,b,c,d已调至重现格点QCD的IMC效应

---

## 🔗 相关公式链

1. **能隙方程**：式(2-55) - 夸克动力学质量 \(M_i = m_i - 4G\phi_i + 2K\phi_j\phi_k\)
2. **稳定平衡条件**：式(2-63) \(\frac{\partial\Omega_{PNJL}}{\partial\phi_u} = \frac{\partial\Omega_{PNJL}}{\partial\phi_d} = \frac{\partial\Omega_{PNJL}}{\partial\phi_s} = \frac{\partial\Omega_{PNJL}}{\partial\Phi} = \frac{\partial\Omega_{PNJL}}{\partial\bar{\Phi}} = 0\)
3. **Polyakov圈有效势**：式(2-50)、(2-51) - 对数形式有效势
4. **重子数密度**：\(\rho_B = (\rho_u + \rho_d + \rho_s)/3\)