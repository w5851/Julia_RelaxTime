## 【提取报告：高雪艳博士论文第2.1节 NJL模型公式】

### 📍 源信息
- **文献**：高雪艳博士论文.pdf
- **位置**：第2章第1节 “2.1 NJL夸克模型”
- **提取深度**：公式 + 直接上下文
- **提取时间**：2024-01-20 15:00:00

---

### 🔢 目标公式组（按原文出现顺序排列）

#### **公式1：QCD对称性群**
**【类型】** 对称性群表示
**【内容】**
\[
S U_{V}(N_{f})\otimes S U_{A}(N_{f})\otimes U_{V}(1)\otimes U_{A}(1) \quad (2-1)
\]
**【上下文关联】**
- **引入文本**：“QCD的拉氏量所具有的对称性为：”
- **解释文本**：（紧接公式后说明） “\(S U_{V}(3)\) 和 \(U_{V}(1)\) 的对称性分别确保了同位旋和重子数守恒。”

---

#### **公式2：两味NJL拉氏量**
**【类型】** 拉格朗日密度
**【内容】**
\[
\mathcal{L} = \bar{q} (i\gamma^{\mu}\partial_{\mu} - \hat{m}_{0} + \gamma^{0}\hat{\mu})q + G\left[(\bar{q} q)^{2} + (\bar{q} i\gamma_{5}\tau q)^{2}\right] \quad (2-3)
\]
**【上下文关联】**
- **引入文本**：“对于两味体系，具有上述对称性的拉氏量可采用如下形式：”
- **解释文本**：（无直接后续解释，转入三味体系）

---

#### **公式3：三味NJL拉氏量（一般形式）**
**【类型】** 拉格朗日密度分解
**【内容】**
\[
\mathcal{L}_{NJL} = \mathcal{L}_0 + \mathcal{L}_4 + \mathcal{L}_6
\]
**【上下文关联】**
- **引入文本**：“对于三味体系，具有上述对称性的拉氏量可采用如下形式[125]：”
- **解释文本**：（后续分段解释各部分） “其中 \(\mathcal{L}_0\) 是自由夸克的拉氏量密度; \(\mathcal{L}_4\) 是 4 夸克相互作用... \(\mathcal{L}_6\) 是味道混合相 6 夸克相互作用...”

---

#### **公式4：夸克场定义**
**【类型】** 场变量定义
**【内容】**
\[
q = \begin{pmatrix} u \\ d \\ s \end{pmatrix} \quad (2-5)
\]
**【上下文关联】**
- **引入文本**：“夸克场为：”
- **解释文本**：（无直接后续解释）

---

#### **公式5：平均场近似下 \(\mathcal{L}_4\) 的化简**
**【类型】** 相互作用项化简
**【内容】**
\[
\begin{aligned} 
\mathcal{L}_4 &= G \left[ \left( \bar{q}\lambda^0 q \right) + \left( \bar{q}\lambda^3 q \right) + \left( \bar{q}\lambda^8 q \right) \right] \\
&= 2G \sum_{f=u,d,s} \left\{ 2 \langle \bar{q}_f q_f \rangle \langle \bar{q}_f q_f \rangle - \langle \bar{q}_f q_f \rangle^2 \right\} 
\end{aligned} \quad (2-6)
\]
**【上下文关联】**
- **引入文本**：“在平均场近似下不考虑 \(\pi\) 介子凝聚 \((\langle \bar{q}i\gamma_5\lambda^a q \rangle = 0)\), 忽略高阶涨落, \(\mathcal{L}_4\) 可化简为：”
- **解释文本**：（后续说明计算中用到的公式(2-7)）

---

#### **公式6：\(\mathcal{L}_6\) 的化简公式**
**【类型】** 行列式展开与平均场近似
**【内容】**
\[
\det(\bar{q}\Lambda q) = \sum_{i,j,k} \varepsilon_{ijk} (\bar{q}_i \Lambda q_i) (\bar{q}_d \Lambda q_j) (\bar{q}_s \Lambda q_k) \quad (2-8)
\]
\[
(\bar{q}_a\Lambda q_b)(\bar{q}_c\Lambda q_d) = \langle \bar{q}_a\Lambda q_b\rangle (\bar{q}_c\Lambda q_d) + (\bar{q}_a\Lambda q_b)\langle \bar{q}_c\Lambda q_d\rangle -\langle \bar{q}_a\Lambda q_b\rangle \langle \bar{q}_c\Lambda q_d\rangle \quad (2-9)
\]
**【上下文关联】**
- **引入文本**：“对于 \(\mathcal{L}_6\) 的化简需要用到以下公式：” 及 “在平均场近似下的公式：”
- **解释文本**：（后续给出味道守恒条件和最终化简结果）

---

#### **公式7：化简后的三味NJL拉氏量密度**
**【类型】** 平均场拉格朗日密度
**【内容】**
\[
\mathcal{L}_{NJL} = \bar{q}\left(i\gamma^\mu \partial_\mu -\hat{M}\right)q - 2G\sum_{i = u,d,s}\phi_i^2 +4K\phi_u\phi_d\phi_s \quad (2-12)
\]
**【上下文关联】**
- **引入文本**：“则三味NJL模型的拉氏量密度可以写为：”
- **解释文本**：（后续定义组分夸克质量 \(\hat{M}\)）

---

#### **公式8：能隙方程（组分夸克质量）**
**【类型】** 自洽方程
**【内容】**
\[
M_{i} = m_{i} - 4G\phi_{i} + 2K\phi_{j}\phi_{k}\quad i\neq j\neq k \quad (2-13)
\]
**【上下文关联】**
- **引入文本**：“其中 \(\hat{M} = \mathrm{diag}(M_u,M_d,M_s)\) 为组分夸克质量,是由于征对称性自发破缺产生的有效质量,其表达式为：”
- **解释文本**：“上式也被称为能隙方程,为计算真空中夸克凝聚,我们采用平均场近似。”

---

#### **公式9：非零化学势下的自由拉氏量**
**【类型】** 拉格朗日密度
**【内容】**
\[
\mathcal{L}_0 = \bar{q}\left(i\gamma^\mu \partial_\mu -\hat{m} +\gamma_0\mu\right)q \quad (2-14)
\]
**【上下文关联】**
- **引入文本**：“对于零化学势和非零化学势的NJL模型其推导过程相似。不失一般性,下面以非零化学势的NJL模型为例,则式(2- 4)中的 \(\mathcal{L}_0\) 变为”
- **解释文本**：“其中 \(\mu\) 为夸克化学势。”

---

#### **公式10：欧氏空间配分函数与作用量**
**【类型】** 配分函数与作用量表达式
**【内容】**
\[
Z = N^{\prime}\int [d\bar{q} ][dq]\exp \left[\int_0^\beta d\tau \int d^3 x\mathcal{L}\right] \quad (2-15)
\]
\[
S_{E} = \int_{0}^{\beta}d\tau \int d^{3}x\bar{q} S_{0}^{-1}q + \beta \left(2G\sum_{i = u,d,s}\phi_{i}^{2} - 4K\phi_{u}\phi_{d}\phi_{s}\right) \quad (2-16)
\]
**【上下文关联】**
- **引入文本**：“在欧氏空间配分函数是虚时作用量,指数函数关于场量的泛函积分：” 及 “作用量：”
- **解释文本**：“\(N^{\prime}\) 为归一化常数。” 及 （定义 \(S_0^{-1}\) 为夸克传播子的倒数）

---

#### **公式11：巨热力学势定义**
**【类型】** 热力学势定义
**【内容】**
\[
\Omega = \frac{-T}{V}\ln Z \quad (2-18)
\]
**【上下文关联】**
- **引入文本**：（紧接作用量后）“则巨热力学势可写为：”
- **解释文本**：（后续转入坐标变换和Matsubara频率）

---

### 🔢 目标公式组（接续上文）

#### **公式12：费米子场的傅里叶展开**
**【类型】** 场量展开式
**【内容】**
\[
q_{\alpha}(\tau ,\pmb {x}) = \frac{1}{\sqrt{V}}\sum_{n,p}e^{i(\pmb {p}\cdot \pmb {x} + \omega_{n}\tau)}\tilde{q}_{\alpha ;n}(\pmb {p}) \quad (2-19)
\]
**【上下文关联】**
- **引入文本**：“为便于计算,将时空坐标换为 \((\omega_{n},\pmb {p})\) ,费米子需要满足反周期边界条件 \(q_{\alpha}(0,\pmb {x}) = - q_{\alpha}(\beta ,\pmb {x})\) ,则有：”
- **解释文本**：其中 \(\omega_{n}\) 为费米子Matsubara频率,且 \(\omega_{n} = (2n + 1)\pi T\) 。

---

#### **公式13：配分函数的行列式表示**
**【类型】** 配分函数推导
**【内容】**
\[
Z_{0} = \left[\prod_{n} \int dq_{\alpha ;n}^{\dagger}(\pmb {p})dq_{\alpha ;n}(\pmb {p})\right]e^{S_{p}} = \operatorname{det}D \quad (2-20)
\]
\[
S_{p} = \sum_{n,p}iq_{\alpha ;n}^{\dagger}(\pmb {p})D_{\alpha \rho}q_{\rho ;n}(\pmb {p}) \quad (2-21)
\]
\[
D_{\alpha \rho} = \beta \left[(-i\omega_{n} + \mu) - \gamma^{0}\gamma \cdot \pmb {p} - M\gamma^{0} + \mu\right] = \beta S^{-1}(i\omega_{n},\pmb {p}) \quad (2-22)
\]
**【上下文关联】**
- **引入文本**：“将 \(q_{\alpha}(\tau ,x)\) 带入到巨热力学势式(2- 18)中可得：” （注：原文为“带入到巨热力学势式(2- 18)”，但实际推导配分函数 \(Z_0\)）
- **解释文本**：（后续利用行列式性质求迹）

---

#### **公式14：求迹与行列式运算**
**【类型】** 迹与行列式关系
**【内容】**
\[
\sum_{n,p}\operatorname{Tr}\ln\left(\frac{1}{T}S^{-1}\left(i\omega_{n},\pmb{p}\right)\right)=\sum_{n,p}\operatorname{Tr}_{c}\operatorname{Tr}_{f}\operatorname{Tr}_{\mathrm{Dirac}}\ln\left(\frac{1}{T}S^{-1}\left(i\omega_{n},\pmb{p}\right)\right) \quad (2-23)
\]
\[
=\sum_{n,p}\operatorname{Tr}_{c}\operatorname{Tr}_{f}\ln\left[\operatorname{Det}_{\mathrm{Dirac}}\left(\frac{1}{T}S^{-1}\left(i\omega_{n},\pmb{p}\right)\right)\right]
\]
**【上下文关联】**
- **引入文本**：“利用 \(\operatorname {Tr}\ln D = \ln (\operatorname {Det}D)\) 可以得到：”
- **解释文本**：\(\operatorname {Tr}_{\mathrm{Dirac}}\) , \(\operatorname {Tr}_{f}\) 和 \(\operatorname {Tr}_{c}\) 分别是对狄拉克空间、味空间和色空间求迹, \(\operatorname {Det}_{\mathrm{Dirac}}\) 表示对狄拉克空间求行列式。

---

#### **公式15：狄拉克行列式的具体形式**
**【类型】** 行列式计算
**【内容】**
\[
\ln \left[\operatorname {Det}_{\mathrm{Dirac}}\left(\frac{1}{T} S^{-1}\left(i\omega_{n},\pmb {p}\right)\right)\right] = \ln \left[\frac{1}{T^{2}}\left(\omega_{n}^{2} + (\hat{E} -\hat{\mu})^{2}\right)\right] + \ln \left[\frac{1}{T^{2}}\left(\omega_{n}^{2} + (\hat{E} +\hat{\mu})^{2}\right)\right] \quad (2-24)
\]
**【上下文关联】**
- **引入文本**：（紧接公式(2-23)）
- **解释文本**：其中有效能量矩阵为 \(\hat{E} = \operatorname {diag}\left(E_{u},E_{d},E_{s}\right)\) ,有效能量的表达式为：\(E_{i} = \sqrt{M_{i}^{2} + \pmb{p}^{2}}\) （式2-25）。

---

#### **公式16：有效能量定义**
**【类型】** 能量-动量关系
**【内容】**
\[
E_{i} = \sqrt{M_{i}^{2} + \pmb{p}^{2}} \quad (2-25)
\]
**【上下文关联】**
- **引入文本**：（紧接公式2-24解释）
- **解释文本**：\(T\) 为温度, \(\pmb{p}\) 为动量。

---

#### **公式17：展开后的求和表达式**
**【类型】** 分味求和
**【内容】**
\[
\begin{aligned}
& \sum_{n}\mathrm{Tr}_{c}\mathrm{Tr}_{f}\ln \left[\mathrm{Det}_{\mathrm{Diac}}\left(\frac{1}{T} S^{-1}\left(i\omega_{n},\pmb {p}\right)\right)\right] \\
& = N_{c}\sum_{n}\left\{\ln \left[\frac{1}{T^{2}}\left(\omega_{n}^{2} + (E_{u} - \mu_{u})^{2}\right)\right] + \ln \left[\frac{1}{T^{2}}\left(\omega_{n}^{2} + (E_{u} + \mu_{u})^{2}\right)\right]\right. \\
& \qquad \left. + \ln \left[\frac{1}{T^{2}}\left(\omega_{n}^{2} + (E_{d} - \mu_{d})^{2}\right)\right] + \ln \left[\frac{1}{T^{2}}\left(\omega_{n}^{2} + (E_{d} + \mu_{d})^{2}\right)\right]\right. \\
& \qquad \left. + \ln \left[\frac{1}{T^{2}}\left(\omega_{n}^{2} + (E_{s} - \mu_{s})^{2}\right)\right] + \ln \left[\frac{1}{T^{2}}\left(\omega_{n}^{2} + (E_{s} + \mu_{s})^{2}\right)\right] \right\}
\end{aligned} \quad (2-26)
\]
**【上下文关联】**
- **引入文本**：（紧接公式2-24，原文中公式编号为(2-25)，但根据上下文应为(2-26)）
- **解释文本**：（后续利用Matsubara频率求和公式进行求和）

---

#### **公式18：Matsubara频率求和公式**
**【类型】** 求和恒等式
**【内容】**
\[
T\sum_{n}\ln\left(\frac{1}{T^{2}}\left(\omega_{n}^{2} + \lambda_{k}^{2}\right)\right) = \lambda_{k} + 2T\ln\left(1 + e^{-\lambda_{k}/T}\right) \quad (2-27)
\]
**【上下文关联】**
- **引入文本**：“利用公式”
- **解释文本**：并对动量积分有：（得到 \(\ln Z_0\) 的表达式）

---

#### **公式19：配分函数对数表达式**
**【类型】** 配分函数推导结果
**【内容】**
\[
\ln Z_{0} = 2N_cV\int \frac{d^3p}{(2\pi)^3} \sum_{i=u,d,s} \left[ E_i + \ln\left(1 + e^{-(E_i-\mu_i)/T}\right) + \ln\left(1 + e^{-(E_i+\mu_i)/T}\right) \right] \quad (2-28)
\]
**【上下文关联】**
- **引入文本**：“并对动量积分有：”
- **解释文本**：将上式带入式（2-18）可以得到NJL模型的巨热力学势为：（引出公式2-29）

---

#### **公式20：NJL模型巨热力学势（原始形式）**
**【类型】** 热力学势积分表达式
**【内容】**
\[
\Omega = -2N_c \int \frac{d^3p}{(2\pi)^3} \sum_{i=u,d,s} \left[ E_i + T\ln\left(1 + e^{-(E_i-\mu_i)/T}\right) + T\ln\left(1 + e^{-(E_i+\mu_i)/T}\right) \right] + 2G \sum_{i=u,d,s} \phi_i^2 - 4K\phi_u\phi_d\phi_s \quad (2-29)
\]
**【上下文关联】**
- **引入文本**：“将上式带入式（2-18）可以得到NJL模型的巨热力学势为：”
- **解释文本**：积分中的第二项和第三项还可以化简，以第二项为例：（开始推导化简形式）

---

#### **公式21：对数项的化简过程（分部积分）**
**【类型】** 积分化简步骤
**【内容】**
\[
\begin{aligned}
& \int \frac{d^3p}{(2\pi)^3} T \ln\left(1 + e^{-(E_i-\mu_i)/T}\right) \\
& \Rightarrow \frac{T}{2\pi^2} \int p^2 dp \ln\left(1 + e^{-(E_i-\mu_i)/T}\right) \quad (\text{球坐标系}) \\
& \Rightarrow \frac{T}{6\pi^2} \left[ p^3 \ln\left(1 + e^{-(E_i-\mu_i)/T}\right) \Big|_0^\infty - \int_0^\infty \frac{d}{dp} \ln\left(1 + e^{-(E_i-\mu_i)/T}\right) \cdot p^3 \, dp \right] \quad (\text{分部积分}) \\
& \Rightarrow \frac{1}{6\pi^2} \int_0^\infty \frac{1}{e^{(E_i-\mu_i)/T} + 1} \frac{p^4}{E_i} \, dp \quad (\text{化简}) \\
& \Rightarrow \int \frac{d^3p}{(2\pi)^3} \frac{p^2}{3E_i} \frac{1}{e^{(E_i-\mu_i)/T} + 1} \quad (\text{直角坐标系})
\end{aligned} \quad (2-30)
\]
**【上下文关联】**
- **引入文本**：“积分中的第二项和第三项还可以化简，以第二项为例：”
- **解释文本**：式（2-29）中的第三项可以用同样的方法化简为：（公式2-31）

---

#### **公式22：第三项的化简结果**
**【类型】** 积分化简结果
**【内容】**
\[
\int \frac{d^3p}{(2\pi)^3} \frac{p^2}{3E_i} \frac{1}{e^{(E_i+\mu_i)/T} + 1} \quad (2-31)
\]
**【上下文关联】**
- **引入文本**：式（2-29）中的第三项可以用同样的方法化简为：
- **解释文本**：将式（2-30）和式（2-31）带入式（2-29）可以得到：（公式2-32）

---

#### **公式23：巨热力学势的化简形式（含分布函数）**
**【类型】** 热力学势最终表达式
**【内容】**
\[
\begin{aligned}
\Omega_{\text{NJL}} = & -2N_c \sum_{i=u,d,s} \int \frac{d^3p}{(2\pi)^3} \frac{p^2}{3E_i} \left[ n(E_i - \mu_i) + \bar{n}(E_i + \mu_i) \right] \\
& -2N_c \sum_{i=u,d,s} \int^\Lambda \frac{d^3p}{(2\pi)^3} E_i \\
& + 2G \left( \phi_u^2 + \phi_d^2 + \phi_s^2 \right) - 4K \phi_u \phi_d \phi_s + C
\end{aligned} \quad (2-32)
\]
**【上下文关联】**
- **引入文本**：将式（2-30）和式（2-31）带入式（2-29）可以得到：
- **解释文本**：其中夸克和反夸克的分布函数分别为：（公式2-33）

---

#### **公式24：夸克与反夸克分布函数**
**【类型】** 统计分布函数
**【内容】**
\[
n(E_i - \mu_i) = \frac{1}{e^{(E_i - \mu_i)/T} + 1}, \quad \bar{n}(E_i + \mu_i) = \frac{1}{e^{(E_i + \mu_i)/T} + 1} \quad (2-33)
\]
**【上下文关联】**
- **引入文本**：其中夸克和反夸克的分布函数分别为：
- **解释文本**：式（2-32）第三项的动量积分是发散的，所以需要使用重整化方法去掉发散部分，本文采用三动量截断减除发散，截断动量为 \(\Lambda\) ，常数 \(C\) 的引入是为了保证系统在真空 \((T = \mu = 0)\) 时的压强为零， \(N_{c} = 3\) 是夸克色自由度的个数。

---

#### **公式25：稳定条件（极值条件）**
**【类型】** 变分方程
**【内容】**
\[
\frac{\partial \Omega_{\text{NJL}}}{\partial \phi_{u}} = \frac{\partial \Omega_{\text{NJL}}}{\partial \phi_{d}} = \frac{\partial \Omega_{\text{NJL}}}{\partial \phi_{s}} = 0 \quad (2-34)
\]
**【上下文关联】**
- **引入文本**：三味 NJL 模型中有三个独立变量 \(\phi_{u}\) 、 \(\phi_{d}\) 和 \(\phi_{s}\) ，对于给定的条件 \((T, \mu)\) 或 \((T, \rho)\) ，它们的值通过热力学势取极值来确定：
- **解释文本**：在手征极限 \((m_{0} = 0)\) 时，夸克手征凝聚 \(\phi\) 是手征对称性的序参量...（后续讨论相变判断标准）

---

### 📖 本节后续物理内容总结
本节后续部分（公式2-19至2-34）完成了NJL模型在有限温度 \(T\) 和化学势 \(\mu\) 下的巨热力学势 \(\Omega_{\text{NJL}}\) 的推导，并给出了稳定条件（能隙方程）。关键步骤包括：
1.  **动量空间表示**：通过傅里叶展开和Matsubara频率求和，将配分函数表示为动量积分。
2.  **巨热力学势表达式**：
    - **原始形式**（式2-29）：包含真空能、夸克/反夸克热激发贡献以及凝聚项。
    - **化简形式**（式2-32）：通过对数项的分部积分，引入分布函数 \(n\) 和 \(\bar{n}\)，使物理意义更清晰（热激发部分表现为粒子数加权积分）。
3.  **发散处理**：真空积分发散，采用三动量截断 \(\Lambda\) 正规化，并引入常数 \(C\) 保证真空压强为零。
4.  **稳定条件**（式2-34）：通过巨热力学势对三个夸克凝聚 \(\phi_u, \phi_d, \phi_s\) 的变分极值，确定系统的基态（即能隙方程，与式2-13在有限温密下的推广对应）。

## 🔧 NJL模型参数表格（补充）

根据您提供的原文片段（应出自高雪艳博士论文后续章节），现提取并整理 **NJL模型的标准参数表**如下：

### ⚙️ **NJL模型参数表（文献[141]拟合结果）**

| 参数符号 | 物理意义 | 数值 | 单位 | 备注 |
|----------|----------|------|------|------|
| \(m_u\) | u夸克流质量 | 5.5 | MeV | 与 \(m_d\) 相同，保持同位旋对称性 |
| \(m_d\) | d夸克流质量 | 5.5 | MeV | 与 \(m_u\) 相同 |
| \(m_s\) | 奇异夸克流质量 | 140.7 | MeV | 拟合K介子质量得到 |
| \(\Lambda\) | 三动量截断 | 602.3 | MeV | 用于正规化发散的真空积分 |
| \(G\Lambda^2\) | 四费米子耦合强度（无量纲组合） | 1.835 | 无量纲 | 拟合π介子质量和衰变常数 \(f_\pi\) |
| \(K\Lambda^5\) | 六费米子耦合强度（无量纲组合） | 12.36 | 无量纲 | 拟合η′介子质量，破缺 \(U_A(1)\) 对称性 |

### 📝 **参数说明与拟合依据**
- **拟合目标**：该组参数通过拟合以下实验/格点数据确定：
  - π介子质量 \(m_\pi = 135.0 \text{ MeV}\)
  - K介子质量 \(m_K = 497.7 \text{ MeV}\)
  - η′介子质量 \(m_{\eta'} = 957.8 \text{ MeV}\)
  - π介子衰变常数 \(f_\pi = 92.4 \text{ MeV}\)
- **参数个数**：虽然列出6个参数，但 \(m_u = m_d\) 视为一个独立参数，故常称为“五个参数”。
- **使用约定**：在后续计算（如能隙方程、热力学量计算）中，这些参数通常保持不变，作为模型基准参数集。

### 🔗 **与之前提取公式的关联**
- 该参数表对应**2.1节NJL模型**中出现的：
  - 流质量 \(m_i\)（式2-13）
  - 耦合常数 \(G, K\)（式2-12）
  - 截断 \(\Lambda\)（式2-32中的积分上限）
- 在数值计算中，使用该参数集可求解：
  - 能隙方程（式2-13）→ 得到组分夸克质量 \(M_i\)
  - 巨热力学势（式2-32）→ 计算热力学量
  - 稳定条件（式2-34）→ 确定基态凝聚 \(\phi_i\)

