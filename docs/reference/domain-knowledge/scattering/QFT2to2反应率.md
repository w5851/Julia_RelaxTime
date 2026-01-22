- QFT 中 **(2→2)** 反应率计算方法。分为三层：从最基础的 S-矩阵 → 微分截面 → 反应率（含热平均）。

---

# ✅ 1. 从 S-矩阵到散射截面

对于  

\[
a(p_1)+b(p_2)\to c(p_3)+d(p_4)
\]

**不含平均与求和的矩阵元记为 \(\mathcal{M}\)**。

## (1) 微分截面（质心系）

最常用的形式：  

\[
\frac{d\sigma}{d\Omega} = \frac{1}{64\pi^2 s}\cdot\frac{|\mathbf p_f|}{|\mathbf p_i|}\cdot\overline{|\mathcal M|^2}.
\]

若粒子质量不等：  

\[
|\mathbf p_i| = \frac{\lambda^{1/2}(s,m_a^2,m_b^2)}{2\sqrt{s}},\qquad
|\mathbf p_f| = \frac{\lambda^{1/2}(s,m_c^2,m_d^2)}{2\sqrt{s}},
\]

其中 Källén 函数：  

\[
\lambda(x,y,z)=x^2+y^2+z^2-2xy-2yz-2zx.
\]

**总截面**：  

\[
\sigma(s)=\int d\Omega\,\frac{d\sigma}{d\Omega}.
\]

---

# ✅ 2. “反应率” 的标准定义（常用在热场论、输运理论）

量纲：**速率/体积（per unit volume per time）**

最常见的定义：  

\[
R_{ab\to cd} = \int \frac{d^3p_1}{(2\pi)^3}\frac{d^3p_2}{(2\pi)^3}
f_a(p_1)f_b(p_2)\, v_{\text{rel}}\, \sigma(s).
\]

其中  
**Møller 相对速度**（洛伦兹不变）：  

\[
v_{\text{rel}}=\frac{\sqrt{(p_1\cdot p_2)^2 - m_a^2 m_b^2}}
{E_1E_2}.
\]

极其重要：  
这个 \(v_{\text{rel}}\) 正是保证整个反应率表达式 Lorentz-invariant 的关键。

---

# ✅ 3. 若系统处于热平衡（Boltzmann/FD/BE），可定义 **热平均反应率**

常用形式：  

\[
\langle\sigma v_{\text{rel}}\rangle
= \frac{
\displaystyle\int \frac{d^3p_1}{(2\pi)^3}\frac{d^3p_2}{(2\pi)^3}
f_a f_b\,
v_{\text{rel}}\,\sigma(s)
}{
\displaystyle\int \frac{d^3p_1}{(2\pi)^3}\frac{d^3p_2}{(2\pi)^3}
f_a f_b}.
\]

如果使用 Boltzmann 近似（如在早期宇宙或夸克散射中常用）可简化为一个单变量积分：  

\[
\langle \sigma v\rangle =
\frac{1}{8m_a^2m_b^2T\,K_2(m_a/T)K_2(m_b/T)}
\int_{s_{\text{min}}}^\infty ds\,
\sigma(s)\,(s-m_a^2-m_b^2)\sqrt{s}\,
K_1\!\left(\frac{\sqrt{s}}{T}\right).
\]

（这个公式在热暗物质、夸克散射等都大量使用。）

---

# ✅ 4. QFT 中更原始、更严谨的形式（Lorentz 不变测度）

反应率可以写成完全 Lorentz 不变的 4 维形式：  

\[
R = \int d\Pi_1 \,d\Pi_2 \,d\Pi_3 \,d\Pi_4\,
(2\pi)^4 \delta^{(4)}(p_1+p_2-p_3-p_4)\,
f_1 f_2\,
|\mathcal M|^2,
\]

其中  

\[
d\Pi_i = \frac{d^3p_i}{(2\pi)^3\,2E_i}.
\]

对最终态若求和则乘以 \((1\pm f_3)(1\pm f_4)\)（玻色增强/泡利阻塞）。

这是你在各向异性分布函数、非平衡动力学（RTA、Boltzmann 方程）里经常会用到的 **最底层公式**。

---

# ✅ 5. 若体系是各向异性的（如你前面讨论的 R/S 各向异性）

则：

* 你不能对 \(\phi^*\) 积分掉（角度依赖不能消去）
* \(\sigma(s,\theta^*,\phi^*)\) 会有额外角度依赖
* 反应率变成 **六维积分**

\[
R = \int \frac{d^3p_1}{(2\pi)^3}\frac{d^3p_2}{(2\pi)^3}
f_1(\mathbf p_1)\, f_2(\mathbf p_2)\,
v_{\text{rel}}\,
\int d\Omega^*\;\frac{d\sigma}{d\Omega^*}(s,\theta^*,\phi^*).
\]

其中 \( d\sigma/d\Omega^* \) 完全由 \( |\mathcal M|^2 \) 决定，包括介子极化函数、传播子角度依赖等。

你在前面问过的点仍然成立：  
**截面的角度依赖来自矩阵元 \(|M|^2\)**，若体系分布各向异性，则最终反应率必然也显式依赖这些角度。
