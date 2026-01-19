# 介子质量计算（RPA 传播子极点）

本文汇总 PNJL/RPA 框架下介子（π、K、η、η′等）质量的计算公式与步骤，并与项目中既有极化函数/传播子实现保持一致。极化函数记号与 [../Polarization_极化函数byB0.md](../Polarization_极化函数byB0.md) 和 [../Propagator_传播子byPolarization.md](../Propagator_传播子byPolarization.md) 对齐。

## 1. 介子传播子公式

### 1.1 π 与 K 介子传播子
在随机相位近似（RPA）下，赝标量（P）或标量（S）介子传播子为：

$$
\mathcal{D}_{\pi} = \frac{2K_{3}^{\pm}}{1 - 4K_{3}^{\pm}\Pi_{u\bar{u}}^{P(S)}(p_{0},\bar{p})},
\quad
\mathcal{D}_{K} = \frac{2K_{4}^{\pm}}{1 - 4K_{4}^{\pm}\Pi_{u\bar{s}}^{P(S)}(p_{0},\bar{p})}
$$

其中：
- $K_{3}^{\pm}, K_{4}^{\pm}$ 为有效耦合常数（参见 [../K_有效耦合常数byA.md](../K_有效耦合常数byA.md)）。
- $\Pi_{ff'}^{P(S)}(p_{0},\mathbf{k})$ 为夸克-反夸克极化函数（参见 [../Polarization_极化函数byB0.md](../Polarization_极化函数byB0.md)）。

### 1.2 η 与 η′ 介子传播子（含混合）
赝标量介子 η、η′（以及标量 σ、σ′）存在单态/八重态混合，其传播子结构为：

$$
\mathcal{D} = 2\frac{\det K}{M_{00}M_{88} - M_{08}^2} 
\left( M_{00}\bar{\psi}\lambda_0\psi \cdot \bar{\psi}'\lambda_0\psi' 
+ M_{08}\bar{\psi}\lambda_0\psi \cdot \bar{\psi}'\lambda_8\psi'
+ M_{80}\bar{\psi}\lambda_8\psi \cdot \bar{\psi}'\lambda_0\psi' 
+ M_{88}\bar{\psi}\lambda_8\psi \cdot \bar{\psi}'\lambda_8\psi' \right)
$$

矩阵元：

$$
\begin{aligned}
M_{00} &= K_0^+ -\frac{4}{3}\det K(\Pi_{u\bar{u}} + 2\Pi_{s\bar{s}}) \\
M_{08} &= K_{08}^+ -\frac{4}{3}\sqrt{2}\det K(\Pi_{u\bar{u}} - \Pi_{s\bar{s}}) \\
M_{88} &= K_8^+ -\frac{4}{3}\det K(2\Pi_{u\bar{u}} + \Pi_{s\bar{s}}) \\
\det K &= K_0^+ K_8^+ -K_{08}^{+2}
\end{aligned}
$$

## 2. 介子质量方程（零动量极点）
介子质量定义为零动量传播子极点。

### 2.1 π 与 K 介子

$$
1 - 4K_{3}^{+}\Pi_{u\bar{u}}^{P}(M_{\pi} + i\frac{\Gamma}{2}, 0) = 0
$$

$$
1 - 4K_{4}^{+}\Pi_{u\bar{s}}^{P}(M_{K} + i\frac{\Gamma}{2}, 0) = 0
$$

其中 $\Gamma$ 为衰变宽度（虚部）。当介子进入共振态（Mott 分解）时，$\Gamma$ 由 0 变为非零。

### 2.2 η 与 η′ 介子

$$
M_{\eta}^{-1}(m_{\eta},\vec{0}) = 0, \quad
M_{\eta'}^{-1}(m_{\eta'},\vec{0}) = 0
$$

$$
\begin{aligned}
M_{\eta}^{-1} &= M_{00} + M_{88} - \sqrt{(M_{00} - M_{88})^{2} + 4M_{08}^{2}} \\
M_{\eta'}^{-1} &= M_{00} + M_{88} + \sqrt{(M_{00} - M_{88})^{2} + 4M_{08}^{2}}
\end{aligned}
$$

## 3. 计算流程（与工程实现对齐）
1. 输入 $T,\mu$ 与模型参数（$G,K,m_{0f},\Lambda$ 等）。
2. 在平均场模块求解动力学夸克质量 $M_f = m_{0f} - 2G\sigma_f$。
3. 调用极化函数模块计算 $\Pi_{ff'}^{P(S)}(p_0,\mathbf{0})$。
4. 求解介子质量：
   - π/K：在复平面搜索极点方程根，实部为 $M$、虚部为 $\Gamma/2$。
   - η/η′：先构造 $M_{00}, M_{08}, M_{88}$，再解极点方程。
5. 输出各介子质量 $M$ 与宽度 $\Gamma$（若适用）。

## 4. 注意事项
- 极化函数含柯西主值积分，会产生虚部；相关细节见 [../Polarization_极化函数byB0.md](../Polarization_极化函数byB0.md)。
- Mott 点满足 $M_M = M_{q_1}+M_{q_2}$，宽度开始非零。
- 标量介子通道与赝标量形式类似，但需使用标量耦合与标量极化函数。

### 4.1 复数能量输入的物理意义与实现现状
介子质量方程中常用 $p_0 = M + i\frac{\Gamma}{2}$ 的复数能量形式，其物理意义与数值处理说明见：
[../../domain-knowledge/relaxtime/介子质量计算中的复数输入与物理意义.md](../../domain-knowledge/relaxtime/%E4%BB%8B%E5%AD%90%E8%B4%A8%E9%87%8F%E8%AE%A1%E7%AE%97%E4%B8%AD%E7%9A%84%E5%A4%8D%E6%95%B0%E8%BE%93%E5%85%A5%E4%B8%8E%E7%89%A9%E7%90%86%E6%84%8F%E4%B9%89.md)。

实现检查结论：当前极化函数接口以实数能量为输入并返回实部/虚部，尚未直接支持复数 $p_0$ 作为参数。具体表现为：
- `PolarizationAniso.polarization_aniso` 的 `k0` 参数为 `Float64`；
- `PolarizationCache.PolarizationKey` 的 `k0` 字段为 `Float64`；
- `OneLoopIntegrals.B0` 的 `λ` 也按实数路径计算。

若需要在复平面上寻根，应在调用层通过“实部/虚部方程”或显式复数封装处理，并同步扩展上述接口类型。