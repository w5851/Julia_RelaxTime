# 夸克物质的输运系数（弛豫时间近似）

在弛豫时间近似下，基于相对论动力学，夸克物质的剪切粘滞系数、体粘滞系数和电导率表达式如下。

本文档聚焦 `eta`、`zeta`、`sigma`。`kappa_BB / kappa_QQ / kappa_SS / lambda` 已单独整理到关联文档中，避免与静态涨落中的 `\kappa\sigma^2` 记号混淆；其中 `kappa` 家族现在已用 Das 2022 的扩散矩阵框架重新锚定。需要特别注意的是，关联文档现已区分标准 `lorenz_number = lambda / (sigma * T)` 与 legacy `lorentz_legacy = lambda / (sigma / T)`，不要把两者混写成同一个量。

## 相关文档

- 弛豫时间的求和与过程键值对应： [RelaxationTime_FromAverageScatteringRate.md](RelaxationTime_FromAverageScatteringRate.md)
- 流算符向量在夸克味空间的对应关系： [FlowOperatorVector_QuarkFlavorSpaceWavefunctionMapping.md](FlowOperatorVector_QuarkFlavorSpaceWavefunctionMapping.md)
- 守恒荷扩散系数与热导率： [KappaSeries_AndLambda_FromLegacyFortran.md](KappaSeries_AndLambda_FromLegacyFortran.md)

## 公式语义、外部依据与可追溯性

### 能量角色约定

本文档采用“各向异性只进入分布函数，不重新定义准粒子色散关系”的约定，并区分：

$$
E_{\mathrm{kin}}=\sqrt{p^2+m^2},
\qquad
E_{\mathrm{dist}}=E_\xi=\sqrt{p^2+m^2+\xi(p\cos\theta)^2}.
$$

- $E_{\mathrm{kin}}$ 是普通在壳准粒子能量，用于 $\eta$、$\sigma$、$\zeta$ 的运动学能量分母，以及 $\zeta$ 的 $x=E\mp\mu$、$dE/dT$、$dE/d\mu_B$ 和 $B$ 核；关联的 $\kappa_{XY}$ 文档将同一能量用于 $p^4/E^2$ 与 Landau-Lifshitz 投影。
- $E_{\mathrm{dist}}$ 是 RS 分布的变形自变量，只通过各输运公式中的 $f_\xi$ 或 $f_\xi(1-f_\xi)$ 进入积分。
- 质量 $m=m(T,\mu_B,\xi)$ 可以来自当前各向异性热力学背景；这里排除的是在普通色散上额外加入 $\xi(p\cos\theta)^2$，不是把质量强制改为各向同性解。

当前 `TransportCoefficients.jl` 已按这一约定拆分 $E_{\mathrm{kin}}$ 与 $E_{\mathrm{dist}}$，并覆盖 $\eta$、$\sigma$、等熵平方核 $\zeta$ 和共享的 $\kappa_{XY}$ 积分状态。旧 production 不因源码修正而改变，其论文输入资格由外部 production registry 单独记录。

### 外部公式与方程映射

| 本文档内容 | 外部来源与方程号 | 映射范围 |
|---|---|---|
| 各向同性 $\eta$；各向异性下普通 $E_{\mathrm{kin}}$ 与 RS 分布的职责分离 | Thakur et al. (2017), Eqs. (18), (21), (25)-(26) | Eq. (18) 给出有限化学势 RTA 剪切黏度；Eqs. (21), (25)-(26) 给出 RS 分布及其小 $\xi$ 展开，并保持普通 $E_f^2=p^2+m_f^2$ 作为运动学分母。它支持本文档的能量职责分离，但其显式结果是小 $\xi$ 展开，不等同于本文档保留完整 RS 分布的全角积分写法。 |
| 各向同性 $\sigma$；各向异性下普通 $E_{\mathrm{kin}}$ 与 RS 分布的职责分离 | Thakur et al. (2017), Eqs. (37)-(38) | Eq. (38) 是各向同性电导率，Eq. (37) 是小 $\xi$ 各向异性结果；二者的运动学分母均为普通 $E_f^2$。 |
| 等熵平方核 $\zeta$ | Albright and Kapusta (2016), Eqs. (111), (138) | Eq. (111) 是经典统计版本；Eq. (138) 是含 $f_a^{eq}(1+d_af_a^{eq})$ 的量子统计版本，与本文档主公式直接对应。该文的 Eqs. (54)-(55) 还给出 $\sigma=s/n_B$、$(\partial\mu_B/\partial T)_\sigma$ 与 $v_n^2$ 的同一套定义。 |
| $\mu_B=0$ 时对 $\zeta$ 的交叉核验 | Mykhaylova and Sasaki (2021), Eq. (10) | 其 $[c_s^2(E^2-T^2\partial\Pi/\partial T^2)-p^2/3]^2$ 与本文档加号导数形式等价。 |
| 备用热力学导数形式 $\zeta$ | Sasaki and Redlich (2009), Eq. (2.19) | 对应本文档“备用公式”，也对应贺伟博博士论文 Eq. (5.5)；它不是当前采用的等熵平方核。博士论文 Eq. (5.6) 是电导率。 |
| 关联文档中的 $\kappa_{XY}$ | Das et al. (2022), Eq. (55) | 给出普通 $E_a$ 的运动学分母和 Landau-Lifshitz 投影结构。完整 RS 角积分是仓库采用的分布替换约定，详见 `KappaSeries_AndLambda_FromLegacyFortran.md`。 |

历史文档和源码曾把等熵平方核称为“公式 A26”。在已核对的 Albright-Kapusta、Mykhaylova-Sasaki、Sasaki-Redlich、Thakur 等文献及贺伟博博士论文中均未找到可验证的 Eq. (A26) 对应关系。因此：

- “A26”只保留为 legacy Fortran/内部别名，不能作为外部引用；
- 公式本身可追溯到 Albright and Kapusta (2016) Eq. (138)，经典统计版本为 Eq. (111)；
- 代码、审计或论文引用应使用作者、年份和精确方程号，不得继续以“A26”代替文献来源。

## 各向同性与各向异性情况

### 各向同性 (ξ = 0)

当系统各向同性时，分布函数 $f(\vec{p}) = f(|\vec{p}|) = f(p)$ 仅依赖于动量大小。此时角度积分可以简化：

$$
\int \frac{d^3p}{(2\pi)^3} = \frac{4\pi}{(2\pi)^3} \int_0^\infty p^2 dp = \frac{1}{2\pi^2} \int_0^\infty p^2 dp
$$

### 各向异性 (ξ ≠ 0, Romatschke-Strickland 形式)

当系统存在动量各向异性时，采用 Romatschke-Strickland (RS) 形式的分布函数：

$$
f_{RS}(\vec{p}) = f_{iso}\left(\sqrt{p^2 + m^2 + \xi (p \cos\theta)^2}\right)
$$

其中 $\theta$ 是动量方向与各向异性方向（通常取为 z 轴）的夹角，$\xi$ 是各向异性参数。

此时角度积分不能简化，需要完整计算：

$$
\int \frac{d^3p}{(2\pi)^3} = \frac{2\pi}{(2\pi)^3} \int_0^\infty p^2 dp \int_{-1}^{1} d(\cos\theta) = \frac{1}{4\pi^2} \int_0^\infty p^2 dp \int_{-1}^{1} d(\cos\theta)
$$

注意：φ 角积分已完成（乘以 2π），但 cosθ 积分必须显式进行。

---

## 1. 剪切粘滞系数 (η)

### 各向同性情况 (ξ = 0)

$$
\eta = \frac{1}{15T} \sum_a \int \frac{d^3p}{(2\pi)^3} \frac{p^4}{E_a^2} \left[ d_a \tau_a f_a^0 (1 - f_a^0) + d_{\bar{a}} \tau_{\bar{a}} f_{\bar{a}}^0 (1 - f_{\bar{a}}^0) \right]
$$

展开为一维积分：

$$
\eta = \frac{1}{15T} \cdot \frac{1}{2\pi^2} \sum_a \int_0^\infty dp \, \frac{p^6}{E_a^2} \left[ d_a \tau_a f_a^0 (1 - f_a^0) + d_{\bar{a}} \tau_{\bar{a}} f_{\bar{a}}^0 (1 - f_{\bar{a}}^0) \right]
$$

**注意**：积分核是 $p^6/E^2$，来自 $p^2 \cdot p^4/E^2$（相空间测度 × 物理因子）。

### 各向异性情况 (ξ ≠ 0)

$$
\eta = \frac{1}{15T} \cdot \frac{1}{4\pi^2} \sum_a \int_0^\infty dp \int_{-1}^{1} d(\cos\theta) \, \frac{p^6}{E_a^2} \left[ d_a \tau_a f_a (1 - f_a) + d_{\bar{a}} \tau_{\bar{a}} f_{\bar{a}} (1 - f_{\bar{a}}) \right]
$$

其中分布函数 $f_a = f_a(p, \cos\theta; \xi)$ 采用 RS 形式。

---

## 2. 电导率 (σ)

### 各向同性情况 (ξ = 0)

$$
\sigma = \frac{1}{3T} \sum_a \int \frac{d^3p}{(2\pi)^3} \frac{p^2 q_a^2}{E_a^2} \left[ d_a \tau_a f_a^0 (1 - f_a^0) + d_{\bar{a}} \tau_{\bar{a}} f_{\bar{a}}^0 (1 - f_{\bar{a}}^0) \right]
$$

展开为一维积分：

$$
\sigma = \frac{1}{3T} \cdot \frac{1}{2\pi^2} \sum_a \int_0^\infty dp \, \frac{p^4 q_a^2}{E_a^2} \left[ d_a \tau_a f_a^0 (1 - f_a^0) + d_{\bar{a}} \tau_{\bar{a}} f_{\bar{a}}^0 (1 - f_{\bar{a}}^0) \right]
$$

**注意**：积分核是 $p^4 q^2/E^2$，来自 $p^2 \cdot p^2 q^2/E^2$。

### 各向异性情况 (ξ ≠ 0)

$$
\sigma = \frac{1}{3T} \cdot \frac{1}{4\pi^2} \sum_a \int_0^\infty dp \int_{-1}^{1} d(\cos\theta) \, \frac{p^4 q_a^2}{E_a^2} \left[ d_a \tau_a f_a (1 - f_a) + d_{\bar{a}} \tau_{\bar{a}} f_{\bar{a}} (1 - f_{\bar{a}}) \right]
$$

---

## 3. 体粘滞系数 (ζ)

体粘滞系数采用等熵声速形式。量子统计主公式直接对应 Albright and Kapusta (2016) Eq. (138)；其经典统计版本见 Eq. (111)。“A26”仅是无法溯源到外部文献的 legacy Fortran/内部别名。

### 主公式

$$
\zeta = \frac{1}{9T} \sum_a \int d\Gamma_a \frac{\tau_a(E_a)}{E_a^2} \left[ p_a^2 + 3 v_n^2 T^2 E_a \frac{\partial}{\partial T} \left( \frac{E_a - \mu_a}{T} \right)_\sigma \right]^2 f_a^{eq} (1 + d_a f_a^{eq})
$$

方括号中第二项前的**加号**是本文档和实现应采用的正确符号。内部核对所用稿件 *Bulk viscosity of quark matter across the QCD phase transitions* 的 PDF 副本在 Eq. (1) 中曾印成减号；作者已于 2026-07-14 确认该处为文稿错误，正确形式为加号。加号形式同时与 Albright and Kapusta Eq. (138)、其 $\mu_B=0$ 极限下的 Mykhaylova and Sasaki Eq. (10)，以及 legacy Fortran 实现一致。

其中：
- $d\Gamma_a = (2s_a + 1) \frac{d^3p_a}{(2\pi)^3}$ 是相空间积分测度（包含自旋简并）
- $d_a = (-1)^{2s_a}$ 是统计因子：
  - **费米子**（$s_a = 1/2$，半整数自旋）：$d_a = (-1)^1 = -1$，因子为 $(1 - f_a^{eq})$，代表**泡利阻塞**
  - **玻色子**（$s_a = 0, 1, ...$，整数自旋）：$d_a = (-1)^0 = +1$，因子为 $(1 + f_a^{eq})$，代表**玻色增强**
- $\sigma = s/n_B$ 是熵密度与重子数密度之比

### 各向同性情况 (ξ = 0)

对于夸克体系，简并度 $g_a = 2 N_c = 6$，统计因子 $d_a = -1$（费米子）：

$$
\zeta = \frac{N_c}{9\pi^2 T} \sum_{a=u,d,s} \int_0^\infty dp \, \frac{p^2}{E_a^2} \left[ \tau_a f_a(1-f_a) \cdot B_a^2 + \tau_{\bar{a}} f_{\bar{a}}(1-f_{\bar{a}}) \cdot B_{\bar{a}}^2 \right]
$$

其中括号项 $B$ 定义为：

$$
B_a = p^2 + 3 v_n^2 T^2 E_a \cdot \frac{\partial}{\partial T}\left(\frac{E_a - \mu_a}{T}\right)_\sigma
$$

$$
B_{\bar{a}} = p^2 + 3 v_n^2 T^2 E_a \cdot \frac{\partial}{\partial T}\left(\frac{E_a + \mu_a}{T}\right)_\sigma
$$

**注意**：
- 对于费米子，$d_a = -1$，所以 $f(1 + d_a f) = f(1 - f)$，代表泡利阻塞效应
- 积分核是 $p^2/E^2 \times f(1-f) \times B^2$，注意 $B$ 是**平方**的

### 各向异性情况 (ξ ≠ 0)

在 Romatschke-Strickland 各向异性下，分布函数替换为各向异性形式，角度积分不能简化：

$$
\zeta = \frac{N_c}{9 \cdot 2\pi^2 T} \sum_{a=u,d,s} \int_0^\infty dp \int_{-1}^{1} d(\cos\theta) \, \frac{p^2}{E_a^2} \left[ \tau_a f_a(1-f_a) \cdot B_a^2 + \tau_{\bar{a}} f_{\bar{a}}(1-f_{\bar{a}}) \cdot B_{\bar{a}}^2 \right]
$$

**系数说明**：各向异性情况下分母为 $9 \cdot 2\pi^2$，因为：
- 方位角积分 $\int_0^{2\pi} d\phi = 2\pi$ 消去一个因子 2
- 简并度 $2N_c$ 中的因子 2 与方位角积分的 $2\pi$ 中的 2 相消

---

### 备用公式：热力学导数形式（已弃用）

> **注意**：以下公式在物理上也是正确的，但在不同的近似条件下推导得到。我们决定采用上述等熵声速形式，以下公式仅供参考。

$$
\zeta = -\frac{1}{3T} \cdot \frac{1}{2\pi^2} \int_0^\infty dp \, p^2 \left\{ \frac{M^2}{E} \left[ g \tau f (1 - f) + \bar{g} \bar{\tau} \bar{f} (1 - \bar{f}) \right] \times \text{bracket} - \frac{M^2}{E} \left[ g \tau f (1 - f) - \bar{g} \bar{\tau} \bar{f} (1 - \bar{f}) \right] \times \left( \frac{\partial P}{\partial n} \right)_\varepsilon \right\}
$$

其中 bracket 定义为：
$$
\text{bracket} = \frac{p^2}{3E} - \left( \frac{\partial P}{\partial \varepsilon} \right)_n \left( E - T \frac{\partial E}{\partial T} - \mu \frac{\partial E}{\partial \mu} \right) + \left( \frac{\partial P}{\partial n} \right)_\varepsilon \frac{\partial E}{\partial \mu}
$$

该公式使用热力学导数 $(∂P/∂ε)_n$ 和 $(∂P/∂n)_ε$，积分核是**线性**的（不是平方），与等熵声速形式在数值上会有差异。

---

## 4. 热力学速度 $v_n^2$

### 定义

$$
v_n^2 = \frac{s \cdot \chi_{\mu\mu} - n_B \cdot \chi_{\mu T}}{T \left( \chi_{TT} \cdot \chi_{\mu\mu} - \chi_{\mu T}^2 \right)}
$$

其中 $\chi_{xy} = \frac{\partial^2 P}{\partial x \partial y}$ 是压力的二阶偏导数：
- $\chi_{TT} = \frac{\partial^2 P}{\partial T^2}$
- $\chi_{\mu\mu} = \frac{\partial^2 P}{\partial \mu_B^2}$
- $\chi_{\mu T} = \frac{\partial^2 P}{\partial \mu_B \partial T}$

### 物理意义

$v_n^2$ 是一个热力学响应系数，出现在熵每重子数 $\sigma = s/n_B$ 固定的条件下的导数中。它在体粘滞系数公式中起到关键作用，反映了系统在压缩或膨胀过程中熵与重子数耦合的热力学行为。

### 相关速度

**固定熵密度下的声速平方**：
$$
v_s^2 = -\frac{n_B \cdot \chi_{TT} - s \cdot \chi_{\mu T}}{\mu_B \left( \chi_{TT} \cdot \chi_{\mu\mu} - \chi_{\mu T}^2 \right)}
$$

**固定熵每重子数下的声速平方**：
$$
v_\sigma^2 = \frac{v_n^2 \cdot T \cdot s + v_s^2 \cdot \mu_B \cdot n_B}{w}
$$

其中 $w = \varepsilon + P$ 是焓密度。

---

## 5. 体粘滞系数的辅助量

### 固定 $\sigma = s/n_B$ 时的导数

**化学势对温度的导数**：
$$
\frac{\partial \mu}{\partial T}\bigg|_\sigma = -\frac{\partial \sigma / \partial T}{\partial \sigma / \partial \mu}
$$

其中：
$$
\frac{\partial \sigma}{\partial T} = \frac{1}{n_B}\frac{\partial s}{\partial T} - \frac{s}{n_B^2}\frac{\partial n_B}{\partial T}
$$
$$
\frac{\partial \sigma}{\partial \mu} = \frac{1}{n_B}\frac{\partial s}{\partial \mu} - \frac{s}{n_B^2}\frac{\partial n_B}{\partial \mu}
$$

### 括号项中的导数

对于夸克（$x = (E-\mu)/T$）：
$$
\frac{\partial x}{\partial T}\bigg|_\sigma = \frac{\partial E}{\partial T} + \left(\frac{\partial E}{\partial \mu} - b_a\right) \cdot \frac{\partial \mu}{\partial T}\bigg|_\sigma
$$

$$
\frac{\partial}{\partial T}\left(\frac{E-\mu}{T}\right)_\sigma = \frac{1}{T}\frac{\partial x}{\partial T}\bigg|_\sigma - \frac{E-\mu}{T^2}
$$

对于反夸克（$x = (E+\mu)/T$）：
$$
\frac{\partial x}{\partial T}\bigg|_\sigma = \frac{\partial E}{\partial T} + \left(\frac{\partial E}{\partial \mu} + b_a\right) \cdot \frac{\partial \mu}{\partial T}\bigg|_\sigma
$$

$$
\frac{\partial}{\partial T}\left(\frac{E+\mu}{T}\right)_\sigma = \frac{1}{T}\frac{\partial x}{\partial T}\bigg|_\sigma - \frac{E+\mu}{T^2}
$$

其中 $b_a = 1/3$ 是夸克的重子数（$\mu_q = \mu_B/3$）。

### 能量导数

$$
\frac{\partial E}{\partial T} = \frac{M}{E} \cdot \frac{\partial M}{\partial T}
$$
$$
\frac{\partial E}{\partial \mu} = \frac{M}{E} \cdot \frac{\partial M}{\partial \mu}
$$

---

## 6. 电导率公式说明

### 与C++/Fortran的对比

Julia和C++/Fortran使用相同的公式结构，在相同的电荷约定下应该给出一致的结果。

**Julia公式**：
$$
\sigma = \frac{1}{3T} \cdot \frac{1}{2\pi^2} \cdot 2N_c \cdot \int dp \, \frac{p^4 q^2}{E^2} \tau f(1-f) = \frac{N_c}{3\pi^2 T} \cdot \int dp \, \frac{p^4 q^2}{E^2} \tau f(1-f)
$$

使用自然单位电荷 $q = q_{frac} \cdot e = q_{frac} \cdot \sqrt{4\pi\alpha}$：
$$
\sigma = \frac{4\alpha N_c}{3\pi T} \cdot \int dp \, \frac{p^4 q_{frac}^2}{E^2} \tau f(1-f)
$$

**C++/Fortran公式**：
$$
\sigma = \frac{4 \alpha N_c}{3\pi T} \cdot \int dp \, \frac{p^4 q_{frac}^2}{E^2} \tau f(1-f)
$$

**比值**：Julia / C++Fortran = 1（一致）

### 电荷约定

Julia使用自然单位制电荷：$q = q_{frac} \cdot e = q_{frac} \cdot \sqrt{4\pi\alpha}$

其中 $q_{frac}$ 是分数电荷（u: 2/3, d: -1/3, s: -1/3）。

---

## 7. 量纲分析与单位转换

### 自然单位制 ($\hbar = c = k_B = 1$)

在自然单位制中，所有物理量都可以用能量的幂次表示。基本量纲关系：

$$
[E] = [p] = [m] = [T] = [\mu] = \text{MeV} = \text{fm}^{-1}
$$

其中 $\hbar c \approx 197.327$ MeV·fm 用于 MeV 和 fm⁻¹ 之间的转换。

### 输运系数的量纲

| 输运系数 | 物理量纲 | 自然单位 | fm 单位 |
|---------|---------|---------|---------|
| 剪切粘滞系数 η | $[E]^3$ | MeV³ | fm⁻³ |
| 体粘滞系数 ζ | $[E]^3$ | MeV³ | fm⁻³ |
| 电导率 σ | $[E]^1$ | MeV | fm⁻¹ |

### 量纲验证

**剪切粘滞系数 η**：
$$
[\eta] = \frac{1}{[T]} \cdot \frac{[p]^6}{[E]^2} \cdot [\tau] = \frac{1}{[E]} \cdot \frac{[E]^6}{[E]^2} \cdot [E]^{-1} = [E]^3 = \text{fm}^{-3}
$$

**电导率 σ**：
$$
[\sigma] = \frac{1}{[T]} \cdot \frac{[p]^4 [q]^2}{[E]^2} \cdot [\tau] = \frac{1}{[E]} \cdot \frac{[E]^4 \cdot 1}{[E]^2} \cdot [E]^{-1} = [E]^1 = \text{fm}^{-1}
$$

注意：在自然单位制中，电荷 $q$ 是无量纲的（$[q] = 1$）。

**体粘滞系数 ζ**：
$$
[\zeta] = \frac{1}{[T]} \cdot \frac{[p]^2}{[E]^2} \cdot [\tau] \cdot [E]^4 = \frac{1}{[E]} \cdot \frac{[E]^2}{[E]^2} \cdot [E]^{-1} \cdot [E]^4 = [E]^3 = \text{fm}^{-3}
$$

### 单位转换

从 MeV 单位转换到 fm 单位：

$$
1 \text{ MeV} = \frac{1}{\hbar c} \text{ fm}^{-1} \approx \frac{1}{197.327} \text{ fm}^{-1}
$$

因此：
- $\eta [\text{MeV}^3] = \eta [\text{fm}^{-3}] \times (\hbar c)^3$
- $\sigma [\text{MeV}] = \sigma [\text{fm}^{-1}] \times (\hbar c)$

### 常用无量纲比值

在物理分析中，常使用以下无量纲比值：

| 比值 | 物理意义 | 量纲验证 |
|------|---------|---------|
| $\eta/s$ | 剪切粘滞系数与熵密度之比 | $[E]^3/[E]^3 = 1$ |
| $\zeta/s$ | 体粘滞系数与熵密度之比 | $[E]^3/[E]^3 = 1$ |
| $\sigma/T$ | 电导率与温度之比 | $[E]/[E] = 1$ |
| $\eta/T^3$ | 剪切粘滞系数的温度标度 | $[E]^3/[E]^3 = 1$ |

---

## 8. 符号说明

| 符号 | 物理意义 | 具体值 (QCD, N_c=3) |
|------|----------|---------------------|
| $a$ | 夸克味（u, d, s） | - |
| $g_a$ | 夸克简并度 | 6（自旋2×颜色3） |
| $d_a$ | 统计因子 | -1（费米子） |
| $e$ | 元电荷（自然单位） | $e = \sqrt{4\pi\alpha} \approx 0.303$ |
| $q_a$ | 夸克电荷（自然单位） | $q_u = \frac{2}{3}e$, $q_d = q_s = -\frac{1}{3}e$ |
| $\alpha$ | 精细结构常数 | $\alpha \approx 1/137$ |
| $\tau_a, \tau_{\bar{a}}$ | 弛豫时间 | 依赖 T, μ, m |
| $E_a$ | 能量 | $\sqrt{p^2 + m_a^2}$ |
| $f_a^0$ | 各向同性平衡态分布函数 | PNJL 分布函数 |
| $f_a$ | 各向异性分布函数 | RS 形式 |
| $\xi$ | 各向异性参数 | $\xi > 0$ 表示沿 z 轴压缩 |
| $M$ | 等效质量 | 依赖 T, μ |
| $P, \varepsilon, n_B$ | 压强、能量密度、重子数密度 | - |
| $s$ | 熵密度 | - |
| $v_n^2$ | 热力学速度（与体粘滞相关） | - |
| $\sigma$ | 熵每重子数 $s/n_B$ | - |
| $b_a$ | 夸克重子数 | $1/3$ |

### 自然单位制中的电荷

在自然单位制（$\hbar = c = 1$）中，精细结构常数定义为：
$$
\alpha = \frac{e^2}{4\pi} \approx \frac{1}{137}
$$

因此元电荷为：
$$
e = \sqrt{4\pi\alpha} \approx 0.303
$$

夸克电荷（自然单位）：
- u夸克：$q_u = \frac{2}{3}e \approx 0.202$
- d夸克：$q_d = -\frac{1}{3}e \approx -0.101$
- s夸克：$q_s = -\frac{1}{3}e \approx -0.101$

---

## 9. 数值实现要点

### 1. 相空间积分测度

- **各向同性**：$\frac{4\pi}{(2\pi)^3} \int p^2 dp = \frac{1}{2\pi^2} \int p^2 dp$
- **各向异性**：$\frac{2\pi}{(2\pi)^3} \int p^2 dp \int d(\cos\theta) = \frac{1}{4\pi^2} \int p^2 dp \int d(\cos\theta)$

### 2. 积分核（含相空间测度 p²）

| 输运系数 | 各向同性积分核 | 各向异性积分核 |
|---------|---------------|---------------|
| η | $p^6/E^2$ | $p^6/E^2$ |
| σ | $p^4 q^2/E^2$ | $p^4 q^2/E^2$ |
| ζ | $p^2/E^2 \times B^2$ | $p^2/E^2 \times B^2$ |

### 3. Gauss-Legendre 积分

使用 Gauss-Legendre 求积时：
- 动量积分：$\int_0^{p_{max}} f(p) dp \approx \sum_i w_i f(p_i)$
- 角度积分：$\int_{-1}^{1} g(\cos\theta) d(\cos\theta) \approx \sum_j w_j g(c_j)$

**重要**：积分核中的 $p^2$ 因子必须显式包含，因为 Gauss-Legendre 权重对应的是 $dp$ 而非 $p^2 dp$。

### 4. 体粘滞系数的二阶导数

计算 $v_n^2$ 需要压力的二阶偏导数 $\chi_{TT}$, $\chi_{\mu\mu}$, $\chi_{\mu T}$。可以使用：
- 数值微分（对一阶导数再微分）
- 自动微分（ImplicitDifferentiation.jl 支持高阶导数）

---

## 10. 关键假设与说明

1. **弛豫时间近似**：假设非平衡态分布函数相对于平衡态的偏离为一小量。
2. **准粒子图像**：强相互作用物质用具有温度/化学势依赖质量的准粒子描述。
3. **弛豫时间 τ**：依赖于热力学参数（T, μ）及散射截面，需要从微观计算获得。
4. **各向异性参数 ξ**：描述系统偏离各向同性的程度，$\xi = 0$ 对应各向同性极限。
5. **体粘滞中的 n**：指净重子数密度 $n_B = (\rho_u + \rho_d + \rho_s)/3$。
6. **等熵声速形式**：体粘滞系数公式采用固定熵每重子数 $\sigma = s/n_B$ 的形式，需要计算相关的热力学导数。
7. **统计因子**：费米子 $d_a = -1$（泡利阻塞），玻色子 $d_a = +1$（玻色增强）。

---

## 11. 参考文献与 Zotero 映射

1. M. Albright and J. I. Kapusta, “Quasiparticle Theory of Transport Coefficients for Hadronic Matter at Finite Temperature and Baryon Density,” *Phys. Rev. C* **93**, 014903 (2016). DOI: [10.1103/PhysRevC.93.014903](https://doi.org/10.1103/PhysRevC.93.014903), arXiv: [1508.02696](https://arxiv.org/abs/1508.02696). Zotero item `NKDVM97F`; BibTeX key `Albright:2015fpa`.
2. V. Mykhaylova and C. Sasaki, “Impact of quark quasiparticles on transport coefficients in hot QCD,” *Phys. Rev. D* **103**, 014007 (2021). DOI: [10.1103/PhysRevD.103.014007](https://doi.org/10.1103/PhysRevD.103.014007), arXiv: [2007.06846](https://arxiv.org/abs/2007.06846). Zotero item `J6TSE7IA`; BibTeX key `Mykhaylova:2020pfk`.
3. L. Thakur, P. K. Srivastava, G. P. Kadam, M. George, and H. Mishra, “Shear viscosity $\eta$ to electrical conductivity $\sigma_{el}$ ratio for an anisotropic QGP,” *Phys. Rev. D* **95**, 096009 (2017). DOI: [10.1103/PhysRevD.95.096009](https://doi.org/10.1103/PhysRevD.95.096009), arXiv: [1703.03142](https://arxiv.org/abs/1703.03142). Zotero item `PGUUMFXD`; BibTeX key `Thakur:2017hfc`.
4. C. Sasaki and K. Redlich, “Bulk viscosity in quasi particle models,” *Phys. Rev. C* **79**, 055207 (2009). DOI: [10.1103/PhysRevC.79.055207](https://doi.org/10.1103/PhysRevC.79.055207), arXiv: [0806.4745](https://arxiv.org/abs/0806.4745). Zotero item `5QRS2YMA`; BibTeX key `Sasaki:2008fg`.
