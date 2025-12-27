# 夸克物质的输运系数（弛豫时间近似）

在弛豫时间近似下，基于相对论动力学，夸克物质的剪切粘滞系数、体粘滞系数和电导率表达式如下。

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

### 各向同性情况 (ξ = 0)

$$
\begin{aligned}
\zeta = & -\frac{1}{3T} \cdot \frac{1}{2\pi^2} \int_0^\infty dp \, p^2 \Bigg\{ \frac{M^2}{E} \left[ g \tau f_a^0 (1 - f_a^0) + \bar{g} \bar{\tau} f_{\bar{a}}^0 (1 - f_{\bar{a}}^0) \right] \\
& \times \left[ \frac{p^2}{3E} - \left( \frac{\partial P}{\partial \varepsilon} \right)_n \left( E - T \frac{\partial E}{\partial T} - \mu \frac{\partial E}{\partial \mu} \right) + \left( \frac{\partial P}{\partial n} \right)_\varepsilon \frac{\partial E}{\partial \mu} \right] \\
& - \frac{M^2}{E} \left[ g \tau f_a^0 (1 - f_a^0) - \bar{g} \bar{\tau} f_{\bar{a}}^0 (1 - f_{\bar{a}}^0) \right] \times \left( \frac{\partial P}{\partial n} \right)_\varepsilon \Bigg\}
\end{aligned}
$$

### 各向异性情况 (ξ ≠ 0)

$$
\begin{aligned}
\zeta = & -\frac{1}{3T} \cdot \frac{1}{4\pi^2} \int_0^\infty dp \int_{-1}^{1} d(\cos\theta) \, p^2 \Bigg\{ \frac{M^2}{E} \left[ g \tau f_a (1 - f_a) + \bar{g} \bar{\tau} f_{\bar{a}} (1 - f_{\bar{a}}) \right] \\
& \times \left[ \frac{p^2}{3E} - \left( \frac{\partial P}{\partial \varepsilon} \right)_n \left( E - T \frac{\partial E}{\partial T} - \mu \frac{\partial E}{\partial \mu} \right) + \left( \frac{\partial P}{\partial n} \right)_\varepsilon \frac{\partial E}{\partial \mu} \right] \\
& - \frac{M^2}{E} \left[ g \tau f_a (1 - f_a) - \bar{g} \bar{\tau} f_{\bar{a}} (1 - f_{\bar{a}}) \right] \times \left( \frac{\partial P}{\partial n} \right)_\varepsilon \Bigg\}
\end{aligned}
$$

---

## 电导率公式说明

### 与C++/Fortran的对比

Julia和C++/Fortran使用相同的公式结构，在相同的电荷约定下应该给出一致的结果。

**Julia公式**：
$$\sigma = \frac{1}{3T} \cdot \frac{1}{2\pi^2} \cdot 2N_c \cdot \int dp \, \frac{p^4 q^2}{E^2} \tau f(1-f) = \frac{N_c}{3\pi^2 T} \cdot \int dp \, \frac{p^4 q^2}{E^2} \tau f(1-f)$$

使用自然单位电荷 $q = q_{frac} \cdot e = q_{frac} \cdot \sqrt{4\pi\alpha}$：
$$\sigma = \frac{4\alpha N_c}{3\pi T} \cdot \int dp \, \frac{p^4 q_{frac}^2}{E^2} \tau f(1-f)$$

**C++/Fortran公式**：
$$\sigma = \frac{4 \alpha N_c}{3\pi T} \cdot \int dp \, \frac{p^4 q_{frac}^2}{E^2} \tau f(1-f)$$

**比值**：Julia / C++Fortran = 1（一致）

### 电荷约定

Julia使用自然单位制电荷：$q = q_{frac} \cdot e = q_{frac} \cdot \sqrt{4\pi\alpha}$

其中 $q_{frac}$ 是分数电荷（u: 2/3, d: -1/3, s: -1/3）。

---

## 符号说明

| 符号 | 物理意义 | 具体值 (QCD, N_c=3) |
|------|----------|---------------------|
| $a$ | 夸克味（u, d, s） | - |
| $d_a, g$ | 夸克简并度 | 6（自旋2×颜色3） |
| $d_{\bar{a}}, \bar{g}$ | 反夸克简并度 | 6 |
| $e$ | 元电荷（自然单位） | $e = \sqrt{4\pi\alpha} \approx 0.303$ |
| $q_a$ | 夸克电荷（自然单位） | $q_u = \frac{2}{3}e$, $q_d = q_s = -\frac{1}{3}e$ |
| $\alpha$ | 精细结构常数 | $\alpha \approx 1/137$ |
| $\tau_a, \tau_{\bar{a}}$ | 弛豫时间 | 依赖 T, μ, m |
| $E_a$ | 能量 | $\sqrt{p^2 + m_a^2}$ |
| $f_a^0$ | 各向同性平衡态分布函数 | PNJL 分布函数 |
| $f_a$ | 各向异性分布函数 | RS 形式 |
| $\xi$ | 各向异性参数 | $\xi > 0$ 表示沿 z 轴压缩 |
| $M$ | 等效质量 | 依赖 T, μ |
| $P, \varepsilon, n$ | 压强、能量密度、粒子数密度 | - |

### 自然单位制中的电荷

在自然单位制（$\hbar = c = 1$）中，精细结构常数定义为：
$$\alpha = \frac{e^2}{4\pi} \approx \frac{1}{137}$$

因此元电荷为：
$$e = \sqrt{4\pi\alpha} \approx 0.303$$

夸克电荷（自然单位）：
- u夸克：$q_u = \frac{2}{3}e \approx 0.202$
- d夸克：$q_d = -\frac{1}{3}e \approx -0.101$
- s夸克：$q_s = -\frac{1}{3}e \approx -0.101$

---

## 数值实现要点

### 1. 相空间积分测度

- **各向同性**：$\frac{4\pi}{(2\pi)^3} \int p^2 dp = \frac{1}{2\pi^2} \int p^2 dp$
- **各向异性**：$\frac{2\pi}{(2\pi)^3} \int p^2 dp \int d(\cos\theta) = \frac{1}{4\pi^2} \int p^2 dp \int d(\cos\theta)$

### 2. 积分核（含相空间测度 p²）

| 输运系数 | 各向同性积分核 | 各向异性积分核 |
|---------|---------------|---------------|
| η | $p^6/E^2$ | $p^6/E^2$ |
| σ | $p^4 q^2/E^2$ | $p^4 q^2/E^2$ |
| ζ | $p^2 \times (\text{复杂表达式})$ | $p^2 \times (\text{复杂表达式})$ |

### 3. Gauss-Legendre 积分

使用 Gauss-Legendre 求积时：
- 动量积分：$\int_0^{p_{max}} f(p) dp \approx \sum_i w_i f(p_i)$
- 角度积分：$\int_{-1}^{1} g(\cos\theta) d(\cos\theta) \approx \sum_j w_j g(c_j)$

**重要**：积分核中的 $p^2$ 因子必须显式包含，因为 Gauss-Legendre 权重对应的是 $dp$ 而非 $p^2 dp$。

---

## 关键假设与说明

1. **弛豫时间近似**：假设非平衡态分布函数相对于平衡态的偏离为一小量。
2. **准粒子图像**：强相互作用物质用具有温度/化学势依赖质量的准粒子描述。
3. **弛豫时间 τ**：依赖于热力学参数（T, μ）及散射截面，需要从微观计算获得。
4. **各向异性参数 ξ**：描述系统偏离各向同性的程度，$\xi = 0$ 对应各向同性极限。
5. **体粘滞中的 n**：指净重子数密度 $n_B = (\rho_u + \rho_d + \rho_s)/3$。
