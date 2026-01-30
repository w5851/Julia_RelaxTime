# 各向异性夸克物质中的有效分布函数

## 1. 各向异性分布函数的引入

### 1.1 基本形式

在重离子碰撞中，由于沿束流方向的快速膨胀，系统会产生动量空间各向异性。采用Romatschke-Strickland形式的各向异性分布函数：

```math
f^{\text{aniso}}(\mathbf{p}) = f^{\text{iso}}\left(\sqrt{\mathbf{p}^2 + \xi(\mathbf{p} \cdot \mathbf{n})^2}\right)
```

其中：
- $\mathbf{n}$：各向异性方向的单位矢量
- $\xi$：各向异性参数，量化动量空间各向异性程度
- $f^{\text{iso}}$：各向同性分布函数

### 1.2 各向异性参数的定义

```math
\xi = \frac{\langle p_\perp^2 \rangle}{2\langle p_\parallel^2 \rangle} - 1
```

其中：
- $p_\perp = |\mathbf{p} - (\mathbf{p} \cdot \mathbf{n})\mathbf{n}|$：垂直于$\mathbf{n}$的动量分量
- $p_\parallel = \mathbf{p} \cdot \mathbf{n}$：平行于$\mathbf{n}$的动量分量

参数范围：
- $-1 < \xi < 0$：动量分布沿各向异性方向收缩
- $\xi > 0$：动量分布沿各向异性方向拉伸

## 2. 弱各向异性近似

### 2.1 分布函数的展开

在弱各向异性极限（$|\xi| \ll 1$）下，将分布函数展开至$\xi$的一阶：

```math
f^{\text{aniso}}(\mathbf{p}) = f^0 - \frac{\xi(\mathbf{n} \cdot \mathbf{p})^2}{2ET} f^0(1 - f^0)
```

其中$f^0$是平衡态的费米-狄拉克分布函数。

### 2.2 坐标系选择

选择各向异性方向沿z轴：
```math
\mathbf{n} = (0, 0, 1)
```

动量球坐标表示：
```math
\mathbf{p} = p(\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)
```

此时各向异性项简化为：
```math
\xi(\mathbf{n} \cdot \mathbf{p})^2 = \xi p^2\cos^2\theta
```

### 2.3 质心系中的分布函数近似（式60）

在计算散射截面的质心系中，各向异性分布函数采用以下近似形式：

```math
f^{\text{aniso}}(p_{\text{cm}},\mu,x) = f^{0}(p_{\text{cm}},\mu) - \frac{p_{\text{cm}}^{2}\xi x^{2}}{2E_{\text{cm}}T} f^{0}(p_{\text{cm}},\mu)(1 - f^{0}(p_{\text{cm}},\mu))
```

其中：
- $p_{\text{cm}}$：质心系中的动量大小
- $E_{\text{cm}} = \sqrt{s}/2$：质心系中的能量，$s$为曼德尔斯坦变量
- $x = \cos\theta$：动量与各向异性方向的夹角余弦
- $\mu$：化学势

**物理意义和推导：**

1. **坐标变换**：将实验室系的分布函数变换到质心系，保持各向异性参数$\xi$不变
2. **能量关系**：在质心系中，$E_{\text{cm}} = \sqrt{s}/2$，其中$s = 4m_q^2 + 4p_{\text{cm}}^2$
3. **角度变量**：$x$表示质心系中粒子动量与各向异性方向$\mathbf{n}$的夹角
4. **小参数展开**：保持$\xi$的一阶项，忽略高阶贡献

**在散射截面计算中的应用：**

这个近似形式用于计算热平均散射截面：
```math
\sigma_{ab\rightarrow cd} = \int_{s_0}^{\infty} ds \int_{t_{\text{min}}}^{t_{\text{max}}} dt \frac{d\sigma_{ab\rightarrow cd}}{dt} \sin^2\theta_p \int_{-1}^{1} dx_3 (1 - f_c^{\text{aniso}}(p_{\text{cm}},\mu,x_3)) \times \int_{-1}^{1} dx_4 (1 - f_d^{\text{aniso}}(p_{\text{cm}},\mu,x_4)) \mathcal{L}(s,\mu,x_1,x_2)
```

其中$(1 - f_{c,d}^{\text{aniso}})$是泡利阻塞因子，考虑了终态费米子的量子统计效应。

## 3. 在各向异性介质中的应用

### 3.1 修正的gap方程

在各向异性介质中，组分夸克质量的gap方程修正为：

```math
0 = \left[ 2N_c N_f G \int_0^\infty \frac{p^2 dp}{\pi^2} \frac{m_q}{E} \left( -f_q^0 - \bar{f}_q^0 + 1 + \frac{p^2 \xi F_p}{6ET} \right) \right] + m_0 - m_q
```

其中$F_p = f_q^0(1 - f_q^0) + \bar{f}_q^0(1 - \bar{f}_q^0)$。

### 3.2 介子偏振函数

介子偏振函数$\Pi_M$在各向异性介质中也需要修正，其中函数$A$和$B_0$分别修正为：

函数$A$的修正：
```math
A = 4 \int \frac{p^2 dp}{E} \left[ f_q^0 + \bar{f}_q^0 - 1 - \frac{\xi p^2 F_p}{6ET} \right]
```

函数$B_0$的修正较为复杂，分为两种情况：
- $\mathbf{k} = 0, k_0 \neq 0$
- $\mathbf{k} \neq 0, k_0 \neq 0$

### 3.3 输运系数

#### 剪切粘滞系数
```math
\eta_a = -\frac{\xi d_a}{180T^2} \int \frac{dp}{\pi^2} \frac{\tau_a p^8}{E_a^3} f_a^0(1 - f_a^0)\left(1 - 2f_a^0 + \frac{T}{E_a}\right) + \frac{d_a}{30T} \int \frac{dp}{\pi^2} \frac{\tau_a p^6}{E_a^2} f_a^0(1 - f_a^0)
```

#### 电导率
```math
\sigma_{\text{el}a} = \frac{e_a^2 d_a \tau_a}{6T} \int \frac{dp}{\pi^2} \frac{p^4}{E_a^2} f_a^0(1 - f_a^0)\left(1 + \frac{\xi}{3}\right) - \frac{\xi e_a^2 d_a \tau_a}{36T^2} \int \frac{dp}{\pi^2} \frac{p^6}{E_a^3} f_a^0(1 - f_a^0)\left(1 - 2f_a^0 + \frac{T}{E_a}\right)
```

## 4. 物理结果

### 4.1 手征相变
- 各向异性增强（$\xi$增大）会阻碍手征对称性恢复
- 临界端点(CEP)随$\xi$增大向更高化学势和更低温度移动
- 对各向异性参数的敏感程度在温度和化学势方向几乎相同

### 4.2 介子性质
- $\pi$介子和$\sigma$介子质量对$\xi$的依赖关系相反
- 介子衰变宽度随$\xi$增大而减小

### 4.3 输运性质
- 标度剪切粘滞系数$\eta/T^3$和标度电导率$\sigma_{el}/T$在临界温度附近呈现凹陷结构
- Seebeck系数$S$为正号，表明主导载流子为上夸克
- 在低温区域，$S$随$\xi$显著增强

## 5. 结论

各向异性分布函数的引入为研究重离子碰撞中产生的非平衡QCD物质提供了有力工具。通过弱各向异性近似，可以系统地研究动量空间各向异性对手征相变、介子性质和输运系量的影响，为理解QCD相图和重离子碰撞实验观测提供理论支持。

质心系中的分布函数近似（式60）特别重要，它使得在计算散射截面和弛豫时间时能够自洽地包含各向异性效应，这是准确计算输运系数的关键步骤。

---

*本文基于Zhang等人在《Nuclear Science and Techniques》2022年发表的论文整理*
对于PNJL模型中的有效分布函数，各向异性修正的一阶展开形式如下：

---

### ✅ **PNJL模型各向异性分布函数的一阶展开**

#### **1. 基本各向异性变换**
按照Romatschke-Strickland方法，各向异性分布函数通过以下变换得到：
\[
f^{\text{aniso}}(\mathbf{p}) = f^{\text{iso}}\left(\sqrt{\mathbf{p}^2 + \xi(\mathbf{n} \cdot \mathbf{p})^2}\right)
\]

#### **2. 一阶展开（弱各向异性近似，|\xi| ≪ 1）**
在弱各向异性极限下，对上述函数在ξ=0处进行泰勒展开，保留到ξ的一阶项：

\[
f^{\text{aniso}}(\mathbf{p}) \approx f^{\text{iso}}(p) + \xi \cdot \frac{\partial f^{\text{aniso}}}{\partial \xi}\bigg|_{\xi=0}
\]

其中：
\[
\frac{\partial f^{\text{aniso}}}{\partial \xi}\bigg|_{\xi=0} = \frac{1}{2} \frac{(\mathbf{n} \cdot \mathbf{p})^2}{E} \cdot \frac{\partial f^{\text{iso}}}{\partial E}
\]

因此，**一阶修正项**为：
\[
\delta f = \frac{\xi}{2} \frac{(\mathbf{n} \cdot \mathbf{p})^2}{E} \cdot \frac{\partial f^{\text{iso}}}{\partial E}
\]

---

### ✅ **应用于PNJL模型的夸克分布函数**

#### **夸克分布函数的一阶各向异性修正：**
\[
f_q^{\text{aniso}} = f_q^+ + \frac{\xi}{2} \frac{(\mathbf{n} \cdot \mathbf{p})^2}{E} \cdot \frac{\partial f_q^+}{\partial E}
\]

#### **反夸克分布函数的一阶各向异性修正：**
\[
f_{\bar{q}}^{\text{aniso}} = f_q^- + \frac{\xi}{2} \frac{(\mathbf{n} \cdot \mathbf{p})^2}{E} \cdot \frac{\partial f_q^-}{\partial E}
\]

---

### ✅ **具体导数计算**

#### **夸克分布函数的导数：**
令 \( x = \beta(E - \mu) \)，则：
\[
\frac{\partial f_q^+}{\partial E} = \beta \cdot \frac{\partial f_q^+}{\partial x}
\]

通过商法则计算：
\[
\frac{\partial f_q^+}{\partial x} = \frac{N'(x)D(x) - N(x)D'(x)}{[D(x)]^2}
\]

其中：
- \( N(x) = \Phi e^{-x} + 2\bar{\Phi} e^{-2x} + e^{-3x} \)
- \( D(x) = 1 + 3\Phi e^{-x} + 3\bar{\Phi} e^{-2x} + e^{-3x} \)
- \( N'(x) = -\Phi e^{-x} - 4\bar{\Phi} e^{-2x} - 3e^{-3x} \)
- \( D'(x) = -3\Phi e^{-x} - 6\bar{\Phi} e^{-2x} - 3e^{-3x} \)

#### **反夸克分布函数的导数：**
令 \( y = \beta(E + \mu) \)，类似计算：
\[
\frac{\partial f_q^-}{\partial E} = \beta \cdot \frac{\partial f_q^-}{\partial y}
\]

---

### ✅ **最终表达式**

#### **各向异性夸克分布函数（一阶修正）：**
\[
f_q^{\text{aniso}} = f_q^+ + \frac{\xi \beta}{2} \frac{(\mathbf{n} \cdot \mathbf{p})^2}{E} \cdot \frac{N'(x)D(x) - N(x)D'(x)}{[D(x)]^2}
\]

#### **各向异性反夸克分布函数（一阶修正）：**
\[
f_{\bar{q}}^{\text{aniso}} = f_q^- + \frac{\xi \beta}{2} \frac{(\mathbf{n} \cdot \mathbf{p})^2}{E} \cdot \frac{C'(y)D(y) - C(y)D'(y)}{[D(y)]^2}
\]

其中 \( C(y) = \bar{\Phi} e^{-y} + 2\Phi e^{-2y} + e^{-3y} \)，\( D(y) \) 类似定义。

---

### 📌 **总结：**
是的，对于PNJL模型，只需要在各向同性分布函数的基础上，**加入一个与 \( \xi(\mathbf{n} \cdot \mathbf{p})^2 \) 成正比的一阶修正项**，该修正项包含分布函数对能量的导数。这就是PNJL模型中引入动量各向异性的**弱各向异性近似方法**。