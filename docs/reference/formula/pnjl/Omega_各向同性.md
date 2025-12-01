# PNJL模型：巨热力学势与稳定平衡态序参量求解

（基于论文第 19–23 页）

---

## 1. PNJL 模型巨热力学势（Grand Potential）

（式 2.26）


```math
\Omega 
= 2G \sum_{i=u,d,s} \phi_i^2
 - 4K\,\phi_u\phi_d\phi_s
 - 2N_c \int_0^\Lambda \frac{d^3p}{(2\pi)^3} \sum_{i=u,d,s} E_i
 - 2T \sum_{i=u,d,s} \int \frac{d^3p}{(2\pi)^3} ( \mathcal{Q}_1 + \mathcal{Q}_2 )
+ U(\Phi,\bar\Phi;T)
```

其中：

```math
E_i = \sqrt{p^2 + M_i^2}
```

NJL 有效质量能隙方程（式 2.17）：


```math
M_i = m_i - 4G\phi_i + 2K\phi_j\phi_k,\quad i\neq j\neq k.
```

---

## 2. Polyakov 圈有效势（对数形式）

（式 2.23–2.25）


```math
\frac{U(\Phi,\bar\Phi;T)}{T^4}
= -\frac{a(T)}{2}\Phi\bar\Phi
  + b(T)\ln\!\left[
    1 - 6\Phi\bar\Phi + 4(\Phi^3+\bar\Phi^3) - 3(\Phi\bar\Phi)^2
  \right]
```

```math
a(T)=a_0 + a_1(T_0/T) + a_2(T_0/T)^2
```

```math
b(T)=b_3 (T_0/T)^3
```

---

## 3. 有效费米对数函数 ( \mathcal{Q}_1, \mathcal{Q}_2 )

（式 2.27–2.28）


```math
\mathcal{Q}_1 = \ln\!\left(
1 + 3\Phi e^{-(E_i-\mu_i)/T}
  + 3\bar\Phi e^{-2(E_i-\mu_i)/T}
  + e^{-3(E_i-\mu_i)/T}
\right)
```

```math
\mathcal{Q}_2 = \ln\!\left(
1 + 3\bar\Phi e^{-(E_i+\mu_i)/T}
  + 3\Phi e^{-2(E_i+\mu_i)/T}
  + e^{-3(E_i+\mu_i)/T}
\right)
```

---

## 4. 稳定平衡态序参量条件（能隙方程）

（式 2.29–2.30）


```math
\frac{\partial\Omega}{\partial \phi_u}
= \frac{\partial\Omega}{\partial \phi_d}
= \frac{\partial\Omega}{\partial \phi_s}
= 0
```

```math
\frac{\partial\Omega}{\partial \Phi}
= \frac{\partial\Omega}{\partial \bar\Phi}
= 0
```

---

## 5. 热力学量（以热力学势 Ω 的偏导定义）

你要求使用 **原始定义**（适用于自动微分），因此以下给出偏导形式。

### 5.1 压强

（式 2.31）


```math
P = -\Omega
```

---

### 5.2 夸克数密度（原始定义：Ω 关于 μ 的偏导）

```math
\rho_i = -\frac{\partial \Omega}{\partial \mu_i}
```

若只需要重子数密度：

```math
\rho_B = \frac{1}{3}(\rho_u+\rho_d+\rho_s)
```

---

### 5.3 熵密度（原始定义：Ω 关于 T 的偏导）

不展开至显式分布函数形式：

```math
s = -\frac{\partial \Omega}{\partial T}
```

---

### 5.4 能量密度（不依赖显式展开）

（式 2.34）

```math
\epsilon = Ts + \sum_i \mu_i \rho_i - P
```

---

# 6. PNJL 模型参数（论文表格整理）

以下参数来自论文表 2.1（NJL 部分）与表 2.2（Polyakov 部分）。

---

## 6.1 NJL 模型参数（表 2.1）

```text
Λ      = 602.3 MeV      （三动量截断）
GΛ²    = 1.835
KΛ⁵    = 12.36
m_u    = 5.5 MeV
m_d    = 5.5 MeV
m_s    = 140.7 MeV
```

---

## 6.2 Polyakov 场参数（表 2.2）

```text
T₀     = 210 MeV        （纯规范理论）
a₀     = 3.51
a₁     = -2.47
a₂     = 15.2
b₃     = -1.75
```

---

