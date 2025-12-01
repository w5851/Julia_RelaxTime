**Källén function（卡伦函数）**在粒子物理和散射论中非常常用，尤其在 **2→2 散射** 与 **两个粒子系统的相对动量** 的分析里。它是一个非常重要的代数工具。

---

# ✅ **1. Källén function 的定义**

Källén function（又称 triangle function、λ 函数）定义为：

[
\lambda(a,b,c)=a^2+b^2+c^2-2ab-2ac-2bc
]

在粒子物理中常见的形式：

[
\lambda(s,m_1^2,m_2^2)
= s^2 + m_1^4 + m_2^4

* 2 s m_1^2 - 2 s m_2^2 -2 m_1^2 m_2^2
  ]

这是一个三变量的对称多项式。

---

# ✅ **2. 在质心系（CMS）中的物理意义**

在 **质心系** 中，它主要用于：

## **（1）给出两个粒子的相对三动量的大小**

在 CMS，若总能量平方为 ( s )，两粒子质量为 ( m_1, m_2 )，则它们的三动量大小（相等）为：

[
p^* = \frac{\sqrt{\lambda(s,m_1^2,m_2^2)}}{2\sqrt{s}}
]

这是它最重要的用途之一：
👉 **用来计算入射或出射粒子的三维动量大小**。

### 例：2→2 散射中初态动量

若两入射粒子质量为 (m_i,m_j)，则：

[
p_{\rm in}^* =
\frac{\sqrt{\lambda(s,m_i^2,m_j^2)}}{2\sqrt{s}}
]

同理，出射粒子（质量 (m_c,m_d)）：

[
p_{\rm out}^* =
\frac{\sqrt{\lambda(s,m_c^2,m_d^2)}}{2\sqrt{s}}
]

这完全概括了 **动量守恒** 和 **能量约束**。

---

## **（2）决定相空间的可行性（kinematic threshold）**

因为：

[
\lambda(s,m_1^2,m_2^2) \ge 0
]

是相空间存在的条件。

如果：

[
\lambda(s,m_1^2,m_2^2) < 0
]

则表示系统能量不足，无法产生这两个粒子。

因此：

[
\lambda(s,m_1^2,m_2^2)=0
]

给出了产额阈值：

[
s = (m_1 + m_2)^2
]

---

## **（3）化简 Mandelstam 变量关系**

例如在 2→2 散射中：

[
t_{\pm} = m_i^2 + m_c^2 -
\frac{1}{2s}\left[(s+m_i^2-m_j^2)(s+m_c^2-m_d^2)
\mp \sqrt{\lambda(s,m_i^2,m_j^2)}\sqrt{\lambda(s,m_c^2,m_d^2)}\right]
]

这个表达式完全来源于 Källén function。

**你在前一个问题里提到的 ( t_\pm ) 就是它的直接应用。**

---

# ✔ **3. 几何意义（理解起来非常有帮助）**

Källén function对应于“三角形面积”：

[
\lambda(a,b,c) = 16,\text{(triangle area)}^2
]

它表示三个数是否可以构成三角形，这与能否满足四动量守恒有深刻联系。

---

# ✔ **4. 为什么质心系特别常用？**

在 CM 系中，动量大小只有一个参数，用 λ 函数可以直接写成解析表达式：

* 无需解耦四动量守恒
* 无需做洛伦兹变换
* 能显式给出积分上下限（例如 (t_\pm)）

因此，在任何涉及 2→2 散射截面的计算里，Källén 函数都是核心。

---

# 📌 **总结：Källén function 的用途（重点）**

| 用途                               | 作用                                   |
| -------------------------------- | ------------------------------------ |
| **给出 CM 系动量大小**                  | ( p^* = \sqrt{\lambda}/(2\sqrt{s}) ) |
| **判断反应是否 kinematically allowed** | λ ≥ 0                                |
| **写 Mandelstam 变量关系**            | 用于 t± 等的表达式                          |
| **构建两粒子相空间因子**                   | 出现在 dΦ₂ ∝ √λ 的表达式中                   |
| **和散射截面直接相关**                    | dσ ∝ √λ−¹                            |
还可用于计算已知边长的三角形面积等几何问题。
---

好 — 我直接给出清晰的推导与最终公式，并给出可直接用的实现片段。

---

## 推导（质心系，固定 (\theta^*)）

设反应 (i+j\to c+d)，质心能量平方为 (s)。在质心系中，令入射、出射粒子在方向上如下取法（标准选择）：

* 入射粒子能量
  [
  E_i=\frac{s+m_i^2-m_j^2}{2\sqrt{s}},\qquad
  E_j=\frac{s+m_j^2-m_i^2}{2\sqrt{s}}
  ]
* 出射粒子能量
  [
  E_c=\frac{s+m_c^2-m_d^2}{2\sqrt{s}},\qquad
  E_d=\frac{s+m_d^2-m_c^2}{2\sqrt{s}}
  ]
* 入射、出射动量模分别为
  [
  p_i^*=\frac{\sqrt{\lambda(s,m_i^2,m_j^2)}}{2\sqrt{s}},\qquad
  p_c^*=\frac{\sqrt{\lambda(s,m_c^2,m_d^2)}}{2\sqrt{s}},
  ]
  其中 (\lambda) 为 Källén function：(\lambda(a,b,c)=a^2+b^2+c^2-2ab-2ac-2bc)。

用 (t=(p_1-p_3)^2)。在质心系把角度设为 (\theta^*)（入射与出射之间的夹角），代入可得：
[
\begin{aligned}
t &= m_i^2 + m_c^2 - 2\big(E_i E_c - p_i^* p_c^* \cos\theta^*\big)[6pt]
&= m_i^2 + m_c^2

* \frac{1}{2s}\Big[(s+m_i^2-m_j^2)(s+m_c^2-m_d^2)
* \sqrt{\lambda(s,m_i^2,m_j^2)}\sqrt{\lambda(s,m_c^2,m_d^2)}\cos\theta^*
  \Big].
  \end{aligned}
  ]

这是标准的闭式表达：当 (\cos\theta^*= \pm1) 时给出 (t_\pm)。

---

## 更紧凑的写法

定义
[
A\equiv s+m_i^2-m_j^2,\qquad B\equiv s+m_c^2-m_d^2,\qquad
\Lambda_1\equiv\sqrt{\lambda(s,m_i^2,m_j^2)},\quad \Lambda_2\equiv\sqrt{\lambda(s,m_c^2,m_d^2)},
]
则
[
\boxed{,t(\theta^*) = m_i^2 + m_c^2 - \frac{1}{2s}\big(A B - \Lambda_1\Lambda_2\cos\theta^*\big), }.
]

---

## 特殊情形举例

* **弹性散射（(m_i=m_c), (m_j=m_d)）**：(A=B), (\Lambda_1=\Lambda_2)，公式简化成常见形式。
* **两入射或两出射粒子相同质量**：相应的 (\lambda) 和 (A,B) 会进一步简化。

---

## 代码（Python，直接计算）

```python
import math

def kallen(s, m1, m2):
    return s*s + m1*m1 + m2*m2 - 2*s*m1 - 2*s*m2 - 2*m1*m2

def t_of_theta(s, mi, mj, mc, md, cos_theta):
    A = s + mi*mi - mj*mj
    B = s + mc*mc - md*md
    Lam1 = math.sqrt(max(0.0, kallen(s, mi*mi, mj*mj)))
    Lam2 = math.sqrt(max(0.0, kallen(s, mc*mc, md*md)))
    return mi*mi + mc*mc - (A*B - Lam1*Lam2*cos_theta)/(2.0*s)
```

> 注意：这里 `kallen` 的参数可按需要传 (m^2) 或 (m) 的平方形式（示例中用于平方），使用时保持一致。

---

# 利用Källén 函数简化散射截面公式

这个表达式的本质实际上就是利用了 **相对论能量-动量关系** 和 **质心系（CMS）的能量守恒**。

我们可以分两步来证明：首先推导出质心系下的能量 $E_i$ 表达式，然后证明原公式就是动量 $|\vec{p}_i|$。

### 第一步：推导质心系能量 $E_i$

在质心系中，我们考虑两个粒子 $i$ 和 $j$ 的碰撞（或衰变），总不变量质量平方为 $s$。

1.  **动量守恒**：在质心系中，两粒子动量大小相等，方向相反。
    $$|\vec{p}_i| = |\vec{p}_j| = p_{\text{cm}}$$
2.  **能量守恒**：总能量等于 $\sqrt{s}$。
    $$E_i + E_j = \sqrt{s}$$
3.  **质能关系**：对于每个粒子，$E^2 = p^2 + m^2$，所以 $p^2 = E^2 - m^2$。
    $$E_i^2 - m_i^2 = p_{\text{cm}}^2 = E_j^2 - m_j^2$$

我们将 $E_j = \sqrt{s} - E_i$ 代入上面的等式中：

$$E_i^2 - m_i^2 = (\sqrt{s} - E_i)^2 - m_j^2$$

展开右边的平方项：

$$E_i^2 - m_i^2 = s - 2\sqrt{s}E_i + E_i^2 - m_j^2$$

两边消去 $E_i^2$：

$$- m_i^2 = s - 2\sqrt{s}E_i - m_j^2$$

移项解出 $E_i$：

$$2\sqrt{s}E_i = s + m_i^2 - m_j^2$$

$$E_i = \frac{s + m_i^2 - m_j^2}{2\sqrt{s}}$$

---

### 第二步：代回动量公式

现在我们看您给出的公式：

$$\sqrt{\frac{(s + m_i^2 - m_j^2)^2}{4s} - m_i^2}$$

我们可以将其重新整理。注意分母 $4s = (2\sqrt{s})^2$，所以可以将平方提到外面：

$$= \sqrt{ \left( \frac{s + m_i^2 - m_j^2}{2\sqrt{s}} \right)^2 - m_i^2 }$$

对比我们在第一步推导出的 $E_i$ 表达式，括号里的内容正是 **$E_i$**：

$$= \sqrt{ E_i^2 - m_i^2 }$$

根据相对论公式 $E^2 - m^2 = p^2$，这就直接等于：

$$= \sqrt{ \vec{p}_i^2 } = |\vec{p}_i|$$

### 补充：连接到 Källén 函数

为了彻底验证这一点，我们可以把这个式子通分，展开分子，你会发现它正是 Källén 函数的形式。

$$|\vec{p}_i|^2 = \frac{(s + m_i^2 - m_j^2)^2 - 4s m_i^2}{4s}$$

展开分子：
$$\text{分子} = (s + m_i^2 - m_j^2)^2 - 4s m_i^2$$
$$= s^2 + m_i^4 + m_j^4 + 2s m_i^2 - 2s m_j^2 - 2 m_i^2 m_j^2 - 4s m_i^2$$
$$= s^2 + m_i^4 + m_j^4 - 2s m_i^2 - 2s m_j^2 - 2 m_i^2 m_j^2$$

这正是 $\lambda(s, m_i^2, m_j^2)$ 的定义。所以：

$$|\vec{p}_i| = \frac{\sqrt{\lambda(s, m_i^2, m_j^2)}}{2\sqrt{s}}$$

### 结论
您给出的那个根号表达式：
1.  物理意义是：**质心系能量的平方减去质量平方（$E_i^2 - m_i^2$）**。
2.  结果就是：**质心系三维动量的大小**。


