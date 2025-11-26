# 散射截面与微分散射截面的关系（修订版）

该文档计算的是 **过程 $i j \to c d$** 的总散射截面，说明如何由微分散射截面积分得到总散射截面，并对各符号给出说明。

## 公式

总散射截面可以写为对微分散射截面的积分：

```math
\sigma(s, T, \mu_q)
=\int_{t_-}^{t_+} \mathrm{d}t\; \frac{\mathrm{d}\sigma}{\mathrm{d}t}\,
\bigl(1 - f(E_c, T, \mu)\bigr)\bigl(1 - f(E_d, T, \mu)\bigr)\,.
\tag{5.13}
```

**注意**：本项目只考虑夸克-夸克散射（费米子），因此统一使用 Pauli blocking 因子 $(1-f)$，不涉及玻色子的 Bose enhancement $(1+f)$。

其中上下限 $t_\pm$（Mandelstam 变量 $t$ 的取值范围）为：

```math
\begin{aligned}
t_\pm &= m_i^2 + m_c^2 
- \frac{1}{2s}\,(s + m_i^2 - m_j^2)(s + m_c^2 - m_d^2)\\[6pt]
&\qquad\pm 2\sqrt{\left(\frac{(s + m_i^2 - m_j^2)^2}{4s} - m_i^2\right)
\left(\frac{(s + m_c^2 - m_d^2)^2}{4s} - m_c^2\right)}\,.
\end{aligned}
\tag{5.14}
```

---

## 符号说明（针对 $i j \to c d$ 过程）

- $i,j$：初态入射粒子（incoming particles）。  
- $c,d$：末态产物粒子（outgoing particles）。  
- $m_i,m_j,m_c,m_d$：分别为粒子 $i,j,c,d$ 的质量。  
- $s$：Mandelstam 变量 $s$，即总能量平方（中心质心系下通常为 $(p_i+p_j)^2$）。
- $t$：Mandelstam 变量 $t$（四动量转移的平方）。
- $\dfrac{\mathrm{d}\sigma}{\mathrm{d}t}$：对 $t$ 的微分散射截面（由矩阵元计算得到）。
- $f(E, T, \mu)$：Fermi-Dirac 分布函数，表示能量为 $E$ 的费米子的占据数，温度为 $T$，化学势为 $\mu$。
- $(1 - f)$：Pauli blocking 因子，费米子不能占据已被占据的态。
- $E_c, E_d$：末态粒子 $c,d$ 的能量（在相应自变量下确定）。

---

