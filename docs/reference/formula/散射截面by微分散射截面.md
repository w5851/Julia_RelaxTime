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

对应的微分截面公式采用 \(s_{ij}^\pm = s - (m_i \pm m_j)^2\)，其分母组合满足 \(s_{ij}^+ s_{ij}^- = \lambda(s,m_i^2,m_j^2)\)，后续可直接以 Kall\'en 函数形式代入。

当分布函数各向异性时，需要显式注明 \(f\) 可能依赖末态的出射角，例如 \(f = f(E, T, \mu, \theta^*)\)，其中 \(\theta^*\) 是质心系下的散射角。

在 PNJL 各向异性模型中，末态分布函数可直接写为

```math
f = f(E, T, \mu, \xi, \cos\theta^*) \equiv \texttt{quark\_distribution\_aniso}\bigl(p, m, \mu, T, \Phi, \bar\Phi, \xi, \cos\theta^*\bigr),
```

其中 \(p = \sqrt{E^2 - m^2}\)。对反粒子需使用相应的反夸克分布；末态粒子 d 的动量方向与 c 相反，因此可用 \(-\cos\theta^*\)。

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

引入 Kall\'en 函数 \(\lambda(x,y,z)=x^2+y^2+z^2-2xy-2xz-2yz\) 可将上式改写为更紧凑的形式：

```math
t_\pm = m_i^2 + m_c^2 - \frac{1}{2s}\Big[(s+m_i^2-m_j^2)(s+m_c^2-m_d^2)
\mp \sqrt{\lambda(s,m_i^2,m_j^2)}\,\sqrt{\lambda(s,m_c^2,m_d^2)}\Big].
```

质心系下的散射角 \(\theta^*\) 与 \(t\) 的对应关系为

```math
t(\theta^*) = m_i^2 + m_c^2 - \frac{1}{2s}\Big[(s+m_i^2-m_j^2)(s+m_c^2-m_d^2)
- \sqrt{\lambda(s,m_i^2,m_j^2)}\,\sqrt{\lambda(s,m_c^2,m_d^2)}\cos\theta^*\Big],
```

其中 \(\theta^*\in[0,\pi]\)，当 \(\cos\theta^*=\pm1\) 时分别给出 \(t_\mp\)。反解可得到

```math
\cos\theta^* = \frac{(s+m_i^2-m_j^2)(s+m_c^2-m_d^2) - 2s\,(m_i^2+m_c^2-t)}
{\sqrt{\lambda(s,m_i^2,m_j^2)}\,\sqrt{\lambda(s,m_c^2,m_d^2)}}.
```

因此在数值积分中可任选 \(t\) 或 \(\theta^*\) 作为积分变量；若相空间或分布函数对角度敏感（如各向异性分布），需在被积函数中保留对 \(\theta^*\) 的依赖。

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

