在球坐标系下，考虑轴对称各向异性分布函数 $f_i(p_i, \cos\theta_i)$ 和 $f_j(p_j, \cos\theta_j)$，其中 $\theta_i$ 和 $\theta_j$ 为动量与某一固定轴（如 $z$ 轴）的夹角。通过对两个方位角 $\phi_i$ 和 $\phi_j$ 的积分进行简化，利用被积函数仅依赖于其差值 $\phi = \phi_i - \phi_j$ 的性质，可得平均散射率的简化表达式为：

$$
\omega_{ij} = \frac{d_q^2}{\rho_i \rho_j (2\pi)^5} \int_0^\Lambda p_i^2 dp_i \int_0^\Lambda p_j^2 dp_j \int_{-1}^1 d\cos\theta_i \int_{-1}^1 d\cos\theta_j \int_0^{2\pi} d\phi \ f_i(p_i,\cos\theta_i) f_j(p_j,\cos\theta_j) \ v_{\text{rel}} \ \sigma_{ij\to cd}(s,T,\mu_q)
$$

其中：
- $d_q = 6$ 为简并度（PNJL 模型）。
- $\Lambda$ 为动量模的截断。
- $s$ 为 Mandelstam 变量，$s = (p_i + p_j)^2 = m_i^2 + m_j^2 + 2E_iE_j - 2p_i p_j \cos\Theta$，其中 $\cos\Theta = \cos\theta_i \cos\theta_j + \sin\theta_i \sin\theta_j \cos\phi$。
- $E_i = \sqrt{p_i^2 + m_i^2}$，$E_j = \sqrt{p_j^2 + m_j^2}$（或根据模型的有效质量）。
- $v_{\text{rel}}$ 为相对速度，通常取 $v_{\text{rel}} = \frac{\sqrt{(p_i \cdot p_j)^2 - m_i^2 m_j^2}}{E_i E_j}$，其中 $p_i \cdot p_j = E_iE_j - p_i p_j \cos\Theta$。
- 数密度 $\rho_i$ 和 $\rho_j$ 由下式给出：
  $$
  \rho_i = d_q \frac{1}{(2\pi)^2} \int_0^\infty p_i^2 dp_i \int_{-1}^1 d\cos\theta_i \ f_i(p_i,\cos\theta_i)
  $$
  类似地计算 $\rho_j$。

**重要说明**：数密度积分应使用半无穷积分范围 $[0, \infty)$，而非有限截断 $[0, \Lambda]$。这是因为数密度是物理可观测量，不应依赖于模型的动量截断参数。在数值实现中，可通过变量替换 $p = \text{scale} \cdot t / (1-t)$ 将半无穷积分映射到有限区间 $[0, 1)$。

该表达式已将两个方位角积分简化为一个方位角 $\phi$ 的积分，保留了必要的相对方向依赖性。
# 利用对称性的简化
基于分布函数的偶函数性质 \( f(\cos\theta) = f(-\cos\theta) \) 以及被积函数在联合变换 \( (\cos\theta_i, \cos\theta_j) \to (-\cos\theta_i, -\cos\theta_j) \) 和 \( \phi \to -\phi \) 下的对称性，平均散射率可简化为：

$$
\omega_{ij} = \frac{d_q^2}{4 \pi^5 \rho_i \rho_j} \int_0^\Lambda p_i^2 dp_i \int_0^\Lambda p_j^2 dp_j \int_{0}^1 d\cos\theta_i \int_{0}^1 d\cos\theta_j \int_{0}^{\pi} d\phi \ f_i(p_i,\cos\theta_i) f_j(p_j,\cos\theta_j) \ v_{\text{rel}} \ \sigma_{ij\to cd}(s,T,\mu_q)
$$

其中：
- \( d_q = 6 \)。
- 截断 \( \Lambda \) 应用于动量模积分。
- \( s = (p_i + p_j)^2 = m_i^2 + m_j^2 + 2E_iE_j - 2p_i p_j \cos\Theta \)，\( \cos\Theta = \cos\theta_i \cos\theta_j + \sin\theta_i \sin\theta_j \cos\phi \)。
- 相对速度 \( v_{\text{rel}} = \frac{\sqrt{(p_i \cdot p_j)^2 - m_i^2 m_j^2}}{E_i E_j} \)，其中 \( p_i \cdot p_j = E_iE_j - p_i p_j \cos\Theta \)。
- 数密度 \( \rho_i \) 和 \( \rho_j \) 简化为：
  $$
  \rho_i = \frac{d_q}{2\pi^2} \int_0^\infty p_i^2 dp_i \int_{0}^1 d\cos\theta_i \ f_i(p_i,\cos\theta_i)
  $$
  同理计算 \( \rho_j \)。

**重要说明**：数密度积分应使用半无穷积分范围 $[0, \infty)$，而非有限截断 $[0, \Lambda]$。

该表达式充分利用了对称性，将对 \( \cos\theta_i \) 和 \( \cos\theta_j \) 的积分区域缩减至 \( [0,1] \)，对方位角 \( \phi \) 的积分缩减至 \( [0,\pi] \)，减少了计算量。