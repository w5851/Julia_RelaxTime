在球坐标系下，考虑轴对称各向异性分布函数 $f_i(p_i, \cos\theta_i)$ 和 $f_j(p_j, \cos\theta_j)$，其中 $\theta_i$ 和 $\theta_j$ 为动量与某一固定轴（如 $z$ 轴）的夹角。通过对两个方位角 $\phi_i$ 和 $\phi_j$ 的积分进行简化，利用被积函数仅依赖于其差值 $\phi = \phi_i - \phi_j$ 的性质，可得平均散射率的简化表达式为：

$$
\omega_{ij} = \frac{d_q^2}{\rho_i \rho_j (2\pi)^5} \int_0^\Lambda p_i^2 dp_i \int_0^\Lambda p_j^2 dp_j \int_{-1}^1 d\cos\theta_i \int_{-1}^1 d\cos\theta_j \int_0^{2\pi} d\phi \ f_i(p_i,\cos\theta_i) f_j(p_j,\cos\theta_j) \ v_{\text{rel}} \ \sigma_{ij\to cd}(s,T,\mu_q)
$$

其中：
- $d_q = 6$ 为简并度（PNJL 模型）。
- $\Lambda$ 为动量模的截断。
- $s$ 为 Mandelstam 变量，$s = (p_i + p_j)^2 = m_i^2 + m_j^2 + 2E_iE_j - 2p_i p_j \cos\Theta$，其中 $\cos\Theta = \cos\theta_i \cos\theta_j + \sin\theta_i \sin\theta_j \cos\phi$。
- **注意**：在计算相对速度 $v_{\text{rel}}$ 时，应使用质心系下的能量 $E_i^*$ 和 $E_j^*$，而非实验室系下的能量。具体计算公式为（参考 PDF 第 4 页公式 5.11）：
  $$
  E_i^* = \frac{s + m_i^2 - m_j^2}{2\sqrt{s}}, \quad
  E_j^* = \frac{s - m_i^2 + m_j^2}{2\sqrt{s}}
  $$
- 相对速度 $v_{\text{rel}}$ 的计算公式为（参考 PDF 第 4 页公式 5.11）：
  $$
  v_{\text{rel}} = \frac{\sqrt{(E_i^* E_j^* - |\mathbf{p}_i^*||\mathbf{p}_j^*|)^2 - (m_i m_j)^2}}{E_i^* E_j^*}
  $$
  其中 $|\mathbf{p}_i^*|$ 和 $|\mathbf{p}_j^*|$ 为质心系下的动量大小，具体表达式见 PDF。
- 数密度 $\rho_i$ 和 $\rho_j$ 由下式给出：
  $$
  \rho_i = d_q \frac{1}{(2\pi)^2} \int_0^\infty p_i^2 dp_i \int_{-1}^1 d\cos\theta_i \ f_i(p_i,\cos\theta_i)
  $$
  类似地计算 $\rho_j$。

**重要说明**：数密度积分应使用半无穷积分范围 $[0, \infty)$，而非有限截断 $[0, \Lambda]$。这是因为数密度是物理可观测量，不应依赖于模型的动量截断参数。在数值实现中，可通过变量替换 $p = \text{scale} \cdot t / (1-t)$ 将半无穷积分映射到有限区间 $[0, 1)$。

该表达式已将两个方位角积分简化为一个方位角 $\phi$ 的积分，保留了必要的相对方向依赖性。


---

**注意**：文档中先前给出的将角度积分简化到 $\cos\theta \in [0,1]$、$\phi \in [0,\pi]$ 的做法不正确。原因在于

$$
\cos\Theta = \cos\theta_i\cos\theta_j + \sin\theta_i\sin\theta_j\cos\phi
$$

并不在联合变换 $(\cos\theta_i,\cos\theta_j,\phi) \to (-\cos\theta_i,-\cos\theta_j,-\phi)$ 下保持不变，因此无法通过对称性将积分区域简单折半。数值实现必须保留原始积分范围：$\cos\theta_i,\cos\theta_j \in [-1,1]$ 和 $\phi \in [0,2\pi]$。

数密度的常数因子可等价写为：

$$
\rho_i = d_q \frac{1}{(2\pi)^2} \int_0^{\infty} p_i^2 dp_i \int_{-1}^1 d\cos\theta_i\ f_i(p_i,\cos\theta_i)
= \frac{d_q}{2\pi^2} \int_0^{\infty} p_i^2 dp_i \int_{0}^1 d\cos\theta_i\ f_i(p_i,\cos\theta_i)
$$

（已移除）此前文档中提到的“3D 各向同性特化/Fortran 风格参考实现”相关内容容易造成混淆，当前实现统一走通用公式路径。

