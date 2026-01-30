# 磁场下的PNJL模型相关部分

## 一、有磁场情况下的PNJL模型

### 1. 磁场下热力学势

\[
\Omega = \sum_{f=u,d,s} \left( \Omega_f^0 + \Omega_f^T \right) + 2G (\phi_u^2 + \phi_d^2 + \phi_s^2) - 4K \phi_u \phi_d \phi_s + U(\Phi, \bar{\Phi}, T)
\]

其中：

- 真空项：
\[
\Omega_f^0 = -N_c \frac{|q_f| eB}{2\pi} \sum_{n=0}^{\infty} \alpha_n \int \frac{dp_z}{2\pi} E_{f,n}
\]

- 温度项：
\[
\Omega_f^T = -T \frac{|q_f| eB}{2\pi} \sum_{n=0}^{\infty} \alpha_n \int \frac{dp_z}{2\pi} \left( Z_f^+ + Z_f^- \right)
\]

### 2. 夸克能级与密度

\[
E_{f,n} = \sqrt{2n|q_f|B + p_z^2 + M_f^2}
\]
\[
\rho_{f,n} = 3 \frac{q_f B}{2\pi} \alpha_n \int \frac{dp_z}{2\pi} \left[ f_{f,n}^+ - f_{f,n}^- \right]
\]

分布函数类似，但含 Landau 能级指标 \(n\)。