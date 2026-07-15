# Kappa 系列与热导率 lambda

## 文档状态

- 来源类型：一手论文 + legacy Fortran 对照
- 理论文献：Arpan Das, Hiranmaya Mishra, Ranjita K. Mohapatra, Phys. Rev. D 106, 014013 (2022)
- 代码来源：`Relaxtime_fortran/codes/main/coefficients_tmu.f90`
- 当前用途：为 Julia 的 `kappa_BB`、`kappa_QQ`、`kappa_SS`、完整扩散矩阵和 `lambda` 提供正式公式锚点，并说明与 legacy 实现的对应关系
- 限制：Das 2022 在强子 HRG 框架下推导扩散矩阵；本仓库落地时使用的是夸克侧 RTA 实现，因此应区分“理论结构一致”和“模型细节完全相同”
- 各向异性能量约定：$E_{\mathrm{kin}}=\sqrt{p^2+M^2}$ 用于运动学分母与 Landau-Lifshitz 投影，$E_{\mathrm{dist}}$ 只作为 RS 分布自变量

---

## 1. 适用范围

本文档记录守恒荷扩散输运家族：

- $\kappa_{BB}$
- $\kappa_{QQ}$
- $\kappa_{SS}$
- $\lambda$
- 以及其下游派生量 $\mathrm{Lorentz}$、$\mathrm{Pr}$、$R_{\eta\sigma}$、$\zeta_\eta$、$RI$

这些量属于 conserved-charge transport family，不属于静态涨落或磁化率。

---

## 2. Das 2022 的理论框架

Das 2022 在多守恒荷系统中给出的广义 Fick 定律是：

$$
\left( \begin{array}{c}
\Delta J_B^i \\
\Delta J_Q^i \\
\Delta J_S^i
\end{array} \right)=
\left( \begin{array}{ccc}
\kappa_{BB} & \kappa_{BQ} & \kappa_{BS} \\
\kappa_{QB} & \kappa_{QQ} & \kappa_{QS} \\
\kappa_{SB} & \kappa_{SQ} & \kappa_{SS}
\end{array} \right)
\left( \begin{array}{c}
D^i \alpha_B \\
D^i \alpha_Q \\
D^i \alpha_S
\end{array} \right)
$$

其中：

- $\alpha_q = \mu_q / T$
- $\omega = \epsilon + P$ 是焓密度
- $D^\mu = \partial^\mu - u^\mu (u \cdot \partial)$ 是正交于流体四速度的导数

这说明：

- 对角元 `BB / QQ / SS` 描述各守恒荷对自身驱动力的扩散响应
- 非对角元 `BQ / BS / QS` 描述不同守恒荷之间的交叉扩散耦合

Julia 当前已公开完整的 $B/Q/S$ 对称扩散矩阵接口；legacy Fortran 只显式实现了三个对角元。

---

## 3. 扩散矩阵元的一般公式

Das 2022 的核心结果是公式 (55)：

$$
\kappa_{qq'}
= \sum_a \int \frac{d^3 p_a}{(2\pi)^3} \frac{p_a^2}{3 E_a^2} \left(q_a - \frac{n_q E_a}{\omega}\right) \tau_a \left(q'_a - \frac{n_{q'} E_a}{\omega}\right) f_a^{(0)}
$$

这里：

- $a$ 是粒子种类求和指标
- $q, q' \in \{B, Q, S\}$ 是守恒荷标签
- $q_a, q'_a$ 是粒子 $a$ 携带的相应守恒荷量子数
- $n_q, n_{q'}$ 是对应守恒荷密度
- $\tau_a$ 是弛豫时间
- $f_a^{(0)}$ 是平衡分布函数

这个结构有三个对本项目最重要的结论：

1. $\kappa$ 系列是输运系数，因为显式依赖 $\tau_a$
2. 投影项 $q_a - n_q E_a / \omega$ 来自 Landau-Lifshitz 条件，不是静态 susceptibilities 的导数结构
3. 理论上完整对象是一个扩散矩阵，而不是只有三个彼此无关的标量

### 3.1 RS 各向异性下的能量职责

本仓库采用“各向异性只进入分布函数，不重新定义准粒子色散关系”的统一约定：

$$
E_{\mathrm{kin},a}=\sqrt{p^2+M_a^2},
\qquad
E_{\mathrm{dist},a}=\sqrt{p^2+M_a^2+\xi(p\cos\theta)^2}.
$$

在完整 RS 角积分下，Das 2022 Eq. (55) 的项目实现目标写成

$$
\kappa_{qq'}^{(\xi)}
= \sum_a \int \frac{d^3 p_a}{(2\pi)^3}
\frac{p_a^2}{3E_{\mathrm{kin},a}^2}
\left(q_a-\frac{n_qE_{\mathrm{kin},a}}{\omega}\right)
\tau_a
\left(q'_a-\frac{n_{q'}E_{\mathrm{kin},a}}{\omega}\right)
f_{a,\xi}.
$$

其中 $f_{a,\xi}$ 以 $E_{\mathrm{dist},a}$ 为 RS 变形自变量。质量 $M_a=M_a(T,\mu_B,\xi)$ 可以来自当前各向异性热力学背景；这里不允许的是把 $\xi(p\cos\theta)^2$ 额外加入准粒子色散。

这项 RS 扩展是本仓库的物理约定：Das 2022 Eq. (55) 提供扩散矩阵、普通 $E_a$ 分母和 Landau-Lifshitz 投影结构，但该文并未直接推导本仓库的完整全角 RS 形式，不应把两者表述为完全相同的公式。

当前 `TransportCoefficients.jl` 已在共享扩散 kernel 中使用 $E_{\mathrm{kin}}$ 计算 $p^4/E^2$ 与 Landau-Lifshitz 投影，并只用 $E_{\mathrm{dist}}$ 生成 RS 占据分布。旧 production 不应被描述为已经采用本节的能量分工，其状态由外部 production registry 记录。

---

## 4. 与 legacy Fortran 的对应关系

在本项目的 legacy Fortran 中，粒子种类求和被具体化为：

- $j \in \{q, \bar q\}$：夸克 / 反夸克分支
- $f \in \{u, d, s\}$：味道

并使用：

$$
E_f(p) = \sqrt{p^2 + M_f^2}
$$

$$
F_{q,f}(p) = f1f(E_f - \mu_f, T, \Phi, \bar{\Phi})
$$

$$
F_{\bar q,f}(p) = f1fb(E_f + \mu_f, T, \Phi, \bar{\Phi})
$$

因此 Das 2022 的一般公式在 legacy 口径下可写成：

$$
\kappa_{q q'}
= \sum_{j,f} \frac{N_c\,\tau_{j f}}{3\pi^2} \int_0^{\infty} dp\, \frac{p^4}{E_f(p)^2} F_{j f}(p) \left(q_{j f} - \frac{n_q E_f(p)}{\epsilon + P}\right) \left(q'_{j f} - \frac{n_{q'} E_f(p)}{\epsilon + P}\right)
$$

其中：

- $N_c = 3$
- $q_{j f}, q'_{j f}$ 在代码中分别对应 `arr_quark_baryon`、`arr_quark_charge`、`arr_quark_strange`
- $F_{j f}(p)$ 对应 `arrFF` 一类占据核

当 $q = q'$ 时，就得到对角元；这正是 legacy 已实现的 $\kappa_{BB}$、$\kappa_{QQ}$、$\kappa_{SS}$。

---

## 5. 三个对角元

### 5.1 $\kappa_{BB}$

$$
\kappa_{BB}
= \sum_{j,f} \frac{N_c\,\tau_{j f}}{3\pi^2} \int_0^{\infty} dp\, \frac{p^4}{E_f(p)^2} F_{j f}(p) \left(B_{j f} - \frac{n_B E_f(p)}{\epsilon + P}\right)^2
$$

它描述的是重子守恒流对 $D^i \alpha_B$ 的扩散响应。

### 5.2 $\kappa_{QQ}$

按 Das 2022 的矩阵公式，对角电荷元应写成：

$$
\kappa_{QQ}
= \sum_{j,f} \frac{N_c\,\tau_{j f}}{3\pi^2} \int_0^{\infty} dp\, \frac{p^4}{E_f(p)^2} F_{j f}(p) \left(Q_{j f} - \frac{n_Q E_f(p)}{\epsilon + P}\right)^2
$$

legacy Fortran 还额外乘了 $4\pi \alpha_{em}$。这应理解为实现层的电流归一化约定，而不是扩散矩阵结构本身的改变。因此 Julia 落地时应把：

- 扩散矩阵本体 $\kappa_{QQ}$
- 与电磁流归一化相关的实现约定

明确区分。

Das 2022 还指出，在 $\mu_B = 0$ 时，$\kappa_{QQ} / T^2 \equiv \sigma_{el} / T$；这给出了 $\kappa_{QQ}$ 与电导率的一个重要检查点，但并不意味着二者在一般化学势下恒等。

### 5.3 $\kappa_{SS}$

$$
\kappa_{SS}
= \sum_{j,f} \frac{N_c\,\tau_{j f}}{3\pi^2} \int_0^{\infty} dp\, \frac{p^4}{E_f(p)^2} F_{j f}(p) \left(S_{j f} - \frac{n_S E_f(p)}{\epsilon + P}\right)^2
$$

它描述奇异守恒荷扩散通道，对 Julia 来说与 $\kappa_{BB}$ 共用同一积分骨架，只需更换投影量。

---

## 6. 非对角元

Das 2022 明确表明，理论上的扩散矩阵还包含：

- $\kappa_{BQ} = \kappa_{QB}$
- $\kappa_{BS} = \kappa_{SB}$
- $\kappa_{QS} = \kappa_{SQ}$

其统一写法就是上一节的一般公式，只需取 $q \neq q'$：

$$
\kappa_{q q'}
\propto
\left(q_{j f} - \frac{n_q E_f}{\epsilon + P}\right)
\left(q'_{j f} - \frac{n_{q'} E_f}{\epsilon + P}\right)
$$

这给出一个重要的实现边界：

- legacy Fortran 只显式输出三个对角元
- Julia 当前已经公开 `kappa_BQ`、`kappa_BS`、`kappa_QS` 及完整对称 `diffusion_matrix`
- 从理论上，非对角元与对角元属于同一扩散矩阵，而不是附会出来的独立公式

因此，Issue #130 的能量职责修正必须同时覆盖对角元和非对角元，不能只修改 legacy 曾输出的三个对角通道。

把一般公式具体写开，可得到当前最值得记录的三个交叉通道：

### 6.1 $\kappa_{BQ}$

$$
\kappa_{BQ}
= \sum_{j,f} \frac{N_c\,\tau_{j f}}{3\pi^2} \int_0^{\infty} dp\, \frac{p^4}{E_f(p)^2} F_{j f}(p) \left(B_{j f} - \frac{n_B E_f(p)}{\epsilon + P}\right) \left(Q_{j f} - \frac{n_Q E_f(p)}{\epsilon + P}\right)
$$

### 6.2 $\kappa_{BS}$

$$
\kappa_{BS}
= \sum_{j,f} \frac{N_c\,\tau_{j f}}{3\pi^2} \int_0^{\infty} dp\, \frac{p^4}{E_f(p)^2} F_{j f}(p) \left(B_{j f} - \frac{n_B E_f(p)}{\epsilon + P}\right) \left(S_{j f} - \frac{n_S E_f(p)}{\epsilon + P}\right)
$$

### 6.3 $\kappa_{QS}$

$$
\kappa_{QS}
= \sum_{j,f} \frac{N_c\,\tau_{j f}}{3\pi^2} \int_0^{\infty} dp\, \frac{p^4}{E_f(p)^2} F_{j f}(p) \left(Q_{j f} - \frac{n_Q E_f(p)}{\epsilon + P}\right) \left(S_{j f} - \frac{n_S E_f(p)}{\epsilon + P}\right)
$$

在夸克侧，这些非对角元不是形式上更复杂的新对象，而是同一个参数化 kernel 在不同守恒荷投影下的混合项。当前 Julia 实现已在共享层完成 $E_{\mathrm{kin}}/E_{\mathrm{dist}}$ 拆分，并由全部 $BQ/BS/QS$ 通道继承。

---

## 7. $\lambda$

legacy Fortran 并未单独构造新的动量积分核，而是由 `kappa_BB` 推出热导率：

$$
\lambda
= \kappa_{BB}
\left(\frac{\epsilon + P}{n_B T}\right)^2
$$

这说明：

- $\lambda$ 与 $\kappa_{BB}$ 是同一输运家族的上下游量
- Julia 实现应先稳定实现 $\kappa_{BB}$
- 在 $n_B \to 0$ 附近需要做数值稳定性处理

Das 2022 本身聚焦扩散矩阵，没有直接把 $\lambda$ 作为主结果展开，因此 $\lambda$ 当前仍主要依赖 legacy 关系式作为工程实现锚点。

不过，结合 [docs/reference/DAS_kappa_formula_temp.md](docs/reference/DAS_kappa_formula_temp.md) 中整理的说明，可以把这条关系的框架语义再说清楚：

- Das 2022 处理的是 Landau-Lifshitz 框架下的守恒荷扩散矩阵，而不是独立热流方程
- 在这一框架下，热流与重子扩散流并不是完全独立的 dissipative 自由度
- 因而 $\lambda$ 可以被视为由 $\kappa_{BB}$ 和热力学量组合出来的派生热输运量，而不是另一个与 $\kappa_{BB}$ 并列、需要单独做新动量积分的基础系数

因此，当前更稳妥的表述是：

- $\kappa_{BB}$ 有 Das 2022 提供的一手扩散矩阵公式支撑
- $\lambda = \kappa_{BB} ((\epsilon+P)/(n_B T))^2$ 没有在 Das 2022 文中直接写出
- 但它与 Landau-Lifshitz 框架下“热流可由重子扩散流和热力学关系导出”的物理图像相容，并与 legacy 实现一致

需要区分两个名为 energy 的对象：$\lambda$ 公式中的 $\epsilon$ 是热力学能量密度，不是扩散积分中的单粒子能量。`lambda_from_kappa_BB` 不直接选择 $E_{\mathrm{kin}}$ 或 $E_{\mathrm{dist}}$，但会继承 $\kappa_{BB}$ 的能量语义；当前实现已由修正后的 $\kappa_{BB}$ 重新计算 $\lambda$。

这也解释了为什么当前 Julia 路线应当是“先补 $\kappa_{BB}$，再补 $\lambda$”，而不是把二者当作两个完全独立的主方程对象。

---

## 8. 派生量

### 8.1 标准 Lorenz 数与 legacy `Lorentz`

当前需要明确区分两种不同口径：

- 标准热电输运文献里的 Lorenz 数通常写作

$$
L = \frac{\lambda}{\sigma T}
$$

- legacy Fortran `coefficients_tmu.f90` 实际输出的 `Lorentz` 则是

$$
\mathrm{Lorentz}_{\mathrm{legacy}} = \frac{\lambda}{\sigma/T} = \frac{\lambda T}{\sigma}
$$

因此，Julia 主线实现应将二者视为两个不同的派生诊断量：

- `lorenz_number = lambda / (sigma * T)`：标准 Lorenz 数
- `lorentz_legacy = lambda / (sigma / T)`：保留 legacy compare 口径

### 8.2 剪切粘滞与电导比

$$
R_{\eta\sigma} = \frac{\eta/s}{\sigma/T}
$$

这里并不是常见教科书中的标准命名量，而更接近一个“粘滞-导电耦合诊断比值”。Julia 侧当前将其视为 legacy-compatible derived diagnostic。

### 8.3 Prandtl 数

$$
\mathrm{Pr} = \frac{\eta\,c_p}{\lambda\,\rho_{\text{mass}}}
$$

legacy 代码中这两个背景量的口径可以进一步明确为：

$$
\rho_{\text{mass}} = \sum_{f \in \{u,d,s\}} \left(n_f + n_{\bar f}\right) M_f
$$

即把夸克与反夸克的数密度分别乘对应有效质量后求和。Fortran `quantity.f90` 中对应的是：

$$
\rho_{\text{mass}} = \sum_f N_q(f) M_f + \sum_f N_{\bar q}(f) M_f
$$

对于 `c_p`，legacy `coefficients_tmu.f90` 直接使用的代码级定义是：

$$
c_p^{\mathrm{legacy}} = T \left( \frac{d s}{d T} - \frac{s}{n_B} \frac{d s}{d \mu_B} \right)
$$

这里保留 legacy 记号 `c_p`，因为其确实被用于 Prandtl 数的背景量；但当前更稳妥的表述是“legacy transport background coefficient”，而不是在没有进一步文献锚点时直接把它强断言为某个标准教材口径。

Julia 当前已实现：

- `rho_mass_from_densities(...)`：按 legacy 质量密度公式计算 `\rho_{\text{mass}}`
- `legacy_transport_c_p(...)`：按 legacy 导数组合计算 `c_p`

workflow 侧现已自动给出 `rho_mass`，而 `c_p` 在 `compute_bulk=true` 时可随导数计算一并获得；也允许外部显式传入两者作为背景量。

### 8.3a Julia 公共 API 边界

当前公共 API 边界建议明确如下：

- 稳定公开的 advanced transport API 位于 `Main.TransportCoefficients` / `Main.RelaxTime`：
	- `kappa_BB`, `kappa_QQ`, `kappa_SS`
	- `kappa_BQ`, `kappa_BS`, `kappa_QS`
	- `diffusion_matrix`
	- `lambda_from_kappa_BB`
	- `lorenz_number`, `lorentz_legacy`, `viscous_conductive_coupling_ratio`, `prandtl_number`, `bulk_to_shear_viscosity_ratio`
	- `rho_mass_from_densities`
- 面向上层统一入口的公开边界位于 workflow：`solve_gap_and_transport(...).transport` 与 `solve_transport_from_equilibrium(...).transport`。
- 当前不再把 `diffusion_matrix` 与非对角元视为“内部预留但未公开”；它们属于稳定的 advanced public API。
- 但也不额外升格为一个新的 `Models` 顶层 facade 契约；主推荐入口仍是 workflow 结果或 `Main.TransportCoefficients` 模块本身。

### 8.4 粘滞比值

$$
\zeta_\eta = \frac{\zeta}{\eta}
$$

### 8.5 $RI$

legacy 中直接定义为：

$$
RI = \frac{4}{\eta / \rho_{\text{mass}}}
$$

当前仍未确认其标准物理名称与一手文献来源，因此建议继续视为低优先级派生量。

---

## 9. 输出口径

legacy 输出文件中对应口径为：

- `coe_diffusion.dat`：

$$
\kappa_{BB}/T^2,
\qquad
\kappa_{QQ}/T^2,
\qquad
\kappa_{SS}/T^2
$$

- `coe_lambda.dat`：

$$
\kappa_{BB}/T^2,
\qquad
\lambda/T^2
$$

因此 Julia 对 legacy 时，首选验证口径仍应为 `kappa_*/T^2` 与 `lambda/T^2`。

---

## 10. 与静态涨落的区别

不要把这些量与下列对象混同：

- 守恒荷广义磁化率 $\chi_{ijk}^{BQS}$
- 累积量 $C_{ijk}^{BQS}$
- 峰度观测量 $\kappa\sigma^2$

更精确地说：

- `kappa_BB / QQ / SS` 是动态输运系数
- $\chi$ 是静态热力学导数
- 二者可能通过 Einstein/Onsager 型关系间接相连，但在本项目中不能直接等同

`kappa` 家族不是静态涨落，最直接的证据就是：

- 显式含有弛豫时间 $\tau$
- 使用了 Landau-Lifshitz 投影项 $q - nE/(\epsilon + P)$
- 理论上属于扩散矩阵，而不是压强对化学势的高阶导数

静态涨落请见：

- [docs/reference/formula/models/pnjl/守恒荷广义磁化率与累积量.md](docs/reference/formula/models/pnjl/守恒荷广义磁化率与累积量.md)

---

## 11. Julia 实现建议

建议按以下顺序实现：

1. 抽象统一的 conserved-charge diffusion kernel
2. 先做对角元 `kappa_BB / QQ / SS`
3. 用稳定热力学量组合实现 `lambda`
4. 再决定是否扩展到非对角元 `BQ / BS / QS`
5. 最后暴露 `Lorentz`、`Pr` 等派生比值

实现时应单独记录：

- 单位链
- $n_B \to 0$ 时的稳定性策略
- 与 legacy 输出文件的归一化口径对应关系
