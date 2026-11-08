# 有限 RS 背景下的凝聚与标量输运约定

2026-11-07 核对 Arena [#313](https://github.com/w5851/Julia_RelaxTime/issues/313)，源码基于 `714079df5629c267062b3bdf962ccadde4fb4a1b`。作者已接受保留现行计算定义，issue 已按[审计结论](https://github.com/w5851/Julia_RelaxTime/issues/313#issuecomment-6706495944)关闭。本页是五项审计的公式与裁决入口；数值算法、基线与 v15 数据保持原值。本文档及配套测试尚待提交集成。

## 当前定义

| 对象 | 现行实现 | 本次裁决 |
| --- | --- | --- |
| 背景场 | `GapSolver.gap_residual` 对 `omega` 作 AD；RS 只进入热项 | 保留 RS 变分态，场变量记为 `phi_Omega` |
| 介子顶点 | `AFieldBuilder` 的完整 `A_aniso`；`calculate_G_from_A=-Nc*M*A/(4pi²)` | 顶点输入记为 `phi_A=Nc*M*A/(4pi²)`，helper 是 `H=-phi_A` |
| η、σ、ζ | 普通在壳核乘 `f_xi(1-f_xi)`；物种平均 τ | 占据加权标量 RTA 量；有限 ξ 完整张量响应未实现 |
| κ | 普通在壳核乘 `f_xi` 和两侧 `q-nE/w`；另有 `1/T²` 归一化 | 保留现有经典形权重和标量投影，不与 η/σ 的统计权重混写 |
| 焓 | workflow 给出 `epsilon=-P+Ts+mu*n`、`P=-Omega`，`w=epsilon+P` | 热力学标量焓；不是从动理学 EMT 得到的方向焓 |
| A 的新约化 | `A_aniso_reduced` 保留变换后的角权重，有限上限另分支处理 | 与同一 A 积分的算法等价；不消除两种凝聚定义之差 |

### 1. 变分凝聚与动理学凝聚

固定 `M,T,mu,Phi,Phibar`，令 `F(E)=f+(E)+f-(E)`，真空积分固定同一物理 cutoff。两种凝聚为

\[
\phi_\Omega=-2N_cM\int^\Lambda_{\bf p}\!\frac1{E}
 +2N_cM\int_{\bf p}\!\frac{F(E_\xi)}{E_\xi},\qquad
\phi_A=-2N_cM\int^\Lambda_{\bf p}\!\frac1{E}
 +2N_cM\int_{\bf p}\!\frac{F(E_\xi)}E.
\]

因此

\[
\phi_\Omega-\phi_A=-N_cM\xi\int_{\bf p}
 \frac{p^2c^2}{E^3}F(E)+O(\xi^2)
=-\frac{N_cM\xi}{3}\int_{\bf p}\frac{p^2}{E^3}F(E)+O(\xi^2).
\]

热件之比的一阶差为 `-xi*<v²>/6`，其中平均权重是 `F(E)/E`。这是热件的相对差，不是完整负凝聚的相对差，也不是最终 τ 的误差估计。在相同积分支持和归一化下 ξ=0 两定义相等；有限热 cutoff 不同则还会留下独立的尾部差。

纯径向函数有精确换元 `q_perp=p_perp, q_L=sqrt(1+xi)p_L`，故 `∫g(E_xi)d³p=(1+xi)^(-1/2)∫g(sqrt(q²+M²))d³q`，适用于全空间或正确变换后的积分域。`PNJLModel` 的 `E_dist=sqrt(q²+M²)` 是同一个分布自变量，不能更名为原物理 p 下的 `E_kin`；本次增加说明注释。带额外 `1/E_kin` 的 A 也可以解析做角积分，但权重不是简单 Jacobian。#313 的“必须直接二维积分”被已完成的 A 约化算法取代。

**传播评估。** 对不混合通道 `D=1-4K_M Pi`，固定背景的顶点扰动满足
`delta D=-4Pi delta K_M`，简单根 `delta z=4Pi delta K_M/(partial_z D)`。混合通道以左右本征向量给出 `delta z=-l^T delta D r/(l^T partial_z D r)`。例如 `K_pi,P=G-K phi_A,s/2`，故改变凝聚约定给出 `delta K=-K delta phi_s/2`。再由 `Gamma_i=∑n_j R_ij` 得

\[
\frac{\delta\tau_i}{\tau_i}=-\frac{\sum_j(\delta n_j R_{ij}+n_j\delta R_{ij})}{\Gamma_i},
\qquad \delta X=\sum_i(\delta\tau_i J_i+\tau_i\delta J_i).
\]

改变 gap 方程本身时还需 `delta X_bg=-J_gap^{-1} delta F_gap`。这些系数依赖根导数、混合间隙、背景 Jacobian 和碰撞通道；近阈值或近临界可放大。因而没有从“凝聚热件约 0.1%”推得“输运约 0.1%”的普适结论。本批选择保留现行模型；备选 gap/顶点自洽模型及其整链传播属于新的模型研究，不伪报为已数值量化。

### 2. 标准核与一个明确的固定 RS 响应问题

在固定质量、固定 τ、固定变形背景和均匀弱电场下，忽略维持背景的源对扰动的变化，设 `f=F(E_xi)`。用 `W=-T F'`，动量导数为

\[
\partial_{p_j} f_\xi=-\frac{W(E_\xi)}{T}
\frac{p_j+\xi p_L n_j}{E_\xi}.
\]

普通 FD 极限 `W=f(1-f)`。由 `delta f=-tau q E_j partial_pj f` 与普通在壳电流速度得到

\[
\sigma_\perp^{\rm src}=\frac{gq^2\tau}{T}\int_{\bf p}
\frac{p^2(1-c^2)}{2EE_\xi}W(E_\xi),\qquad
\sigma_\parallel^{\rm src}=\frac{gq^2\tau}{T}\int_{\bf p}
\frac{p^2c^2(1+\xi)}{EE_\xi}W(E_\xi).
\]

这是已经说明源、当前速度与碰撞近似的一个响应问题；固定 RS 分布不是无源 Boltzmann 碰撞项的普遍定态，故不称为唯一 RS 电导。普通 FD 下，取迹平均与当前 `sigma_occ` 的几何差为

\[
\frac{{\rm tr}\sigma^{\rm src}}3-\sigma_{\rm occ}
=\frac{gq^2\tau\xi}{9T}\int_{\bf p}
\frac{p^2}{E^2}\left(1-\frac{p^2}{2E^2}\right)f_0(1-f_0)+O(\xi^2).
\]

无质量 FD、常 τ、固定 T,mu 时，归一到 ξ=0：

| 量 | 一阶系数 |
| --- | ---: |
| 当前 `sigma_occ` | -1/2 |
| `sigma_perp_src` | -2/5 |
| `sigma_parallel_src` | -1/5 |
| `tr(sigma_src)/3` | -1/3 |
| 迹平均减当前量 | +1/6 |
| 当前 `eta_occ` | -5/6 |

这些是固定背景解析校验，不是 v15 完整计算的百分比修正。

**PNJL 额外区别。** #313 写下的 FD 导数恒等式不直接适用于颜色平均的 PNJL 占据。令 `z=exp(-(E-mu)/T)`、`D=1+3Phi z+3Phibar z²+z³`，`f=(Phi z+2Phibar z²+z³)/D`，则

\[
W_\Phi=-T\partial_E f=\frac{\Phi z+4\bar\Phi z^2+3z^3}{D}-3f^2.
\]

`Phi=Phibar=1` 时 `W=f(1-f)`；二者为零时 `W=3f(1-f)`。差别在 ξ=0 也存在。本次只澄清现有统计闭合的定义；完整颜色响应及碰撞算符需共同选择，不能只换这个权重就宣称全部响应已修正。

### 3. 张量与标量焓

轴对称电导/扩散一般为 `sigma_perp(delta_ij-n_i n_j)+sigma_parallel n_i n_j`。剪切粘滞是四阶张量，不是 issue 简写的二阶 `eta_xx,eta_zz`。当前标量投影及热力学焓适用于已定义的占据加权量。重做有限 RS 输运需给出背景源、frame、守恒匹配和张量基，不能用 `(2P_T+P_L)/3` 替换 P 就视为完成。

### 4. 平衡恒等式的解析反例

无质量经典气体，`a=sqrt(1+xi c²)`，角平均记作 `< >`：

\[
n_\xi/n_0=\langle a^{-3}\rangle=(1+\xi)^{-1/2},\quad
\bar P_\xi/P_0=\langle a^{-4}\rangle=1-2\xi/3+O(\xi^2),\quad
\bar P_\xi/(n_\xi T)=1-\xi/6+O(\xi^2).
\]

方向压强一阶为 `P_L/P0=1-6xi/5`、`P_T/P0=1-2xi/5`。这是实际在壳 EMT 的压强；变分热势的标量 `P_Omega=n_xi T` 在该特例仍成立。两者区别使各向同性 Euler/LL 的最终积分恒等式不能直接搬到 RS 态。Das 2022 本身采用经典占据；其 Eq. (40) 是经过平衡热力学与流体关系重写的源项，不应把 issue 的近似编号当作原文逐式对应。

### 5. 标量剪切提取

令 `Q_ij=p_i p_j-delta_ij p²/3`，标准各向同性标量定义相当于对五个正交归一、无迹对称基作平均：`sum_A(e_A:Q)²=2p⁴/3`，由此得到 `p⁴/(15TE²)`。这在任何占据上定义一个标量，但不证明背景的真实响应张量也由它完全决定。

固定 FD 背景的简单占据张量投影在无质量常 τ 极限的 `xy`、`xz`、轴向无迹分量一阶系数分别为 `-5/14,-15/14,-55/42`，按 `(2,2,1)/5` 加权给出 `-5/6`。这展示了提取方案的区别。Thakur 2017 Eq. (25) 含其推导的额外项；issue 中括号的压缩重写不足以作为系数证据。本文通过显式源和投影定义比较，不指定一篇文献的标量为唯一通用系数。

## 验证与任务裁决

在 quote-remote 的独立目录 `D:/CodexCompute/arena313-conventions-20261107`，Julia 1.12.5 完成已有能量公式映射 3/3 与新增解析极限/PNJL 导数检查 29/29。新增检查覆盖 Jacobian、经典反例、FD 电导分量、剪切投影、PNJL 导数以及有质量背景下两凝聚差的斜率与 A 积分归一化；不调用 gap 或输运生产扫描。中心差分使用 h=2e-5，另以确切解析系数与 AD 比较，未放宽基线容差。

五项均已有当前源码、公式和计算定义的裁决。保留现行混合模型与标量输出，论文按定义表述；不要求重跑全网格。尚未完成的是备选自洽凝聚/完整颜色张量响应的模型研究与整链数值传播，这些不作为现行算法 bug，也不假称已由本次 32 个检查解决。它们的重启条件是作者确需将有限 ξ 输入直接嵌入张量演化或比较另一套 gap 闭合。

文献：Romatschke--Strickland, PRD 68, 036004 (2003), hep-ph/0304092；Albright--Kapusta, PRC 93, 014903 (2016), arXiv:1508.02696；Thakur et al., PRD 95, 096009 (2017), arXiv:1703.03142；Das et al., PRD 106, 014013 (2022), arXiv:2109.01543。
