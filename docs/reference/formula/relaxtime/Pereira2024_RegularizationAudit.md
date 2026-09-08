# Pereira 2024：双线域、配套单线项与 retarded 对照

## 范围与结论层级

本文件说明 charged GBU 研究候选的极化函数依据与独立验证方法，不选定最终
论文 regulator，不授权 production。运行和失败证据见
[strict-audit 第22.16节](../../../dev/active/2026-08-30_charged-rpa-bu-strict-audit.md)。
完整背景与部分产额契约见 [研究方法 v1](ChargedGBU_ResearchMethod_v1.md)。

直接公式来源是 Pereira、Costa、Providência、Moreira，
“New approach to the 3-momentum regularization of the in-medium one- and
two-fermion line integrals with application to cross sections in the NJL model”，
[Phys. Rev. C 109, 025206 (2024)](https://doi.org/10.1103/PhysRevC.109.025206)。
作者提供的57页PDF页序与论文页码025206-N一致；运行manifest记录原PDF的SHA-256。

| 来源位置 | 文献事实 | 本项目如何使用 |
|---|---|---|
| 第6--7页 Eq.(26)--(28) | 极化函数的双线项与配套单线项继承同一正则化域 | 独立核验 P/S 投影和 C(q)，不调用平衡态 A 代替 |
| 第9--10页 Sec.II C、Fig.3 | 两条夸克线各自满足三动量截断，形成双球交集 | 第三个条件保留；q=0 两球重合，q>=2L 交集消失 |
| 第11--12页 Eq.(44)--(56) | 区分质量移位与外频率移位，并分开 pair/scattering | 采用外部 k0+i0 的 retarded 处方 |
| 第13页 Eq.(60)--(62) | 同域单线积分、平移后的透镜域与真空解析式 | 直接计算两个单线项，并以独立真空原函数核验 |
| 第32--33页 C6--C19 | 实部 PV；外频率移位改变两种 cut 的虚部符号 | pair 与 Landau 分项核验，不照搬 mass-i0 虚部 |
| 第35--36页 C37、D2--D10 | E、epsilon 坐标及允许域 | 实轴 cut 的独立坐标 oracle；q=0 单独处理 Jacobian |

论文研究的是 NJL 极化函数与散射，不能由此声称已完成 PNJL-BQS 的 GBU 密度、
热力学反馈或冻结线实验匹配。数值占据函数替换与热项延伸分别属于项目推广。
Sec.II C 第9页明确选择把截断施加于发散与收敛的全部项；热尾延伸不是
该文原样处方。原文重核及其对完整物理谱解释的影响见
[热谱适用性审查](ChargedGBU_ThermalAdmissibility.md)。

## 1. 同域 A/B 分解没有失效

令第一条线动量为 p，第二条为 p-q，x 为 p 与 q 的夹角。每个分量的共同域为

```math
\mathcal D_L(q)=\{0\le p<L,\;-1\le x\le1,\;
p^2+q^2-2pqx<L^2\}.
```

在同一域内定义

```math
\mathcal A_i(q;L)=-2\int_{\mathcal D_L(q)}p^2\,dp\,dx\,
\frac{1-f_i(E_i)-\bar f_i(E_i)}{E_i},
\qquad C(q;L)=\mathcal A_1(q;L)+\mathcal A_2(q;L).
```

第二项中的 E2 使用第二条线的动量；通过 p->q-p 交换可写为论文 Eq.(28)/(60)
的同域形式。分母消去一条线不等于解除那条线原来的截断条件。
P/S 通道分别令 M_P=m1-m2、M_S=m1+m2，则

```math
\Pi_X(z,q)=\frac{N_c}{8\pi^2}
\left[(z^2-q^2-M_X^2)B_0(z,q)-C(q)\right],
\quad z=k_0+\mu_1-\mu_2+i\gamma,\quad \gamma>0.
```

这与论文 Eq.(26) 的符号/归一化相同；RPA 通道耦合仍遵守本项目自己的
1-4K Pi 规范，不能从不同文献的 flavor 归一化直接移植一个系数。

平衡态 A_i 是单球积分，而此处的 mathcal A_i 是双球交集积分。只有 q=0、
相同占据函数、相同真空/热截断时才恢复 C(0)=A1+A2。有限 q 一般不等。
因此问题不是“A/B 分解不能用”，而是不能混用不同积分域的 A 和 B。
q>=2L 时整个该分量均为零，包括两个配套单线项；不能只令 B0 为零。
若真空用 Lambda、热项用更大的 L_th，q>=2Lambda 仅使真空分量消失，
不自动使完整泡消失。

同域代数还给出 B0(z,q)=C(q)/z^2+O(z^-3)，故上述 Pi 的高能常数项
精确消去。这是同谱无自由减法常数的依据；若只把 C(q) 换成平衡态常数，
Pi 在无穷远会残留 Nc[C(q)-A1-A2]/(8pi^2)。这不是允许任意调节的 subtraction，
也不是采用不同 regulator 后仍必须与本处方相等的要求。合成占据测试独立
检查 B0 的 1/z 系数为零、1/z^2 系数等于配套 contact。

## 2. 独立参考核

analysis-only `causal_gbu_pereira_reference.jl` 不复用现有 (p,x) 域划分、
四留数求和、cut 区间或 Cauchy 重建。仅复用 Gauss-Legendre 节点生成器。
复频率直接在柱坐标透镜内求积：

```math
t=p_z-q/2,\quad y=p_\perp^2,\qquad
|t|<L-q/2,\quad 0<y<L^2-(|t|+q/2)^2,\quad
p^2dp\,dx=\tfrac12\,dt\,dy.
```

这里是**将同一双球交集平移到中心**，不是既有 centered-sphere 敏感性探针
所用的单个相对动量球；前者是换坐标，后者是换模型处方。

以 n_i、bar n_i 表示粒子和反粒子占据，S=E1+E2、epsilon=E2-E1，参考核使用

```math
B_0=\int_{\mathcal D_L}\frac{p^2dp\,dx}{E_1E_2}
\left[
\frac{1-\bar n_1-n_2}{S+z}
+\frac{1-n_1-\bar n_2}{S-z}
-\frac{n_1-n_2}{\epsilon+z}
-\frac{\bar n_1-\bar n_2}{\epsilon-z}
\right].
```

前两项为 pair，后两项为 scattering。真空/热拆分是线性的：真空占据为零，
热修正去掉上述常数1，不能把热修正当作另一份完整泡叠加。
两个 mathcal A_i 独立直接积分，不从 B0 留数矩反推。真空 A 另用径向原函数
的闭式计算核验，测试同时检查解析透镜体积。

实轴使用另一组坐标 E=(E1+E2)/2、epsilon=E2-E1。pair 固定
E=-eta*lambda/2，scattering 固定 epsilon=-eta*lambda；允许域来自
Appendix D 的质量壳与 cutoff 条件。符号按 C18/C19 的 k0+i0 选择。
这里论文的 eta=+/-1 是求和指标，绝非有限展宽 gamma。
外部密度频率仍为 k0=lambda-(mu1-mu2)，Bose 核不再减一次化学势差。

## 3. PNJL 占据与热尾必须分账

| 变体 | 占据函数 | 真空截断 | 热截断 | 声明边界 |
|---|---|---|---|---|
| fermi_hard | Phi=PhiBar=1，普通 Fermi | Lambda | Lambda | Pereira 原始域/统计函数，质量和化学势仍冻结，不是重新求出的 NJL 平衡 |
| pnjl_hard | 当前 PNJL 多项式占据 | Lambda | Lambda | 只改变占据函数 |
| pnjl_tail8/12/16/24 | 同一 PNJL 占据 | Lambda | 8/12/16/24 fm^-1 | 单独评估热项延伸及尾部求积 |

所有变体使用相同质量、化学势、温度和耦合；不调用平衡求解器。
Phi=PhiBar=1 恢复 Fermi，Phi=PhiBar=0 恢复相应三夸克指数形式；
两条占据差的代数不依赖于这些函数是否为普通 Fermi。
这支持作为固定背景泡的明确推广，但**不证明**其与含介子反馈的 PNJL
完整热力学泛函、Ward 恒等式及稳定性全部闭合。
thermal-extended 的 Pauli 修正超出真空谱上界后可有符号，不能裁零。

## 4. 检查范围、运行和停止条件

`audit_causal_gbu_pereira.jl` 默认读取保留的 T170/muB240 BQS 输入，不读取
temp7 曲线充当数值真值。四通道、六 q、六变体分别记录复频 Pi/B0、同域
contact、真空单线解析式、双分辨率变化及正负实轴 pair/Landau cut。
数值为自然单位；Pi/C 为 fm^-2，B0 无量纲，频率/q/L 为 fm^-1。
每个量分别使用声明单位下的绝对门槛，不能把跨单位最大值解读成物理量。

输出目录由 `GBU_PEREIRA_OUTPUT` 指定且必须不存在；`GBU_PEREIRA_PDF`
可指定原始论文以绑定其哈希。384/192节点受控运行示例：

```powershell
$env:GBU_PEREIRA_OUTPUT='D:\w\jrt-ord\data\outputs\results\relaxtime\analysis\charged_rpa_phase_backend\pereira_new_diagnostic'
julia --project=. -e 'include("scripts/analysis/relaxtime/audit_causal_gbu_pereira.jl"); CausalGBUPereiraAudit.main(nodes=384)'
```

主 manifest 绑定源码快照、输入、配置和输出，并保留失败案例。进程运行中
不能修改被快照绑定的源码/配置。固定门槛下未收敛时保留结果并报错，不能
复用原目录覆盖失败。

另在 K+、q=1.4 的代表点用 Poisson 测试函数核验

```math
\int dx\,\frac{a}{\pi[(x-x_c)^2+a^2]}\operatorname{Im}B_0^R(x)
=\operatorname{Im}B_0(x_c+ia),\qquad a>0.
```

利用同一解析函数的 Poisson 半群，将有限 gamma 加到 a 上，检验 gamma->0
时这一加权积分趋近 PV+retarded 的结果。实轴积分由独立 (E,epsilon) cut
计算，右侧由柱坐标复泡计算，不要求阈值点值与有限 gamma 逐点相等。
这是 B0 的一个弱极限探针，**不是**非线性 GBU 相位密度的 eta 收敛认证，
也不能代替全谱根计数、UHP 稳定性、Mott 或冻结线验证。

本轮不调整内部域，不改 PNJLCore/旧生产默认/baseline，不计算新的介子密度。
确认公式一致只回答“是否正确实现该处方”；最终采用哪种正则化及其模型误差
仍需作者评审，不能从 parity 通过自动推导“物理结果应当唯一”。
