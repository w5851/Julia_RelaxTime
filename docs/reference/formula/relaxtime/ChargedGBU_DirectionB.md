# Charged GBU 方向B：保留热延伸的静态基础与响应匹配

## 授权、问题与阶段结论

作者2026-09-06要求先尝试方向B，再视结果考虑同硬截断方向A。
作者随后明确：热积分的模型目标是无穷，本文历史探针中的有限Lth只是数值
上限，不是新模型参数。当前无限目标和验收入口见
[无限热积分](ChargedGBU_InfiniteThermal.md)；下文保留各轮历史范围与证据。
本阶段不修改 PNJLCore、生产默认或方法TOML，不重求背景、不计算介子密度；
属于当前 charged track 的 required_follow_up。

**B未被排除，但尚未形成可采用的新物理响应。** 第6节给出明确的固定投影
减法泛函，其二阶核对应现有泡，不是消除负热尾的修正。
当前热延伸势的质量/化学势导数
可与单线凝聚和夸克净密度一致；这不推出它的任意动态延拓都具有物理谱。
上一轮的负热尾反例针对现有完整泡的全频物理谱解释，不证明上游热延伸平衡
本身不合法。B失败的操作性含义是无法满足本项目目标观察量的响应匹配与验收，
不是必须先证明所有热延伸模型不成立才能研究A。

后续B2/B3的[辅助响应与观察量审查](ChargedGBU_ObservableClosure.md)明确区分
Gaussian行列式、标量GBU恒等式与完整Phi驻点来源；B仍为主候选，尚未完成密度验收。
证据与执行见 [strict-audit 第22.19节起](../../../dev/active/2026-08-30_charged-rpa-bu-strict-audit.md)；
前序负谱判定及其必要条件边界见 [热谱适用性](ChargedGBU_ThermalAdmissibility.md)。

## 1. 针对性文献复核，不构成系统综述

INSPIRE API查询 `NJL and regularization and thermal`、`Nambu Jona Lasinio
and dispersion relations`、`NJL and medium separation scheme`，再按arXiv号
核对3条记录；公开PDF已下载并绑定哈希。内置web网关404且无INSPIRE MCP工具，
使用公开API和arXiv HTTP，不访问机构授权资源。

| 原文 | 直接支持的事实 | 不能据此宣称 |
|---|---|---|
| [Maslov/Blaschke, arXiv:2301.09882](https://arxiv.org/abs/2301.09882)，第3页Eq.(4)--(7) | PNJL真空三动量截断、热积分全动量，并由同一势求平均场 | 其两味mu=0动态处方直接适用charged BQS |
| [Lopes等, PRD103,076023](https://doi.org/10.1103/PhysRevD.103.076023)，arXiv:2102.02844，Sec.II A，第3页 | SU(3)有限同位旋NJL比较TRS/MSS；TRS热积分到无穷，MSS将介质依赖从发散积分分离 | MSS就是当前真空/热拆分，或已给出本项目GBU谱 |
| [Xue/Yu/Wang, Chinese Physics C (2022)](https://doi.org/10.1088/1674-1137/ac2ed3)，arXiv:2105.14323v2，Eq.(7)、结果/结论 | 比较热项截断选择；报告其两味mu=0模型在高温下未截热项可能出现M低于裸质量，而全热硬截断影响高温热力学量 | 一篇研究已普遍证明B不能用，或其现象已在本项目冻结线上出现 |
| [Broniowski等, arXiv:hep-ph/9509363](https://arxiv.org/abs/hep-ph/9509363)，Secs.9--10 | 讨论有限regulator带来的额外奇点、围道项；Pauli-Villars的高能负cut与proper-time复极点是不同问题 | 任意正则化后的辅助响应都应被当作全频正物理谱，或所有负谱方案都不能计算任何关联贡献 |

最后一项说明需要区分“受正热迹约束的完整物理关联函数”和“含regulator
自由度的辅助响应”。要将后者用于GBU，须检查从同一配分函数导出实轴公式时
有没有额外极点/围道/减法项；不能只保留一条相移曲线。
Maslov Appendix A4给定的实部谱上界与正文不截虚部必须单独核对，不能直接
拼接为本项目的因果候选。当前保留原PDF，文本位置已复核；
第3页渲染后图像工具仍不支持image input，未声称视觉验收。

## 2. 项目上游究竟定义了什么

实际入口是 `src/models/Models.jl` 引入的 `pnjl_physics/PNJLModel.jl`、
`pnjl_physics/PNJLCore.jl` 和 `omega.jl`，不是旧的同名兼容core。
在xi=0、给定Phi/PhiBar下，每味夸克动能势为

```math
\Omega_{q,i}=
-\frac{N_c}{\pi^2}\int_0^\Lambda p^2 E_i\,dp
-\frac{T}{\pi^2}\int_0^{L_{\rm th}}p^2
[\log Z_\Phi((E_i-\mu_i)/T)+\log Z_{\bar\Phi}((E_i+\mu_i)/T)]\,dp ,
```

其中Nc=3，`Z_Phi(x)=1+3Phi*exp(-x)+3PhiBar*exp(-2x)+exp(-3x)`。
再加手征势与Polyakov势组成上游Omega；Lth数值默认10 fm^-1，
理论目标为Lth趋无穷，此处有限上限只记录实现。所有质量依赖、占据和
数值上限必须一致求导，并对目标观察量控制相应尾误差。

由 `d log Z_Phi/dx=-3 f_Phi` 得到

```math
\frac{\partial\Omega_{q,i}}{\partial M_i}
=\frac{N_c M_i A_i}{4\pi^2},\qquad
-\frac{\partial\Omega_{q,i}}{\partial\mu_i}
=\frac{N_c}{\pi^2}\int_0^{L_{\rm th}}p^2(n_i-\bar n_i)\,dp .
```

这不要求真空和热使用同一上限。由同一势定义静态导数在数学上成立，
并不能据负热尾否定这些恒等式。也不能由恒等式通过就声称背景是全局最小值。
本轮测试质量而非凝聚坐标的动能势偏导，未检查全部Polyakov/约束Hessian。
这些静态导数也不是 charged pseudoscalar 有限k0/q响应的直接认证。

冻结T170/muB240 BQS下，两种占据、3种热上限、2种节点、3味，共36条：
模型有效入口的AD导数与独立积分通过；所有变体的mass/mu均不重求。
改变Phi到1仅是Fermi控制，不是得到新的NJL自洽平衡。

## 3. 数学上为什么静态匹配还不够

一个恒定场Omega并不能唯一确定所有非均匀源的生成泛函。
待推导的B候选可从以下**符号性结构**出发：

```math
\Gamma_B[\mathcal M]=\Gamma_{\rm local}[\mathcal M]
+\Gamma_{\rm vac}^{\rm reg}[\mathcal M]
+\Delta\Gamma_{\rm med}^{\infty}[\mathcal M] .
```

这里Mcal包括标量与pseudoscalar外源。介质差在同一外源处减真空，
不是把质量固定到真空数值后再作差。必须明确reg算符及其对外源的依赖，
才能对同一泛函求二阶变分；不能只指定两个积分上限。
若有场依赖regulator，变分其边界/顶角产生的项也不能遗漏。

有限的梯度/频率匹配项可能不改变恒定场Omega，却改变Pi(k0,q)。
因此“保持上游不动”没有自动固定所有下游系数。任何新匹配条件都要列来源、
单位、适用范围和自由参数，不能用ratio或负密度来调参。

### 合成反例：正谱和两个匹配系数仍不保证稳定

仅用合成谱 `rho(t)=c*t*theta(t-S)`，c>0、S>0；t和s均为频率平方。
两次减法补项为

```math
\Delta\Pi(s)=\frac{s^2}{\pi}\int_S^\infty\frac{c\,t}{t^2(t-s)}\,dt
=-\frac{c\,s}{\pi}\log(1-s/S).
```

它满足DeltaPi(0)=DeltaPi'(0)=0，s>S时ImDeltaPi=c*s>0。
然而沿s=-Q²，

```math
\Delta\Pi(-Q^2)=\frac{cQ^2}{\pi}\log(1+Q^2/S)\to+\infty .
```

于是最简单的contact-RPA `F=1-4K DeltaPi`，K>0，必在某Q>0穿零，
对应z=iQ的UHP极点。这是合成反例：**不是在当前项目泡中发现了该极点，
也不排除有适当高能定义的其他B候选。**

它只排除“恢复正尾+保持常数/斜率，就已经完成所有自洽性”的论证。
不能任意将本补项加到真实charged核；charged静态k0=0对应
lambda=mu1-mu2，而非lambda=0，同样需要正确的匹配点。

## 4. 具体推进顺序和停止条件

| 阶段 | 下一项交付 | 通过条件 |
|---|---|---|
| B0 静态基础 | 当前Omega、A、密度偏导映射及本阶段36项证据 | 同一热处方下恒等式可复核；不把它当动态/完整驻点认证 |
| B1 生成泛函与动态响应 | 明确一个具名regulator/介质差候选，对charged源求二阶变分 | q0静态关系、flavor反射、同域contact、retarded实虚部及所需减法相互一致；未定系数全部标出 |
| B2 解析结构与有效域 | 区分物理cut和regulator cut/极点，建立近轴及高能验收 | 不用只改实部/只改虚部修谱；若有额外围道项，从公式中完整保留 |
| B3 GBU观察量闭合 | 从该响应的关联泛函核对GBU权重、计数和边界 | 说明是否仍是作者要求的部分粒子产额；不能自行换成另一种有符号相关量 |
| B4 受控数值 | 小矩阵根continuation、Mott/Levinson、PV/eta积分及外层误差 | 原门槛验收并分账regulator依赖，再接统一研究生产脚本 |

B1优先保留现有上游，不自动改热上限、PNJLCore、参数或相移规范。
不要求一次试探就提供全频UV完备理论；但若只能在有限有效域成立，就必须
明确域外处理和密度误差/匹配依据，不能悄悄将外频率截断当作修复。
从Gaussian行列式得到普通BU，并不自动推导GBU的sin(2delta)修正；
该权重与无介子反馈的部分产额约定需要在B3单独闭合。

触发作者复审/A对照的条件包括：在目标域内找不到配套的动态响应、
需要改变既定观察量或上游模型、关键匹配参数无依据，或B的实施成本明显超出
当前论文范围。**不要求先普遍证明“上游不能热延伸”。**
此时可以选择A的完整一致模型、保留B为局限明确的诊断，或调整论文范围，
由作者决定；不按实验更接近自动择模。

## 5. 入口与测试

`causal_gbu_direction_b.jl` 提供analysis-only静态势和合成减法反例，
`test_causal_gbu_direction_b.jl` 使用纯代数/合成输入。
`audit_causal_gbu_direction_b.jl` 只读冻结背景，调用Models的动能势求AD导数，
输出 `static_identities.csv` 与 `synthetic_subtraction.json`。
`GBU_DIRECTION_B_OUTPUT` 必须指向不存在的diagnostic目录，manifest绑定
源码、输入、原文、参数及输出。`full_stationarity_certified`、
`direction_b_response_selected`、生产授权均false。
这不是新增稳定public provider，旧密度入口和baseline完全不变。

## 6. B1：固定动量投影的同源减法泛函

### 6.1 文献事实与项目定义

Pereira 2024 Sec.II B、Eq.(14)--(17)从fermion determinant展开二次作用量；
Sec.II C将每条内部线限制在球内，并对发散与收敛项均用同一硬截断。
这支持逐线投影结构，但没有提供以下热延伸处方的GBU验收。
下面是本项目构造的明确泛函，不冒充Pereira原样模型。

在有限盒中，定义仅作用于空间动量的正交投影：

```math
P_L=\theta(L-|\hat{\boldsymbol p}|),\quad
D_L[J]=\left.P_LD[J]P_L\right|_{\operatorname{ran}P_L},\quad
\Gamma_L[J]=-\operatorname{Tr}_{\operatorname{ran}P_L}\log D_L[J].
```

先限制fermion Hilbert空间再取行列式，不等于任意插入Tr(P_L log D)，
也不是只限制一条线。D[J]对标量/pseudoscalar味非对角源J线性；
P_L不依赖J、质量或Polyakov场。只对J变分，不沿背景解求全导数。
势的绝对归一化采用上游的normal-order约定：去掉填满Dirac海产生的、
与J无关的mu线性项。8x8块对应加回2Nc(mu1+mu2)，绝不能把这一项遗漏后
再声称净夸克密度一致。它不贡献本轮J曲率；第2节上游密度比较使用的是
已经正确归一化的势，不使用这里未归一化的中间特征值和。

具名B候选（只用于公式审查）为

```math
\Gamma_{q,B}[J]=\Gamma_\Lambda^0[J]
 +\left(\Gamma_{L_{\rm th}}^{T,\mu,\Phi}[J]-\Gamma_{L_{\rm th}}^0[J]\right).
```

上标0是同源、同质量、同协变外频率的无热占据参考，不是重求真空质量，
也不能把源先置零再减常数。charged参考必须保留z=lambda=k0+mu1-mu2；
不能在k0固定时擅自把mu差置零。恒定对角场时参考退化为常规真空能，
该组合返回第2节的真空Lambda/热Lth动能势。无限热尾仍需独立收敛。
对周期时变源，参考二次核是在相同外部Matsubara频率上取真空解析核，
再统一做retarded延拓，不是在热圈上把T简单置零而保持所有求和索引不变。

本轮具体实现限制|mu_i|<M_i、T>0、xi=0。此时零温参考负能带为s=-1，
与费米海onset分开。静态源矩阵用H_i=alpha.p+beta.M_i-mu_i的负能
特征值和定义参考。若出现零温费米面、时空依赖Polyakov场、局域规范外源
或源依赖regulator，须重新推导，不能沿用此处证明。

2026-09-07补充：[原连续核与源参考延续](ChargedGBU_EndpointClosure.md)第4节
在Dirac参考能隙保持开启时，以固定rank-four海延续此局部二阶证明。
旧oracle及其限制保留；不将跨单味onset本身判为上游失效。

### 6.2 两线域、Dirac迹与contact

在均匀J=0背景D0与P_L对易，G_L=(D_L[0])^-1。由于delta P_L=0，

```math
\delta^2\Gamma_L=\operatorname{Tr}
  (G_LP_L\delta D P_LG_LP_L\delta D P_L).
```

J连接(p,flavor1)与(p-q,flavor2)，每条线得到一个投影，形成双球交集。
源变分不移动动量边界，所以**此具名定义下没有遗漏delta P边界顶角**。
q导数、场依赖截断和局域规范Ward恒等式属于另一问题，本式不认证它们。

令r=p-q，E_i=sqrt(p_i²+M_i²)，L_i^s=(1+s h_i/E_i)/2，
h_i=alpha.p_i+beta.M_i，V_P=i beta gamma5、V_S=beta，则

```math
\mathcal T_{st}^{P}=\operatorname{tr}_D(L_1^sV_PL_2^tV_P)
 =1-st\frac{\boldsymbol p\cdot\boldsymbol r+M_1M_2}{E_1E_2},\quad
\mathcal T_{st}^{S}=1-st\frac{\boldsymbol p\cdot\boldsymbol r-M_1M_2}{E_1E_2}.
```

定义u_st=sE1-tE2，d_P²=(M1-M2)²、d_S²=(M1+M2)²，得到

```math
\mathcal T_{st}^{P/S}=-\frac{st}{2E_1E_2}(u_{st}^2-q^2-d_{P/S}^2),\quad
\Pi_L^R=N_c\int_{\mathcal D_L}\frac{d^3p}{(2\pi)^3}
  \sum_{s,t}\mathcal T_{st}\frac{N_1^s-N_2^t}{z-u_{st}}.
```

真空N^-=1、N^+=0；热差N^-=-nbar、N^+=n；full N^-=1-nbar、N^+=n。
按tr_flavor(T²)=1实源规范，fermion Hessian为Gamma_q''=-Pi，不是+Pi。
项目D=2K/(1-4KPi)的charged-ladder规范不变。
独立Dirac迹与现核逐留数相同。把u²分子改写为外部z²分子时，同域一阶矩给出

```math
C_L(q)=-2\int_{\mathcal D_L}p^2dp\,dx
\left[\frac{c-n_1-\bar n_1}{E_1}+\frac{c-n_2-\bar n_2}{E_2}\right],\quad
\Pi_L=\frac{N_c}{8\pi^2}[(z^2-q^2-d^2)B_{0,L}-C_L(q)].
```

真空/热/full取(c,n,nbar)=(1,0,0)/(0,n,nbar)/(1,n,nbar)。
C(q)是源泡的代数分解所需项，不是新增可调counterterm；q0才退化为A1+A2。
仅改实contact无法改掉cut符号。

### 6.3 PNJL与独立静态二阶检查

固定Phi/PhiBar可形式化为复化color holonomy，其特征值ell_c满足
sum ell=3Phi、sum ell^-1=3PhiBar、prod ell=1，因而

```math
\prod_c(1+\ell_c e^{-x})=1+3\Phi e^{-x}+3\bar\Phi e^{-2x}+e^{-3x},\quad
\frac13\sum_c\frac{\ell_c e^{-x}}{1+\ell_c e^{-x}}=f_\Phi(x).
```

这解释固定Polyakov背景下占据替换的形式行列式来源；Phi/PhiBar为不同
独立实平均场时，一般不是unitary SU(3)矩阵，**不是正密度矩阵或Lehmann
正性证明**。未积分Polyakov场响应，也未推导PNJL+KMT全部驻点Hessian。

oracle独立构造8x8 Hermitian静态源矩阵
diag(h1-mu1,h2-mu2)+J*[0 V;V 0]/sqrt(2)，对特征值e求
-T*sum(log Z_Phi(e/T))；参考为Nc*sum(e<0 ? e : 0)。
源无关的normal-order常数抵消。此特征值二阶差分不消费B0或相移，应与
-Pi(k0=0,q)相配；同时检验源步长和内部求积节点。
同质量同mu的q0 pseudoscalar极限另检查Pi(0,0)=-(dOmega_i/dM_i)/M_i，
包含简并态的导数极限。这不是介子数密度的相移积分分部检查。

### 6.4 结论与后续门槛

该定义提供**现有核的固定背景fermion泛函来源**，没有构造新正谱核。
双球域、C(q)、PNJL占据和热延伸可同时出现在这一减法定义中，
不能仅因它们并存就称为代码漏项。但

```math
\rho_B=\rho_\Lambda^0+\rho_{L_{\rm th}}^{T,\mu,\Phi}-\rho_{L_{\rm th}}^0
```

是有符号组合；Lambda<p<Lth的q0 pair只有热抑制差项，负热尾仍在，
即使Fermi控制也一样。写出生成泛函不把determinant比值变成正热迹，
更不自动推导GBU的delta-sin(2delta)/2部分产额。

下一阶段B2/B3核对减法辅助核的目标域解析结构、额外围道项及GBU关联泛函。
若只能得到有符号关联压力，不能静默重命名为作者要求的介子部分粒子产额。
若需要新UV完成、源依赖投影或匹配参数，先提交作者评审，不按谱/ratio调参。

入口为causal_gbu_source_functional.jl及audit_causal_gbu_source_functional.jl；
GBU_SOURCE_OUTPUT要求新diagnostic目录。测试包括投影行列式、单线反例、
PNJL多项式、P/S Dirac迹及静态源曲率。实际背景只读既有BQS输入，动态
探针使用external k0，不求平衡或密度。结果见strict-audit第22.20节；
物理响应选择、全驻点、全UHP计数、GBU观察量推导和production授权仍false。
