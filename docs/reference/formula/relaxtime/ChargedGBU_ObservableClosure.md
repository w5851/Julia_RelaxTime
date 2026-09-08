# Charged GBU 方向B：辅助响应与部分产额的闭合条件

## 范围与当前判断

本文接续[方向B源泛函](ChargedGBU_DirectionB.md)第6节。保留固定quark-only
BQS背景、PNJL占据、双球域及真空Lambda/热Lth；不更换物理处方或观察量。
这是B2/B3的公式与有限探针审查，不是完整稳定性、Mott或密度验收。
这里的有限热Lth按作者后续定义是数值上限，模型热积分到无穷；当前实现和
加权验收见[无限热积分](ChargedGBU_InfiniteThermal.md)。固定背景部分产额的
本轮准入条件不要求先完成含介子反馈的完整Phi驻点。

**B仍是可继续研究的主候选，但不能由B1直接宣布GBU产额已从完整驻点导出。**
负热尾不阻止标量GBU恒等式成立，也不单独否定一个辅助减法行列式；
它阻止的是未经说明就把该响应当作全频正物理谱。GBU本身不是正性修复。

## 1. 文献支持分三层

本次是针对性原文核验，非系统综述。INSPIRE公开API检索
`generalized Beth Uhlenbeck Phi derivable`、
`title "Generalized Beth-Uhlenbeck Entropy"`和`arxiv:1912.13162`。
专用MCP缺失时使用公开API及arXiv，不访问机构资源；PDF哈希绑定在本轮manifest。

| 原文 | 直接支持 | 不支持的扩大解释 |
|---|---|---|
| [Blaschke等，arXiv:1305.3907](https://arxiv.org/abs/1305.3907)，Sec.2.2、Eq.(86)--(90) | 固定平均场背景计算Gaussian关联；logdet导出相移导数形式 | 题名里的generalized已经等于本项目的正弦权重 |
| [Blaschke等2020，arXiv:1912.13162](https://arxiv.org/abs/1912.13162)，第7页Eq.(16)、其后Fig.3/4说明 | 明确以delta-sin(2delta)/2计算部分密度比；是目标GBU部分产额的文献先例 | 该文已认证本项目charged BQS双球热延伸或完整反馈驻点 |
| [Blaschke/Röpke/Baym，arXiv:2512.03876v1](https://arxiv.org/abs/2512.03876)，Secs.II A/B、Eq.(1)--(13)、(17)--(24) | 同一两loop Phi泛函的驻点、相关抵消与GBU权重结构 | 固定均场后任意加上GBU密度仍是完整守恒热力学导数 |
| [Blaschke/Dubinin/Turko，arXiv:1612.09556](https://arxiv.org/abs/1612.09556)，Eq.(1)、(4)--(8) | Phi框架中的正弦权重，同时调整夸克准粒子部分、包含碰撞展宽 | 此模型上游等同当前quark-only均场，或直接证明本截断方案 |

2512.03876的个别中间式涉及不同自能/传播子符号与疑似笔误，故本项目不
直接复制其导数等式作为证明；下节从实际D、Pi独立推导。PDF文本与页码已核对，
渲染页无法经当前image-input工具查看，未声称图形/排版验收。

## 2. 项目标量恒等式：需要的是Sigma_M=2Pi

给定实常数K>0，固定背景、频率求导不改变K：

```math
F=1-4K\Pi,\qquad D=\frac{2K}{F},\qquad
D^{-1}=\frac1{2K}-\Sigma_M,\qquad \Sigma_M=2\Pi.
```

Pi、Sigma_M单位fm^-2，K、D为fm^2；F、delta、W无量纲。
此处Sigma_M特指**介子传播子D的自能**；不是历史反馈讨论中的
“介子诱导的夸克自能”Sigma_q^(M)，两者不可混淆。本轮没有加入后者。
实轴非奇点取局部一致分支delta=-arg F，可直接验证：

```math
\operatorname{Im}D=(\operatorname{Im}\Sigma_M)|D|^2,\quad
(\operatorname{Im}\Sigma_M)\operatorname{Re}D=\tfrac12\sin(2\delta),
\quad W=\delta-\tfrac12\sin(2\delta),\quad
\frac{dW}{d\omega}=2\sin^2\delta\,\frac{d\delta}{d\omega}.
```

这些是代数恒等式，不要求ImPi非负。注意正弦项既不是Im(Sigma_M D)，
也不是首阶Born相位Im(4KPi)；两种替代均有纯代数反例。
它们同样不选全局分支、不确定解析gap根数。
每个经独立认证的普通正权束缚态，其delta跳变为+pi时，Delta W=pi；
GBU不会消掉稳定粒子计数。反之，连续区delta'<0时，W'仍非正。

小角度有W=(2/3)delta^3+O(delta^5)，因此GBU压低小角谱尾，但不改变其符号。
代码用小角展开避免浮点消减；不裁相位、不裁谱。

## 3. 减法、RPA与围道不能交换顺序

B1源泛函的二阶核线性相加：Pi_B=Pi_vac,Lambda+Pi_full,Lth-Pi_vac,Lth。
随后只构造一个F_B=1-4KPi_B。一般有

```math
\log(1-4K[\Pi_1+\Pi_2-\Pi_3])
\ne\log(1-4K\Pi_1)+\log(1-4K\Pi_2)-\log(1-4K\Pi_3).
```

因此不能用“三个各自重求和的介子气体相减”解释本核的产额；那会改变模型。
也不能由Gaussian Tr log F的相移项直接推出W替换；完整GBU来源还需要
自能项、Phi项及其驻点抵消。B1的fermion源Hessian不是该两loop认证。

对于**当前有限Lth且同域contact已匹配的核**，可写

```math
\Pi_B(z,q)=\frac1\pi\int_{-S_q}^{S_q}
\frac{\rho_B(x,q)}{x-z}\,dx,\qquad \rho_B=\operatorname{Im}\Pi_B^R.
```

rho_B可正可负。这一定义在开UHP解析；硬截断边缘可在实轴产生cut/端点
奇性，但不会凭空在Pi的UHP产生孤立极点。**F_B的零点仍可能形成D的UHP
极点**，因此log F_B的围道变形仍需检查零点、实轴根缩进及分支绕数。
不能由Pi解析直接宣称完整RPA稳定。

若支持位于[-S,S]且总变差int|rho|<=M已获统一上界，|z|>=R>S时

```math
|\Pi_B(z)|\le\frac{M}{\pi(R-S)},\qquad
4KM<\pi(R-S)\ \Longrightarrow\ F_B(z)\ne0.
```

这给出未来外圆闭合的充分条件，对有符号谱同样成立。当前仅以合成有符号
离散谱验证该界；真实核的采样/求积误差估计不是认证M，不得置certified=true。
Lth趋于无穷时须另证收敛/衰减，不能直接沿用紧支撑证明。

后续[端点闭合审查](ChargedGBU_EndpointClosure.md)给出原投影核的解析总变差界，
同时发现q0有限热硬端点会强制产生额外实零点；外圆排除不等于圆内计数完成。
此前谱插值的全gap检查不能直接作为原连续核的全谱认证。

### 可解析反例的用途

合成F(z)=(z²-E²)/(z²-P²)，E,P>0，无UHP零点/极点。其Gaussian Matsubara和
与带正确真空项的实轴结果一致，即
(E-P)/2+T log[(1-exp(-E/T))/(1-exp(-P/T))]。
E>P时辅助Pi的正频谱权可为负，带符号计数g(E)-g(P)<0。
这说明“负谱就必定无合法围道表示”不成立，但它也**不是正介子气体**。
若将E²改为负数，立刻出现UHP零点，普通实轴公式不再可无条件使用。
两者均为测试模型，不是项目发现的极点/密度。

## 4. 固定背景的部分产额不等于全部化学势导数

沿用已选观察量，而不将其升级为全热力学数密度。对每个charged通道，
在已确定物理分支的前提下，形式记账源nu定义

```math
P_{c,\mathrm{tag}}(\nu)=d_c\int_0^\infty\frac{q^2dq}{2\pi^2}
\int_0^\infty\frac{d\omega}{\pi}
g(\omega-\nu)W(\delta_c(\omega,q)),\qquad
n_{c,\mathrm{partial}}=\left.\partial_\nu P_{c,\mathrm{tag}}\right|_{\nu\to0^-,\delta\ \mathrm{fixed}}.
```

这里g(x)=1/(exp(x/T)-1)，omega=k0是external能量。
核使用lambda=omega+mu1-mu2；如果改用lambda积分，下限也须平移，
不能在g(k0)里重复扣除mu差。nu只是固定profile的形式偏导记号，最后取0，
**不是新拟合的介子化学势，不改变BQS/背景，也不是新增生产入口**。
对从omega=0开始的连续谱，取nu<=0的左导数，或先在a>0有限窗口求导再
检查a趋零的极限；不宣称nu>0时跨Bose极点的积分存在。合成测试使用a=0.03。
对实际BQS化学势求总导数时，profile、质量、K和边界都可能变化；
不能借助quark-only均场驻点省掉关联核的所有隐含导数。

在光滑有限窗口[a,b]（略去q测度），

```math
\int_a^b\frac{d\omega}{\pi}\,\beta g(1+g)W
=\int_a^b\frac{d\omega}{\pi}\,g\,\partial_\omega W
-\frac{[gW]_a^b}{\pi}.
```

有跳变时右侧导数须含离散Delta W delta_Dirac；连续谱积分加独立根贡献，
不是让unwrap数束缚态。若profile依赖某参数mu，还会出现
int g*(2sin²delta)*partial_mu(delta)/pi及移动边界项。
新增合成测试显式对照“固定profile偏导”和“同时改变profile的导数”，
验证二者不同及其差额，不宣称该合成参数就是实际BQS总导数。

GBU部分产额的固定背景近似有2020文献先例；完整守恒热力学则需要同一
Phi泛函和相应夸克自能/驻点。当前没有也不新增介子反馈，因此论文不得
声称后者已闭合。也不能以这一额外要求直接否定已约定的固定背景研究目标。

## 5. 红外、尾项与下一步门槛

若Bose端点处物理分支delta=c*omega+O(omega²)，则
beta*g(1+g)*W约为(2Tc³/3)*omega，GBU红外可积；
这只是充分条件。若端点保留非零pi平台，不能用同一结论掩盖发散。

若某尾区全域Re F>=a>0、|Im F|<=b，且**已确认近零无绕数分支**，
则|delta|<=atan(b/a)，|W|<=2atan(b/a)^3/3。
从omega0到无穷的每q壳层部分产额绝对界（d_c=1）为
q²/(2pi³)*[2atan(b/a)^3/3]*g(omega0)。
这不删除尾部；没有统一界或分支证据时，数值只作为估计，认证为false。

下一步仍优先B：用原连续cut核完成目标区域的RPA零点/外圆与近轴控制；
再把独立gap根、Mott/Levinson、Bose端点和负尾带误差界地接入部分产额。
完整研究生产扫描在上述验收之后，通过统一入口执行，不手动加密旧诊断。
若只能通过更改观察量、补无依据UV参数或加入反馈来闭合，提交作者评审，
不自动切换A或修改PNJLCore。

## 6. 实现与证据边界

analysis-only模块causal_gbu_observable_closure.jl提供标量恒等式、非交换反例、
固定profile导数、条件尾界和可解析辅助模型；全为纯代数/合成输入单元测试。
audit_causal_gbu_observable_closure.jl使用GBU_OBSERVABLE_OUTPUT指定不存在的目录，
读取T170/muB240保留BQS，使用独立cut+直接PV积分，不消费折叠相位。
实际结果仅验证局部主值相位的恒等式，不认证该主值为密度分支。
结果、失败与哈希记录见strict-audit第22.21节；src、旧默认及baseline不变。
