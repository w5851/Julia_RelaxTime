# Charged GBU：原连续核的热端点、计数与源参考延续

## 范围和结论边界

本文接续[观察量闭合条件](ChargedGBU_ObservableClosure.md)，属于方向B研究，
不是新production路线授权。以下端点与能带结论是项目独立数学推导；
不声称是Pereira或2020论文已给出的PNJL热延伸定理。
固定quark-only BQS、PNJL占据、逐线双球域、同域contact和既定GBU部分产额不变。

作者后续明确有限Lth是数值处方，模型热积分到无穷。本文有限端点定理仍成立，
但其整数根数不是无限目标的生产条件；应控制加权产额误差。无限域实现与当前
验收见[无限热积分](ChargedGBU_InfiniteThermal.md)。本文末尾清单保留前轮状态。

**有限热硬端点的原连续核，与端点线性平滑后的谱插值核，不一定有相同的根数。**
这是对此前“全gap数值扫描”的实质性限制，不是发现了新的正常介子，
也不直接证明B路线不稳定或总密度为负。源参考跨夸克起始面的覆盖问题则另行处理。

## 1. 一个不依赖谱插值的外圆界

内部频率为z=lambda=k0+mu1-mu2。有限Lth时支撑包含于[-S,S]，
S=sqrt(Lth²+M1²)+sqrt(Lth²+M2²)。Dirac投影迹满足T_st>=0、sum T_st=4，
full和vacuum占据差的绝对值不大于1。由三角不等式：

```math
\frac1\pi\int|\rho_B(x,q)|dx\le M_q
=\frac{N_c}{2\pi^3}\left[V_\Lambda(q)+2V_{L_{\rm th}}(q)\right],
\qquad V_L(q)=\begin{cases}
\pi(4L+q)(2L-q)^2/12,&0\le q<2L,\\0,&q\ge2L.
\end{cases}
```

Mq单位fm^-3。取R=S+max(1 fm^-1,8K Mq)，有|4K Pi|<=1/2，故此圆外F不为零。
external k0平面相应圆心平移为-(mu1-mu2)。这是实数算术下的解析界，
并非浮点有向舍入认证；它只排除外圆之外，**不排除圆内或近实轴的零点**。

## 2. q=0有限热端点必有额外实零点

对有限T>0、Lth>Lambda，在正频热支撑上端S的左侧，真空已结束，pair权为
-(n1+nbar2)，故rho(S-)=c<0；其值不是恰在S处按代码定义赋的零。
令a=S-h，h>0，显式分出这一端点常数：

```math
\Pi(z)=R_h(z)+\frac c\pi\log\frac{S-z}{a-z},\qquad
R_h(S)\ \text{有限}.
```

z从实轴S右侧趋近时Pi趋+infinity，F=1-4KPi趋-infinity；远处F趋1。
因此热cut上方**至少一个正external频率的实零点**必存在（假设S>mu1-mu2）。
不声称唯一；一般简单负到正穿越的F'>0，对项目谱规范留数-2K/F'<0。
它是截断辅助零点，不等同正常正谱束缚态，也不是UHP不稳定极点。

其距端点的渐近估计为

```math
\log\frac{r-S}{h}\simeq\frac{\pi[1-4K R_h(S)]}{4Kc}.
```

在保留的T170 MeV、muB240 MeV BQS背景，四通道、Lth=10/20/24 fm^-1均
验证了负端点、有限余项与外部正F。Lth=10时log10[(r-S)/(fm^-1)]约
-6174至-10185；20/24时绝对值更大。这是渐近距离估计，**不是直接求根的残差**。
普通Float64扫描再加密不能分辨这样的根。线性谱插值将端点接到零，
消除了对数发散，因此此前插值计数不能直接传递给原核。

q=0本身没有q²相空间权，不能据此归因旧的有限q负密度。
每个外根有Bose权上界g(S-mu1+mu2)，但尚未证明根数上界或有限q范围，
所以此单根上界不是总密度误差界。

## 3. 正常束缚态与辅助零点必须分账

对retarded边值，一个简单实零点局部F(z)=F'(r)(z-r)+...；无论F'(r)正负，
按解析延续得到的delta=-arg F在穿过该零点时均有+pi跳变。
GBU权W=delta-sin(2delta)/2也跳+pi。因此“只数正谱正常根”与“数log F的
所有零点”不是同一件事；后者的补偿也不能让unwrap悄悄代劳。

可解析合成对照：rho(x)=-c sign(x)，a<|x|<S，c>0。

```math
F(z)=1+\frac{4Kc}{\pi}\log\frac{S^2-z^2}{a^2-z^2},\quad
t=e^{-\pi/(4Kc)},\quad r^2=\frac{S^2-ta^2}{1-t}>S^2.
```

各对数取从UHP延续的支，不能任意合并后丢掉绕数。正频连续区W从0降至-pi，
外根带来+pi，净计数为0。只数正留数根会漏掉一单位；这个合成例沿纯虚轴F>1，
额外实根不能被误称为UHP不稳定的证据。它也不是项目的密度结果。

有限q>0时，固定上端附近的壳上能量区间随S-lambda收缩，通常
rho(S-epsilon)=O(epsilon/q)趋零。q0对数端点定理不能直接套到所有有限q。
q趋零、Lth趋无穷、eta趋零的次序必须记录；热尾指数小不等于整数计数连续。
正常低能Mott的根数1→0与阈值相位pi→0，只验证正常支局部补偿，
仍不能代替包括这些辅助结构在内的全谱Levinson检查。

## 4. 跨夸克起始面的固定Dirac海参考

此前[源泛函](ChargedGBU_DirectionB.md)第6节oracle用|mu_i|<M_i，
故“grand-energy<0”与负Dirac能带等价。保留BQS在Mott附近部分夸克已
mu_i>M_i，该oracle不再覆盖；这不意味着PNJL有限密度背景非法。

原定义的参考是**无占据Dirac海**，不是另一个填充至mu的零温介质。
延续时应保留s=-1的四个Dirac态，vacuum占据仍为N_s=(1-s)/2；
full仍用PNJL，thermal=full-vacuum。对两味静态8x8源块，将vacuum参考
定义为最低四条能带的特征值和，而不是所有负grand-energy特征值之和。
其全动量充分能隙为

```math
g_{\rm ref}=\min\{2M_1,2M_2,\sqrt{q^2+(M_1+M_2)^2}-|\mu_1-\mu_2|\}>0.
```

源顶角除sqrt(2)，由特征值扰动界，sqrt(2)|J|<gref保证[-|J|,|J|]内
参考能带不交叉。重求导给出同一Dirac投影核，故**不用改变现有Pi或上游**。
源无关normal-order常数仍加回2Nc(mu1+mu2)。同能带简并用PNJL导数极限；
异能带隙关闭则拒绝套用该证明，不能仅删除旧输入检查。

这是均匀背景局部二阶源Hessian的延续，静态8x8有限J只是曲率oracle，
不是任意非均匀有限源的完整action，也不认证PNJL/KMT全部驻点或GBU产额。
合成P/S、普通及越过onset情况均与直接源曲率、现核对照。
真实证据只使用保留背景，见strict-audit第22.22节。

## 5. 数值合同与剩余验收

causal_gbu_endpoint_closure.jl给出端点定理的诊断量、外圆界及可解析toy。
causal_gbu_continuation_gate.jl对原cut做cosine映射求积，保留各panel及跳变。
遇到Float64不能表示的panel内部节点，使用宽精度坐标重算同一个Pereira cut，
不夹节点到端点，不删小panel。Gauss权仍是Float64，故这不是任意精度认证。
门槛是两阶数变化的数值估计，不是原核统一误差上界。

causal_gbu_source_continuation.jl是独立源参考延续oracle；旧oracle契约不变。
三个audit入口使用GBU_ENDPOINT_OUTPUT、GBU_CONTINUATION_OUTPUT、
GBU_SOURCE_CONTINUATION_OUTPUT指定不存在的新诊断目录。源码快照、输入与CSV
均由manifest绑定；首次失败保留，不覆盖重跑。

剩余工作应按依赖验收，而不能以“已运行”代替“已通过”：

1. 原连续核圆内/近轴结构及所有解析gap根，显式处理硬热端点辅助零点。
2. 正常Mott局部补偿与全谱Levinson分账；单独检查辅助根/连续项的抵消。
3. 带上述计数的GBU加权积分，再做PV/eta分布极限、IR端点、q与omega尾界及节点收敛。
4. 前三项通过后接统一研究生产入口，保留失败点；不手动加密旧冻结线诊断。

可继续评估的同处方数值路线包括端点对数解析分离，或对原本意图的无限热积分
显式计算/控制尾部，再检验加权极限；不得静默删掉有限L的辅助贡献。
如果需要换观察量或引入新的UV定义/参数，则回到作者评审。
