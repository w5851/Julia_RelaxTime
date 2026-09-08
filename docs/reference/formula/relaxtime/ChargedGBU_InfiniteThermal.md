# Charged GBU：无限热积分目标与数值验收

## 定义与授权边界

作者明确：理论热积分遍及全内部动量空间；有限 Lth 是数值积分上限，
不是新模型参数。真空 Lambda 则仍是模型参数。本文接续
[方向B源定义](ChargedGBU_DirectionB.md)与[端点分析](ChargedGBU_EndpointClosure.md)，
不改变 PNJLCore、quark-only BQS 约束或生产默认。

计算量固定为 quark-only 背景上的 GBU 介子部分产额，不要求将介子反馈加入
上游驻点作为此次准入条件，也不将它改称实验最终强子产额。
此前有限 Lth 端点定理继续有效，但不能据辅助实根的存在直接否定无限目标，
也不要求不同数值上限的整数根数一致。

## 核与数值表示

目标是同源减法的二阶响应：

```math
\Pi_B^R(z,q)=\Pi_{\mathrm{vac},\Lambda}^R(z,q)
 +\Delta\Pi_{\mathrm{med},\infty}^R(z,q),\qquad
z=\omega+\mu_1-\mu_2+i0^+ .
```

真空项同时限制两条内部线的动量；热差在全空间积分，保持 PNJL 占据、
contact 与 retarded 处方一致。新增的径向 Dirac 投影积分与原 cut 的
Cauchy 积分是两种独立数值表示，不是两种物理处方。

对半无限能量区间使用坐标映射

```math
x=a+T_*\frac{t}{1-t},\quad dx=\frac{T_*dt}{(1-t)^2},\quad0<t<1.
```

有限但很长的 Landau 区间也使用同一映射，终点取
`tmax=(b-a)/(Tstar+b-a)`，避免线性 Gauss 网格漏掉热衰减边界层。
无穷远不是另一个有限硬端点；有限求积阶数仍需验收。

快速表示从谱中解析分离 unitary 阈值的平方根项，剩余谱分段插值，
真空硬边左右极限分别保留。窗口函数和其 Cauchy 积分为

```math
b(u)=A\sqrt u(1-u/H),\quad0<u<H,
\qquad J(s)=\int_0^H\frac{\sqrt u(1-u/H)}{u-s}\,du .
```

窗口在 H 连续归零；该解析分离不删除任何谱。
数值分区以外的热尾仍通过映射积分到无穷，不将其强制设为零。
真空跳点的 PV 点值不定义，必须使用单侧极限或复轴值。

## 根与产额的分账

无限目标不再具有有限 Lth 人造的最外解析 gap。正常正、负解析 gap
位于 `hypot(q,M1-M2)<abs(lambda)<hypot(q,M1+M2)`；不能把热尾指数小的
区域称作严格解析 gap。独立实根搜索不是全谱完备性的证明，仍须用
复围道及 cut 内非零检查补足，不能让密度 unwrap 承担根数认证。

对已经通过正常拓扑检查的通道和 q，设正频正常根为 omega_j，unitary
阈值为 omega_U，W=delta-sin(2delta)/2，g=1/(exp(omega/T)-1)。分账写成

```math
n_M(q)=\frac{q^2}{2\pi^2}\left[
 \sum_jg(\omega_j)
 +\int_{\mathrm{cuts}}\frac{-g'(\omega)}\pi W(\omega)\,d\omega
 -n_Bg(\omega_U)\right].
```

此式假定 Landau 两端的 W 边界为零、unitary 下端 W=nB*pi、高端趋零。
若实际拓扑不满足这些条件，当前受限候选拒绝计算，而不是把异常相位改零。
`-nB*g(omega_U)` 是连续谱的边界项，不是删除束缚态；它与离散根分账后
恢复相移导数表达式。连续谱贡献可以为负，不能按其符号裁零。

## 工程入口与证据规则

- `causal_gbu_infinite_thermal.jl`：原 cut 与独立无限径向积分。
- `causal_gbu_infinite_profile.jl`：解析阈值、保留跳变的快速表示。
- `causal_gbu_infinite_yield.jl`：独立正常 gap 根和受限拓扑的 GBU 分账。
- `audit_causal_gbu_infinite_thermal.jl`：保留背景、源码快照和所有失败。
- `causal_gbu_infinite_acceptance.jl`：实轴Nyquist、根邻域围道及有限eta带权极限。
- `causal_gbu_infinite_qgate.jl`：从原核阈值零点定位Mott动量，再划分外层积分。
- `audit_causal_gbu_infinite_production_gate.jl`：自动外层分段、热上限差与尾区检查。
- `causal_gbu_infinite_limits.jl`：与根数独立的阈值极限及近端点围道窗口。
- `audit_causal_gbu_infinite_gate_recovery.jl`：带旧证据哈希约束的失败点自动恢复。
- `audit_causal_gbu_infinite_readiness.jl`：五阶段证据完整性/数值验收归并，不授权默认晋升。

上述模块保留 analysis 单一源码；第22.24节新增的显式研究生产适配器
`Models.run_charged_gbu_freezeout_scan` 调用同一核，并逐背景执行自动门禁。
这不等于替换旧 production provider 或晋升全部 analysis API。
原核对照门槛仍为1e-6 fm^-2，产额相对误差目标仍为1%。局部探针通过不等于
全谱证明。真实结果、首次失败及结论在 strict-audit 第22.23–22.24节记录，
不写入 baseline。

有限eta比较同时保留积分上下端的 `g*W/pi` 项。先在固定正的IR端点令
eta趋零，再缩小IR端点，不能要求有限eta在Bose原点或极点逐点等于PV。
外部omega数值上限随q提高，以满足同一相位门槛；这不是更改内部真空Lambda。
正常Mott动量处分割q区间，可以使离散束缚态及连续谱边界项分别收敛，
而不只依赖两者误差碰巧抵消后的总数。

### 阈值极限与根围道的数值合同

近Mott处不能用固定的 `omega_U+h` 代表严格单侧极限。若
`rho(U+u)=A*sqrt(u)+o(sqrt(u))`、K>0、F(U)为非零实数，则由
`Im F=-4K rho`独立得到：F(U)>0时delta(U+)=0；F(U)<0时
delta(U+)=sign(A)*pi。该极限不使用根数作为输入，随后才与独立根计数比较。
F(U)落在原Pi误差预算经4K换算的范围内时拒绝认证；不放宽门槛。

根围道的左右边须处在解析gap内且包含已找到的根。固定1e-7边距可能将
近Mott根排除，得到数值很稳定但对象错误的零计数。当前边距按根到端点距离
缩小，并自适应提高围道节点，仍保留原边界非零及相邻辐角步长门槛。
这排除已知漏根方式，不是任意近端点根的区间算术完备证明。

### 有限Lth对照的解释边界

新路线直接实现无限热目标，因此不需要把有限Lth辅助实根当作物理介子。
`finite_L10_density`是相同正常分支/有限频率窗的数值对照；它没有重新证明
有限硬端点原核的所有辅助根及其补偿。其加权差不能改称旧有限Lth全谱算法的
严格误差上界，也不能据此授权忽略旧算法的辅助根。
对于生产候选，应使用无限热积分表示及其自动拒绝条件；保留有限Lth结果仅作
热尾敏感性诊断。q尾部连续区间衰减检验与omega条件界，也应分别标注为数值
收敛证据及依赖相位分支条件的界，不冒充全频解析定理。

## 当前决策

strict-audit第22.23节的自动归并5/5通过，支持无限热B路线在保留背景上的
研究级GBU部分产额数值可行。旧有限硬端点的全谱辅助根误差仍未获得认证，
不能用正常分支的小差值替代该证明。显式研究生产入口现由统一脚本对每个新
背景运行失败封闭的门禁；完整冻结线状态以该次运行的 manifest 为准，
不把代表点验收扩大为全部背景或默认生产授权。
