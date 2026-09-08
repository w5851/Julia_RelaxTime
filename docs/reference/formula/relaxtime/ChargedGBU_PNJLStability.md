# Charged GBU：PNJL 一致性、实极点权重与稳定性检查

## 定义、来源与授权

本文件是 [研究方法 v1](ChargedGBU_ResearchMethod_v1.md) 的后续审计，不是
production 处方选择。双线域与 retarded 的文献依据见
[Pereira 公式审查](Pereira2024_RegularizationAudit.md)；以下谱稳定性充分条件
是本项目的数学推导，不声称该论文已证明 PNJL-BQS 的完整稳定性。
实际运行和失败见 [strict-audit 第22.17节](../../../dev/active/2026-08-30_charged-rpa-bu-strict-audit.md)。

主候选仍使用 PNJL 占据、双线交集、同域两个单线项与 B0。普通 Fermi
仅作统计函数控制；真空/热分量均同时存在于 A 和 B，不能把 A 统称为真空项。
全硬截断和热项延伸是两种不同模型处方，固定同一上游背景并不使二者都成为
该背景的严格热力学二阶导数。任何数值通过均不授权切换默认、修改 PNJLCore、
重拟合参数、更新 baseline 或加入介子反馈。

## 1. PNJL 单线与冻结背景检查

独立检查占据补关系 f(x;Phi,PhiBar)+f(-x;PhiBar,Phi)=1、单调性及
Phi=PhiBar=1 的 Fermi 极限。q=0、同一热上限下检查 C(0)=A1+A2。
这只是同域代数；不等于当前冻结点已满足任意热处方下的 gap 方程。

对三个味重算

```math
\phi_i=\frac{N_c m_i A_i}{4\pi^2},\qquad
M_i=m_{0i}-4G\phi_i+2K_{\rm KMT}\phi_j\phi_k,
\quad K_{12}=G-\frac{K_{\rm KMT}}2\phi_s,\quad
K_{45}=G-\frac{K_{\rm KMT}}2\phi_d .
```

来源是项目 PNJLCore 的 calculate_mass_vec 与 MesonInteractionKernel 的 P
通道定义。分别记录全硬、上游默认热上限、热延伸24 fm^-1 的质量残差和
耦合变化；不把重算值写回冻结背景。未重算 Polyakov/守恒荷方程，不宣称
完整驻点或 Ward 恒等式认证。合成手征极限检查只验证等质量 q0 的
Pi(0)=-Nc*A/(4pi^2) 与一致耦合下 F(0)=0，不冒充真实 BQS 的 Goldstone 模。

## 2. 实 gap 极点不是由 unwrap 计数

令 F=1-4K_a*Pi、D=2K_a/F，内部坐标 lambda=k0+mu1-mu2。
在已认证开解析 gap 的简单零点 k0b 附近，

```math
D^R(k_0)\simeq\frac{2K_a}{F'(k_{0b})}
\frac1{k_0-k_{0b}+i0},\qquad
\frac{\operatorname{Im}D^R}{\pi}
=-\frac{2K_a}{F'(k_{0b})}\delta(k_0-k_{0b}).
```

这里沿用项目 ImD/pi 的谱符号，不能混用另一套传播子负号规范。
记录 W_b=-2K_a/F'(k0b) 并检查 k0b*W_b>0；正频根的负斜率对应正权重。
导数由同一线性谱的解析 Cauchy 导数计算，并以 gap 内有限差分对照。
有限差分起始步长为 min(1e-5,距 gap 端点/20) fm^-1，随后减半；必须
同时满足两次差分变化与解析导数误差均<1e-5（按 max(1,abs(F')) 归一化）。
只调整验证步长，不移动根、不改谱；最多12次，未通过则保留失败。
正负频根不对单个电荷密度重复相加。W_b 是传播子谱权重；
GBU 离散粒子计数仍为每个合格正频模一个 pi 跳变，不把 W_b 乘到粒子数上。
没有实 gap 根不排除 cut 内共振；简单极点检查不适用于精确 Mott 端点。

## 3. 被动性充分条件与 UHP 检查的边界

为避免与通常含1/pi的谱函数命名混淆，本节定义
rho_Pi(v)=ImPi^R(v)，其中 v 是外部 k0 坐标。对当前无任意减法的同一谱，
假设 K_a>0、rho_Pi(v)/v>=0、F(0)>0，并且零频极限与积分存在，则

```math
\frac{F(w)-F(0)}{w}
=-4K_a\int\frac{dv}{\pi}\,
\frac{\rho_\Pi(v)/v}{v-w},\qquad \operatorname{Im}w>0.
```

右侧虚部非正；若 F(w)=0，左侧 -F(0)/w 的虚部严格为正，矛盾。
因此满足这些条件的谱没有 UHP 零点。这是充分条件，不是必要条件：
**有符号热延伸谱不满足条件，不等于已发现不稳定极点。**

脚本逐线性单元检查 v*rho_Pi 的二次函数端点及内部极值，不只检查节点。
端点使用原始谱值，避免斜率-截距重建把精确零变成舍入负值；不修改谱。
sign_sufficient_condition_float64 是 Float64 条件检查，不是区间算术证明；
numerical_passivity 使用明确1e-10容差，不能据此宣布数学定理已适用于原连续核。

另对 UHP 矩形计算两种分辨率绕数。F 在 UHP 解析，围道不计求积节点为物理极点。
底边分别为 Im(lambda)=0.03、0.01、0.003 fm^-1，实部范围与上边由同一谱的
无限尾排除半径 R 决定。保留原门槛：相角步长<pi/4、边界 min|F|>1e-8、
两种绕数到同一非负整数的偏差<1e-7。对合并后的所有 anchor 区间一起细分，
不能只增加均匀网格而留下陡变区间不动。

这些围道仍留下 0<Im(lambda)<0.003 的近轴条带；它们检查插值谱，不证明
原连续核全 UHP 无零点。高能绝对谱矩界用于排除半径 R 外，不代替近轴检查。
实极点符号、全实 gap 计数、UHP 检查、Mott/Levinson 和密度积分分别保留状态。

rho_Pi 的被动性与 GBU 连续产额的符号不是一回事。后者含
d[delta-sin(2delta)/2]，即使传播子谱满足上述条件，连续项也不必非负。
不把本项充分条件检查静默新增为不可豁免的 production 门槛。

2026-09-06 的[热谱解释后续审查](ChargedGBU_ThermalAdmissibility.md)进一步区分：
上述充分条件失败不证明不稳定；但若要把函数解释为正热密度矩阵下的完整
物理关联函数，Lehmann/KMS 的谱符号是另一个必要条件。当前真空谱外的负热尾
可由 q0 解析式独立确定，触发的是该物理解释的作者评审，不是发现了 UHP 极点。

## 4. 单球对照的准确含义

分析探针 first_line_sphere_no_shift_of_boundary 精确定义为
|p1|<L，第二条线动量 p2=p1-q 不另设 cutoff。两个 A 和 B0 仍在该同域积分；
它不是旧 shifted/log B0 的复写，也不假定 Lambda>>q 的误差已经足够小。

角向 (p1,x) 与径向 (p1,p2) 两坐标互验；后者的区间是
|p1-q|<p2<p1+q，测度为 p1*p2/q。q=0 恢复双线域。
有限 q 下 A1 与平衡 A1 相同，A2 一般仍随 q 变化，因为被积函数依赖 p2。
交换两味并反射 lambda 后，第一条线单球变成另一个球，故可能失去
Pi_12(z,q)=conj(Pi_21(-conj(z),q))。数值两坐标一致并不修复这一处方差别。
交换对称平均或 centered 单球均是另一个明确处方，不在本轮自动采用。

可同时保留多个具名研究路线，但不能把“历史代码”“删第三条件的单球探针”
和“经审查的大截断近似”当成同一路线。先固定公式、坐标、适用范围并独立验收，
再考虑统一生产入口中的显式选择；禁止按结果是否好看自动 fallback。

## 5. 可重复入口与下一步

`audit_causal_gbu_pnjl_stability.jl` 仅读取冻结的 T170/muB240 BQS 输入。
`GBU_HEALTH_Q` 默认为0,1.4,3.2,6.2 fm^-1；
`GBU_HEALTH_MESH` 控制谱网格，`GBU_HEALTH_NODES` 控制直接 loop/单线求积，
两者均比较一半分辨率。`GBU_HEALTH_OUTPUT` 必须指向不存在的新目录。
输出 background_gap、q0_contacts、routing、summary、poles、contours、
mesh_checks CSV；evaluation 异常另存 failures。manifest 绑定源码、输入和输出。

`GBU_HEALTH_CHANNELS` 和 `GBU_HEALTH_VARIANTS` 可缩小复核范围，默认仍为
四通道、两处方；manifest 记录实际覆盖，子集通过不冒充全矩阵通过。
独立报告 `summarize_causal_gbu_pnjl_stability.jl` 使用 `GBU_HEALTH_INPUT/REPORT`，
检查源码快照、背景、输出哈希、行覆盖和记录中原门槛的归约。报告验证已有证据，
不重新认证原连续核；实 gap 的底层排除界由已测试 gap 审计器负责。

进程 exit0 仅表示诊断完成，不能读取为物理全部通过。两个处方的通过布尔、
各项未决风险、solver_called=false、density_computed=false、
production_authorized=false 均显式保留。测试使用合成占据、谱、极点及单球域。

随后先评审热延伸谱及有限 q regulator 的物理定义，解析近轴稳定性与连续核
误差，再做精确 Mott/Levinson、非线性 GBU 的分布/积分极限及外层积分联合验收。
原冻结线图继续作 diagnostic，不重新手动加密；完整曲线由后续统一研究入口触发。
