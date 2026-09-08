# CausalSpectralBubble

固定背景上的 ordered P/S 完整谱泡，研究诊断入口，不替换 `ChargedRPAProvider`
或 `MesonDensity` 默认。公式和 regulator 见
[charged 路线 §5.2.1](../../reference/formula/relaxtime/ChargedRPA_BU_ProductionRoute.md)。

## 构建与求值

```julia
g = RelaxTime.CausalSpectralBubble.build_spectral_bubble(
    q, m1, mu1, m2, mu2, T;
    channel=:P, Phi=0.3, PhiBar=0.4,
    vacuum_cutoff_inv_fm=Lambda, thermal_cutoff_inv_fm=16.0,
    momentum_nodes=128, angle_nodes=64)
r = RelaxTime.CausalSpectralBubble.spectral_bubble(g, k0; eta_inv_fm=0.01)
```

能量/动量参数均为 `fm^-1`，Phi/PhiBar 无量纲且在 `[0,1]`，q 非负，质量、T
和两个 cutoff 必须为正；热 cutoff 不小于真空 cutoff。q=0 角积分精确为 2，
其余 q 按两截断球交集求积。真空与热项始终分别使用完整的两线谱和顶角。

`r.value`、`contact_inv_fm2` 为 `fm^-2`，`B0` 无量纲，`derivative` 为 `fm^-1`。
`contact_identity_residual` 检查同域 C/B0 重构，绝不消费外部 A。`eta=0` 只允许
在正/负频率的保守开 gap 内直接求积；cut、端点和退化静态点抛出 `ArgumentError`。
有限 eta 的求积节点不是物理极点。

## 实轴重建

`spectral_bubble_cut(g,k0;energy_nodes=64)` 返回同域 Pi 的 pair、Landau、总虚部，
以及无量纲 `B0_imaginary`。零频率/零动量退化点保留未认证标记。

```julia
p = RelaxTime.CausalSpectralBubble.build_bubble_dispersion(g; segment_nodes=128)
value = RelaxTime.CausalSpectralBubble.cauchy_transform(p, k0 + mu1 - mu2)
```

注意 `cauchy_transform` 的变量是内部 **lambda**，不是 k0。实轴返回与分段线性
谱完全配套的 PV 实部和虚部；上半平面也可调用，不接受下半平面。
`segment_nodes>=8`，在运动学/cutoff 区间内使用余弦聚簇节点。q=0 cutoff 跳变
被一个网格单元近似，必须检验积分收敛；这不是带自动误差控制的求积器。

谱节点还包括同一两线截断域的完整 pair/Landau 极值，不仅是未截断阈值。
端点为几何意义的零测集；cut 层在计算前排除它们，支撑内部任意小的非零权重
仍保留。外部审计同时检查“插值谱超出几何支撑”和“插值遗漏几何支撑”，不能
用删除微小谱权重创建可认证的解析 gap。

`PiecewiseSpectralFunction(x,rho)` 供纯合成测试，要求节点有限、严格递增、
长度相同且至少 3，紧支撑两端 rho 为零；允许内部 rho 为负。

线性单元使用 `log1p` 与 `H(t)=1-log(1+t)/t` 的小参数展开，其中
`t=(b-a)/(a-z)`。这与原解析积分代数等价，但避免近重合节点的大斜率与远场
距离相乘后消减；不合并或删除谱节点，不改变 PV/retarded 分支。
合成窄单元及远场测试用256-bit独立积分表达式核验 Float64 结果。

## 授权边界

所有数值为 diagnostic；`production_authorized=false`。因果表示、截断选择、
节点收敛、物理束缚态完备性、Levinson/Mott 和最终密度授权是不同验收项。
不含 PNJL 平衡调用、介子反馈、phase fold、Landau 删除或负密度裁零。
测试位于 `tests/unit/relaxtime/test_causal_spectral_bubble.jl`；冻结输入入口为
`scripts/analysis/relaxtime/audit_causal_charged_bu.jl`。

后者使用 `CAUSAL_BU_INPUT/OUTPUT` 指定冻结输入和非覆盖输出目录，
`CAUSAL_BU_Q_NODES/Q_MAX` 开启 q 积分（不设则只计算 `CAUSAL_BU_Q` 壳层）；
`CAUSAL_BU_MESH/W_NODES/W_MAX/THERMAL_MAX/LOWERS` 控制谱、频率和端点。
`audit_causal_charged_eta.jl` 单独检验 q=1 的有限 eta GBU 积分趋近 PV 加离散根，
不是阈值点值检验。运行证据与有效配置见 strict-audit 第21节。

后续固定背景的 v1 合同、同背景 temp7/外推比较、独立围道计数及稀疏冻结线
使用 [研究方法文件](../../reference/formula/relaxtime/ChargedGBU_ResearchMethod_v1.md)
和 analysis-only `audit_causal_gbu_*.jl` 入口。执行和失败记录见 strict-audit 第22节。
其中 `build_bubble_dispersion(...; energy_nodes=64)` 控制 cut 的一维求积，
不是 grid 的二维 `momentum_nodes/angle_nodes`；两者的加密证据必须区分。

`audit_causal_gbu_full_closure.jl` 在一个冻结背景上比较两种 loop 坐标，校验
同域 contact 与 Cauchy 表示，并独立检查全部有符号实轴 gap、端点条带和无限尾。
`GBU_FULL_Q=density_shells` 复用24个密度q节点并加入6个诊断q值，不计算新的
背景或密度；显式列表必须包含0。`GBU_FULL_LOOP_NODES/MESH` 分别控制直接loop
与插值谱，`GBU_FULL_OUTPUT` 必须是新目录。报告入口
`summarize_causal_gbu_full_closure.jl` 使用 `GBU_FULL_INPUT/REPORT` 核对快照、
输入、输出哈希、行覆盖和原门槛。两者均不启动production或选定regulator。
“全部实gap”仅指各紧支撑插值谱；cut内零点、UHP不稳定模、精确Mott端点和
未插值的连续核完备性不在该声明内。

## Pereira 独立公式对照

`scripts/analysis/relaxtime/causal_gbu_pereira_reference.jl` 用同一双球交集的
柱坐标直接计算 pair/scattering、两个单线项和 Pi，并以独立 (E,epsilon)
域计算 retarded cut，不复用本模块的四留数、角向域或 Cauchy 实现。
这属于 analysis oracle，不是新增生产 provider。

`audit_causal_gbu_pereira.jl` 冻结已有 BQS 背景，先比较 Fermi/PNJL
同硬截断，再比较 PNJL 热尾8/12/16/24 fm^-1；不求新平衡或密度。
`GBU_PEREIRA_OUTPUT` 要求新目录，`GBU_PEREIRA_PDF` 可绑定原论文哈希；
节点可通过 `CausalGBUPereiraAudit.main(nodes=384)` 显式指定，粗网格取一半。
完整公式、来源页码、单位及弱极限边界见
[Pereira 审查](../../reference/formula/relaxtime/Pereira2024_RegularizationAudit.md)。
测试为 `tests/unit/relaxtime/test_causal_gbu_pereira_reference.jl`。

## PNJL 与极点稳定性分析

`causal_gbu_pnjl_stability.jl` 提供 analysis-only 的同谱解析导数、实极点
权重、逐单元谱符号和 UHP 矩形检查，以及第一条线单球的两坐标探针。
它不改变本模块或 provider 的默认行为。`audit_causal_gbu_pnjl_stability.jl`
用冻结背景比较全硬和热延伸、核验 q0 contact 和三个味 gap 恒等式，
不重求背景、不积分密度。`GBU_HEALTH_Q/MESH/NODES/OUTPUT` 分别控制
q列表、谱节点、loop节点和非覆盖输出目录。

频率/q/热上限为 fm^-1；Pi/contact 为 fm^-2，极点斜率为 fm，
传播子离散权重为 fm，v*ImPi 为 fm^-3，逆传播子与绕数无量纲。
`sign_sufficient_condition_float64` 不是原连续核的严格 UHP 证明；
`numerical_passivity` 也不等同于密度正性门禁。矩形不覆盖近轴条带，
`full_continuum_stability_certified` 与 production 授权保持 false。
完整定义与限制见 [公式说明](../../reference/formula/relaxtime/ChargedGBU_PNJLStability.md)。

`GBU_HEALTH_CHANNELS/VARIANTS` 可显式缩小复核范围，不改变默认四通道/两处方。
导数检查在 gap 内减半步长，必须同时满足两次差分变化及解析导数对照；
原门槛不变。`summarize_causal_gbu_pnjl_stability.jl` 通过
`GBU_HEALTH_INPUT/REPORT` 独立核验留存哈希、覆盖和门槛归约，拒绝覆盖报告。

## 热谱解释与原连续核探针

`causal_gbu_thermal_admissibility.jl` 是 analysis-only oracle：
q0 壳上 pair 权重、独立原 cut、无谱插值的 PV/UHP Cauchy 求积及 RPA 谱符号。
求积采用解析常数减法，误差字段是总绝对误差的估计，不是严格/统一上界。
输入断点必须有限且严格递增；不接受下半平面或外部支撑端点的实轴探针；
调用者须避免内部硬跳变点。未收敛保留 false，不改谱或放宽目标。

`audit_causal_gbu_thermal_admissibility.jl` 的 `GBU_THERMAL_OUTPUT` 指定新目录；
`settings()` 固定覆盖和门槛。`summarize_causal_gbu_thermal_admissibility.jl`
通过 `GBU_THERMAL_INPUT/REPORT` 核验源码/输入/PDF/输出哈希、行覆盖与布尔归约。
同一核 Pi/contact 为 fm^-2，D 为 fm^2，k0/lambda/q/eta 为 fm^-1。
`continuous_UHP_count_certified`、密度计算和 production 授权均为 false。

负热尾在独立 q0/有限q cut 中同样存在。当前完整物理谱解释需要评审，
不能读作“数值通过即可用于粒子产额”。具体推导和候选边界见
[热谱适用性](../../reference/formula/relaxtime/ChargedGBU_ThermalAdmissibility.md)。

方向B前置分析见
[静态基础与响应匹配](../../reference/formula/relaxtime/ChargedGBU_DirectionB.md)。
`audit_causal_gbu_direction_b.jl` 以 `GBU_DIRECTION_B_OUTPUT` 指定新目录；
只检验有效Models入口的动能势偏导与独立积分，不重求平衡。附带的线性谱
两次减法反例是合成模型，绝不是本模块检出的极点或可直接采用的修正。
`is_stationarity_test`、`is_dynamic_hessian_test`、`is_project_pole` 和
production授权均false；本API的有效实现未改变。

## 方向B源泛函核验（analysis-only）

`causal_gbu_source_functional.jl` 明确固定动量投影后的fermion行列式，
并以“真空Lambda + full Lth - 同源真空Lth”定义二阶响应。
独立Dirac矩阵投影与静态8x8源行列式曲率用于检查本模块的核，
不是替换本模块的新物理provider。推导见
[方向B第6节](../../reference/formula/relaxtime/ChargedGBU_DirectionB.md)。

`audit_causal_gbu_source_functional.jl` 使用 `GBU_SOURCE_OUTPUT` 指定新目录。
`main(nodes=128)` 比较两种内部求积节点；静态源步长为
0.004/0.002/0.001 fm^-1，另使用32/64内部节点。动态频率是external k0，
Pi与源曲率均为fm^-2；局部与积分级合成测试覆盖P/S Dirac迹和源归一化。
源曲率检查不是相移导数/分部积分，也不是完整PNJL/KMT驻点认证。

固定投影不随源变动，故此定义没有额外delta P顶角，但减法核负热尾仍在。
`physical_response_selected`、`positive_physical_trace_certified`、
`GBU_observable_derived`、全UHP计数和production授权均false。
不修改src实现、旧默认或baseline，不求新背景或密度。

## 方向B的GBU观察量审查（analysis-only）

`causal_gbu_observable_closure.jl` 按本项目归一化Sigma_M=2Pi检验
W=delta-Im(Sigma_M)*Re(D)及其导数。它不选择全局相位、不认证根或密度。
`tagged_pressure`/`partial_density_terms`只对合成固定profile展示部分记账偏导
与全导数、端点项的区别；不是新增fugacity、生产密度或热力学provider。
尾项界需要统一界和近零分支同时认证；外圆零点排除需要真正的谱总变差上界。
默认认证均false。有限采样误差估计不能提升为这些上界。

`audit_causal_gbu_observable_closure.jl`以`GBU_OBSERVABLE_OUTPUT`指定新目录，
读取既有BQS背景，比较无插值PV/原cut和局部GBU代数，并展示先减法后RPA
与各自重求和后相减的区别。实际探针不求相移导数、新根或密度。
公式与文献边界见[观察量闭合条件](../../reference/formula/relaxtime/ChargedGBU_ObservableClosure.md)。

## 原连续核端点与源参考延续（analysis-only）

新增causal_gbu_endpoint_closure.jl、causal_gbu_continuation_gate.jl和
causal_gbu_source_continuation.jl，均不是本生产核的新provider或默认开关。
三个audit入口的GBU_ENDPOINT_OUTPUT、GBU_CONTINUATION_OUTPUT、
GBU_SOURCE_CONTINUATION_OUTPUT要求新输出目录，绑定源码与保留BQS背景。
前两者不积分密度，后者只检验静态源曲率与动态Pi。

端点距离字段是log或log10[(r-S)/(fm^-1)]的渐近估计，不是求根残差；
单根Bose上界不等于总密度界。外圆界为实数算术证明，未作有向舍入认证。
continuation的两阶求积差不是统一误差界；原低能Mott检查不包括辅助UV根。
宽坐标只处理无法由Float64表示的panel内部节点，不删panel或改cut。

source_continuation参考保持最低四个Dirac海能带，要求reference_gap>0及
sqrt(2)*source_step<reference_gap；不再要求每味|mu|<M。旧源oracle契约不变。
温度/能隙/频率/源步长为fm^-1，Pi及曲率为fm^-2；拒绝隙关闭或非法输入。
完整公式与状态见[端点闭合说明](../../reference/formula/relaxtime/ChargedGBU_EndpointClosure.md)。
