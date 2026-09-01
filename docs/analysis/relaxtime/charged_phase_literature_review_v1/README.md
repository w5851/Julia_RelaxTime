# charged-RPA/BU 严格路线定向文献综述

**状态**：研究背景与公式路线审计；不是 production 授权、数值基线或实验拟合结论。
**检索截止**：2026-09-01。
**仓库状态**：`codex/charged-rpa-bu-strict-phase-backend`，HEAD `81fcc26d8a7c372db44403cfb0be20869f195001`。
**问题范围**：2+1 味 PNJL/NJL 中的 (K^+/pi^+)、(K^-/pi^-)，不等质量有序 charged RPA，相移的 (i0^+)/PV 极限，束缚态计数，Levinson/Mott 闭合，标准 BU/GBU 密度，完整 KMT 耦合以及冻结线/守恒荷约束。

## 1. 先给结论

### 1.1 没有一篇论文能直接给出本项目的“唯一现成生产流程”

公开文献可以拼出一条自洽的模块化路线，但没有找到一篇同时满足以下全部条件的已发表工作：

1. 完整三味 KMT (mathcal K_{ab})；
2. BQS 或 (ho_Q/ho_B,ho_S) 约束下的上游解；
3. (uar s)、(sar u) 有序 charged kaon RPA；
4. 严格实轴 (i0^+)/PV 相移；
5. q-resolved 标准 BU 与 GBU 密度；
6. 沿同一化学冻结线输出两个 charged ratio；
7. 介子 (Omega_M) 反馈回完整守恒荷和平均场方程。

因此本项目需要做的是**选择一个来源充分、约定明确、从上游到下游闭合的组合路线**，而不是声称文献已经唯一决定所有实现细节。

### 1.2 相移的物理约束是唯一的，数值分支实现不是唯一的

文献共同支持以下物理要求：

- 相移来自 retarded/advanced 传播子或相应 (S)-matrix 的复对数；
- 相位分支必须连续，并与高能端和阈值端点一致；
- 束缚态、共振和散射连续谱必须一起考虑；
- Mott 转变时束缚态从阈下进入连续谱，阈值相位应发生相应跳变；
- Levinson 型端点关系必须闭合，通常写成
  [
  delta(omega_{m th})-delta(infty)=pi N_{m bound}
  ]
  （半束缚态需另加约定）。

但以下实现选择并没有在文献中被唯一规定：`unwrap` 的具体算法、有限网格上的 `high_energy_zero` anchor、(eta) 的数值大小、异常低能模的删除方式、有限 (omega_{max}) 如何代表无穷远，以及 `arg(propagator)` 与 `arg(inverse propagator)` 的代码层约定。它们必须由本项目显式记录并经过端点、Levinson、Mott、节点和截断测试。

### 1.3 标准 BU、Breit--Wigner、稳定粒子和 GBU 是不同近似对象

稳定粒子极限只保留一个窄壳准粒子；reduced BW 用一个极点/宽度近似；标准 BU 使用完整相移导数；GBU 的

[
F_{m GBU}(delta)=delta-	frac12sin(2delta)
\]

来自特定的 (Phi)-derivable/广义 BU 推导，用于处理连续谱与重复计数，不能被当作标准 BU 的数值稳定化因子。四种结果不应被要求收敛到同一数值；应分别报告其公式、适用范围和失败状态。

### 1.4 完整 KMT 交叉耦合主要直接进入中性混合扇区

Rehberg 1996 以及 Mei 等人的现代实现给出

[
K_{03},K_{30},K_{38},K_{83}propto sigma_u-sigma_d,
]

它们直接改变 (pi^0)-(eta)-(eta') 的 (0/3/8) RPA 矩阵。带电 kaon 的非对称映射是 (K^+leftrightarrow K_{45})、(K^-leftrightarrow K_{45}) 的共轭通道，(K^0,ar K^0) 对应 (K_{67}) 及其共轭。因而在当前没有中性混合反馈的 charged-only 后端中，(K_{03}/K_{38}) 不应被直接归因于 (K^pm/pi^pm) 密度变化；它们只有通过重新求解凝聚、质量和极化函数产生间接影响。

### 1.5 反馈是另一项热力学闭合问题

大量 PNJL K/π 工作先求 quark-only MFA，再以后处理方式计算介子极点、相移和密度。Gaussian/one-meson-loop/(Phi)-derivable 工作表明，把 (Omega_M) 放回驻点方程需要同时处理自能、极化、支撑、正则化、重复计数、Ward/Goldstone 关系和额外玻色 cutoff。当前 partial feedback 只能作为守恒荷敏感性诊断，不能标记为完整平衡。

## 2. 综述方法与边界

### 2.1 检索问题

本综述只筛选能回答下列问题的文章：

1. 如何从 PNJL/NJL 的夸克传播子和 RPA 传播子定义相移；
2. 如何处理实轴分支、(i0^+)、PV、束缚态、连续谱和高能尾；
3. 如何用 Levinson/Mott 检查相移是否物理闭合；
4. 如何将相移放入标准 BU 或 GBU 压强/密度；
5. 不等质量 (uar s/sar u) charged kaon 如何排序；
6. 完整 KMT 耦合、显式 (mu_I)、BQS 和介子反馈在流程中处于哪一层；
7. 化学冻结线和实验 (K/pi) 激发函数能提供什么数量级/趋势约束。

这不是按所有数据库逐条穷尽的系统综述，因此“所有文章”在这里具体化为：核心模型、相移/BU、KMT/同位旋、反馈和冻结线五个主题中，能够提供可复核公式或流程的代表性已发表工作；同一工作的 arXiv、期刊和会议版本去重。

### 2.2 数据库与检索记录

- 公开全文优先：arXiv、APS/SCOAP3、CPC/IHEP、INSPIRE 元数据、出版社 DOI 页面。
- 本轮在用户明确确认具有授权后，使用 Chrome 的 XJTU 网关逐篇核对少数确需出版社原文的条目；账号、密码、MFA、cookie 和浏览器存储不进入仓库。核心方法仍优先引用公开全文，机构路由只补充正式版题名/页码/PDF provenance。
- INSPIRE MCP 检索工具在本运行环境中未提供，因此使用公开 arXiv/出版社/本地已下载全文作为替代，并在 `manifest.json` 中记录该限制。
- 可复查查询族见 `tables/search_log.csv`；候选逐篇筛选见 `tables/source_screening.csv`；主题证据矩阵见 `tables/evidence_matrix.csv`。

### 2.3 纳入/排除规则

**纳入**：至少有一项可复核的 MFA/RPA 公式、相移/BU 方程、KMT/同位旋混合、Mott/Levinson 或介子反馈方程，并能说明适用边界。
**排除**：只有实验数据而没有方法信息；纯综述/重复会议摘要；只研究与 charged K/π、相移或反馈无关的 PNJL 现象；只给数值图且不给可重建公式/流程。磁场论文只有在明确给出完整 (0/3/8) 耦合或矩阵 RPA 时保留为结构证据。

## 3. 文献给出的可复用公式闭合

### 3.1 从微观拉氏量到 RPA

共同的结构是

[
mathcal L_{m NJL/KMT}
longrightarrow 	ext{MFA 背景 }(M_f,Phi,sigma_f,mu_f)
longrightarrow Pi^R_{ij}(omega,q)
longrightarrow D^R_{ij}(omega,q)
longrightarrow delta_{ij}(omega,q).
]

无混合的标量通道可写成某种归一化下的

[
D^{-1,R}_{ij}=1-2K_{ij}Pi^R_{ij},
]

或项目旧接口使用的等价因子约定。矩阵通道则必须使用矩阵逆和 determinant/eigenphase；例如中性 (0/3/8) 扇区的物理相位应从

[
delta_{m mat}=-operatorname{Im}lndet D^{-1,R}
]

或等价的 retarded/advanced 比值构造，而不能在强混合时任意挑一个固定分量的相位。

### 3.2 相移定义与分支

Blaschke 等的 GBU 论文给出的代表性定义是

[
delta_X(omega,q)=-operatorname{Im}ln[eta^2S_X^{-1}(omega-mu_X+ieta,q)].
]

PNJL K/π 论文用等价的传播子/逆传播子 arctan 写法。关键不是某一行 `atan`，而是：

1. 选择 retarded 或 advanced 边界；
2. 统一传播子还是逆传播子；
3. 在 (omega) 上连续跟踪复对数分支；
4. 让相位满足阈值束缚态计数和高能端点。

本项目的 `delta=-arg(Delta^R)`、反向 unwrap 和 `high_energy_zero` 是一种**数值实现约定**，不是文献中唯一指定的算法。文献能约束目标端点，不能替本项目选择有限网格上的 anchor。

### 3.3 (i0^+) 与 Cauchy 主值

严格实轴边界的数学对象是

[
Pi^R(omega)=operatorname{PV}Pi_{m real}(omega)-i,operatorname{Im}Pi(omega),
]

其中实部在奇点处使用 Cauchy 主值，虚部由解析 continuation/切割结构给出。项目 `OneLoopIntegrals.B0` 的“PV 实部 + 解析 imaginary cut”最接近 (eta	o0^+) 的严格实现；`B0_retarded` 直接在 (omega+ieta) 上做有限 (eta) 复积分，只能作为 retarded probe 或交叉检查，不能自动等同于严格极限。有限 (eta) 还会移动窄峰、平滑阈值并改变 branch tracking，因此必须记录 (eta)、PV gap、端点和网格。

### 3.4 BU、密度和 GBU

标准 BU 的典型结构为

[
p_X=-d_XT!int!rac{d^3q}{(2pi)^3}
int_0^infty!rac{domega}{2pi}
Bigl[ln(1-e^{-(omega-mu_X)/T})+ln(1-e^{-(omega+mu_X)/T})Bigr]
partial_omegadelta_X.
]

对 (mu_X) 求导得到带正、负电荷的密度；本项目的单电荷 `domega/pi` 版本必须把“只算一个 charged species”与粒子/反粒子是否显式相加写清楚。GBU 在特定 (Phi)-derivable 推导中把相位导数替换为

[
partial_omega!left[delta_X-	frac12sin(2delta_X)ight],
]

以避免部分连续谱重复计数。它不是所有相移密度都应默认使用的普适修正。

### 3.5 束缚态、Levinson 与 Mott

标准检查为：

1. 独立确定每个 (q) 的阈下 pole/bound-state 数 (N_{m bound}(q))；
2. 独立确定连续谱阈值 (omega_{m th}(q))；
3. 检查 (delta(omega_{m th})-delta(infty)=pi N_{m bound})；
4. 在 Mott 点检查束缚态从阈下进入连续谱，阈值相位从 (pi) 变为 0（依赖整体相位 convention）；
5. 检查相位导数的 continuum 贡献能否补偿 pole/resonance 贡献，使高温/高能相关项消失。

当前 `count_subthreshold_roots` 只基于有限实轴网格的 sign change；复杂 pole、双根、阈值切线、有限 (eta) 和多个 (q) 分支都可能漏计。因此 `q -> 0` 或 sign-change 计数只能作为诊断输入，不是生产级束缚态计数。

## 4. 按主题筛选出的有效信息

### 4.1 相移/Levinson/BU 主线

- **Hüfner--Klevansky--Zhuang--Voss (1994)**：把 (1/N_c) 介子与夸克散射贡献写成束缚态和散射相移的 BU 形式，是“相移进入热力学”的早期场论实现。
- **Wergieluk et al. (2013)**：PNJL pion Mott 解离；明确指出 BW 谱展宽单独不能满足 Levinson，必须加入 scattering continuum。
- **Blaschke et al. (2014)**：从玻色化后的 Gaussian fluctuation 得到 meson/diquark 相移和 GBU；明确给出相位、连续谱/共振分解、阈值 bound-state step 和 Levinson 闭合，同时提醒该 phase 是谱参数而非可直接观测的 on-shell scattering phase。
- **Dubinin et al. (2015, 2017)**：把该路线扩展到 2+1 味 pion/kaon；给出 unequal-mass charged kaon、Mott 以及 anomalous low-energy (K^+) 模式，说明完整 continuum 对 K/π ratio 的影响。
- **Maslov--Blaschke (2023)**：保留全动量 dependence，包括 spacelike Landau damping，并展示介子压强对动量 cutoff 的敏感性；这证明“相移公式写对”仍不等于“数值积分已收敛”。

**对本项目的直接启示**：先把 phase object、branch、measure、bound count、threshold、high-energy tail 和 Mott gate 分开测试，再比较 stable/BW/BU/GBU；不能用一张复现图或一个正数 ratio 代替这些门禁。

### 4.2 不等质量 charged kaon 与 anomalous mode

[
K^+=u\bar s,quad K^-=s\bar u,quad K^0=d\bar s,quad \bar K^0=s\bar d.
]

带电通道的 bubble 必须按有序的 flavor pair 处理，化学势差进入 (Pi_{ij}(omega,q))。当 (M_s\neq M_{u,d}) 或 (mu_s\neqmu_{u,d}) 时，粒子/反粒子通道不再是简单相同副本，K+ 可出现低能 anomalous mode。文献讨论了其对 ratio 的敏感性，但没有给出一个可直接移植的“自动删除 anomalous branch”标准算法；因此本项目应保留 full、显式分支和删除分支三种诊断，而非静默裁剪。

### 4.3 完整 KMT 与同位旋不对称

Rehberg et al. (1996) 给出了三味 KMT 收缩后的完整有效耦合；Mei et al. (2023) 在 (sigma_u\neqsigma_d) 的现代算例中显式构造了 (3\times3) 的 (0/3/8) pseudoscalar RPA 矩阵。Xia--He--Zhuang (2013) 则展示了显式 (mu_I) 和 pion condensate 时，KMT 会使 strange sector 与 light sector 耦合，并需要扩大 MFA/BSE 变量。

可复用的层次关系是：

1. 完整微观 KMT 拉氏量决定 MFA 和二次涨落的相互作用核；
2. 若背景保持对角且无 pion/kaon condensate，charged 通道仍是有序单通道，(K_{03}/K_{38}) 不直接出现在 charged propagator；
3. 强中性混合或凝聚出现时，必须使用矩阵 determinant/eigenphase，而不能沿用 2×2 或单分量 phase；
4. 把完整核接入 RPA 与把介子压强反馈回 (Omega) 是两个独立扩展。

### 4.4 (mu_u,mu_d,mu_s)、BQS、(mu_I) 与冻结线

文献存在三种不可直接互换的输入方式：

- 现象学情景：(mu_u=mu_d)，另取 (mu_s=0.2mu) 或 (0.55mu)；
- 显式同位旋：例如 (mu_u=mu_B/3+mu_I/2, mu_d=mu_B/3-mu_I/2)，但不同论文对 (mu_I) 归一化不同；
- 守恒荷求解：固定 ((T,mu_B))，解 (mu_Q,mu_S) 使 B/Q/S 密度约束满足。

Friesen 2019、Blaschke 2020 的 K/π 工作主要属于第一类，并把介子谱/密度作为 quark-only 背景的后处理；不能把 (mu_s=0.55mu_u) 当成 BQS 平衡方程的推论。实验层面 NA49 的 (K^+/pi^+) 激发函数显示低 SPS 能区的非单调 horn、较高能区趋于较平坦；(K^-/pi^-) 通常随能量上升。该信息适合作为趋势和数量级 sanity check，不是当前模型的逐点拟合目标。化学冻结线参数化（如 Cleymans 等）也只能固定一条实验拟合轨迹，不能自动决定上游 PNJL 的 flavor chemical-potential closure。

### 4.5 介子反馈与“完整平衡”

Blaschke 2014 的 Gaussian BU、Oertel--Buballa--Wambach 的 (1/N_c)/meson-loop 比较、Pereira--Costa 2020 的 one-meson-loop gap equation、Blaschke--Ebert 2017 的 variational backreaction 共同说明：

- 介子反馈需要把 (P_M) 对 (M_f,sigma_f,Phi,mu_f) 的导数带回所有驻点方程；
- meson-loop momentum cutoff 和真空 subtraction 是新参数/新门禁；
- 只把已算出的 (Omega_M) 数值加回而不改 quark self-energy 会丢失热力学自洽；
- (Phi)-derivable/one-loop 方案必须审计 double counting 以及 Goldstone/Ward identity。

所以当前 partial feedback 应保留为 diagnostic；它不能用来证明 (ho_s+ho_K=0) 已经实现了完整热力学平衡。

## 5. 与当前仓库实现的逐模块映射

| 代码模块 | 当前行为 | 文献支持的部分 | 仍缺的生产门禁 |
|---|---|---|---|
| `src/relaxtime/OneLoopIntegrals.jl` | `B0` 为 PV 实部加解析切；`B0_retarded` 为有限 (eta>0) 复积分 | 与严格 retarded 边界的分解相符 | PV gap、cut、(eta	o0) 和节点的联合收敛 |
| `src/relaxtime/ChargedRPAProvider.jl` | `ordered_retarded`/legacy ordered bubble；有序单通道 | 对应 unequal-mass charged (q\bar q') RPA | 传播子归一化、pole/threshold 和 ordered/legacy 固定点闭合 |
| `src/relaxtime/ChargedPhaseBackend.jl` | `-arg(\Delta^R)`、unwrap、tail/Levinson/Mott gate | 实现了文献约束的接口骨架 | 物理 bound count、统一高能 anchor、真实 Mott 点和 branch continuation |
| `src/relaxtime/BUPhaseGates.jl` | sign-change 阈下根计数和端点 phase gate | 可作快速 diagnostic | 独立 pole solver 或可靠实轴束缚态判定 |
| `src/relaxtime/MesonDensity.jl` | stable/BW/BU/GBU 多算法，当前仍有低节点/截断诊断 | 与文献的四种近似对象相符 | `d\omega/\pi` 语义、Bose support、负密度状态、q/ω/截断收敛 |
| `src/relaxtime/EffectiveCouplings.jl` | charged 兼容层和完整 kernel diagnostic | (K_{45}/K_{67}) 的非对称映射 | 真实完整 RPA/BU 归一化及 neutral matrix 后端 |
| `src/models/solver/spec/ConstraintModes.jl` | 固定 BQS quark-only 平衡 | 与守恒荷输入层相符 | 介子 (Omega_M) 反馈和全 hadronic BQS |

### 5.1 当前严格分支的状态判定

当前真实 ordered profile 已能输出有限正数的 charged strict-BU 诊断，但已有记录显示阈下 root/Levinson gate 失败、部分 high-energy tail 不稳定，`q -> 0` 计数仍是临时输入。因此本分支是“严格公式后端探索性实现”，不是 production candidate。已有负密度审计还显示，未统一的 high-energy anchor 可在不同 q 壳层施加不同符号平移；不能简单全局翻转相位或把负值裁剪为零。

## 6. 筛选后可执行的路线建议

### 6.1 推荐的文献闭合顺序

1. 固定一种相位 convention：retarded inverse propagator 的复对数，明确 sign、branch 和端点。
2. 以 `B0` 的 PV+analytic-cut 作为 (eta	o0^+) 主线；`B0_retarded` 只做有限-(eta) 交叉检查。
3. 对每个 (q) 独立确定 ordered (u\bar s/s\bar u) 的阈值、pole 和 (N_{m bound})。
4. 以端点/Levinson/Mott gate 为必要条件，再做 (omega,q,eta,omega_{max}) 收敛。
5. 单独比较 stable、reduced BW、q-pole BW、标准 BU 和 GBU；不强求算法数值相同。
6. 完整 (K_{45}) charged kernel 通过上述门禁后，才将其作为 production candidate；(K_{03}/K_{38}) 仍作为 neutral-mixing 后续支线。
7. 最后才沿 `baseline_freezeout` 输出 (K^+/\pi^+)、(K^-/\pi^-)，并与 NA49/STAR 的趋势和数量级比较。

### 6.2 暂不做的事项

- 不把 `high_energy_zero` 当作已证明的 (delta(\infty)=0)；
- 不以 `count_subthreshold_roots` 替代物理束缚态计数；
- 不在当前 strict branch 中启用完整 (Omega_M) 反馈；
- 不修改 `PNJLCore`、transport 或旧 production 默认；
- 机构订阅路径只在用户明确确认授权且公开全文不足、或确需出版社版本核对时使用；遇到 CAPTCHA、MFA 或其他受保护流程则停止，不绕过。

## 7. 证据强度与待作者确认项

证据等级定义：

- **A：直接全文公式**：可逐式核对并映射；
- **B：全文方法流程**：流程清楚但仍有数值约定缺口；
- **C：摘要/项目映射**：只支持方向或适用边界，不足以授权实现。

尚需作者确认/后续实验的关键项：

1. 项目最终采用 `delta=-arg(D^{-1,R})` 还是严格 retarded/advanced log ratio；
2. (omega	oinfty) 的数值锚定如何从统一尾部测试导出；
3. 阈值半束缚态如何计数；
4. anomalous (K^+) 低能模是保留、显式拆分还是做物理投影；
5. GBU (sin^2delta) 权重的 (Phi)-derivable 适用条件；
6. 生产结果是否只采用 quark-only BQS，还是另立完整反馈项目；
7. 实验比较采用哪些 acceptance/rapidity 定义，避免把 mid-rapidity yield 与模型全动量密度直接等同。

## 8. 文件索引

- `manifest.json`：检索日期、仓库 SHA、数据库/查询族、边界和来源清单。
- `tables/search_log.csv`：可重放的检索族、目的、公开全文状态和限制。
- `tables/source_screening.csv`：候选逐篇纳入/排除、路线、可用信息、局限和来源。
- `tables/evidence_matrix.csv`：主题到公式/流程/代码门禁的映射。
- `tables/claim_ledger.csv`：可写入项目文档的主张、证据、强度和禁止外推。
- `tables/publisher_pdf_provenance.csv`：逐篇记录出版社页、授权路径、PDF 文件校验、版权边界和失败/错配原因；PDF 本体保留在 `D:/Desktop/paper/dev/outputs/gateway_fetch`，不纳入 Git。

## 9. 出版社原文核对与 PDF provenance（2026-09-01）

本轮只处理与当前公式路线最直接相关、且公开全文不足以完成版本核对的少数条目。完整字段见
`tables/publisher_pdf_provenance.csv`。

| source | 结果 | 可用于项目的范围 |
|---|---|---|
| S04 Rehberg et al. 1996, PRC 53, 410 | APS 正式题名、卷页和 PDF 均通过；PDF 20 页，`pdftotext` 可读，SHA-256 已记录 | 可作为完整 KMT/带电通道 legacy 公式的出版社版本证据；PDF 本体不入 Git |
| S02 Dashen–Ma–Bernstein 1969, Phys. Rev. 187, 345 | APS 文章页明确显示 XJTU 授权；渲染后确认 PDF 下半页包含正确标题、摘要和正文，页 1 顶部残留另一篇液态合金电阻率文章的版面片段，页 2–26 为目标文章 | 可用于 S-matrix/Levinson 公式核对，但引用页 1 时避开顶部残留并保留该版面 caveat |
| S03 Hüfner et al. 1994, Ann. Phys. 234, 225 | ScienceDirect 文章页显示 XJTU 机构标识，但明确写明 “Xi'an Jiaotong University does not subscribe to this content”；此前入口还触发 CAPTCHA | 版权全文无授权可用；不绕过 CAPTCHA，继续使用已有公开/本地证据 |

版权边界：项目仓库只保存 DOI、出版社 URL、授权路径、文件大小、页数、SHA-256、可读性探针和失败原因；受版权保护的 PDF 保留在项目外的 `D:/Desktop/paper/dev/outputs/gateway_fetch`，不提交、不推送。

本文件与既有的 `docs/analysis/relaxtime/meson_isospin_workflow_literature_survey_2026-08-28.md` 互补：后者覆盖化学势、KMT 与冻结线背景；本文件集中审计 strict charged phase/BU 的公式闭合和生产门禁。

## 10. S02 重新核对后的影响（2026-09-01）

本轮对 Dashen--Ma--Bernstein 的 APS PDF 做了逐页文本和渲染复核。PDF 共 26 页：
页 1 顶部有另一篇文章的残留版面，但页 1 下半页已经包含目标论文的标题、作者、
摘要和 Introduction，页 2--26 为目标论文正文。因此此前仅依据第一页首行作出的
`content_mismatch` 判断已撤回；当前 provenance 状态为
`publisher_pdf_verified_first_page_overlay`。该版面问题不改变正文公式的可读性，但
引用页 1 时必须避开顶部残留。

S02 对前轮结论的实质影响如下：

1. **形式对象需要分层。** 论文的核心公式使用连通、on-shell 的 S 矩阵，并在能量上
   计算 `S^{-1}` 的双向导数（论文 Eq. (4.32)、(5.10)）。它没有直接规定
   retarded/advanced 传播子的复对数、`unwrap` 或 `high_energy_zero`。因此前轮关于
   “连续分支和端点条件不可省略”的结论保持不变，但应把“DMB 直接证明 propagator
   phase”改写为“DMB 给出 S-matrix 祖先，propagator-to-S-matrix 映射由后续场论文献
   另行确定”。
2. **测度归一化得到更强的可核对依据。** 对单一弹性本征道写
   `S(E)=exp(2 i delta(E))` 时，DMB 的 `(4 pi i)^{-1} S^{-1}<->d_E S`
   化为 `(1/pi) d_E delta`。所以项目使用 `domega/pi` 的前提是被积相位确实是
   `delta`；若使用完整 `arg(S)=2 delta`，等价写法才是 `domega/(2 pi)`。这并不
   自动确定 `-arg(Delta^R)` 在项目中对应哪一个相位归一化，仍需 phase-normalization
   gate。
3. **Levinson/Mott 的证据强度需要加条件。** DMB 明确把束缚态作为离散谱处理，并说明
   窄共振在相位快速变化且热权重近似不变时可近似成粒子；同时其相对论性 Levinson
   推广是基于非相对论结果的合理猜测，而非从当时相对论量子力学第一性推导。因此
   项目的 Levinson/Mott gate 是必要的模型一致性门禁，但不能表述成对所有相对论
   PNJL 相移都无条件成立的定理。
4. **化学势层得到独立支持。** 论文在守恒律部分以 (B)、第三同位旋分量 (I) 和
   奇异度 (S) 为独立守恒量配置 Lagrange multipliers。这支持在 BU/S-matrix 后端
   显式记录 μ_B、μ_I/μ_Q、μ_S；但它不提供 PNJL 平衡方程、ρ_Q/ρ_B 闭合或
   `mu_s=0.55 mu_u` 的推导。
5. **适用域再次得到限制。** DMB 的 S-matrix 展开物理图像是稀薄气体；高密度集体
   运动和相变超出其直接适用范围。这进一步支持把当前 strict BU/GBU 保持为诊断/候选
   production 路线，而不是把 S02 当作完整 Ω_M 反馈的依据。

因此，前轮六条核心结论总体不变；需要更新的只有三处措辞：S02 是 S-matrix 与测度
归一化的直接来源，不是 retarded propagator phase 算法的唯一来源；相对论性 Levinson
应标为 conditional；`domega/pi` 与 `domega/(2pi)` 的选择必须先锁定相位变量的
归一化。PV/​i0、ordered unequal-mass charged bubble、KMT (K_{45}/K_{67})、GBU
和 Ω_M 反馈的判断仍分别由其他来源支撑，S02 不改变它们。
