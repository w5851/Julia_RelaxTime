# Andronic 函数族：文献使用与参数包络

研究日期：2026-10-18（Asia/Shanghai）；2026-10-22 阶段审阅后整理进入 PR。状态：`review`、`diagnostic_only`；
作者已要求将三篇国外团队正文与原作者 2026 v2 综述的适用参数纳入新版逐坐标探索域。
作者随后要求与多项式族同构的数值与图像产出；59 组筛选、三组 standard/fine 与两套 PNG 已完成，作者审阅待定。
研究续接对话：`01a144e3-f91e-7230-b41d-04ca2f50ae52`；本次收尾对话来源：`01a14d43-2154-7be2-9f66-72f1d33e52d5`。
默认比较、文献筛选、多项式与 Andronic 族作为同一研究链收纳，
[阶段结论与 PR 证据入口](README.md#2026-10-22-阶段审阅结论)。后续问题待组会讨论后独立确定。

## 结论

Andronic 2018 原文完整参数组的冻结轨迹和 K+/pi+ 已完成 standard/fine 计算，
后续直接复用。此前“先算四条完整原文曲线”的建议没有清楚区分已有基准和新增候选，
在此纠正。现有[多项式探索](parameter_envelope_20261018.md)优化了另一函数族；
本次另完成该族的有限样本探索。最佳成员 `andronic_local_T_lim_GeV_minus` 的 fine 九点 RMSE
为 0.104031，较 2018 默认的 0.154025 降低 32.5%；其峰仍为约 14.5–15 GeV 的宽平台，
没有重现实验约 7–8 GeV 的 horn。此前多项式最优的同口径 RMSE 为 0.069104；
该比较只涉及两批已采样候选，不证明任一函数族的全局极值。

作者先选取 Andronic 2009、Andronic 2018、Yu 2025 Set II（强子与轻核联合）三组中心值，
随后明确要求整理三篇国外团队正文及原作者新综述的参数，并据此更新包络。
当前采用[扩展参数域 v2](../../../../config/physics/freezeout/parameter_domains/andronic_literature_envelope_v2.toml)：
T_lim=153.5-164 MeV、A=2.20-2.82、w=0.40-0.48、D=1050-1356.8 MeV、E=0.286-0.326 GeV^-1。
Yu 2025 Set I（强子、不含轻核）继续排除。
允许跨来源组合，每个候选在所有能量上使用同一套五参数；不由各能量点分别选参数。
该域是本阶段的联合探索范围，不是统计联合置信域。

本次新增计算覆盖这个域内其他完整参数组的轨迹与 ratio，用相同九个 Au+Au 点比较，
并对优选成员核查 standard/fine。原 Andronic 2018 基准直接复用。
独立冻结点资料可作后续解释依据，不作为本轮启动计算的新增前置门槛。

近期国外机构团队仍采用此形式，已核实 2023 年 Elfner–Müller 综述、2024 年 Du 的
流体力学比较，以及 2022 年 Thakur–Hirono 的温度分式应用。四篇的正文参数均已整理到下表：
Du 的 D 下限和原作者 2026 v2 的 A、w 扩大范围；Elfner–Müller 复用 2018，
Thakur–Hirono 只提供温度三参数，均不改变现有极值。

## 统一公式与单位

令 x 为 sqrt(sNN) 按 GeV 计的数值，统一写成

```math
T(x)=\frac{T_{\rm lim}}{1+\exp(A-\ln x/w)},\qquad
\mu_B(x)=\frac{D}{1+Ex}.
```

表内 T_lim、D 使用 MeV，A、w 无量纲；若直接以有量纲的 sqrt(sNN) 写分母，
E 的单位为 GeV^-1。存入项目 GeV 配置时将 T_lim、D 除以 1000，E 数值不变。
指数是 `A-ln(x)/w`，不能改成 `(A-ln(x))/w`。

现有 [Andronic 配置](../../../../config/physics/freezeout/andronic2018.toml)
与该公式一致。A、w 控制温度上升形状；T_lim 是高能极限，
`x_half=exp(A*w)` 为温度达到 T_lim/2 的位置。D、E 控制沿曲线的能量映射，
消去能量后也影响 Andronic 的 T(muB) 几何形状。

## 逐篇中心参数与本阶段选择

| 来源/场景 | T_lim / MeV | A | w | D / MeV | E / GeV^-1 | 本阶段入域 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Andronic 2009 | 164 | 2.60 | 0.45 | 1303 | 0.286 | 是 |
| Andronic 2018 | 158.4 | 2.60 | 0.45 | 1307.5 | 0.288 | 是；复用已有 ratio |
| Yu 2025，Set I，强子且不含轻核 | 161.8 | 0.52 | 0.79 | 1442.2 | 0.329 | 否；作者明确排除 |
| Yu 2025，Set II，强子与轻核联合 | 153.5 | 2.82 | 0.40 | 1356.8 | 0.326 | 是 |
| Elfner–Müller 2023，J. Phys. G | 158.4 | 2.60 | 0.45 | 1307.5 | 0.288 | 是；沿用 2018，非独立新拟合 |
| Du 2024，PRC | 156 | 2.60 | 0.45 | 1050 | 0.288 | 是；流体比较的尺度调整 |
| Thakur–Hirono 2022，JHEP | 158 | 2.60 | 0.45 | 未提供 | 未提供 | 仅温度三参数；该节固定 muB=0 |
| Andronic et al. 2026 v2，新综述 | 156.6 | 2.20 | 0.48 | 1307.5 | 0.288 | 是；原作者低能形状更新 |

完整组共六条来源记录，其中 Elfner–Müller 与 2018 数值相同，共五个不同的完整中心组；
另有一条温度部分来源。取界只使用各来源实际提供的系数，不能把 JHEP 的固定 muB=0
补成 D=0 的完整五参组，也不能把它的缺失 D/E 填成其他来源的数值。

来源定位：

- Andronic 2009，[0901.2909](https://arxiv.org/abs/0901.2909)，
  Acta Phys. Pol. B 40, 1005，出版版 p.1008 Eq.(1) 及前页脚注；
  同源 [0812.1186](https://arxiv.org/abs/0812.1186)，Phys. Lett. B 673, 142，
  Eq.(1) 与 muB 脚注。两篇不能算作两次独立参数提取。
- Andronic 2018，[1710.09425v2](https://arxiv.org/abs/1710.09425v2)，
  Nature 561, 321，[DOI](https://doi.org/10.1038/s41586-018-0491-6)，Fig.3 caption。
- Yu et al.，[2403.16548v3](https://arxiv.org/abs/2403.16548v3)，
  Nucl. Sci. Tech. 36, 65 (2025)，
  [出版信息](https://doi.org/10.1007/s41365-025-01661-z)。
  本次使用开放预印本 PDF pp.6-7 Eq.(5)、Table 4；已逐页渲染核对。
  Set I 包含 pi、K、p、Lambda、Xi 及相应反粒子；Set II 另含 d、反 d、t。

### 报告误差与来源问题

- Andronic 2018 的 T_lim 误差为 1.4 MeV，来自最高五个能量点的拟合。
  同文 ALICE 单能量拟合的 156.5 +/- 1.5 MeV 是不同结果，不能替代此系数。
- Yu Table 4：Set I 给 T_lim=161.8 +/- 3.6 MeV、D=1442.2 +/- 306.5 MeV；
  Set II 给 T_lim=153.5 +/- 0.6 MeV、D=1356.8 +/- 104.3 MeV。
  该文 p.7 正文却把 Set II 温度误差写成 1.4 MeV，与表不一致。
  当前仅采用无歧义的中心组；误差包络必须保留此冲突，不能静默选取。
- Yu Table 4 的 Andronic 参照行写 D=1307、E=0.286，与 Andronic 2018 原文
  的 1307.5、0.288 不同。Andronic 默认值以原文为准。
- 上述位置没有提供可直接使用的完整五参协方差；部分 A、w、E 只给中心值。
  不据此补造独立高斯误差或 68%/95% 联合置信带。
- Andronic 2006 的温度 Eq.(6) 使用不同函数形式，不能把其
  T_lim=161 +/- 4 MeV 当成同一 logistic 形状的一整组参数。
  其 muB Eq.(7) 单独给 D=1303 +/- 120 MeV、E=0.286 +/- 0.049 GeV^-1，
  属于旧拟合的证据，不能当作 2018 新拟合误差。
  见 [nucl-th/0511071](https://arxiv.org/abs/nucl-th/0511071)。

## 使用这套形式的文献多吗

定向检索已核实下列使用实例，2009 两篇同源工作合并列示。
其中只使用温度分式的论文明确标出。这不是穷尽计量，也没有用综述总被引次数
代替公式采用次数。

| 文献 | 实际使用及定位 |
| --- | --- |
| Andronic 2009，上述两篇 | 提出/更新 logistic 温度式，用于粒子产额和 K/pi horn |
| Andronic 2018，Nature | Fig.3 caption 更新高能温度与 muB 映射 |
| [Mukherjee et al. 2017, PRC 96, 025205](https://arxiv.org/html/1611.10144) | Eqs.(16)-(17)，用同形轨迹研究重子数累积量；165、120 MeV 是不同冻结场景，后者接近动理学冻结 |
| [Lafferty & Rothkopf 2020, PRD 101, 056010](https://doi.org/10.1103/PhysRevD.101.056010) | 出版 PDF p.16 Eqs.(64)-(65)，采用 T_lim=158 MeV、A=2.60、w=0.45、D=1307.5 MeV、E=0.288，研究有限密度重夸克偶素；德国海德堡/挪威 Stavanger 团队，年代较早 |
| [Wen & Fu 2021, CPC 45, 044112](https://doi.org/10.1088/1674-1137/abe199) | 出版版 Eqs.(34)-(37)，采用 2018 轨迹研究守恒荷关联，再进行有效理论尺度重标定 |
| [Thakur & Hirono 2022, JHEP 02, 207](https://doi.org/10.1007/JHEP02(2022)207) | 出版 PDF p.18（印刷页码 17）Eq.(4.7)，预印本 Eq.(55)；采用 T_lim=158 MeV、A=2.60、w=0.45 的温度式；该节明确 muB=0，未提供 D/E 映射 |
| [Zhao et al. 2023, CPC 47, 053103](https://doi.org/10.1088/1674-1137/acbf2a) | Eq.(28)，非广延 PNJL 模型中采用 164、2.60、0.45、1303、0.286 |
| [Elfner & Müller 2023, J. Phys. G 50, 103001](https://doi.org/10.1088/1361-6471/ace824) | 出版 PDF p.35 Eq.(37)、Fig.18；对应预印本 Sec.5.1 Eq.(40)，完整列出并使用 2018 五参形式；德国 GSI/Frankfurt 与美国 Duke 团队 |
| [Du 2024, PRC 110, 014904](https://doi.org/10.1103/PhysRevC.110.014904) | Fig.10 讨论，保留形状、改为 T_lim=156 MeV、D=1050 MeV，与流体力学冻结结果比较 |
| Yu et al. 2025，上述来源 | Eq.(5)、Table 4，对不同粒子集合重新拟合同一五参形式 |
| [Andronic et al. 2026, 2606.14799v2](https://arxiv.org/abs/2606.14799v2) | PDF p.9 Eq.(6)，原作者团队更新为 T_lim=156.6 MeV、A=2.20、w=0.48、D=1307.5 MeV、E=0.288；百科综述章节预印本 |
| [Unveiling QCD Criticality with Cross-Rapidity Net-Baryon Cumulants, 2609.01859v1](https://arxiv.org/html/2609.01859v1) | 补充材料 Eq.(S47)，继续使用 2018 默认形式；截至检索日为预印本 |

已核实的多数应用直接沿用固定形状；Yu 提供了明确的形状重拟合，2026 年原作者
综述另给出了低能形状更新。
频繁采用说明可比性和传播范围，不能独立证明某个五维参数区域合理。

### 三篇国外团队正文的参数与使用背景

- **完整形式的综述说明：Elfner–Müller 2023。**
  *The exploration of hot and dense nuclear matter: introduction to relativistic heavy-ion physics*。
  [出版 PDF](../../../../data/inputs/meson_freezeout_andronic_literature_20261018/ElfnerMuller2023_JPhysG50_103001.pdf)
  p.35 Eq.(37)、Fig.18 给出整个 2018 轨迹：158.4、2.60、0.45、1307.5、0.288。
  作者机构为德国 GSI/Frankfurt/FIAS/HFHF 与美国 Duke。
  对应[预印本正文 Eq.(40)](https://arxiv.org/html/2210.12056v1#S5.SS1)。
  该文沿用原参数，没有重新拟合五个参数；入域不改变逐项上下限，也不增加一次独立参数提取。
- **与当前数值研究最接近的应用：Du 2024。**
  *Bulk medium properties of heavy-ion collisions from the beam energy scan with a multistage hydrodynamic model*。
  作者机构为加拿大 McGill。[出版 PDF](../../../../data/inputs/meson_freezeout_andronic_literature_20261018/Du2024_PRC110_014904.pdf)
  p.10（014904-10）Fig.10 讨论中的无编号公式，明确从 2018 的 T_lim=158.4、D=1307.5
  改为 T_lim=156、D=1050 MeV，A=2.60、w=0.45、E=0.288 不变。
  目的为对照热模型与流体冻结，并讨论快度窗口平均。该组按作者本次决定纳入探索域，
  将 D 下限降至 1050 MeV；其性质仍是比较尺度调整，不能称为新的强子产额五参联合拟合或统计误差。
- **温度函数的跨模型应用：Thakur–Hirono 2022。**
  *Spectral functions of heavy quarkonia in a bulk-viscous quark gluon plasma*。
  作者机构为韩国 APCTP/POSTECH、日本 RIKEN。
  [出版 PDF](../../../../data/inputs/meson_freezeout_andronic_literature_20261018/ThakurHirono2022_JHEP02_207.pdf)
  p.18（印刷页码 17）Eq.(4.7)，对应[预印本 Eq.(55)](https://arxiv.org/html/2111.08225v2#S4.SS1)，
  将 T_lim=158 MeV、A=2.60、w=0.45 的温度式用于 psi'/J/psi 的能量依赖。
  紧随正文明确固定 muB=0；没有提供 D/E，故只让温度三参数参与取界，它们均在既有范围内。

这些是领域主流期刊的实际使用例，不能一概称作 Nature/PRL 级的新参数拟合。
本次没有核实到同时满足“近期国外团队、Nature/PRL 级、完整五参新拟合”的新代表。
2018 Nature 仍是当前原始基准的直接来源；2026 预印本不冒充已经发表的高影响力期刊论文。

### 2026 年原作者团队的更新：纳入 v2

*Experimental exploration of the QCD phase diagram*，Andronic、Braun-Munzinger、Redlich、Stachel，
德国 Münster/GSI/Heidelberg、波兰 Wrocław 为主，部分作者兼有 CCNU 机构。
采用[本地开放 PDF](../../../../data/inputs/meson_freezeout_andronic_literature_20261018/Andronic2026_arxiv2606.14799v2.pdf)
的 **v2（2026-09-02）p.9 Eq.(6)**：

```text
T_lim=156.6 MeV, A=2.20, w=0.48, D=1307.5 MeV, E=0.288 GeV^-1
```

作者说明低能形状更新与奇异数正则处理的重新评估有关，并引用 STAR 3 GeV 奇异粒子数据，
而温度平台固定为其 LHC 拟合温度。这里是作者给出的更新参数，不能自行解释成五参联合拟合误差。
本轮下载、提取并渲染核对的是 **v2 的 156.6 +/- 1.7 MeV**；
引用这组更新参数时须同时注明版本与 Eq.(6) 的位置。

A=2.20、w=0.48 超出原三组区间；按作者最新要求纳入 v2，分别更新 A 下限与 w 上限。
T_lim 的 1.7 MeV 误差仅记录，不叠加到中心值包络。
这里使用的是新综述的 v2，不能混用早期版本或检索片段中的 158.4 MeV。
近期国际合作的 cross-rapidity 2026 预印本仍用完整 2018 参数，没有提供新的五参拟合。

## 更新后的联合探索域 v2

在原三组的基础上加入上述四篇的适用中心系数，逐坐标取最小/最大值。
完整组的五参数均参与；JHEP 只贡献温度三参数：

| 参数 | 原三组 v1 | 当前 v2 | 项目配置名 | 当前下限 / 上限来源 |
| --- | ---: | ---: | --- | --- |
| T_lim / MeV | 153.5-164 | 153.5-164 | T_lim_GeV = [0.1535, 0.164] | Yu Set II / Andronic 2009 |
| A | 2.60-2.82 | 2.20-2.82 | log_offset | Andronic 2026 v2 / Yu Set II |
| w | 0.40-0.45 | 0.40-0.48 | log_width | Yu Set II / Andronic 2026 v2 |
| D / MeV | 1303-1356.8 | 1050-1356.8 | d_GeV = [1.050, 1.3568] | Du 2024 / Yu Set II |
| E / GeV^-1 | 0.286-0.326 | 0.286-0.326 | e_GeV_inv1 | Andronic 2009 / Yu Set II |

五个闭区间的笛卡尔积就是本阶段联合探索域。按作者决定允许跨文献配对，
不叠加表列误差，不要求仍然落在原文中心曲线的逐能量上下界内。
每个候选始终是一条由固定五参数定义的完整曲线。

解析代入核对：7.7 GeV 时，原三组盒的范围为 T=130.109-151.594 MeV、
muB=371.204-423.709 MeV；扩展 v2 为 T=123.920-155.472 MeV、muB=299.128-423.709 MeV。
这些是固定能量下的解析坐标范围，不是 ratio 或密度包络。
原三组配置 v1 及其核验记录保持字节不变，原始四组来源清单也保留；当前后续候选以 v2 为准。

以下来源不参与本阶段 v2 取界：

- Yu Set I：作者明确排除；保留文献事实，不参与上下限生成。

- Yu Set III：T_lim=141.7、A=8.00、w=0.21、D=829.7、E=0.184。
  它只拟合轻核；7.7 GeV 缺少足够数据，不能用来支持共同 K/pi 化学冻结线的低温下界。
- Wen & Fu 的整体尺度重标定是有效理论和 QCD 的标尺匹配，不是实验冻结参数误差。
- 从其他函数族的格点涨落冻结温度借用 147/150 MeV 作为 logistic 的 T_lim，
  需要重新约束低能形状；它属于跨方法探索，不能改标为 Andronic 原文置信区间。
- Lafferty–Rothkopf 2020：保留为较早的补充采用例，本次按作者指定的四篇新核实来源更新。
  其参数已处于 v2 范围内。

## Andronic 族同构产出

状态：`review`。采用[固定采样设计](../../../../config/physics/freezeout/parameter_domains/andronic_exploration_v1.toml)，
不改写已确认的 v2 文献域。首轮为四个尚无密度结果的不同完整中心组、32 个分层内部样本、
10 个坐标边界探针，共 46 组；2018 默认参数及 Elfner 的相同组复用已有基准。
JHEP 的温度部分来源参与文献取界，不构造完整五参密度候选。

初筛网格、固定九个 Au+Au 主评分点、RMSE 最优值 1.10 倍的近最优阈值、完整曲线取带、
优选三条曲线的 standard/fine 能量网格及 q0/BQS 数值口径均沿用多项式研究。
首轮后只做一次预定的五轴正负 5% 步长细化，再尝试最多三组域内的保几何能量缩放。
Andronic 保持 T(muB) 的缩放须使用 `k=x_old/x_new, E'=k*E, A'=A-ln(k)/w`，
T_lim、w、D 不变；不能沿用多项式只改 D 的规则。

首轮[设计与参数](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_design/design.json)，
[数值精化计划](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_execution/numerical_plan.json)。
文献域核验时的旧阅读摘要已按其原 hash 保存为
[历史快照索引](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_execution/previous_research_notes.json)。

### 筛选、失败与数值核验

46 组首轮、10 组局部细化、3 组可行能量映射，共 59 组、1003 个 screening 点。
967 点成功，36 点的 K+ 均报告 `ArgumentError: q0 reference onset is not positive real`。
43 组具备完整九点评分，42 组在全部 17 个能量上成功并进入灰带；其中 10 组满足
`RMSE <= 1.10 * sampled_minimum`，进入蓝色斜线带。成员在整个能区固定。
失败点及被排除的完整曲线名单见[最终筛选](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_all_screen_assessment_r2/assessment.json)。

三条优选曲线共完成 standard 66 点、fine 120 点，全部成功、无端点警告，排序一致。
66 个共同有效坐标上的最大相对差为 `5.93004e-5`（0.005930%）。
screening 带仍采用 screening 分辨率，未将其余成员全部重算到 fine；最佳实线采用 fine。
有限窗口、qmax、背景 BQS 和密度内核保持本批约定，standard/fine 比较不证明端点或尾部极限收敛。

| 优选成员 | screening 九点 RMSE | fine 九点 RMSE |
| --- | ---: | ---: |
| `andronic_local_T_lim_GeV_minus` | 0.103764 | 0.104031 |
| `andronic_local_d_GeV_minus` | 0.105177 | 0.105443 |
| `andronic_local_log_offset_plus` | 0.105295 | 0.105561 |

完整中心组的 screening 九点 RMSE 分别为 Andronic 2009 的 0.184728、Yu Set II 的
0.131722、Du 2024 的 0.118745、Andronic 2026 v2 的 0.146592。
这些是当前 PNJL+GBU、BQS 与有限窗口下的诊断分数，不是对来源论文拟合优劣的重新裁决。
保几何映射共接受 3 组、另有 22 次提案越界而排除；接受组未超越局部最优。

自动汇总首次因实验表 hash 不同而停止：旧表与 2026-10-18 勘误表仅有 `NA49_17.27`
中心度 `0-7.2% -> 0-5%` 的元数据差异，所有比值、误差和九点选择相同。
保留旧筛选与失败日志，改用勘误表生成 `all_screen_assessment_r2` 和最终 summary；
评分、排序与带成员不变。详见[恢复记录](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_execution/postprocessing_recovery.json)。

### 最佳参数及两族同口径比较

| 参数 | 本轮最佳值 |
| --- | ---: |
| T_lim / MeV | 154.5086863749 |
| A | 2.7184177421 |
| w | 0.4715028151 |
| D / MeV | 1118.4781730770 |
| E / GeV^-1 | 0.3104982416 |

精确五参数与其他优选成员见[参数表](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_summary/selected_parameters.csv)，
10 组近最优参数见[近最优表](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_summary/near_best_parameters.csv)。

| fine 指标 | Andronic 2018 默认 | 本轮 Andronic 最优 | 此前多项式最优 |
| --- | ---: | ---: | ---: |
| 固定九点 RMSE | 0.154025 | 0.104031 | 0.069104 |
| `R(7.7)` | 0.389881 | 0.263439 | 0.273625 |
| `R(200)` | 0.274506 | 0.260274 | 0.230995 |
| 全曲线采样最大值的能量 / GeV | 9.25 | 15.0 | 7.5 |
| 全曲线采样最大比值 | 0.395486 | 0.307018 | 0.273656 |

本轮 Andronic 的 RMSE 较此前最优默认 FRG mu-independent 改善 18.8%，
但仍比已计算的多项式最优高 50.5%。新汇总重新读取旧多项式 fine 点、复算同一九点
RMSE，核对实验表 hash、数值档位、背景条件及共有物理源码；原多项式结果与图像保持原字节。
此比较没有要求两个参数域具有相同的物理覆盖，也没有证明全局优劣。

Andronic 最优的 14.5、15 GeV 比值分别为 0.307017682、0.307018141，差约 `4.59e-7`，
小于已观测的 standard/fine 最大绝对差。因此报告为约 14.5–15 GeV 的宽平台，
不赋予连续峰位精度。固定 horn 诊断窗口止于 14.5 GeV，summary 的 `sampled_horn`
相应标为边界最大值；上表的 15 GeV 来自 4.87–200 GeV 全部 fine 支点。
实验 7.7 GeV Au+Au 中心值为 0.2227、200 GeV 为 0.159，峰形和高能幅度仍有偏差。
低能 4.87 GeV 不参与九点主评分，最佳曲线该点约为 0.1594，低于观测；未将其隐藏或补入调优权重。

### 两套图与解析投影

[T–muB 彩色图](../../../../data/outputs/figures/relaxtime/meson_density/andronic_envelope_20261018_Tmu_png_review_r1/chemical_freezeout_T_muB_band.png)
· [灰度图](../../../../data/outputs/figures/relaxtime/meson_density/andronic_envelope_20261018_Tmu_png_review_r1/chemical_freezeout_T_muB_band_gray.png)
· [manifest](../../../../data/outputs/figures/relaxtime/meson_density/andronic_envelope_20261018_Tmu_png_review_r1/plot_manifest.json)。

[K+/pi+ 彩色图](../../../../data/outputs/figures/relaxtime/meson_density/andronic_envelope_20261018_png_review_r1/kplus_piplus_parameter_envelope.png)
· [灰度图](../../../../data/outputs/figures/relaxtime/meson_density/andronic_envelope_20261018_png_review_r1/kplus_piplus_parameter_envelope_gray.png)
· [manifest](../../../../data/outputs/figures/relaxtime/meson_density/andronic_envelope_20261018_png_review_r1/plot_manifest.json)。

T–muB 图在固定 muB 处，对 ratio 图同一 42 组完整成员和 10 组近最优成员取解析温度
最小/最大值；近最优参数对应 10 条不同的解析几何曲线。7.7 GeV 圆/方标记分别定位
最佳成员与 2018 默认；它们是参数化坐标。Andronic 消去能量后由
`T_lim,w,D,A+ln(E)/w` 决定几何，能量缩放可以保持该组合。
muB=0–500 MeV 的部分解析延伸超出各组原密度能区，未被当作新密度结果。
156 个既有热力学坐标与解析式的最大差为 `5.69e-14 MeV`，详见
[解析表](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_Tmu_prepared/preparation.json)。

两套图均为 600 dpi PNG，使用与多项式图相同的布局、带语义和原 2018 参照；
原参照在 17.27 GeV 的失败缺口保留。四张彩色/灰度图均实际看图复核，图例、线型、
斜线与空心标记可辨，无裁切或数据遮挡；机器检查通过。图注、源码和审查记录保留在
各自同名 results case。它们仍为 `png_review`，作者审阅 `pending`，未授予论文使用资格。

[完整数值汇总](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_summary/summary.json)
· [分辨率对照](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_summary/resolution_comparison.json)
· [最终核验](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_execution/final_validation.json)。

本轮复用同一 K/pi 资料校准参数，不能称为独立预测验证。包络不具有统计置信度；
没有加入最终态衰变、流场与快度接受度模型。当前产出支持本次有限探索和作者审阅；后续问题与判据待组会讨论后独立确定，当前不安排新增研究。

## 本地证据与复核范围

- [Yu 开放 PDF](../../../../data/inputs/meson_freezeout_andronic_literature_20261018/Yu2025_arxiv2403.16548v3.pdf)
  和 [来源清单](../../../../data/inputs/meson_freezeout_andronic_literature_20261018/source_manifest.json)。
- [本轮近期来源清单](../../../../data/inputs/meson_freezeout_andronic_literature_20261018/recent_sources_manifest.json)
  记录五份开放 PDF 的字节/hash、正文公式位置、INSPIRE 元数据和检索范围，其中四份支撑本次更新、
  一份是此前补充的 Lafferty 2020。
- 当前[扩展参数域 v2](../../../../config/physics/freezeout/parameter_domains/andronic_literature_envelope_v2.toml)
  及[v2 选择与配置核验](andronic_domain_validation_v2_20261018.json)。
- 历史[三组域 v1](../../../../config/physics/freezeout/parameter_domains/andronic_literature_envelope_v1.toml)
  及[当时的核验](andronic_domain_validation_20261018.json)保留，仅证明旧范围，不作为 v2 已通过的依据。
- Yu 的公式、参数表和 Eq.(5) 括号已进行 PDF 页级视觉复核；原表与正文误差不一致保留。
- 原四组及三组域的解析记录保留；v2 从所选完整/部分来源分别取界并重新复核解析坐标。
- 2026 v2 p.9 的 Eq.(6)、邻近参数文字和 Du 出版 PDF p.10 已渲染核对；
  Elfner 出版 PDF p.35 Eq.(37) 与 JHEP 出版 PDF p.18 Eq.(4.7) 及 muB=0 的邻近正文也已渲染核对。
  Lafferty 出版 PDF 的公式已提取核对。
- v2 的系数取界、单位转换、既有 2018 profile、角点/中心组解析坐标、PDF 完整性、
  本地引用与仓库文档检查的实际结果见 v2 核验记录；历史 315 点结果仅属于 v1。
- 文献核验没有改变密度核、回归基线或正式产物；后续已授权的数值探索状态见上节。
