# 化学冻结线的跨文献参数包络

状态：`review`，`diagnostic_only`；数值与 PNG 审阅产物已完成，作者审阅待定。
2026-10-22 阶段审阅将默认比较、文献筛选和两族探索作为同一研究链整理进入 PR，
[阶段结论与证据入口](README.md#2026-10-22-阶段审阅结论)。
作者于 2026-10-18 同意开始带状化，并允许不同
文献之间组合参数。本轮先用文献区间的逐坐标闭包构造探索区域，再比较完整曲线。
这不是任何一篇文献的联合置信域，也不是物理上冻结温度的空间分布。

59 组参数初筛后，三组候选完成 standard/fine 核验。样本最优 `hornmap_001`
的九点 RMSE 为 0.069104，比同能区的 Andronic 2018 降低 55.1%，比此前最优
默认 FRG 2024 mu-independent 降低 46.1%。其加密采样峰位为 7.5 GeV，
但高能幅度与 horn 的尖锐程度仍未匹配实验，且混合后的 muB 映射需要独立物理约束。

[彩色对照图](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_envelope_20261018_png_review_r2/kplus_piplus_parameter_envelope.png)
· [灰度图](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_envelope_20261018_png_review_r2/kplus_piplus_parameter_envelope_gray.png)
· [图像 manifest](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_envelope_20261018_png_review_r2/plot_manifest.json)
· [完整指标](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_summary/summary.json)。

对应的 [T–muB 参数带图](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_envelope_20261018_Tmu_png_review_r2/chemical_freezeout_T_muB_band.png)
· [灰度图](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_envelope_20261018_Tmu_png_review_r2/chemical_freezeout_T_muB_band_gray.png)
· [manifest](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_envelope_20261018_Tmu_png_review_r2/plot_manifest.json)。

## T–muB 图与优化函数族的范围

本轮“最优混合参数”仅在多项式族
`T=a-b*muB^2-c*muB^4, muB=d/(1+e*sqrt(sNN))` 中得到。
**本批多项式图中的 Andronic 2018 使用原文默认参数；本批没有对其 logistic 温度函数族进行参数带扫描
或寻找最优组合。** 因此这批结果比较的是调优后的多项式与默认 Andronic，
不能据此判断两个函数族分别优化后谁更优。作者随后要求的
[Andronic 族同构研究](andronic_family_20261018.md)单独维护，原多项式数值和图像保持原样。

T–muB 图在固定 muB 处对相同的 53 组参数取温度最小/最大值，灰色为完整探索带，
蓝色斜线为原来四组近最优参数的带；蓝色实线为 `hornmap_001`，橙色虚线为
Andronic 默认。没有重新筛选成员，也没有把 ratio 的逐点上下边界直接变换成
某两条冻结线。两种参数带均不具有统计置信度。

多项式的几何形状只依赖 a、b、c；d、e 决定能量沿曲线的位置。
`local_a_minus` 与 `hornmap_001` 因前三个系数相同而完全重合，所以四组近最优
参数对应三条不同的 T(muB) 曲线，却仍能给出四条不同的 ratio(sqrt(sNN)) 曲线。
图中圆/方标记为各自参数化在 7.7 GeV 的位置：最佳多项式为
`(muB,T)=(235.481,138.217) MeV`，Andronic 为 `(406.359,138.428) MeV`。
标记是参数化坐标，不是新的实验冻结点。

显示范围为 `muB=0–500 MeV`，直接代入冻结参数求值；Andronic 通过
`sqrt(sNN)=(d/muB-1)/e` 消去能量，在 muB=0 使用解析极限 T_lim。
其中部分线段超出各组原先 4.87–200 GeV 的密度计算能区，属于解析延伸，
不代表这些点已有 q0 密度结果。Andronic 的原 K+ 积分失败也不改变解析冻结线
本身的定义；旧 ratio 图的缺口与失败记录保持原样。

参数表先单独冻结，再由绘图阶段读取。新的 4 项投影测试通过，156 个已有热力学
坐标与解析式的最大差为 `5.69e-14 MeV`。图例位于同一子图内，彩色/灰度视觉与
manifest 检查通过。派生表与方法见
[preparation.json](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_Tmu_prepared/preparation.json)，
图像核验见
[review.json](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_Tmu_png_review_r2/review.json)。
此图没有新调用密度求解器；全部原始密度、ratio、选择与评分保持不变。

## 文献约束与构造选择

统一记 `x=sqrt(sNN)`（GeV），`muB=d/(1+e*x)`，`T=a-b*muB^2-c*muB^4`。
T、muB 为 GeV，系数单位依次为 GeV、GeV^-1、GeV^-3、GeV、GeV^-1。
归一化四次式严格换算为 `a=T0,b=kappa2/T0,c=kappa4/T0^3`。

原文数值、误差、定位和 PDF hash 保存在
[约束配置](../../../../config/physics/freezeout/parameter_domains/literature_envelope_v1.toml)。
本轮复核了以下来源：

- Cleymans 2006：[hep-ph/0511094v2](https://arxiv.org/abs/hep-ph/0511094v2)，
  PDF p.3 Eqs.(1)-(2)，五参中心值及逐项误差。
- Poberezhnyuk 2019：[1906.01954v2](https://arxiv.org/abs/1906.01954v2)，
  PDF p.9 Table II，IHRG 与 QvdW 两行均保留；它们代表不同理论提取设置。
- Shao 2024：[PRD 109, 034019](https://doi.org/10.1103/PhysRevD.109.034019)，
  出版版 p.8 Table IV，两场景的中心值。源文献未施加奇异数中性，7.7 GeV 未拟合。
- Lu 2026：[2603.09336v1](https://arxiv.org/abs/2603.09336v1)，
  PDF p.3 Table I，mean、lower、upper 三个成组参数。原文的组间相关性记录保留；
  跨组拆分是本轮经作者同意的额外放宽，不重新标作文献置信带。
- Bazavov 2016：[1509.05786v2](https://arxiv.org/abs/1509.05786v2)，PDF p.4 Table I：
  STAR0.8、STAR2.0、PHENIX/STAR2.0 分别给 `Tf0=145(2),147(2),155(4) MeV`。
  这确实是从高能涨落提取的零化学势冻结温度，而非 crossover 温度。它们仅补充 a
  的范围，没有提供完整的 c、d、e。源文献讨论的接受度及净质子/净重子转换系统误差
  仍然适用；不同列不是三次独立测量。负曲率估计未移入现有 `b,c>=0` 的单调模型。
- Borsanyi 2014：[1403.4576v2](https://arxiv.org/abs/1403.4576v2)，PDF p.4 Fig.4 与正文：
  `Tf=144+-6 MeV` 是 `sqrt(s)>=39 GeV` 的共同轻味涨落冻结温度范围，不能直接当成
  a 的 `138-150 MeV` 区间。此项用于解释、未扩展 a 的下界。

这些来源也不保证 K 与 pi 在共同温度冻结。本轮将其用于共同介子冻结轨迹的作用，
属于模型探索。DSE 的任意有限尺寸低温分支、味分离温度和纯 crossover 结果未混入
五参上下限。检索沿用已有定向文献池，并非穷尽全部冻结线文献。

| 参数 | 文献区间的逐坐标包络 | 首轮内部采样区间 |
| --- | ---: | ---: |
| a / GeV | 0.143–0.168 | 0.1455–0.1655 |
| b / GeV^-1 | 0–0.155 | 0.0155–0.1395 |
| c / GeV^-3 | 0.032–1.09734738 | 0.13853474–0.99081264 |
| d / GeV | 0.70896–1.336 | 0.771664–1.273296 |
| e / GeV^-1 | 0.147–0.281 | 0.1604–0.2676 |

首轮在每个区间两端各留 10% 余量；这是研究设计选择，不是文献给出的概率区间。
包络采用固定 GeV 单位的 `(a,b,c,d,e)` 坐标；若改用 `(T0,kappa2,kappa4,d,e)`
重新取逐坐标包络，会得到不同区域。本轮放宽并不声称这种构造具有坐标不变性。
保留 8 个原文参数组作为参照，加入 6 个 `a=147/150 MeV` 的明确混合场景与
32 个确定性分层提案接受样本。种子为 20261018。只要求当前模型的符号条件以及
4.87–200 GeV 上温度为正；不同来源的相关性不会被偷偷当成独立高斯误差。

## 计算与评价合同

已冻结的[初始设计](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_design/design.json)
包含 46 条曲线、全部配置、来源 hash 和被排除的提案。采样/分析入口为
[freezeout_parameter_envelope.py](../../../../scripts/analysis/relaxtime/freezeout_parameter_envelope.py)。
每组五参在整个能区固定，不逐能量重新拟合。曲线温度与化学势分别与原文公式核验。

数值入口仍为 [q0 比较工作流](../../../guides/scripts/meson_freezeout_comparison.md)，
使用 `--parameter-policy literature_envelope_exploration` 标识此次参数探索。
quark-only BQS、`nS=0,nQ/nB=0.4,xi=0`、外部 omega Bose 权重、有限积分窗口及
密度核不改。有限窗口积分与原有数值失败均保留；结果没有加入衰变或实验接受度。

首轮在 17 个能量直接计算，主评分事先固定为 7.7–200 GeV 的 9 个 Au+Au 观测点。
低能 4.87 GeV、Pb+Pb 的 horn 点及更密的理论支点均不增加主评分权重；另报此前
11.5–200 GeV 的 8 点指标。这个新的 9 点结果不能直接与旧 14 点 RMSE 比大小。
候选缺失任一主评分点即不参与排名，不通过缩小公共点集改善分数。

首次筛选后只进行一轮预设的轴向局部探索：每个参数正负移动完整包络宽度的 5%，
仍受完整包络及正温度限制。此步骤可以超出初始内部采样区间，因此也将单独标记。
随后对较优参数作 standard/fine 复核并增加 horn 附近的能量支点。
加密最多选择三条：首先保留九点 RMSE 最优者，再保留 RMSE 在其 1.10 倍内且
4.87→7.7→11.5 GeV 先升后降的最佳三点形状候选，其余名额按九点 RMSE 补足。
这一三点形状诊断与密集网格上的峰位分别报告，不将其当作实验峰位精确测量。
首轮发现幅度最优样本的采样峰位后移至11.5 GeV，因此增加一个同范围的能量映射
诊断：固定 `(a,b,c)`，用原文中心值 `e=0.273 GeV^-1`，令
`d=muB_at_sampled_peak*(1+e*7.7)`，将该热力学点映射到7.7 GeV。只接收仍在原包络
内的组合，并直接重新计算q0密度，不将重标能量或插值当作数值结果。由于峰位最初
来自有限网格，加密后的实际峰位仍需核查。这是在看到峰位偏移后增加的诊断，
不改主评分或参数范围，也不把拟合后的结果称为独立预测。
“最优”仅表示已计算集合中的最优值，不保证五维连续区域的全局最优。

绘制全部完整曲线的采样包络，以及主 RMSE 不超过样本最小值 1.10 倍的集合。
后者是预先声明的描述性筛选，不能称为 68%/95% 置信区间。两类包络均保持一个
固定的完整曲线成员集合；缺失点不通过更换逐点成员来掩盖。上下边缘可由不同
参数组给出，逐点覆盖实验不等于某一参数组同时描述整个实验能区。

## 同一观测集合上的结果

下表全部使用同一九个 Au+Au 点的等权 RMSE；旧方案从各自冻结 fine CSV 重新
评分，没有拿历史十四点或八点分数直接比较。R 为 K+/pi+，无量纲。

| 方案 | 九点 RMSE | R(7.7) | R(200) |
| --- | ---: | ---: | ---: |
| Cleymans 2006 默认 | 0.182908 | 0.415393 | 0.298910 |
| Andronic 2018 默认 | 0.154025 | 0.389881 | 0.274506 |
| FRG 2024 mu-independent 默认 | 0.128127 | 0.341983 | 0.260447 |
| local_a_minus | 0.073414 | 0.264836 | 0.232425 |
| hornmap_003 | 0.075232 | 0.280740 | 0.235143 |
| **hornmap_001** | **0.069104** | **0.273625** | **0.230995** |
| 实验中心值 | — | 0.222700 | 0.159000 |

FRG 在 7.7 GeV 的行内结果是本项目对其公式的延伸，原文没有拟合该点。
IHRG、另一 FRG 分支和 functional QCD 的相同九点结果也保存在完整指标中。
最佳样本在九个评分点均偏高；7.7 与 200 GeV 分别高出实验中心值约 22.9% 和
45.3%。不能将 RMSE 的降低解释为已经与数据相符。

加密后的三组候选峰位分别为 `hornmap_001: 7.5`、`hornmap_003: 8.0` 和
`local_a_minus: 11.0 GeV`。最优样本的峰高为 0.273656，直接相邻支点为
7.0 与 7.62 GeV；这是离散采样峰，不是连续峰位的置信区间。
Andronic 默认的原加密网格峰位为 9.25 GeV，峰高 0.395486。

峰位接近尚未解决 horn 形状。最优曲线在 4.87→7.7→11.5 GeV 的比值变化为
`+0.01012,-0.00519`，实验中心值为 `+0.02070,-0.02090`，尤其峰后下降仍然过缓。
`hornmap_003` 的三点差分形状误差略小，但整体幅度评分较差，因此保留为形状
候选，没有替换主 RMSE 最优者。三点差分只作描述性诊断，未构造带实验协方差的拟合。

## 参数组合及物理解释边界

最优五参为：

| a / GeV | b / GeV^-1 | c / GeV^-3 | d / GeV | e / GeV^-1 |
| ---: | ---: | ---: | ---: | ---: |
| 0.1456442197 | 0.1036414873 | 0.5463741671 | 0.7304851723 | 0.2730000000 |

显示的小数位用于复算，不表示拟合精度；完整值见
[选中参数表](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_summary/selected_parameters.csv)。
`local_a_minus` 与 `hornmap_001` 的 a、b、c 相同，能量映射不同；初筛时前者
11.5 GeV 与后者 7.7 GeV 对应相同热力学点，直接计算的比值也相同。
这支持映射可以移动峰位，但不会单独制造更尖锐的热力学峰。

降温与映射的作用需要分开陈述。最优组在 200 GeV 给出 `T=145.626 MeV`，
Andronic 为 158.384 MeV；而在 7.7 GeV，两者温度近似相同（138.217 与
138.428 MeV），muB 却分别为 235.481 与 406.359 MeV。
因此不能把整个改善都归因于降低零化学势温度。

逐坐标包络允许产生超出所有原始曲线的轨迹。本轮纳入五参池的八组原文中心参数
在 7.7 GeV 给出 `muB=332.548–495.269 MeV`，最优混合组的 235.481 MeV 已在
该范围之外；在 200 GeV 也出现类似情况。相关数值与来源 hash 见
[轨迹检查](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_execution/trajectory_review.json)。
这符合本轮明确放宽参数相关性的规则，却不等于文献联合支持这条轨迹。
这里的原文中心曲线范围仅作诊断，不是新增排除门槛；低能外推也不能当作该能量
的独立实验约束。本轮没有据此改动评分或隐藏最优样本。

因此当前建议保留“最优曲线 + 探索带”的呈现，并说明其混合参数性质。
当前没有独立的 muB(s) 或其他粒子产额约束，因此不能将混合最优轨迹称为文献共同支持的
论文采用冻结线；只针对同一 K+/pi+ 数据选参仍属校准。后续问题待组会讨论后独立确定。
局域密度比与中快度最终产额的差异仍涉及冻结面、衰变和接受度。仅凭 y≈0 不能
断言实验比值应当更低，原讨论保留在[默认参数报告](README.md#快度解释与后续参数带)。

## 包络、批次与失败记录

初始 46 组、轴向邻域 10 组与能量映射 3 组合计 59 组，各有 17 个直接计算支点。
1003 个 screening 点中 991 个成功，245 个点记录端点警告。12 个失败均保留
`q0 reference onset is not positive real` 的 K+ 门禁原因；背景求解本身没有失败。
五组缺少主评分点而不参与排名，六组曲线不完整而不进入整条曲线包络：
`envelope_002,004,017,023,024,027`。其中 004 仅缺低能辅助支点，因此主排名
与完整曲线成员数量不同。详情见
[初筛评估](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_all_screen_assessment/assessment.json)。

灰色带采用固定的 53 条完整初筛曲线。蓝色斜线带采用四条近最优完整曲线：
`hornmap_001,local_a_minus,hornmap_002,hornmap_003`；其五参及初筛分数见
[近最优参数表](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_summary/near_best_parameters.csv)。
这四个组合不是允许任意再混合的独立参数区间，更不是统计置信带。
上下边界仅在 17 个直接支点上取极值，绘制直线连接；没有生成额外插值点或平滑曲线。

选中的三条曲线完成 66 个 standard 点和 120 个 fine 点，全部成功且端点警告为零。
本轮共 1189 条数值计算记录，1177 条成功；这个数量包含各档位的重复坐标。
原始批次为
[standard](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_standard/manifest.json)、
[fine_0](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_fine_0/manifest.json) 和
[fine_1](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_fine_1/manifest.json)。
完整初筛带没有统一重算 fine，图上的实线最优曲线与阴影带的分辨率口径分别保留。

初筛时的相对 profile 路径快照错误与加密驱动的入口语法错误，均在生成数值点前
修正。旧目录、计划、源码快照与失败日志保留；本次成功加密使用
[numerical_plan_r2.json](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_execution/numerical_plan_r2.json)。
图像 r1 缺少三条 series 的支撑/缺失值元数据，未通过 manifest 校验；r2 补齐
元数据并通过，彩色/灰度 PNG 字节与已检查的 r1 相同。没有覆盖中间记录。

## 验证与交付范围

66 个共同有效坐标的 standard/fine 最大相对差为 `4.85456e-5`，即 0.004855%，
以 standard 为分母；最大绝对比值差为 `1.07162e-5`。三组候选排序保持不变，
没有失败状态改变。细节见
[分辨率比较](../../../../data/outputs/results/relaxtime/meson_density/freezeout_envelope_20261018_summary/resolution_comparison.json)。
独立复算九点评分与密度商一致，fine 背景残差最大 `8.322e-10`，小于原 `1e-7`
门槛。与旧默认比较的共同物理源码 hash 不变；runner 差异限于参数来源标识、
相对路径快照修复和已显式传能量的 CLI 默认表。没有修改密度核或放宽容差。
这只支持固定窗口内的分辨率稳健性，不替代 qmax、omega 下限或完整有限 q 验收。

Julia 定向集成 71 个断言通过，包含真实 q0 运行、恢复身份、BQS/频率合同和相对
profile 快照；Python 两组共 17 项通过，覆盖参数换算、范围、映射、固定观测点、
整曲线缺失与包络、候选选择和旧比较语义。此次没有改变密度核，复用此前的核
回归证据，没有为参数扫描重复运行全套测试。当前检查记录见
[验证记录](parameter_envelope_validation_20261018.json)。

新版实验 CSV 仅修订 `NA49_17.27.centrality: 0-7.2% → 0-5%`，22 行其余字段
逐项相同，所有数值不变；来源与父版本 hash 见
[2026-10-18 实验输入](../../../../data/inputs/meson_freezeout_experiment_20261018/README.md)。
旧实验表、默认计算和图像的冻结 hash 均保留；没有将这次勘误扩大为全部元数据审定。

最终 PNG 使用同一子图内图例，彩色与灰度均完成代理视觉检查，manifest 校验通过。
图例不遮挡曲线或实验点，保留 Andronic 17.27 GeV 的失败缺口。图宽为 6.75 英寸、
600 dpi，建议不缩小到 6.296 英寸以下；不适合直接缩为单栏。
交付状态为 `png_review`，作者审阅待定，没有生成正式矢量图或写入论文图层。

原始计算记录仍绑定 `D:/w/meson-freezeout-defaults` 与基线
`698feaf8fada670fde81c93bab6093b50093935f`。本次 PR 在独立工作树
`D:/w/meson-freezeout-review`、`codex/meson-freezeout-review` 整理；
[证据恢复说明](../../../../data/outputs/results/relaxtime/meson_density/freezeout_review_20261022/README.md)
保留原始点、源码快照和全部失败记录，原工作树的其他研究保持原状。
