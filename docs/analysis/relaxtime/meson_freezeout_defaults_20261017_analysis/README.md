# 化学冻结线与 K+/pi+：研究链阶段审阅

研究阶段：`review`、`diagnostic_only`；默认比较、文献筛选、两族探索和 standard/fine、PNG 已完成。
2026-10-22 收尾审阅：将上述阶段作为同一研究链整理进入 PR，图像作者审阅仍为 `pending`。
本页保留默认参数阶段的证据和评分口径。Ledger item：
`meson-freezeout-default-comparison`。

同一研究链依次包含默认比较、[近期文献筛选](literature_candidates_20261017.md)、
[多项式五参数探索](parameter_envelope_20261018.md)和 [Andronic 族探索](andronic_family_20261018.md)。
两族统一使用固定九个 Au+Au 点评价；本页的历史十四点、文献筛选的八点指标各自保留。

作者随后要求将三篇国外团队正文和原作者 2026 v2 综述纳入
[Andronic 参数探索域 v2](andronic_family_20261018.md)。原 Andronic 2009、2018 与 Yu Set II 保留，
新增 Du 2024 与 Andronic 2026 的适用参数；Elfner–Müller 沿用 2018，JHEP 2022 仅贡献温度三参数。
当前 T_lim=153.5-164 MeV、A=2.20-2.82、w=0.40-0.48、D=1050-1356.8 MeV、E=0.286-0.326 GeV^-1，
继续排除 Yu Set I。Andronic 同构产出已完成：59 组筛选、42 组完整带、10 组近最优带，
三组 standard/fine 与 T–muB、K+/pi+ 两套彩色/灰度 PNG。最佳九点 RMSE 为 0.104031，
较 2018 默认降低 32.5%，宽峰仍在约 14.5–15 GeV；此前多项式同口径最优为 0.069104。
方法、精确参数及全部图像入口见上述 Andronic 研究页；两族当前均为 PNG 审阅阶段。

## 2026-10-22 阶段审阅结论

| 固定九点 fine 指标 | Andronic 2018 默认 | 多项式最优 `hornmap_001` | Andronic 最优 |
| --- | ---: | ---: | ---: |
| RMSE | 0.154025 | 0.069104 | 0.104031 |
| 直接采样最大值能量 / GeV | 9.25 | 7.5 | 约 14.5–15 的宽平台 |

两族均覆盖 59 组候选；多项式有 53 组完整曲线、4 组近最优，Andronic 有 42 组完整曲线、
10 组近最优。每族三条优选完成 66 个 standard 与 120 个 fine 点，排序稳定；共同坐标
最大相对差分别为 0.004855% 和 0.005930%。筛选中的失败记录和不完整曲线保留。

这支持“本批有限采样的多项式候选匹配较好”的阶段判断，不证明函数族全局优劣。
多项式峰位接近实验，但幅度和尖锐程度仍有偏差，最优混合 muB(s) 超出原始中心曲线范围；
Andronic 的宽平台后移。两族均未完整重现实验 horn。采样带不是统计置信带，同一 K/pi
数据用于选参后不能再称独立预测验证；standard/fine 仅验证当前有限窗口内的分辨率稳健性。

当前 PR 保留实现、参数和冻结证据，图像仍为 `png_review`、作者审阅待定。
[证据恢复说明](../../../../data/outputs/results/relaxtime/meson_density/freezeout_review_20261022/README.md)
定位数值及检索 ZIP、原始逐点记录和源码快照，
[本次 PR 验证](pr_validation_20261022.json) 记录在独立工作树的复核。
后续物理约束、连续优化或论文采用问题待组会讨论后独立确定，当前不预设后续任务。

## 结论与选择范围

**初轮两条完整能区方案中，Andronic 2018 较优。** 在加密扫描中，14 个 Au+Au
点上的等权 RMSE 从 Cleymans 的 0.185087 降至 0.143604，改善约 22.4%；标准档
得到相同排序。两条线仍明显高于数据，Andronic 的直接采样峰位为 9.25 GeV，
Cleymans 为 9.5 GeV；均未匹配实验约 7–8 GeV 的 horn。选择较优默认方案不等于
解决 horn 形状、峰位和绝对幅度的差异。

后续问答中的[近期文献与默认参数补筛](literature_candidates_20261017.md)已找到
IHRG 2019、FRG 2024 和 functional QCD 2026 的更好候选，并完成标准/加密核验。
其中 FRG 的 mu-independent 分支在共同 8 个 Au+Au 点上降低 RMSE 约 15.2%，
但源文献的守恒荷条件、低能拟合范围及新曲线 horn 形状仍需分别处理。8 点结果
不能替代本页历史 14 点指标；后续参数包络已单独按固定九点评价。

[最终彩色对照图](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_defaults_20261017_fine_png_review_r2/kplus_piplus_freezeout_defaults.png)
· [灰度图](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_defaults_20261017_fine_png_review_r2/kplus_piplus_freezeout_defaults_gray.png)
· [图像 manifest](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_defaults_20261017_fine_png_review_r2/plot_manifest.json)
· [完整比较指标](../../../../data/outputs/results/relaxtime/meson_density/freezeout_defaults_20261017_fine_png_review_r2/comparison_metrics.json)。

本轮比较两条同时给出 `Tf(sqrt(sNN))` 和 `muB(sqrt(sNN))` 的完整化学冻结线。
先前文献筛选中的 HotQCD 2019、Borsányi 2020 给出 crossover，遵照作者澄清不纳入
本轮化学冻结线候选。涨落冻结分析给出的约 147–150 MeV 温度、曲率界限仍可约束
后续研究，但并不单独给出覆盖 horn 至 LHC 的完整默认能量参数化；本轮没有自行
拼接它们或补造束流能量映射。

因此，本轮“相对较优”严格指 Cleymans 与 Andronic 这两条已计算的完整经验线。
Bazavov 2016 的低温化学冻结分支尚未计算，不能由当前排序判定其优劣；它仍可在
有文献支持的高能区作为后续比较对象。

## 文献公式与实现

以 `z=sqrt(sNN)/GeV` 为碰撞能量的无量纲数值，以下公式中的 `T,muB` 使用 GeV。

| 方案 | 温度 | 化学势 |
| --- | --- | --- |
| Cleymans 2006 | `T=0.166-0.139*muB^2-0.053*muB^4` | `muB=1.308/(1+0.273*z)` |
| Andronic 2018 | `T=0.1584/(1+exp(2.60-log(z)/0.45))` | `muB=1.3075/(1+0.288*z)` |

来源：Cleymans et al., PRC 73, 034905 (2006),
[hep-ph/0511094v2](https://arxiv.org/abs/hep-ph/0511094)，PDF p.3 Eqs.(1)–(2)；
Andronic et al., Nature 561, 321 (2018),
[1710.09425v2](https://arxiv.org/abs/1710.09425)，PDF p.6 Fig.3 图注。
Andronic 指数的除法只作用于 `log(z)`；158.4 MeV 是跨能量参数化的极限温度，
不同于同文 ALICE 单点拟合的 156.5 MeV。该参数化不是 Au+Au 专属拟合。

冻结线由 [FreezeoutProfiles.jl](../../../../src/models/scans/FreezeoutProfiles.jl)
按类型分派，多项式旧构造器和数值保留；Andronic 使用独立的
[TOML](../../../../config/physics/freezeout/andronic2018.toml)。未知 family、缺少
新 family 系数、坏 TOML 或不存在的命名 profile 明确报错，不静默使用另一条线。
完整 API 与单位见 [FreezeoutScan.md](../../../api/models/scans/FreezeoutScan.md)。

## 计算口径

每个冻结点直接求解现有 quark-only `FixedMuBConservedCharges` 背景，然后计算
`pi_plus`、`K_plus` 密度。`nS=0`、`nQ/nB=0.4`、`xi=0`，无介子反馈；背景节点 48、
残差上限 `1e-7`。背景沿能量降序延续，各参数化重置 seed。

使用 PR #325 已有的 `q0_lambda_reference`：
`lambda=omega+(mu1-mu2)`，`lambda0=sqrt(lambda^2-q^2)` 只在 `lambda>=q` 时使用，
其他区域相位为零；Bose 权重使用外部 `omega`。这是有限 Bose 窗口，
`omega_lower=1e-5 fm^-1`，`qmax=8 fm^-1`。没有使用旧等高网格插值。

| 档位 | 能量数/方案 | q nodes | mesh | cut/tail nodes | omega nodes |
| --- | ---: | ---: | ---: | ---: | ---: |
| standard | 48 | 16 | 128 | 64/64 | 96 |
| fine | 50 | 32 | 256 | 96/96 | 128 |

fine 除提高积分分辨率，还在 9.5 GeV 两侧直接增加 9.25、9.75 GeV。
完整能量、物理配置、运行身份、Julia 版本和源码 hash 在各批次 `run.json`；
实际源码与配置留在 `source_snapshot/`。`manifest.json` 记录每个点和 CSV 的 hash。
失败点保留 NaN/null，不填零，不跨缺口连接。

入口：[运行与绘图说明](../../../guides/scripts/meson_freezeout_comparison.md)。
`Models.run_charged_gbu_freezeout_scan` 的直接有限 q 生产默认、密度核和 BQS 求解器
均未由本轮修改。本轮结果只能支持默认方案筛选及相同有限窗口内的分辨率核验，
不声称零端点极限、尾部收敛或完整生产验收。

## 实验输入与评分

[输入表及来源](../../../../data/inputs/meson_freezeout_experiment_20261017/README.md)
保留用户投影片所引综述的实验组成，并优先恢复原始发布表：

- 14 个主比较 Au+Au 点：AGS 5 点、STAR BES 5 点、14.5 GeV 1 点及高能 STAR 3 点。
- NA49 Pb+Pb 5 点单独显示，作为 horn 辅助；不并入 Au+Au 主 RMSE。
- ALICE 2 点是 `(K+ + K-)/(pi+ + pi-)`，只作高能辅助展示。
- HADES 2.4 GeV 是读图点，只展示；像素读图范围不当作实验误差棒。

主指标是实验能量处的等权 RMSE，所有值来自同能量的直接计算。模型曲线加密不改变
实验权重；失败时两种方案使用同一公共有效集合，并记录被排除观测。误差棒保留
各表的不对称误差，合并可用的统计与系统误差。STAR 14.5 GeV 的 `26.7/141` 为
发布产额的派生比值，误差传播明示零协方差假设。评分是描述性方案选择，未把共享
系统误差当作独立误差做正式 chi-square 拟合或参数置信区间。

重要来源核对：STAR 2009 的 HEPData Fig.30a 实际只有负电荷比值，改用原论文
Table IX 正电荷列；NA49 Table 10 的低能 AGS 点实际为 Au+Au，不能照搬表的全局
Pb+Pb qualifier。表中 2.35 GeV 维持发布值，不擅自换为 2.7 GeV。

问答复核补充：发现 CSV 将 `NA49_17.27` 的中心度误写为 0–7.2%，原始 NA49 2002
Table I 为 0–5%；旧 SPS 点的快度字段也需要按原始实验逐能量回溯。详见
[实验输入勘误](../../../../data/inputs/meson_freezeout_experiment_20261017/README.md#元数据勘误与核验边界2026-10-17-问答复核)。
发布比值、已绘图数值与 Au+Au 主评分不因该中心度勘误改变。冻结输入保留原字节，
当前报告不将“22 行可复现”扩大为“全部实验元数据已经逐项审定”。

## 加密结果与失败记录

| 指标 | Cleymans 2006 | Andronic 2018 |
| --- | ---: | ---: |
| Au+Au 14 点 RMSE | 0.185087 | 0.143604 |
| Au+Au 4–20 GeV RMSE | 0.198035 | 0.161149 |
| Au+Au >=19.6 GeV RMSE | 0.168749 | 0.138705 |
| `R(7.7)`；STAR = 0.2227 | 0.415393 | 0.389881 |
| `R(200)`；STAR = 0.159 | 0.298910 | 0.274506 |
| 采样峰能量 (GeV) | 9.5 | 9.25 |
| 采样峰高 | 0.419792 | 0.395486 |

standard 共 96 点，95 点成功。Andronic 在 17.27 GeV 的 K+ 计算触发
`ArgumentError: q0 reference onset is not positive real`，发生在 q shell 积分前；
该点 π+ 成功。`T=154.6924 MeV`、`muB=218.8739 MeV`。原失败点及原因保留在
`points/andronic2018/energy_019.json`。它属于 NA49 辅助观测，不影响 14 点 Au+Au
主评分；两方案的 NA49 辅助评分使用共同有效的 4 点。没有放宽门限或补线。

fine 共 100 点，99 点成功；同一个 17.27 GeV 点重复相同失败。加密本身没有解决
该门禁问题，也不能据此判定某种相变或凝聚。最终图在此断线。
standard/fine 的端点警告点数分别为 56/61；有限窗口保留了这些诊断信息，
没有将警告重新标为完整端点可积性验收。

95 个共同有效坐标上的比值差异如下，最大差异均在 2.35 GeV：

| 方案 | 共同有效点 | 最大绝对差 | 最大相对差（以 standard 为分母） |
| --- | ---: | ---: | ---: |
| Cleymans | 48 | `2.9184e-5` | 0.01413% |
| Andronic | 47 | `2.3158e-5` | 0.01580% |

详见 [分辨率比较](../../../../data/outputs/results/relaxtime/meson_density/freezeout_defaults_20261017_fine_png_review_r2/resolution_comparison.json)。
差异远小于当前模型与实验的幅度差，但只检验已声明的采样加密；qmax 和有限窗口下限
未变。Andronic 峰旁直接点为 9.0、9.5 GeV，Cleymans 为 9.25、9.75 GeV；这些是
采样邻点，不是连续峰位的统计置信区间。NA49 的发布点峰在 7.62 GeV、比值 0.253；
STAR Au+Au 的本次发布点峰在 7.7 GeV、比值 0.2227。两种实验体系分别保留。

原始批次：[standard](../../../../data/outputs/results/relaxtime/meson_density/freezeout_defaults_20261017_standard/manifest.json)
和 [fine](../../../../data/outputs/results/relaxtime/meson_density/freezeout_defaults_20261017_fine/manifest.json)。
最终图只读取 fine 冻结数据，图例全部置于同一 panel 内；字形、600 dpi、四边向内
刻度、遮挡、裁切、颜色对比及灰度可辨性已经检查。物理宽度为 6.75 英寸，最小
建议插入宽度约 6.296 英寸，不能直接缩到单栏。PNG 与图注的完整 hash 由 manifest
保存，图注和派生比较表放在同名 results case。

之前的 standard 绘图尝试及 fine 首版保留；最终使用 `fine_png_review_r2`。R2 补齐
显式图例位置和 `prd_review_v1` 声明，未改变 fine 数值。首个 standard case 曾因
manifest 缺少三个验证字段退出，不能把该中间尝试当最终交付。

## 快度解释与后续参数带

本轮比较的是均匀介质的局域、指定介子部分产额数密度与实验的中快度最终产额。
完整局域动量积分不等于已经计算整场碰撞的全相空间最终产额；映射仍依赖冻结面、
流场、纵向热力学变化及共振衰变。仅凭 `y≈0` 不能推断实验比值必然更低。
例如 NA49 同实验 30A GeV 的中快度中心值 `21.2/83.0≈0.2554` 高于全相空间
`52.9/241≈0.2195`，见 PRC 77, 024903 (2008),
[0710.0118v2](https://arxiv.org/abs/0710.0118)，Table I。
因此接受度可列为待检验机制，不能作为本轮整体偏高的已证实解释。

极限温度、低能形状及束流化学势映射对幅度和峰位的后续探索见
[参数包络研究](parameter_envelope_20261018.md)。文献未给出的协方差不能补造；
参数包络与统计置信带分开。参数调优使用的 K+/pi+ 数据属于校准数据，不能同时
作为独立预测验证。本页的默认参数结果没有使用后续调优值。

## 工作树与验证

PR 整理工作树为 `D:/w/meson-freezeout-review`，分支 `codex/meson-freezeout-review`，
基线为 `698feaf8fada670fde81c93bab6093b50093935f`。
原计算工作树 `D:/w/meson-freezeout-defaults` 的输运优化、SciPlot 待办及其他修改保留；
根工作树未由本次收尾改写。PR 仅纳入本研究链的实现和证据。

历史默认阶段完成 Julia 613 个断言与 Python 7 项验证；多项式阶段完成集成 71 个断言和
Python 17 项，Andronic 阶段完成 Python 25 项。原记录分别见
[默认验证](validation.json)、[多项式验证](parameter_envelope_validation_20261018.json)、
[Andronic 最终核验](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_execution/final_validation.json)。
这些是当时运行事实；本次独立工作树检查、证据 hash 和评分复算见
[PR 验证记录](pr_validation_20261022.json)。
