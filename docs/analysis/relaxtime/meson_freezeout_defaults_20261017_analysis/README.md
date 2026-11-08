# 化学冻结线与 K+/pi+：研究结果与论文用法

当前状态：默认参数比较、两族参数带和 M4 已于 2026-11-07 获作者接受；更新实验图件、三张矢量 PDF 与论文接入已完成。作者于 2026-11-08 确定采用简短的快度讨论，并授权本批计算任务收尾；后续文字审阅由论文库承接。
[PR #327](https://github.com/w5851/Julia_RelaxTime/pull/327) 已于 2026-10-22 合并实现、参数和冻结证据。
论文按下述决定使用这些结果；原产物中的 `png_review`、`diagnostic_only` 保留生成时的身份，当前接受状态由本页记录。
本页保留默认参数阶段的证据和评分口径。Ledger item：
`meson-freezeout-default-comparison`。

首轮[冻结公式与快度口径](../meson_freezeout_experiment_20261017_analysis/README.md)
保留文献参数约定和 NA49 中快度/全相空间对照。

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
方法、精确参数及全部图像入口见上述 Andronic 研究页；两族结果均已获作者接受。

## 2026-11-07 作者接受与论文用法

组会反馈支持优先使用文献原始、未调优的成套参数，将 K+/pi+ 作为偏定性的模型结果讨论。
论文重点比较峰出现的碰撞能量及其随冻结路径的变化，并如实报告幅度和峰形的差异。
此前两族带状调优结果一并接受，作为参数灵敏度与路径依赖的分析依据；参数带仍是有限候选包络。
现有 Andronic 2018 默认参照和 Cleymans 2006 对照均保留，作者尚未指定只保留其中一组。

已有采样峰位为：Cleymans 2006 约 9.5 GeV、Andronic 2018 约 9.25 GeV；
NA49、STAR 发布点中的最大值分别在 7.62、7.7 GeV。原始参数比调优后 Andronic 的
14.5–15 GeV 宽平台更接近这些实验峰点；多项式最优 `hornmap_001` 的峰位约 7.5 GeV，
也接近实验。因此论文选择原始参数的理由是保留已发表的完整参数组合并作定性比较，
而非所有调优结果的峰位都较差。离散实验峰点也不等于已精确确定的连续峰位置。

原始参数下的幅度仍明显偏高，当前模型尚无最终态衰变、流场与快度接受度的完整映射。
峰位是本论文强调的比较量，实验的幅度、峰形和系统大小依赖同样有物理信息。
后续实验参考优先核对 Au+Au、Pb+Pb 等重碰撞体系的同电荷、中快度 K+/pi+；
最新数据的筛选与采用建议见[实验输入入口](../../../../data/inputs/meson_freezeout_experiment_20261107/README.md)。

### 论文接入位置与范围

作者于 2026-11-07 明确：本研究接入 `D:/Desktop/paper/My_Paper/pnjl_meson_mott_xi` 介子期刊稿。
本批已重写标题、摘要、方法、结果和结论，保留 pi/K 谱图并接入三张冻结线图、九点评分比较和重核实验讨论。
`sections/density_method.tex` 与 `sections/freezeout_results.tex` 说明 `xi=0`、quark-only BQS 背景及 q0 类时外推，
与历史各向异性谱扫描分别呈现。核对源码后，旧“完整复极点”叙述改为实轴环函数与质量宽度联立处方。
谱与阈值变化、密度积分、冻结路径及实验比值之间的联系按已有证据建立；当前密度结果不提供 ratio 的各向异性扫描。

M4 的二维等高图已获作者接受，用于观察 ratio 趋势并指导已完成的参数调优，见[归档记录](../../../dev/archived/2026-09-30_charged-gbu-phase-contours.md)。
作者随后授权先处理 [Arena Issue #313](https://github.com/w5851/Julia_RelaxTime/issues/313)。
五项约定的[本地审计](../../../reference/formula/relaxtime/transport/RS_Response_Convention.md)与 quote-remote 32/32 检查已完成，
未修改数值算法、基线或 v15；作者接受保留现行定义，GitHub issue 已于 2026-11-07 按[审计结论](https://github.com/w5851/Julia_RelaxTime/issues/313#issuecomment-6706495944)关闭。完整颜色张量响应与备选 gap 闭合属于可选模型研究。
作者于 2026-11-08 接受三张更新图及实验范围，矢量版已完成并接入两份稿件。计算侧交付已收口，后续论文审阅留在论文库；PALC 性能诊断及全温区 PR3 保持未排期研究。

### 本批显示更新与论文草稿（2026-11-08）

| 内容 | 已交付矢量 PDF |
| --- | --- |
| 文献默认参数 | [PDF](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_experiment_20261108_vector/kplus_piplus_freezeout_defaults.pdf) |
| 多项式参数带 | [PDF](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_experiment_20261108_vector/kplus_piplus_polynomial_envelope.pdf) |
| Andronic 参数带 | [PDF](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_experiment_20261108_vector/kplus_piplus_andronic_envelope.pdf) |

三图均复用已冻结数值，11-08 显示更新使用 39 行实验表：除 STAR 54.4 GeV 和 NA61/SHINE
7.62 GeV 全相空间预备点外，加入 HEPData 25202 的全部 11 点、BRAHMS 62.4 GeV 中快度点和
BES-II 14.6、19.6、27 GeV 的三个派生比值。全相空间点用空心标记，BES-II 单列 `derived, prelim.`，
其比值误差按公开产额误差和零协方差假设传播。参数带图扩展实验显示至 2.35–200 GeV，含 36 点，
默认图含全部 39 点；理论支持与缺口保持原样。本批先前的未接受中间图已原位替换。
当前三个 PNG case 是本次已接受的图面依据，矢量包直接引用，未重复复制 PNG。
原九点/十四点评分、数值源、带成员和失败点均经哈希或字段复核保持一致。
近优带使用密斜线以改善灰度区分，三图的输入/输出 manifest 校验通过。
三张 PDF 均为 6.75 × 6 英寸的单页矢量图，嵌入非 Type 3 字体、无栅格包装；重新渲染 PNG 与作者接受文件逐字节一致。
[唯一交付 manifest](../../../../data/outputs/figures/relaxtime/meson_density/freezeout_experiment_20261108_vector/plot_manifest.json)绑定原接受文件、冻结输入、独立导出器与 PDF。
定向绘图及交付测试为 105 passed。没有重新进行密度扫描，精度结论仍来自原 standard/fine 对照。
交付 manifest 记录渲染源码 hash；结果目录不重复保存源码。
新增实验及 BRAHMS、STAR BES-II、NA49 全相空间候选的使用建议统一记在[实验输入入口](../../../../data/inputs/meson_freezeout_experiment_20261107/README.md)。

介子稿 `D:/Desktop/paper/My_Paper/pnjl_meson_mott_xi/main.pdf` 已生成 13 页；
学位论文 `D:/Desktop/paper/My_Paper/xjtu_master_thesis/` 使用学校 2025 模板，按同组硕士论文风格重组为六章两附录，快度讨论精简后为 114 页，
删除工程过程章节，集中于 PNJL/RS、相移关联、冻结路径、v15 输运及 Mott 响应的物理解释；正文引用 237 篇去重文献。
新增已接受的二维 ratio 及两族 T–muB 路径图，并展开守恒荷、角动量阈值与可观测响应的推导；快度影响按作者决定作简短定性讨论。
预览为 `build/main.pdf`，同组篇幅比较与最新编译核验记录集中在该目录 README。
本批图面、实验范围已于 2026-11-08 获作者接受并完成矢量交付，计算侧任务归档；2026-11-07 对原有数值、参数带和 M4 的接受保持有效。稿件后续审阅由论文库管理。
本次不把固定窗口的数值稳定性推及零频极限，也不把局域关联密度直接等同于接受度内的最终态产额。

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

2026-10-22 的交付记录为 `png_review`；作者接受和论文用途已于 2026-11-07 在上节确定。
[证据恢复说明](../../../../data/outputs/results/relaxtime/meson_density/freezeout_review_20261022/README.md)
定位数值及检索 ZIP、原始逐点记录和源码快照，
[本次 PR 验证](pr_validation_20261022.json) 记录在独立工作树的复核。
本批不再追加参数优化；后续工作是重核实验数据核对和按上述用途组织论文分析。

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

## 快度解释与参数带范围（2026-11-08 核对）

图例 `4π` 表示全相空间总产额，通常包括对未测快度和横动量区域的实验外推。
中快度是碰撞质心系中 `y=0` 附近的区域，实验表中的 `YRAP` 是粒子快度。

作者确定将接受度作为可能影响简短讨论。[BRAHMS HEPData 104454](https://www.hepdata.net/record/104454)
在 62.4 GeV、0–10% Au+Au 中显示离开中快度后 `K+/pi+` 有偏大的趋势。
前向区域可能保留更多净重子，相应的化学组成、较大的奇异化学势及 kaon–超子伴随产生
可能提高正 kaon 相对 pion 的丰度。这为中快度比值可能较小提供了一种物理解释，
纵向组成与接受度也可能影响本文局域热比值的实验比较。两份稿件已采用这一简短表述，
本批快度讨论与计算到此收口。

默认两条理论曲线各含 50 个能量点，覆盖 **2.35–5020 GeV**。
两族参数包络各使用 **4.87–200 GeV 的 17 个初筛能量**，其原计算域由
`literature_envelope_v1.toml` 和 `andronic_exploration_v1.toml` 指定。
近优候选使用 **7.7–200 GeV 的九个 Au+Au 点**选参，4.87 GeV 及 Pb+Pb horn 支点补充形状诊断。
本次实验显示扩展至 2.35–200 GeV，未重新选参。较窄的参数带来自此前计算域，
若扩展到更低能或 LHC，需要沿每条候选路径补算密度，单独延长解析冻结线不能提供比值带。

极限温度、低能形状及束流化学势映射对幅度和峰位的后续探索见
[参数包络研究](parameter_envelope_20261018.md)。文献未给出的协方差不能补造；
参数包络与统计置信带分开。参数调优使用的 K+/pi+ 数据属于校准数据，不能同时
作为独立预测验证。本页的默认参数结果没有使用后续调优值。

## 工作树与验证

历史 PR 整理工作树为 `D:/w/meson-freezeout-review`，计算基线为
`698feaf8fada670fde81c93bab6093b50093935f`；PR #327 合并提交为
`9c727d106abb9e931ba22dd26c206dc38d2efdff`。该 worktree 和分支已经回收。
原计算工作树 `D:/w/meson-freezeout-defaults` 也已在输运重算交付后回收；
其冻结线证据由[五个归档](../../../../data/outputs/results/relaxtime/meson_density/freezeout_review_20261022/README.md)保留。
当前研究使用本页及 ledger 的接受状态，无需恢复旧工作树即可阅读结果。

历史默认阶段完成 Julia 613 个断言与 Python 7 项验证；多项式阶段完成集成 71 个断言和
Python 17 项，Andronic 阶段完成 Python 25 项。原记录分别见
[默认验证](validation.json)、[多项式验证](parameter_envelope_validation_20261018.json)、
[Andronic 最终核验](../../../../data/outputs/results/relaxtime/meson_density/andronic_envelope_20261018_execution/final_validation.json)。
这些是当时运行事实；本次独立工作树检查、证据 hash 和评分复算见
[PR 验证记录](pr_validation_20261022.json)。
