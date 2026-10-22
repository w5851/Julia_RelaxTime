# 近期化学冻结线：文献默认参数的补充筛选

日期：2026-10-17。状态：本次定向文献筛选与限定能区数值核验完成，
`diagnostic_only`。只用文献中心参数，没有拟合、参数带或正式产物晋升。

## 当前结论与比较范围

**存在比 Andronic 2018 更接近当前实验输入的默认参数化。** 在共同的
11.5–200 GeV、8 个 Au+Au 点上，FRG 2024 的 mu-independent 分支、IHRG 2019、
functional QCD 2026 的 mean 分支分别将等权 RMSE 降低约 15.2%、13.5%、9.5%。
这不是原先 14 点主评分的替代结果，也不表示 horn 的峰位或尖锐度已经改善。

每种新方案在 7.7、11.5、14.5、19.6、27、39、62.4、130、200 GeV 直接计算，
标准和加密各 36 点，全部成功。由于 FRG 原文明确未将 7.7 GeV 纳入冻结线拟合，
统一评分预先限定为其余 8 点；7.7 GeV 的计算单列为 horn 附近的诊断。
以下全部为 fine 结果，旧两条线复用已冻结的同档、同能量计算。

| 参数化 | 8 点 RMSE | 相对 Andronic 改善 | R(7.7 GeV) | R(200 GeV) |
| --- | ---: | ---: | ---: | ---: |
| Cleymans 2006 | 0.181647 | -19.27% | 0.415393 | 0.298910 |
| Andronic 2018 | 0.152300 | 0 | 0.389881 | 0.274506 |
| IHRG 2019 | 0.131766 | 13.48% | 0.421767 | 0.254386 |
| FRG 2024，mu-dependent | 0.138199 | 9.26% | 0.398670，拟合范围外 | 0.263496 |
| FRG 2024，mu-independent | **0.129190** | **15.17%** | 0.341983，拟合范围外 | 0.260447 |
| functional QCD 2026，mean | 0.137768 | 9.54% | 0.362050 | 0.270226 |
| 实验 | — | — | 0.2227 | 0.159 |

IHRG 的高能平台改善最明显，却在 7.7 GeV 比 Andronic 更高。因此只看零化学势
温度不足以筛选整个能区。FRG 的 mu-independent 分支数值上最有希望，但其原文
守恒荷条件及低能适用范围需要单独处理。2026 方案实际纳入 BES-II 低能涨落数据，
适合作为 horn 附近完整轨迹比较的优先对象；它仍是预印本。

这些方案在已比较点上仍整体偏高。9 个离散能量不足以确定新曲线的连续峰位或
尖锐程度，不能由较低的 R(7.7) 直接宣称解决 horn。后续默认方案的完整能区选择
应保留上述候选，再决定参数带的起点，而非继续仅在旧两方案之间选择。

## 可直接复用的原文参数

### Poberezhnyuk et al. 2019：IHRG 分支

来源：Poberezhnyuk et al., *Phys. Rev. C* **100**, 054904 (2019)，
[1906.01954v2](https://arxiv.org/abs/1906.01954v2)，Eq.(5)、PDF pp.8–9 Table II。
这是该文的 Ideal-HRG 行，不是另一组 QvdW 参数。

令 z 为 sqrt(sNN) 按 GeV 计的数值，T、muB 使用 GeV：

```math
\mu_B(z)=\frac{1.310}{1+0.278z},\qquad
T_f=0.152-0.026\mu_B^2-0.219\mu_B^4.
```

原文五个系数均有误差；本轮逐项采用中心值。其形式和现有多项式配置兼容，
[ihrg2019.toml](../../../../config/physics/freezeout/ihrg2019.toml) 可直接加载。
在 7.7 GeV，它给出 `(T,muB)=(140.847,417.118) MeV`；Andronic 为
`(138.428,406.359) MeV`，解释了“高能极限更低”不等于 horn 附近也更冷。
原文拟合涉及 HRG 谱、衰变与 gamma_s 等设置，本次只输入其 T、muB 轨迹，
不把原文 HRG 的产额拟合质量移植为本项目 GBU 的正确性证明。

### Shao、Fu、Liu 2024：两组 FRG 化学冻结线

来源：*Phys. Rev. D* **109**, 034019 (2024)，
[DOI](https://doi.org/10.1103/PhysRevD.109.034019)，
[2305.16763v2](https://arxiv.org/abs/2305.16763v2)。Eqs.(29)–(30)、Table IV
位于预印本 PDF p.7、正式出版版 p.8；两版参数已对照一致并视觉复核。

下式的 T0、muB、d 使用 MeV，z 仍按 GeV 计：

```math
\mu_B(z)=\frac{d}{1+ez},\qquad
T_f=T_0\left[1-\kappa_2(\mu_B/T_0)^2-\kappa_4(\mu_B/T_0)^4\right].
```

| 原文场景 | T0 / MeV | kappa2 | kappa4 | d / MeV | e / GeV^-1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| mu-dependent gluon potential | 154.8 | 4.05e-7 | 0.0022 | 708.96 | 0.147 |
| mu-independent gluon potential | 153.9 | 7.35e-7 | 0.0040 | 751.31 | 0.159 |

这里的 kappa2、kappa4 对应原文 a、b，不能误当作直接乘 muB^2、muB^4 的系数。
写入现有多项式配置时，以 GeV 表示 T0，严格使用 `b=kappa2/T0`、
`c=kappa4/T0^3`。两个场景均保留并计算，没有在两行之间拼接系数。
配置为 [mu-dependent](../../../../config/physics/freezeout/frg2024_mudep.toml)
和 [mu-independent](../../../../config/physics/freezeout/frg2024_muindep.toml)。

原文从净质子涨落与重子数易感率比较提取化学冻结条件，并非将 crossover 曲线
直接当冻结线。但原文 Sec.III（出版版 p.3）明示 `muI=muS=0`，未施加奇异数中性；
7.7 GeV 点因误差较大未用于参数拟合。当前计算保留项目的 `nS=0,nQ/nB=0.4`
BQS 背景，所以这里检验的是文献 T、muB 轨迹在本模型中的作用，未复现原文的
全部热力学条件。源文献的这些限制不能因 ratio 较低而省略。

### Lu et al. 2026：functional QCD 的 mean 分支

来源：[2603.09336v1](https://arxiv.org/abs/2603.09336v1)，
*Extracting freeze-out conditions in beam energy scan via functional QCD*，
2026-03-10 预印本；Eqs.(3)–(5)、PDF p.3 Table I。
作者包括 Lu、Fischer、Gao、Liu、Pawlowski。

采用与上节相同的归一化四次式，原文 mean 行为：

| T0 / MeV | kappa2 | kappa4 | d / MeV | e / GeV^-1 |
| ---: | ---: | ---: | ---: | ---: |
| 157 | 0.0153 | 0.00096 | 953.7 | 0.2026 |

[functional_qcd2026_mean.toml](../../../../config/physics/freezeout/functional_qcd2026_mean.toml)
逐项采用这些值。它在 7.7 GeV 给出 `(138.697,372.536) MeV`；温度与 Andronic
接近，muB 较低。这是新的束流能量映射，而非简单下调高能温度。

作者使用 BES-I 的 200–27 GeV、BES-II 的 27–7.7 GeV 涨落提取点，并加入
SPS 6.3 GeV 的产额冻结点约束低能展开。他们令 `Tf(0)=T_chiral(0)=157 MeV`
且 `kappa2_f=kappa2_chiral=0.0153`，再约束冻结线的高阶项和能量映射。
这是原文明确的物理假设，不是本项目额外作出的等同，也不能称为格点 QCD
独立预测的低温冻结线。

Table I 还给出成组的 upper/lower 参数：`(d,e,kappa4)=(1027.0,0.2143,0.00073)`
和 `(913.9,0.1977,0.00147)`。原文说明温度带与化学势带相关，后续可以采用成组
变体；本轮只用 mean，不抽取最低温端，不独立混配边界，不生成参数带。

## 其他近期候选及不优先采用的原因

- Lu et al., PRD **105**, 034012 (2022),
  [2109.09912v1](https://arxiv.org/abs/2109.09912v1)，PDF p.5 Table II：
  DSE 的有限尺寸 `L=3.5,2.8 fm` 分别给出 `T0=150.3,144.6 MeV`，还有更低的
  分支。它们是带尺寸条件的不同场景，不是一个普适的低温默认值；部分低能点
  在 Table I 中没有可用提取结果。选择 L 必须有与当前体系对应的物理依据。
- Flor、Olinger、Bellwied, PLB **814**, 136098 (2021),
  [2009.14781](https://arxiv.org/abs/2009.14781)，Table 2：公布了 11.5 GeV–5.02 TeV
  的逐能量拟合点。其 `pi,K,p` 子集在 200 GeV 为 `T=153.9 MeV`，19.6 GeV 为
  `149.0 MeV`，是可进一步比较的输入；该子集本身包含 K，不能简单说它完全不涉及
  kaon。文中的 `TL=150.2±2.6 MeV`、`TS=165.1±2.7 MeV` 来自味依赖曲线分析，
  不能只拿 TL 替换全能区所有物种的共同温度。连续插值方案也需单独说明。
- Lysenko et al., PRC **111**, 054903 (2025),
  [2408.06473](https://arxiv.org/abs/2408.06473)：用 `E/N=0.951 GeV` 和 HRG 的
  `nS=0,nQ/nB=0.4` 隐式确定冻结线，PDF p.3 明示零化学势温度为 **160 MeV**。
  不是一个把高能平台进一步降到 150 MeV 的方案。
- Singh et al. [2307.05930](https://arxiv.org/abs/2307.05930) 的排斥 HRG、
  Sarkar [2511.15707](https://arxiv.org/abs/2511.15707) 的 PCA/Bayesian 分析并不
  给出普遍更低的 Au+Au 冻结温度，因此未列为本次数值优先项。

这些文献事实与本项目的效果分开：后三类没有在本轮运行 GBU，不能据此断言
它们在本项目一定较差。HotQCD 2019、Borsanyi 2020 的 crossover 仍不作为
直接化学冻结参数化加入比较。

## 数值、来源与验证

- [参数来源、原始系数和坐标核对](literature_candidate_parameters.json)。
- [standard manifest](../../../../data/outputs/results/relaxtime/meson_density/freezeout_literature_candidates_20261017_standard/manifest.json)
  与 [fine manifest](../../../../data/outputs/results/relaxtime/meson_density/freezeout_literature_candidates_20261017_fine/manifest.json)。
- [比较指标与分辨率差异](../../../../data/outputs/results/relaxtime/meson_density/freezeout_literature_candidates_20261017_assessment/assessment.json)、
  [逐点残差](../../../../data/outputs/results/relaxtime/meson_density/freezeout_literature_candidates_20261017_assessment/residuals.csv)、
  [验证记录](../../../../data/outputs/results/relaxtime/meson_density/freezeout_literature_candidates_20261017_assessment/verification.json)。
- 运行入口沿用 [q0 比较脚本](../../../guides/scripts/meson_freezeout_comparison.md)；
  限定能区汇总由 [assess_freezeout_literature_candidates.py](../../../../scripts/analysis/relaxtime/assess_freezeout_literature_candidates.py)
  生成，并在 assessment 目录保存同字节源码快照。

背景、积分路线和有限窗口与首轮完全相同：quark-only BQS、`xi=0`、无介子反馈、
`q0_lambda_reference`、外部 omega 的 Bose 权重、`omega_lower=1e-5 fm^-1`、
`qmax=8 fm^-1`。每个点直接求解，未对旧等高图插值。
主评分只读取既有 Au+Au 输入的 8 个点，不受 NA49 中心度元数据勘误影响。

standard/fine 的最大相对比值差为 `3.535835e-5`，即 **0.003536%**，排序一致。
本次没有失败点；端点诊断警告分别出现在 11、12 个点，仍保留为有限窗口结果。
这一检查不认证零端点极限、q 尾部或正式生产收敛。

核对了 40 个源码/配置快照 hash，全部 72 点的 Julia 输出坐标与原文归一化
MeV 公式独立求值一致，最大差小于 `6e-14 MeV`；6 组 RMSE 也以高精度小数
独立重算。没有修改 `src/` 或密度核，因此未重复运行此前 613 个断言；旧数值
源码与背景配置的共用 hash 均核对一致。旧 runner 的差异仍仅是先前已记录的
CLI 默认能量表修正，两批均显式给定能量，旧快照保留。

本节保留当时的四个物理 profile、限定能区汇总及研究记录。随后已完成两族参数探索、
horn 附近加密及 PNG；它们作为同一研究链进入本次 PR，见[阶段结论](README.md#2026-10-22-阶段审阅结论)。
本节八点评分及首轮冻结输入保持原字节。

## 文献管理

定向 INSPIRE 检索覆盖 2019–2026：标题查询返回 39 条、近期涨落查询 16 条，
广义参数化查询读取最新 60/121 条。它不是穷尽性系统综述。
新文献 S17–S24 及 FRG 的正式出版版 S18b 已保存到：

`D:/Desktop/PNJL论文/01_论文与综述/相图与守恒荷涨落/化学冻结线与Kpi_20261017/`

检索原始响应、PDF hash、版本与索引在：

`D:/Desktop/PNJL论文/05_阅读笔记与复现计划/化学冻结线与Kpi_20261017/recent_default_candidates_20261017/`

Zotero 本地 API 未响应；随后使用 `mode=ro`、`query_only` 查标题和八个 arXiv ID，
未检出这批新候选。未改写 Zotero 原库或 `paper/bib`。

本次 PR 保存[仓库内文献副本与检索记录](../../../../data/inputs/meson_freezeout_polynomial_literature_20261017/README.md)，
并通过[证据归档清单](../../../../data/outputs/results/relaxtime/meson_density/freezeout_review_20261022/archive_manifest.json)
核对原始字节。上述绝对路径保留为历史来源定位。
