# 重核碰撞 K+/pi+ 实验参考更新

核对日期：2026-11-08。本版用于既有冻结线曲线的论文对照，继承
[2026-10-18 版](../meson_freezeout_experiment_20261018/README.md)的 22 行，并新增
STAR 54.4 GeV 的正式发表点及 NA61/SHINE 7.62 GeV Pb+Pb 全相空间预备点。
2026-11-08 按作者要求进一步加入全相空间扫描 11 点、BRAHMS 中快度 1 点及 BES-II 派生比值 3 点，共 39 行。原有行的全部字段保持一致，包括 NA49 17.27 GeV
已修正的 0–5% 中心度；旧版关于接受度、派生误差和 HADES 读图点的说明继续适用。

## 可直接补入图中的正式新点

STAR, *Identified charged hadron production in Au+Au collisions at 54.4 GeV with the STAR detector*,
Phys. Rev. C **113**, 054907 (2026)，[论文](https://arxiv.org/abs/2512.06415)、
[DOI](https://doi.org/10.1103/4161-dflc)、[HEPData Table 14](https://doi.org/10.17182/hepdata.168419.v1/t14)。

| 字段 | 发布值/口径 |
| --- | --- |
| 碰撞体系、能量 | Au+Au，54.4 GeV |
| 中心度、快度 | 0–5%，\(|y|<0.1\)；论文 Fig. 9 |
| \(K^+/\pi^+\) | 0.15879 |
| 统计误差、系统误差 | 0.00005、0.01715 |
| 绘图合成误差 | 0.017150072886142495，按论文将两项平方相加开根 |
| 产额定义 | 中快度积分产额之比；pion 谱已校正弱衰变和 muon 污染 |

直接采用已发布的 ratio 及其误差；论文对 ratio 的系统误差已处理共同 tracking efficiency
成分，不从两个产额重新作无协方差传播。Table 14 的正电荷结果与同论文
[Table 10](https://doi.org/10.17182/hepdata.168419.v1/t10)中最中心的 \(\langle N_{part}\rangle=346\)
行逐项一致。Table 14 还给出 \(K^-/\pi^-=0.13201\)，留在原始 JSON 中，本正电荷 CSV 不混入该值。

该点补充高能区对照，对约 7–8 GeV 的峰位没有直接约束。新行 `primary_auau=False`、
`auxiliary_horn=False` 表示它是本次显示补充：现有模型没有按本次要求在 54.4 GeV 新算点，
既有十四点/九点评分继续使用各自冻结观测集合。这两个字段不是实验质量等级。
绘图入口 `scripts/analysis/relaxtime/plot_freezeout_comparison.py` 的 `--experiment` 可指向本 CSV；
本批三张 `20261107_png_review` 图及介子稿、学位稿中的副本使用本表。按作者要求原位替换此次未接受的中间图，不另存被替换版本。
论文接入目标为 `D:/Desktop/paper/My_Paper/pnjl_meson_mott_xi`；结果用途和章节重组范围见
[介子研究入口](../../../docs/analysis/relaxtime/meson_freezeout_defaults_20261017_analysis/README.md#论文接入位置与范围)。

## 本批加入的 NA61/SHINE 峰区参考点

来源为 [arXiv:2606.00597v1](https://arxiv.org/abs/2606.00597)，§4、Fig. 7（PDF 第 11 页），
作者为 A. Bazgir，代表 NA61/SHINE。采用 0–7.2% 中心 Pb+Pb、30A GeV/c，
\(\sqrt{s_{NN}}=7.62\) GeV；使用外推到全相空间的平均多重数比
\(\langle K^+\rangle/\langle\pi^+\rangle\)，其中 \(\pi^+\) 由测得的 \(\pi^-\) 和 NA49 电荷比估计。

| 量 | 源矢量图读数（近似） |
| --- | ---: |
| ratio | 0.23625 |
| 统计误差 | 0.00043 |
| 系统误差 | 0.01895 |
| 绘图误差 | 0.0189548780（两者平方和开根） |

未取得官方数表，数值来自 arXiv 源附件 `hornPlus4Pi_corr.pdf` 的矢量坐标。
横轴按 log10、纵轴按线性刻度标定，并用未参与拟合的 100 GeV 和 0.1 刻度复核。
同一测量在 [arXiv:2608.04758v1](https://arxiv.org/abs/2608.04758) Fig. 1 左图中的独立图文件
给出中心值差小于 \(10^{-6}\)；其图注明确竖线为统计误差、色带为系统误差。
两篇会议稿展示同一个预备结果，不能计为两个实验点。图形读数精度也不构成额外实验误差。

三图使用独立空心绿色菱形，图例注明 `NA61/SHINE (4pi, prelim.)`。该点与中快度测量的
接受度不同；`primary_auau=False`、`auxiliary_horn=False`，不加入原九点/十四点评分。
它支持在已知峰能区增加独立重体系参照，单个点本身不能确定峰形、峰宽或更精确峰位。

只保留 [17.8 kB 源图](raw/na61_2026_pbpb_4pi_fig7.pdf)和[坐标、误差及复核记录](raw/na61_2026_pbpb_4pi.json)，
不保留整份下载包、渲染截图或重复过程日志。

## 峰附近与系统大小的新资料

| 来源 | 已核对内容 | 建议用途 |
| --- | --- | --- |
| NA61/SHINE 2026，[arXiv:2606.00597](https://arxiv.org/abs/2606.00597)，§4、Fig. 7 | Pb+Pb，30A GeV/c，约 7.6 GeV；全相空间预备结果，定义和读数见上节 | 已加入三图，保持全相空间与预备结果标签 |
| NA61/SHINE 2026，[arXiv:2608.04758](https://arxiv.org/abs/2608.04758)，§2、Fig. 1 | 汇总 Pb+Pb/Au+Au 的 horn 与 p+p、Be+Be、Ar+Sc、Xe+La 的系统大小差异；左图是全相空间比值，其他面板含中快度量；Xe+La 使用 \(y=0.4–0.6\) | 支持将 Au+Au、Pb+Pb 作为峰位主比较；较小体系用于单独的系统大小讨论 |
| NA61/SHINE 2024，[arXiv:2308.16683](https://arxiv.org/abs/2308.16683)，EPJC **84**, 416，[HEPData](https://www.hepdata.net/record/ins2692441) | 0–10% Ar+Sc，13A–150A GeV/c，已正式发表多粒子谱和 K/pi 比值 | 有正式数表，可作为系统大小对照；本次未并入重体系评分 |
| Xe+La 2025 学位论文，[DOI:10.17181/f9r5b-na547](https://doi.org/10.17181/f9r5b-na547) | 六个能量的 charged hadron 研究，INSPIRE 3133203 将其标为 thesis | 补充来源；纳入图前仍需逐表核对快度、中心度和发表状态 |

这里的“轻/重”指入射核构成的碰撞体系，不是末态氘核、氦核等轻核产额。
现有证据支持重体系中 horn 更明显的选图理由；Xe+La 虽比 Ar+Sc 大，在上述汇总中也没有
出现相同的尖峰。因此不把质量数或“重核”标签直接当作 horn 必然出现的判据。
K/pi 可讨论奇异产生相对 pion 产额的变化，联系奇异度/熵的代理量时仍需说明模型假设。

## 进一步找到的可比资料

### 2026-11-08 纳入的能量点

- [HEPData 25202 / Table 9](https://www.hepdata.net/record/25202) 的正电荷全相空间比值全部 11 点已纳入：
  AGS Au+Au 五点（2.35、3.32、3.83、4.29、4.87 GeV），NA49 Pb+Pb 五点
  （6.27、7.62、8.77、12.32、17.27 GeV），BRAHMS Au+Au 一点（200 GeV）。
  最后一项的原始来源为 PRL **94**, 162301 (2005)，[nucl-ex/0403050](https://arxiv.org/abs/nucl-ex/0403050)，
  NA49 原文 Fig.4 说明了这一归属。沿用原表已发布比值和误差，NA49 的非对称系统误差分别合成。
  全相空间和中快度用空心/实心标记区分，保留实际实验身份。
- [HEPData 104454 / Table 4A](https://www.hepdata.net/record/104454) 加入 BRAHMS 62.4 GeV 的 `YRAP=0` 一点，
  比值为 0.159，误差为统计 0.0034 和系统 0.012 的平方和开根。其余七行是同一能量下的粒子快度扫描。
  该点为 0–10% 中心 Au+Au，图例与 AGS/STAR 分开。
- STAR BES-II 的 [EPJ Web Conf. **296**, 14009 (2024)](https://doi.org/10.1051/epjconf/202429614009)
  Fig.1（PDF 第2页）采用 0–5% Au+Au 的 `y=0` 图示 bin。提取实心蓝圆 `pi+` 和实心红方 `K+`，
  后者图示值和误差均除以 2.5。横条约为 y∈[-0.05,0.05]，这是图形范围，正文未另给同样精确的切选声明。
  54.4 GeV 继续使用已纳入的正式发表值。

| 能量 / GeV | pi+ 的 dN/dy | K+ 的 dN/dy | 派生 K+/pi+ | 零协方差传播误差 |
| --- | ---: | ---: | ---: | ---: |
| 14.6 | 144.11 | 24.83 | 0.17230 | 0.01332 |
| 19.6 | 157.29 | 26.91 | 0.17107 | 0.01240 |
| 27 | 167.27 | 31.70 | 0.18950 | 0.01544 |

这些 BES-II 数值是预备产额的读图派生量。图中原误差为统计与系统合成误差，采用
`sigma_R/R=sqrt((sigma_K/K)^2+(sigma_pi/pi)^2)`；K 与 pi 的协方差未发布，因此图中误差不是合作组发布的 ratio 误差。
对任意相关系数 rho，传播式中增加 `-2*rho*(sigma_K/K)*(sigma_pi/pi)`。
由取点与刻度误差给出的比值读图区间分别为 [0.16993,0.17469]、[0.16890,0.17325]、[0.18740,0.19162]，
与上述实验产额传播误差分别记录。图例使用 `BES-II (derived, prelim.)`。

读图保留[原图](raw/star_besii2024_fig1.png)及[标定、误差和公式](raw/star_besii2024_midrapidity.json)；
直接表格保留[全相空间原表](raw/na49_4pi_energy_hepdata.json)和[BRAHMS 原表](raw/brahms_62p4_rapidity_hepdata.json)。
取点叠加与辅助模型审查放在临时目录，最终数据经源图、独立刻度与倍率复核。
全部新增点为显示比较，不改变原有选参集合。参数带图扩展显示至 2.35 GeV，以完整呈现全相空间扫描；理论曲线仅显示原有计算支持。

| 数据 | 实际内容与可比性 | 建议 |
| --- | --- | --- |
| BRAHMS，Au+Au 62.4 GeV，PLB **687**, 36 (2010)，[arXiv:0911.2586](https://arxiv.org/abs/0911.2586)，[HEPData Table 4A](https://doi.org/10.17182/hepdata.89451.v1/t36) | 0–10% 中心，\(y=0\) 的正电荷比为 \(0.159\pm0.0034\,\mathrm{(stat.)}\pm0.012\,\mathrm{(sys.)}\)；另外提供多个快度点。这是相对现用 STAR 的独立实验测量 | y=0 一点已于 11-08 加入；其他快度用于独立的快度依赖研究 |
| STAR BES-II，EPJ Web Conf. **296**, 14009 (2024)，[DOI](https://doi.org/10.1051/epjconf/202429614009)，Fig. 1 | 0–5% Au+Au 的 \(\pi^\pm,K^\pm,p,\bar p\) 预备 \(dN/dy\)，图中有 14.6、19.6、27 GeV，快度扩展到 0.9。论文另讨论 54.4 GeV | 已加入三个 y=0 派生比值，零协方差假设和读图精度见上文；其余快度可用于纵向分布讨论 |
| AGS/NA49/BRAHMS 全相空间能量扫描，[HEPData Table 9](https://doi.org/10.17182/hepdata.25169.v1/t9) | 共 11 点，覆盖 2.35–200 GeV；其中五个 NA49 点的比值为 0.214、0.219、0.204、0.178、0.174 | 全部已加入。NA49 部分与新 NA61 的接受度接近，也与现有 NA49 中快度点来自同批测量；两种接受度适合并列比较，统计分析时须考虑其关联 |
| NA49，PRC **86**, 054903 (2012)，[arXiv:1207.0348](https://arxiv.org/abs/1207.0348) | 40A、158A GeV 的 charged pion/kaon 系统大小与中心度研究，覆盖 Pb+Pb 和较小系统 | 可研究中心度和系统大小对比；不构成新的密集峰区能量扫描 |

近期 STAR fixed-target [2026 会议稿](https://doi.org/10.1051/epjconf/202636411001)虽覆盖
3.2–6.2 GeV，主要报告 \(K_S^0,\Lambda\) 及其比值；正式 [JHEP10(2024)139](https://doi.org/10.1007/JHEP10(2024)139)
同样针对 \(K_S^0,\Lambda\)。这些结果不能直接替代正电荷 K/pi。
PHENIX [PRC **88**, 024906](https://doi.org/10.1103/PhysRevC.88.024906) 的公开 HEPData
主要给出 \(p_T\) 微分谱和比值，不能把其中任意 \(p_T\) bin 当作积分 K/pi 点。
缺少尖峰的 Ar+Sc、Xe+La 仍是模型系统大小解释的有用约束；不能只选择出现尖峰的数据。

## 检索范围与来源

本次通过 INSPIRE REST 按最新顺序检索 STAR identified/kaon/pion production 与 NA61/SHINE
2024–2026 更新，并核对上述开放论文及 HEPData。STAR 的正式新点由查询
`collaboration STAR and date > 2023 and (title identified or title pion or title kaon) and (title production or title spectra)`
定位。2026 年 STAR BES-II 综述 [arXiv:2607.29616](https://arxiv.org/abs/2607.29616)
与 bulk/strangeness 综述 [arXiv:2607.11640](https://arxiv.org/abs/2607.11640)中的新奇异粒子、流和关联测量，
不能替代这里需要的同电荷 K/pi 产额表。本次核对的资料中未找到可替换既有峰区扫描的
BES-II 正式同口径 K+/pi+ 全能区表；已找到的预备产额列于上表，这一结论仅限本次检索范围。

检索于 2026-11-07 完成，11-08 追加数据表、快度及图形倍率核验；使用 INSPIRE REST 的 STAR identified/FXT/BES-II、
BRAHMS kaon/pion Au+Au、PHENIX identified Au+Au、NA49 energy/system-size 条目，
再核正文与 HEPData。关键查询为
`collaboration STAR and (title production or title bulk) and (title BES-II or title beam energy or title 3 GeV)`、
`collaboration BRAHMS and (title production or title yields or title spectra) and (title identified or title pion or title kaon) and title Au+Au`。
保留 INSPIRE 2805139、836865、1120799；剔除流、关联和模型计算，合并同一测量的重复报道。

- 输入：`experimental_ratios.csv`；已有 22 行字段一致，新增实验的身份、接受度和误差定义见各节。
- 原表：[raw/star2026_54p4_kpi_hepdata.json](raw/star2026_54p4_kpi_hepdata.json)，仅 1545 bytes。
- 获取入口：[HEPData record 168419 v1](https://www.hepdata.net/record/ins3090662?format=json)，table id `1918269`；公开状态 `finished`，表审查 `passed`。
- 原表 SHA-256：`49090f30337c4074a8077764a03e3977535a6206732e879264321f84ca6fc580`。
- 父 CSV SHA-256：`8fc1b0846f975bf30b4f6bfdb585a054b4480ba784dc6536344b108969ae6c50`。
