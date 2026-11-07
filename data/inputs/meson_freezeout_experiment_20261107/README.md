# 重核碰撞 K+/pi+ 实验参考更新

核对日期：2026-11-07。本版用于既有冻结线曲线的论文对照，继承
[2026-10-18 版](../meson_freezeout_experiment_20261018/README.md)的 22 行，并新增
STAR 54.4 GeV 的正式发表点，共 23 行。原有行的全部字段保持一致，包括 NA49 17.27 GeV
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
采用时输出到新的图层目录。当前已准备输入，尚未替换论文图像。
论文接入目标为 `D:/Desktop/paper/My_Paper/pnjl_meson_mott_xi`；结果用途和章节重组范围见
[介子研究入口](../../../docs/analysis/relaxtime/meson_freezeout_defaults_20261017_analysis/README.md#论文接入位置与范围)。

## 峰附近与系统大小的新资料

| 来源 | 已核对内容 | 建议用途 |
| --- | --- | --- |
| NA61/SHINE 2026，[arXiv:2606.00597](https://arxiv.org/abs/2606.00597)，§3、Fig. 7 | Pb+Pb，30A GeV/c，约 7.6 GeV；预备结果与 NA49 horn 点相容。使用全相空间平均多重数，\(\pi^+\) 由 \(\pi^-\) 和 NA49 的电荷比估计 | 可在讨论中引用支持峰附近趋势；与现有中快度点另列。本次未获得可直接接入的正式同口径数表 |
| NA61/SHINE 2026，[arXiv:2608.04758](https://arxiv.org/abs/2608.04758)，§2、Fig. 1 | 汇总 Pb+Pb/Au+Au 的 horn 与 p+p、Be+Be、Ar+Sc、Xe+La 的系统大小差异；左图是全相空间比值，其他面板含中快度量；Xe+La 使用 \(y=0.4–0.6\) | 支持将 Au+Au、Pb+Pb 作为峰位主比较；较小体系用于单独的系统大小讨论 |
| NA61/SHINE 2024，[arXiv:2308.16683](https://arxiv.org/abs/2308.16683)，EPJC **84**, 416，[HEPData](https://www.hepdata.net/record/ins2692441) | 0–10% Ar+Sc，13A–150A GeV/c，已正式发表多粒子谱和 K/pi 比值 | 有正式数表，可作为系统大小对照；本次未并入重体系评分 |
| Xe+La 2025 学位论文，[DOI:10.17181/f9r5b-na547](https://doi.org/10.17181/f9r5b-na547) | 六个能量的 charged hadron 研究，INSPIRE 3133203 将其标为 thesis | 补充来源；纳入图前仍需逐表核对快度、中心度和发表状态 |

这里的“轻/重”指入射核构成的碰撞体系，不是末态氘核、氦核等轻核产额。
现有证据支持重体系中 horn 更明显的选图理由；Xe+La 虽比 Ar+Sc 大，在上述汇总中也没有
出现相同的尖峰。因此不把质量数或“重核”标签直接当作 horn 必然出现的判据。
K/pi 可讨论奇异产生相对 pion 产额的变化，联系奇异度/熵的代理量时仍需说明模型假设。

## 检索范围与来源

本次通过 INSPIRE REST 按最新顺序检索 STAR identified/kaon/pion production 与 NA61/SHINE
2024–2026 更新，并核对上述开放论文及 HEPData。STAR 的正式新点由查询
`collaboration STAR and date > 2023 and (title identified or title pion or title kaon) and (title production or title spectra)`
定位。2026 年 STAR BES-II 综述 [arXiv:2607.29616](https://arxiv.org/abs/2607.29616)
与 bulk/strangeness 综述 [arXiv:2607.11640](https://arxiv.org/abs/2607.11640)中的新奇异粒子、流和关联测量，
不能替代这里需要的同电荷 K/pi 产额表。本次核对的资料中未找到可替换既有峰区扫描的
BES-II 正式同口径 K+/pi+ 全能区表；这一结论仅限本次检索范围。

- 输入：`experimental_ratios.csv`；已有 22 行字段一致，新 STAR 值和分量误差来自下列原表。
- 原表：[raw/star2026_54p4_kpi_hepdata.json](raw/star2026_54p4_kpi_hepdata.json)，仅 1545 bytes。
- 获取入口：[HEPData record 168419 v1](https://www.hepdata.net/record/ins3090662?format=json)，table id `1918269`；公开状态 `finished`，表审查 `passed`。
- 原表 SHA-256：`49090f30337c4074a8077764a03e3977535a6206732e879264321f84ca6fc580`。
- 父 CSV SHA-256：`8fc1b0846f975bf30b4f6bfdb585a054b4480ba784dc6536344b108969ae6c50`。
