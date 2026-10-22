# 近期化学冻结线默认参数补筛

日期：2026-10-17。范围为定向检索，不声称穷尽全部近期文献。

- 检索原始响应、实际返回数量和选定条目见 `screening_index.json`。
- PDF 放在既有“化学冻结线与Kpi_20261017”论文目录；每份 `.source.json` 记录路径、来源和 SHA-256。
- 本次8篇候选文献以及 FRG 的正式出版版本均已下载。未改写 Zotero 原库或 paper/bib。
- Zotero 本地 API 未响应，随后用只读数据库查标题及8个 arXiv ID，未检出这些新候选。

| ID | 文献 | 本次用途 |
| --- | --- | --- |
| S17 | Lu et al., PRD 105, 034012 (2022), 2109.09912 | DSE 及有限尺寸低温分支；Table II，不能任取尺寸把低温当普适冻结条件 |
| S18 / S18b | Shao, Fu, Liu, PRD 109, 034019 (2024), 2305.16763 | 两组完整FRG冻结线；预印本p.7及出版版p.8 Table IV参数一致 |
| S19 | Flor, Olinger, Bellwied, PLB 814, 136098 (2021), 2009.14781 | 味依赖冻结；Table 2给离散点，不将约150 MeV当全能区统一参数 |
| S20 | Lysenko et al., PRC 111, 054903 (2025), 2408.06473 | E/N=0.951 GeV隐式冻结判据；其muB=0温度160 MeV，不能假定更低 |
| S21 | Lu et al., 2603.09336v1 (2026预印本) | BES-I/II泛函QCD冻结线；Table I有mean和成组上下界，本轮只用mean |
| S22 | Poberezhnyuk et al., PRC 100, 054904 (2019), 1906.01954 | Table II IHRG分支完整参数化，T0=152 MeV |
| S23 | Singh, Sarkar, Biswas, 2307.05930 (2023预印本) | 排斥HRG；温度/化学势并非普遍降低，不作优先候选 |
| S24 | Sarkar, 2511.15707 (2025预印本) | PCA/Bayesian逐能量提取；Au+Au温度偏高，未提供本轮优先的低温完整曲线 |

实际q=0外推、参数转换、比较能区及后续结论统一保存在：
[D:/w/meson-freezeout-defaults/docs/analysis/relaxtime/meson_freezeout_defaults_20261017_analysis/literature_candidates_20261017.md](D:/w/meson-freezeout-defaults/docs/analysis/relaxtime/meson_freezeout_defaults_20261017_analysis/literature_candidates_20261017.md)。

FRG原文没有施加奇异数中性，7.7 GeV未用于冻结线拟合；2026方案在零化学势处采用作者的手征线约束。这些限制随参数化保留。
