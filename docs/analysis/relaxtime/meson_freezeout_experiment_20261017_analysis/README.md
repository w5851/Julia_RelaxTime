# 化学冻结条件与实验快度口径

本页保留首轮文献公式、原生网格和 NA49 对照的解释。默认参数比较、两族探索
及当前审阅状态统一见[研究链主报告](../meson_freezeout_defaults_20261017_analysis/README.md)。
首轮数值及文献定位保留在 [evidence.json](evidence.json)，原字节不变。

## 文献与公式来源

通过 INSPIRE 公开只读 API 定向检索，核对开放 arXiv 全文；查询记录、筛选范围、PDF 版本和 SHA-256 保留在文献索引。本次是围绕实现问题的首轮筛选，不声称穷尽最新文献。

- PDF：`D:/Desktop/PNJL论文/01_论文与综述/相图与守恒荷涨落/化学冻结线与Kpi_20261017/`
- 索引、`sources.json`、`search_log.json`、可选导入的 RIS：`D:/Desktop/PNJL论文/05_阅读笔记与复现计划/化学冻结线与Kpi_20261017/`
- 本机 Zotero 已有 S04（item `2RTQWGXG`）和 S07（item `4QU4JZ5G`）。S04 附件与新下载 PDF 的 SHA-256 一致；S07 直接复用附件。

| ID | 文献与原文入口 | 核对位置及用途 |
| --- | --- | --- |
| S01 | Cleymans et al., PRC **73**, 034905 (2006), [hep-ph/0511094v2](https://arxiv.org/abs/hep-ph/0511094), [DOI](https://doi.org/10.1103/PhysRevC.73.034905) | PDF p.3, Eqs.(1)–(2)：现有冻结线及五个参数的误差 |
| S02 | Andronic et al., Nature **561**, 321 (2018), [1710.09425v2](https://arxiv.org/abs/1710.09425), [DOI](https://doi.org/10.1038/s41586-018-0491-6) | PDF p.6 Fig.3 图注：完整能量参数化；p.4：ALICE 中快度拟合 |
| S03 | Bazavov et al. (HotQCD), PLB **795**, 15 (2019), [1812.08235v2](https://arxiv.org/abs/1812.08235), [DOI](https://doi.org/10.1016/j.physletb.2019.05.013) | 摘要及 Sec.4.3：`nS=0, nQ=0.4 nB` 的交叉转变线 |
| S04 | Borsányi et al., PRL **125**, 052001 (2020), [2002.02821v1](https://arxiv.org/abs/2002.02821), [DOI](https://doi.org/10.1103/PhysRevLett.125.052001) | Eqs.(1),(3),(6)；PDF p.11 Table III：有相关误差的公开参考带 |
| S05 | Bazavov et al., PRD **93**, 014512 (2016), [1509.05786v2](https://arxiv.org/abs/1509.05786), [DOI](https://doi.org/10.1103/PhysRevD.93.014512) | PDF p.4 Table I、Eq.(8)：`Tf,0=147(2) MeV` 的条件及曲率上界 |
| S06 | Borsányi et al., PRL **113**, 052301 (2014), [1403.4576v2](https://arxiv.org/abs/1403.4576), [DOI](https://doi.org/10.1103/PhysRevLett.113.052301) | PDF p.3：`Tf ≲148 MeV`，允许两倍标准差时上界为 151 MeV；不是一条完整冻结线 |
| S07 | Borsányi et al., PRL **111**, 062005 (2013), [1305.5161v1](https://arxiv.org/abs/1305.5161), [DOI](https://doi.org/10.1103/PhysRevLett.111.062005) | PDF p.5：早期方法与 `Tf ≲157 MeV` 上界；Zotero 既有文献 |
| S08 | Bazavov et al., PRL **109**, 192302 (2012), [1208.1220v1](https://arxiv.org/abs/1208.1220), [DOI](https://doi.org/10.1103/PhysRevLett.109.192302) | 守恒荷涨落作为温度、化学势探针的方法来源 |
| S09 | STAR, PRC **96**, 044904 (2017), [1701.07065v2](https://arxiv.org/abs/1701.07065), [DOI](https://doi.org/10.1103/PhysRevC.96.044904) | Au+Au BES；`|y|<0.1`，横动量外推、衰变修正和 Fig.24 的 K/pi |
| S10 | NA49, PRC **77**, 024903 (2008), [0710.0118v2](https://arxiv.org/abs/0710.0118), [DOI](https://doi.org/10.1103/PhysRevC.77.024903) | PDF p.4 Table I：Pb+Pb 中快度与全相空间的同实验比较 |
| S11 | Andronov, Kuich, Gaździcki, Universe **9**, 106 (2023), [2205.06726v5](https://arxiv.org/abs/2205.06726), [DOI](https://doi.org/10.3390/universe9020106) | PDF p.5 Fig.1：用户投影片所引综述的 Au+Au/Pb+Pb 分组及 `y≈0` 定义 |
| S12 | Bluhm & Nahrgang, EPJC **79**, 155 (2019), [1806.04499v1](https://arxiv.org/abs/1806.04499), [DOI](https://doi.org/10.1140/epjc/s10052-019-6661-3) | Fig.2、Sec.IV：奇异观测量给出的冻结温度可高于轻味涨落结果 |

## 可比较的冻结条件

### 3.1 当前基线：Cleymans 2006

S01 的式(1)–(2)，以下 `T,muB` 使用 GeV，`z=sqrt(sNN)/GeV` 为无量纲数值：

\[
T_f=a-b\mu_B^2-c\mu_B^4,\qquad
\mu_B=\frac{d}{1+e\sqrt{s_{NN}}}.
\]

`a=0.166±0.002 GeV`，`b=0.139±0.016 GeV^-1`，`c=0.053±0.021 GeV^-3`；`d=1.308±0.028 GeV`，`e=0.273±0.008 GeV^-1`。因此分母必须把能量按 GeV 代入。原文拟合没有纳入 200 GeV 点，属于较早的数据汇总。

本仓库 [default.toml](../../../../config/physics/freezeout/default.toml) 和 [cleymans2006.toml](../../../../config/physics/freezeout/cleymans2006.toml) 当前使用这些中心值。它是必要对照，不能为了缩小比值而悄悄覆盖。

### 3.2 优先经验候选：Andronic 2018

S02 Fig.3 图注给出：

\[
T_f(z)=\frac{158.4\ \mathrm{MeV}}
 {1+\exp[2.60-\ln(z)/0.45]},\qquad
\mu_B(z)=\frac{1307.5\ \mathrm{MeV}}{1+0.288z}.
\]

极限温度的不确定度为 `1.4 MeV`，由最高五个能量点确定。注意指数是 `2.60-ln(z)/0.45`，不是 `(2.60-ln(z))/0.45`。同文的 ALICE 单能量拟合为 `156.5±1.5 MeV`，与跨能量参数化的 `158.4±1.4 MeV` 不是同一个拟合量。该参数化汇总中快度与部分全相空间的中心碰撞拟合，不能直接称为 Au+Au 专属参数。

这是优先候选，因为它同时提供温度与束流能量映射，且高能温度低于当前基线。原文没有在该图注给出所有形状参数的协方差，不能凭空为它们添加高斯误差。

### 3.3 首轮保留的 crossover 参考：HotQCD 2019（本轮不纳入）

S03 明确满足 `nS=0, nQ=0.4 nB`。以 `x=muB/T0`：

\[
T_c(\mu_B)=T_0[1-\kappa_2x^2-\kappa_4x^4],\quad
T_0=156.5\pm1.5\ \mathrm{MeV},
\quad\kappa_2=0.012(4),\quad\kappa_4=0.000(4).
\]

这是与本项目守恒荷约束最接近的格点参考，原文讨论的可靠比较区主要为 `muB ≲300 MeV`。不要误用该文 `muQ=muS=0` 的另一组曲率。它约束手征交叉转变；将其用于化学冻结还需明示 `Tf=Tc` 或另外的冻结判据。

### 3.4 首轮保留的 crossover 参考带：Borsányi 2020（本轮不纳入）

S04 原式(1)使用 **`muB/Tc(muB)`**：

\[
\frac{T_c(\mu_B)}{T_0}
=1-\kappa_2\left(\frac{\mu_B}{T_c(\mu_B)}\right)^2
-\kappa_4\left(\frac{\mu_B}{T_c(\mu_B)}\right)^4+\cdots,
\]

`kappa2=0.0153±0.0018`，`kappa4=0.00032±0.00067`；两者误差相关。式(6)报 `Tc(LT=4,muB=0)=158.0±0.6 MeV`。与 S03 的分母约定不同；重写为 `muB/T0` 的四阶展开时，必须用 `kappa4_explicit=kappa4+2*kappa2^2=0.00078818`，而且这只保证到四阶的代数等价。

首轮画参考带更适合直接使用补充材料 **Table III 的完整解析延拓输出**。例如 `(muB,Tc,error)` 为 `(0,158.08,0.47)`、`(200,153.77,1.04)`、`(300,146.16,4.62)` MeV；完整 16 行已保存在 [evidence.json](evidence.json)。表中的 `sigma` 列是过渡宽度，**不是 `Tc` 的误差**；不能把约 15 MeV 的过渡宽度当成可自由调节的冻结温度误差。

原文明确 Table III 的点间误差相关。可画所报逐点误差带，但没有联合协方差时不能把各点独立随机抽样当成物理冻结曲线。已核对 `nS=0`；本轮原文中未找到与 S03 同样明确的 `nQ/nB=0.4` 声明，因此不宣称其电荷约束已经与本项目逐项等同。使用范围仍以 `0–300 MeV` 为限。

### 3.5 低温研究分支：涨落提取的冻结条件

S05 Table I 的 STAR2.0 列报 `Tf,0=147(2) MeV`、未作所讨论修正时 `kappa2_f=-0.001(12)`；式(8)给出的是 **上界 `kappa2_f<0.011`**，不是最佳拟合中心值。拟合使用 `sqrt(sNN)≥27` 或 `39 GeV`；质子横动量范围及净质子到净重子的替代会影响结论。S06 的 `Tf≲148 MeV` 也是受当时涨落误差限制的上界，不能当成误差很小的全能区常数。

可研究低温锚点与小曲率是否改善高能幅度，但必须标为有观测量限制的候选。若采用 `Tf,0∈[145,149] MeV, kappa2_f∈[0,0.011]` 做包络，正曲率区间是额外的单调下降先验，不是原文联合置信区间。不能将这个包络直接外推至 horn 能区并标成“格点预测”。

还需保留 S12 的反向证据：净 kaon 涨落与奇异反重子/重子产额比给出的高能冻结温度高于轻味涨落提取。低温轻味方案是否适用于 K+/pi+ 的共同冻结，需要检验，不能仅挑选能压低比值的文献。

## 4. 对既有结果的复查与 horn 区间

原始远程来源为 [Actions 37172320742](https://github.com/w5851/Julia_RelaxTime/actions/runs/37172320742)，源码 `c22f5e3e267bf2dce2431cf901a9dc3508d27600`。本次恢复 8 个分片，共 1813 个不同坐标，1806 个 `screened`、7 个 `gate_failed`；manifest 所列 11 个源码文件均与该提交的 blob 哈希相符。

路线为 `density_route=q0_lambda_reference`、`density_prescription=finite_window`，`omega_lower_inv_fm=1e-5`，`q_nodes=8`，`omega_nodes=64`，`mesh=64`，`diagnostic_only=true`，`production_default=false`。输入背景是保存 seed 的代数恢复。本次只读原始数值，没有新背景求解、插值补点或新密度积分。

| T / MeV | muB / MeV | 原生网格 K+/pi+ |
| --- | --- | --- |
| 140 | 0 | 0.210137 |
| 145 | 0 | 0.225227 |
| 150 | 0 | 0.240112 |
| 155 | 0 | 0.254765 |
| 160 | 0 | 0.269175 |
| 165 | 0 | 0.283367 |

这些点均为 `screened` 且 K+ 端点警告为 false；仍不升级为生产验收点。表中 165 MeV 是原生网格点，不是基线 166 MeV 的重算。源分片、关键字段、源码和文献表格保留信息见 [evidence.json](evidence.json)。

降低至 150 MeV 把比值降至约 0.24，但不意味着已经匹配实验。S11 Fig.1 的高能重离子平台大致为 0.15–0.17；`0.2` 不能作为所有能量的统一目标。正式拟合须使用原始实验表、误差、中心度和碰撞体系。

由 S01/S02 的公式直接计算，`sqrt(sNN)=7.7 GeV` 时：

| 参数化 | muB / MeV | Tf / MeV |
| --- | --- | --- |
| Cleymans 2006 | 421.650 | 139.612 |
| Andronic 2018 | 406.359 | 138.428 |

因此高能极限温度降低约 8 MeV，并不意味着 horn 附近也降低同样幅度。horn 大约位于 `sqrt(sNN)=7–8 GeV`、`muB≈0.4 GeV`，已超出上述格点参考的 `muB≤0.3 GeV` 范围。低化学势涨落曲率也不能直接控制这一段。局部极值还取决于冻结轨迹方向：沿当前 BQS 约束消去其他化学势后，`dR/dz=(∂R/∂T)Tf'(z)+(∂R/∂muB)muB'(z)`。高度、峰位和尖锐程度应分别报告。

## 实验快度口径

### 比较对象

- S11 的纵轴是中快度的正电荷比值 `(dN_K+/dy)/(dN_pi+/dy)`，不是总 K/pi 或电荷平均比值。图中 AGS/RHIC 的 Au+Au 与 SPS/LHC 的 Pb+Pb 分开标注；最尖的 horn 证据包含 NA49 Pb+Pb，不能统称 Au+Au。
- S09 的 Au+Au BES 数据采用 `|y|<0.1`。公布的 `dN/dy` 已进行探测效率修正并外推未测横动量区域；不是只数探测器接受到的部分原始粒子。pion 弱衰变贡献有专门修正，强共振衰变口径需保留。
- 本项目 [ChargedGBUResearchWorkflow](../../../../src/models/workflow_apps/ChargedGBUResearchWorkflow.jl) 的合同是 `fixed_quark_only_gbu_partial_yield`。q0 参考的 `shell` 在上述源码提交中使用 `q^2/(2pi^2)` 及频率/相移积分，再积分局域动量；没有碰撞冻结超曲面、纵向流或快度选择，也没有完整强子共振衰变链。`partial_yield` 指模型中的产额定义，不是实验快度截取。

窄态极限下，均匀静止介质有 `n_i=g_i/(2pi^2)∫q^2 f_i dq`；三动量积分完整并不等于已经计算整场碰撞。映射到最终观测量需冻结面、流场和相应的衰变/接受度处理。共享温度与化学势的 boost-invariant 中央区域，在适用的热发射近似下可以满足 `dN_i/dy≈(dV_eff/dy)n_i`，公共体积在比值中消去；中快度热模型拟合因此并非天然不可比。

### 方向没有通用结论：原始实验反例

S10 Table I 给出同一 Pb+Pb 实验的全相空间和中快度产额。用表中中心值直接相除：

| 固定靶束流能量 | 中快度 K+/pi+ | 全相空间 <K+>/<pi+> |
| --- | --- | --- |
| 20A GeV | 16.4 / 72.9 = 0.22497 | 40.7 / 190 = 0.21421 |
| 30A GeV（约 sqrt(sNN)=7.6 GeV） | 21.2 / 83.0 = 0.25542 | 52.9 / 241 = 0.21950 |

这里只比较原文中心值，没有忽略相关误差后做显著性或拟合判断。原文还说明正 pion 全产额由负 pion 及电荷比估算。这个 Pb+Pb 例子不能量化 Au+Au 的接受度效应，但足以否定“局部测量必然给出更小比值”的通用解释。

一个简单解析例子也说明方向取决于种类：若两种粒子的 `dN/dy` 是中心单高斯，则 `R(y=0)=R_4pi * sigma_pi/sigma_K`。kaon 的快度分布更窄时，中快度比值更高。NA49 实际使用双高斯；不能拿单高斯例子代替对其数据的拟合。

### 可用于论文讨论的表述边界

> 本文计算的是均匀介质在给定热力学条件下的电荷分辨介子部分产额数密度比；实验则给出特定中心度及中快度窗口内、经效率与横动量外推处理的最终强子产额比。二者的定量对应还依赖冻结条件、纵向热力学分布、流场及共振衰变口径。因而这里的比较应视为热力学趋势检验；仅凭实验标注 y≈0，尚不能确定差异的方向或将整体偏高归因于快度选择。

这一段是讨论候选，不是已经验证的机制归因。后续应分别检验：冻结温度和 `muB(sqrt(sNN))`；当前 quark-only BQS 对 `muQ,muS` 的影响；最终 pion/kaon 的衰变贡献；必要时再做纵向发射和接受度映射。加入衰变时注意避免与已有相移谱重复计数，不能只给分母乘一个经验常数。
