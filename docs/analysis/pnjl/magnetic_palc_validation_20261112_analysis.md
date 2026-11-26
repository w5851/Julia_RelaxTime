# 磁场 PNJL 正确性核验与 PALC 性能基准（2026-11-12–13）

更新：2026-11-13。**正确性核验已获作者接受，同机性能基准已完成**。独立文献、公式与实际入口的核对支持当前所检磁场计算流程；原学位论文的转折偏差保留既有结论，按作者意见结束追查。新增基准表明主要成本位于 Landau 热积分和磁真空项的求导，下一步集中优化这两部分。性能结果见[下节](#性能基准2026-11-13)。

## 范围、代码与独立性

- 工作树：`D:/w/e6bf/Julia_RelaxTime`；分支：`codex/magnetic-palc-validation`。修复基于 `dc183dfa561dc79ef1f8594948b1a1e705679996`，目前源码、测试和说明尚未提交/合入。最终源码/参考脚本及Project/Manifest共14个文件与远端哈希一致，记录在 `literature_validation_closure.json.code_context`。
- 数值计算：quote-remote 的 `D:/CodexCompute/magnetic-palc-20261112-v1`，Julia 1.12.5、单线程、根 Project/Manifest；本地编辑、提取 PDF 文字、传输和汇总。
- 当前物理范围：三味均场、正磁场、`xi=0`、相同夸克化学势；有限温磁场入口、MFIR 默认路线，以及显式 smooth 文献路线。精确零温使用独立解析积分，作为有限温核的物理交叉证据。
- `MagneticReference.jl` 和文献脚本不调用 `Models`、`PNJLCore` 或生产积分生成参考结果。独立核以 loggamma 凝聚式、Euler–Maclaurin 磁真空势、精确填充海减除及分段热修正计算；生产核使用 Abel–Plana 加小参数递推、asinh 分段和 ForwardDiff。共享的只是通用数值库及文献物理模型。
- 本轮按用户指示转向新文献交叉验证。原作者未披露的程序、截断和数值设置不再作为待追查任务，也不推测作者实际使用了哪一种变体。

## 文献选择与口径

按“同一理论公式、公开参数、可直接比较的正文/表格数值”选取目标，本轮新增对照不依赖读图。来源 DOI、版本与 PDF 哈希保存在结果目录的 `literature_sources.json`。INSPIRE REST 用于来源定位；后续宽泛检索曾有 TLS 错误，选定的开放 PDF 均成功获取。本轮是定向交叉验证，不是穷尽性文献综述。

| 来源 | 采用的目标与明确约定 | 独立性 |
| --- | --- | --- |
| [Fu, PRD 88, 014009 (2013)](https://doi.org/10.1103/PhysRevD.88.014009)，arXiv:1306.5804v1，式5–11、Sec.III | 三味 smooth 真空 `fΛ²,N=10`；多项式 Polyakov 势，`T0=270 MeV`；`mu_f=0`；`|dMu/dT|` 峰位 | 独立研究组；直接检验 smooth 有限温路线 |
| [Grunfeld 等, PRD 90, 044024 (2014)](https://doi.org/10.1103/PhysRevD.90.044024)，arXiv:1402.4731v1，式3–7、Table I、Sec.III | 冷三味 NJL，MFIR，常数 G；Set 1、`eB=.01 GeV²` 的三次压力相等跃迁及两组真空质量 | 另一独立研究组；覆盖窄跃迁、Landau 占据与真空核 |
| [Costa 等, PRD 92, 036012 (2015)](https://doi.org/10.1103/PhysRevD.92.036012)，arXiv:1508.07870，Sec.III 的双跃迁段落 | `eB=.09 GeV²,T=0,GV=0`、常数 G，两次相共存 | 与 Ferreira 系列有共同作者，作为补充，不能另算完全独立实现 |
| [Ferreira、Costa、Providência, PRD 97, 014014 (2018)](https://doi.org/10.1103/PhysRevD.97.014014)，已发表开放版，Tables I–II | MFIR，对数 Polyakov 势，`T0=270 MeV`，常数 G 与 G(B)；选取八个轻/奇异味 CEP | 检验 IMC、有限温密度和多临界点结构 |
| [Shao、He、Gao, PRD 100, 014020 (2019)](https://doi.org/10.1103/PhysRevD.100.014020)，arXiv:1907.02201；高雪艳2021第五章 | smooth，`T0=210 MeV`；早期选定的五个图5-3转折与强场真空质量 | 同一课题组/同一研究链；保留原比较，不能当新增独立支持 |

Ferreira 2014（PRD 90, 016012，arXiv:1406.3608）已筛选阅读；由于参数印刷口径问题，本轮未把它列作新增数值验收来源。

标准参数为 `Λ=602.3 MeV, ml=5.5 MeV, ms=140.7 MeV, GΛ²=1.835, KΛ⁵=12.36`。Grunfeld Set 1 按原文改为 `Λ=631.4 MeV, ms=135.7 MeV, KΛ⁵=9.29`。Fu 的 Polyakov 参数为 `(a0,a1,a2,a3,b3,b4)=(6.75,-1.95,2.625,-7.44,.75,7.5)`，场强单位使用其明确给出的 `m_pi=135 MeV`。没有通过调参拟合目标值。

### 公式与实现

令 `b_f=|Q_f|eB`，电荷比例 `(2/3,1/3,1/3)`；内部采用 fm 自然单位。

| 理论关系 | 核对结果 |
| --- | --- |
| `E_fn=sqrt(Mf²+pz²+2n b_f)`，`alpha_n=2-delta_n0` | 能谱、最低层简并度、Landau 测度正确 |
| `Ωvac_f=-Nc b_f/(2π²) Σ alpha_n ∫₀∞dpz fΛ² E_fn` | smooth 真空的偶积分因子、颜色因子和截断平方正确；自定义 Λ 传播已修复 |
| `ΩT_f=-T b_f/(2π²) Σ alpha_n ∫₀∞dpz(log Z+ + log Z-)` | Polyakov 多项式含色迹，不再额外乘 Nc；粒子/反粒子的 Φ 交换正确 |
| `rho_f=-∂Ω/∂mu_f` | 分布积分与势导数一致，单位 `fm^-3`；`rhoB=sum(rho_f)/3` |
| `Mf=mf-4G(B)phi_f+2K phi_j phi_k`；`Ωchi=2GΣphi_f²-4K phi_u phi_d phi_s` | 质量、KMT 符号和手征势一致 |
| `Ωmag=-Nc b_f²/(2π²)[ζ′(-1,x)-(x²-x)log(x)/2+x²/4]`，`x=Mf²/(2b_f)` | MFIR 势与独立 Euler–Maclaurin 式、loggamma 导数一致 |
| `∇phi Ω=JMᵀ g`，`g_f=∂(Ωvac+ΩT)/∂Mf-phi_f` | 链式法则已核对；新增独立凝聚自洽接纳条件，排除 `JM` 奇异伪驻点 |

## 新文献定量结果

### Fu：smooth 有限温转变

| eB | 文中峰位 T (MeV) | 独立复现 T (MeV) |
| --- | ---: | ---: |
| 0 | 224 | 224.420018 |
| `20 m_pi²=.3645 GeV²` | 240 | 239.627969 |

分段阶数/层数/温度差分依次 `(32,96,.5 MeV)`、`(48,160,.25)`、`(64,240,.125)`；最后一次峰位变化分别 `.00101/.00145 MeV`。有限场峰态代入生产 smooth 核，`p_num=128→256→512` 的凝聚差为 `3.16e-6→6.69e-11→3.87e-15 fm^-3`。**smooth 真空仍须显式升阶**，不能把 MFIR 热积分的改善解释为它已在所有默认预算下收敛。此结果较好复现文献明示的磁催化峰位移动。

### Grunfeld / Costa：零温多次相共存

| 文献与条件 | 文献位置 (MeV) | 独立复现 (MeV) |
| --- | ---: | ---: |
| Grunfeld Set 1，`.01 GeV²`，第一跃迁 muq | 335.3 | 335.378588 |
| 同上，第二跃迁 muq | 342.3 | 342.156881 |
| 同上，第三跃迁 muq | 345.4 | 345.219425 |
| Costa，`.09 GeV²`，第一跃迁 muB | 1061 | 1061.492866 |
| 同上，第二跃迁 muB | 1083 | 1083.613434 |

独立公式直接积分 T=0 填充海，联立两相三味质量 gap 与压力相等求位置。Grunfeld 两组 B=0 真空 Mu 为 `335.510963/367.648250 MeV`，对应原文约 `336/368 MeV`。

最初固定 Mu 的目录漏掉第二个窄跃迁；局部反解补回它，两种起点给出的共存 muq 相差约 `6e-10 MeV`。`literature_cold_refined.json` 中较远路径保存有162次失败，其整段折线不用于证明完备性。

本轮又对五组共存相检查了凝聚方向 Hessian：三种差分步长下最小特征值均为正，最小约 `.0470`（自然单位）。在各自 `muq±.01 MeV` 重新解两相，`Ω高密-Ω低密` 从正变负；已知候选中没有更低势根。这支持相共存与稳定性交换的解释。个别补充种子失败已保留，未声称穷尽未知分支。Grunfeld 最大 muq 差 `.181 MeV`、Costa 最大 muB 差 `.614 MeV`；这些是小的剩余文献差异，不能全部当作打印取整误差。

### Ferreira：八个有限温 CEP

定位方法为固定一个凝聚量反解平滑曲线 `mu(phi_k,T)`，求 `dmu/dphi_k=d²mu/dphi_k²=0`；轻味使用 k=u，奇异味使用 k=s。

下表三元组为 `(T/MeV, muB/MeV, rhoB/rho0)`，场强单位 `GeV²`。

| 目标 | 文献 | 独立复现 |
| --- | --- | --- |
| 轻味，0.3，常数 G | (192, 674, 3.54) | (191.877, 675.179, 3.54243) |
| 轻味，0.6，常数 G | (214, 692, 7.22) | (213.630, 694.275, 7.23219) |
| 轻味，0.3，G(B) | (177, 627, 2.63) | (177.289, 625.354, 2.63829) |
| 轻味，0.6，G(B) | (171, 535, 3.90) | (171.500, 533.409, 3.89752) |
| 奇异味，0.3，常数 G，第一端点 | (62, 1330, 6.67) | (62.126, 1330.474, 6.66357) |
| 奇异味，0.6，常数 G | (124, 1234, 12.10) | (123.683, 1234.337, 12.07254) |
| 奇异味，0.3，G(B)，第一端点 | (48, 1193, 5.40) | (47.729, 1193.578, 5.41541) |
| 奇异味，0.3，G(B)，第二端点 | (18, 1539, 10.75) | (18.132, 1539.446, 10.75569) |

控制组 `(积分阶数,n_max,delta_phi)` 为 `(32,120,.012)`、`(48,180,.006)`、`(64,240,.003)`。最后一步最大位移 `ΔT=.00533 MeV, ΔmuB=.00324 MeV`。相对文献，最大 T 差约 `.50 MeV`、muB 差 `2.28 MeV`、密度相对差小于 `.32%`。数值细化远小于剩余文献偏差，因此结论是定量近似交叉支持，未宣称全部差异由取整解释。

选择一个轻味 CEP 和最低温奇异味 CEP，再比较 `T_CEP±.5 MeV` 的邻域：下方各出现两次 `dmu/dphi` 变号，上方均为零次，符合回弯对在端点合并的结构。不是只在取整表值代入 gap 后判断。

八个 CEP 态代入生产核后，凝聚差最大 `6.74e-12 fm^-3`，密度差最大 `6.87e-15 fm^-3`。这里只复现选定的八点，没有声称两张表的所有端点已重算。

## 已修复的问题与数值证据

| 问题 | 本轮修复 | 直接证据 |
| --- | --- | --- |
| 低温整区间 Gauss 热积分漏分辨费米面，制造细回弯 | 每层用 `t=asinh(pz/m_n)`，按 `E=|mu|+(-12,-6,0,6,12,40)T` 分段；势与密度共用节点；节点保留 AD | 最差代表态旧 p128 密度误差约 `.0882 fm^-3`；新每段 p64 误差约 `6.3e-16 fm^-3` |
| 小质量下默认 zeta256 的窄层误差 | 精确递推 `ζ′(-1,x)=ζ′(-1,x+1)-x log x` | `x=7.60e-5` 的凝聚导数差由 `1.60e-5` 降至约 `2.94e-12 fm^-3`；`1e-8≤x≤10` 抽样对 loggamma 核验 |
| 仅势梯度接纳会保留非自洽伪根 | 新增 `magnetic_condensate_residual`；与五场 residual 一并检查 | 旧反例 `g≈5.18`、势梯度约 `7.5e-11`、`sigma_min(JM)≈1.4e-11`；当前方程重新求出的伪驻点被公共入口拒绝 |
| legacy 真空使用全局 Λ，而非模型 Λ | 真空调用传入 `params.Λ_inv_fm/N_color`，自动 pz 上限也使用模型 Λ | 800 MeV profile 测试；修复前同态真空差 `47.21 fm^-4` |
| legacy 自动层数只控制热尾 | smooth 路线要求显式、已检查收敛的 `n_max` | 缺省求值明确报错；热尾预算不再冒充真空求和预算 |
| 低层 T=0 原先返回零热密度 | 要求有限正 T，精确零温留给解析填充海 | T=0、负层数错误路径及无质量最低 Landau 层解析密度单测 |

热项 `p_num` 现在表示**每段节点数**；legacy 真空仍为整区间节点数。默认 MFIR 配置值、物理参数、生产多初值策略未更换。新实现对热层保留40T尾，自动层数仍按配置30T估计；它们是可收敛的数值控制，需结合实际消费者检查。

### 完整低温 PALC 路径

实际生产方程加单高端锚点 PALC，`T=5 MeV, muB=500–1700 MeV`；`p_num=128,n_max=64,pz_max=40,zeta_num=256`，`ds=.002,dsmax=.006,Newton tol=1e-10,linesearch=true`。参考路径独立采用凝聚方程及更细步长；逐点固定态核验用参考阶数96、nmax128。

| eB (GeV²) | 生产路径点数 | 回弯数（与参考相同） | 最大独立凝聚 gap (`fm^-3`) | 预选根最大状态距离 |
| --- | ---: | ---: | ---: | ---: |
| .05 | 448 | 16 | 1.10e-10 | 4.92e-9 |
| .155 | 530 | 14 | 1.24e-10 | 7.72e-9 |
| .4 | 628 | 8 | 1.54e-10 | 4.02e-9 |

共1606个保存顶点全部核对独立凝聚 gap、势梯度和密度；密度差最大 `1.13e-14 fm^-3`。原始 PALC 顶点在低密真空端有约 `2.39e-10` 的微小负 Polyakov 场，原 `physical=false` 标签保留；投影到零边界后的独立残差仍不超过 `1.54e-10`、密度变化不超过 `7e-25 fm^-3`。这是原始续接顶点的边界数值误差，不能把它们全部称作公共 solver 已接纳的状态。

在 `muB=[872.08,900,1000,1050,1100,1150,1200,1500] MeV`，三场强的根数依次为：

- `.05`：`1,1,1,3,3,1,1,1`；
- `.155`：`3,3,7,5,1,1,1,1`；
- `.4`：`3,3,3,1,3,3,1,1`。

共52个目标根与独立路径交叉种子重新校正的根双向一一匹配，独立校正无失败。后两场强另做每段p64全路径，根与p128最大差 `4.37e-14`。`.155,muB=872.08` 原旧整区间p512的五根现为三根，多出的数值回弯对消失。`.05` 完整生产路径的旧缺口已补齐。本轮未把诊断计时用作性能结论，也未把单锚点 PALC 设为默认扫描策略。

### 默认自动层数与实际求解

先前34个参考根被实际 solver 接纳，独立 residual 最大 `6.85e-12`；由于初值已很准，这只检验接纳与核值。

为覆盖真正求解过程，又选14个低温根（包含七根截面）及三个240MeV高温状态，扰动全部初值后使用公共入口与自动 `n_max`。17次均在2–3次迭代后收敛，独立残差最大 `3.15e-11`、密度差最大 `2.93e-15 fm^-3`；实际层数低温3–13、高温98–389。高温外部状态的197.33单位仅用作初值，最终按项目197.3269804单位重新求解，不把其状态漂移当成来源误差。

## 旧外部锚点与测试

外部 `pnjl_mag@e1fc81...` 的原作者数值和SHA保留。本轮尝试获取该固定版本原文件返回404，因此没有宣称重跑原作者程序；来源依据是既有固定SHA源码审计/外部回放。采用重建的旧离散规则及独立公式隔离误差：

| 目标 | 本轮归因 | 验收处理 |
| --- | --- | --- |
| 固定态 p16 Ω | 重建旧值与原值差 `9.56e-9`；相对独立收敛值误差 `1.2905e-5 fm^-4`；旧p32差 `4.25e-9`，p64差 `5.17e-12` | 原 expected/rtol/atol 保留，另加来源离散预算 `1.5e-5 fm^-4`；同时测试当前积分升阶一致性 |
| T240、eB=.8 状态 | 原态与独立解距离 `2.9288e-6`；仅p升阶仍差 `2.9064e-6`，仅zeta升阶仍差 `1.0567e-7`；两项升阶后差 `1.85e-11` | 原态、Ω及 `state_atol=1e-7` 保留，另加来源状态分辨率 `3.0e-6` |
| T240、eB=.2/.4 | 原态距独立解 `1.53e-9/1.20e-8` | 原状态容差足够，来源追加预算为零 |

这两个已解释漂移不通过替换外部数值消除。新列及解释已写入原 target/provenance 文件，历史同节点观察保留。来源离散预算用于原外部粗结果的比较；当前内核精度另由解析极限、公式导数及细化检查约束。

远端测试先整组 **220/220** 通过；随后增强伪驻点的正质量/Polyakov合法域断言，单独重跑model **55/55** 通过。7个相关文件共 **222项**（integrals23、thermodynamics51、model55、scan21、smoke30、external kernel5、external anchors37）；另开nightly磁场固定态范围 **60/60** 通过。没有修改内部 smoke/nightly 基线，也没有运行整个仓库或 GitHub CI。

新增伪驻点断言最初直接套用旧种子而失败：当前积分下种子势梯度为 `7.64e-6`，种子并非新方程根。测试已改为先解当前驻点方程，再检验“小势梯度、大凝聚gap、公共入口拒绝”，最终通过。早期测试包装器的 `Constants_PNJL` 加载错误亦已修复；上面的最终通过记录均来自实际成功运行。

## 原学位论文比较的保留结论

高雪艳2021论文PDF哈希 `421109ff06e8000d4bdd904eb27aaa1081f8692338480e1db0db92887086a904`；第五章图5-3为PDF83/印刷69。早期已用原生矢量折线校准五个预选转折：额外x刻度最大残差约 `.104 MeV`，半线宽约 `1.09 MeV`。采用文献 smooth 公式与已核对参数，结果为：

| 图5-3目标 | 文献矢量 muB (MeV) | 收敛 smooth muB (MeV) | 差 (MeV) |
| --- | ---: | ---: | ---: |
| 第一低密极小 | 1124.974 | 1154.031 | +29.058 |
| 随后的极大 | 1164.767 | 1192.103 | +27.336 |
| 第二低密极小 | 1149.844 | 1176.525 | +26.681 |
| 随后的极大 | 1169.430 | 1194.540 | +25.110 |
| 主回弯极小 | 1051.606 | 1066.651 | +15.045 |

独立阶数48→96影响低于 `1.2e-12 MeV`；nmax32→64影响至多 `.00498 MeV`，64→96至多 `1.06e-5 MeV`；继续到128的固定态凝聚变化低于 `4.14e-10 fm^-3`。生产smooth整区间p2048与这些转折的独立核最大密度差 `6.56e-9`、凝聚差 `2.45e-9 fm^-3`。因此该15–29MeV差异远大于本轮数值变化及读图分辨率。

上游 Ferreira 原文为 `a=.0108805`，学位稿/期刊对应位置印作 `.108805`；十倍变体把五点进一步向高muB移 `36–64 MeV`、强场质量变为约1171/1136MeV，不能复现正文431/400MeV。`KΛ⁵=12.56` 相对标准/原始磁场期刊的12.36也将点移得更高。完整smooth强场质量431.259/399.591MeV与正文相符。有限Landau截断的早期受控变体表明它可能产生明显位移，但没有同时匹配五点，也不提供作者实际算法的信息。

新文献的独立支持，加上两套数值公式对同一模型的收敛一致，显著提高了我们对当前所检计算流程的置信度。剩余图5-3偏差仍记为**按文献明示公式与可核对参数存在未解释差异**；不会把论文错误或未披露的某项实现当成已证实事实。

## 复现与后续

关键脚本在 `scripts/analysis/pnjl_magnetic/`，关键结果在
`data/outputs/results/pnjl/magnetic_palc_validation/20261112/`：

- `crosscheck_literature.jl` → `literature_cold.json`、`literature_cold_refined.json`、`literature_cep_preview.json`、`literature_cep.json`。
- `crosscheck_fu2013.jl` → `literature_fu2013.json`。
- `check_corrected_entrypoint.jl` → `corrected_entrypoint.json`，含修复核值、34根接纳、17个扰动初值及自动层数检查。
- `compare_single_anchor_palc.jl` → `corrected_palc_p64.json` / `corrected_palc_p128.json`；独立参考为 `resolved_mfir_o48_fine_v2.json`。
- `close_literature_validation.jl` → `literature_validation_closure.json`，保存1606顶点审计、双向根匹配、零温稳定性和CEP邻域。
- `audit_external_anchor_drift.jl` → `external_anchor_drift.json`，保留来源离散的独立归因矩阵。
- 原论文公式/转折审计：`thesis_equation_audit.json`、`thesis_landau_cutoff_audit.json`、`kernel_closure_audit.json`、`thesis_fig5_3_vector.json`；原先粗积分曲线与图仍保留其历史诊断意义，本轮没有重画或晋升新图。

典型命令在远端隔离根目录运行，输入文件复制到 `outputs/`：

```text
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/crosscheck_literature.jl cold outputs/literature_cold.json
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/crosscheck_literature.jl cold-refine outputs/literature_cold_refined.json
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/crosscheck_literature.jl cep-preview outputs/literature_cep_preview.json
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/crosscheck_literature.jl cep outputs/literature_cep.json
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/crosscheck_fu2013.jl outputs/literature_fu2013.json
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/check_corrected_entrypoint.jl outputs outputs/corrected_entrypoint.json
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/check_corrected_entrypoint.jl outputs outputs/corrected_entrypoint.json --consumer
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/compare_single_anchor_palc.jl outputs/corrected_palc_p128.json --T=5 --p-num=128 --pz-max=40 --n-max=64 --ds=.002 --dsmax=.006 --max-steps=1800 --anchor-side=high --eB=0.05,0.155,0.4 --targets=872.08,900,1000,1050,1100,1150,1200,1500 --jacobian=forward --linesearch=true --catalog=false --newton-tol=1e-10 --newton-iterations=28
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/close_literature_validation.jl outputs outputs/literature_validation_closure.json
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/audit_external_anchor_drift.jl outputs/external_anchor_drift.json
julia --startup-file=no --project=. scripts/analysis/pnjl_magnetic/run_focused_tests.jl
```

nightly 复核在 PowerShell 设置 `$env:REGRESSION_MAGNETIC_SCOPE="nightly"` 后，运行同一 runner 并传入 `tests/regression/pnjl/test_magnetic_fixedpoint_regression.jl`。

API 的每段节点语义、正温域、显式smooth层数和凝聚残差已同步到 magnetic 主题文档。上轮的正确性工作为下面的性能比较提供固定目标；本次交付为基准与优化点分析。

本任务使用此既有分析入口维护状态；primary track 和并行 PR3 台账沿用各自任务的状态。磁场修复及本次基准仍在本隔离分支，后续源码提交/合入与算法优化继续推进。

## 性能基准（2026-11-13）

### 口径与复现

- 主机为 quote-remote，Intel Core i5-14400F；Julia 1.12.5、Julia/BLAS 各一线程，固定 CPU0，按顺序运行两种模式。
- 磁场使用上述已核验的工作树源码；无磁场使用当前 main 的提交 `44fe531d534ef89a055c0f226896db24c1de152d`，含 PR #346 的三次预测、复用校正器与 residual-first Newton。两个独立源码目录共用相同 Project/Manifest。
- 入口为 [benchmark_magnetic_palc.jl](../../../benchmark/pnjl/benchmark_magnetic_palc.jl)。每案例先跑一整轮预热，再跑三轮并轮换案例顺序；下表为后三轮中位数。无磁场每轮含五次路径求解；磁场每轮一条完整路径。批前 GC、核验和写盘在计时外，计时内部的 GC 保留。
- 结果只新增 [magnetic.json](../../../data/outputs/results/pnjl/magnetic_palc_validation/20261113_performance/magnetic.json) 和 [nonmagnetic.json](../../../data/outputs/results/pnjl/magnetic_palc_validation/20261113_performance/nonmagnetic.json)，合计约102 KiB；包含各轮时间、分配量、调用数、有效预算、根匹配、微基准、采样摘要和源码哈希。两种模式均 `complete=true`、`source_unchanged=true`，本地源码与运行哈希已核对。

磁场求解时间包含单初值高 μB 锚点、双向 PALC、八个目标化学势的全部交点校正和这些根的 Ω/密度/残差重建。另列所有 PALC 顶点的诊断输出重建时间。无磁场采用原有 `solve_path + selected_path_results`，包含全部候选求解与选中态的完整热力学重建。两组窗口、积分预算、输出字段及物理解集不同，因此总时间用于了解现有链路成本；同温同化学势的残差/Jacobian 参考另列于后。

计数时发现 BifurcationKit 在 `bothside=true` 下会 `deepcopy` 迭代器，原计数器只看见正向一侧。基准进程对两个诊断计数器共享复制后的计数对象，完整记录双向工作；数值状态和求解策略仍分别复制。这修正了工作量统计，原总耗时本来已覆盖两个方向。早期漏计结果未作为最终输出。微基准中一次 `SVector` 被传给要求可变向量的 FD 辅助函数也已修正，最终两种模式均实际运行成功。

```powershell
# 在 quote-remote 中先固定父 PowerShell 的 CPU，再顺序运行两个源码目录。
[Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity = [IntPtr]1
$env:OPENBLAS_NUM_THREADS = '1'
$env:PALC_BENCH_CPU_AFFINITY = 'CPU0; inherited parent affinity=1'

# 磁场目录；MAGNETIC_REFERENCE_RESULT 指向已有 p128 验证结果。
$env:MAGNETIC_REFERENCE_RESULT = 'D:\CodexCompute\magnetic-palc-20261112-v1\outputs\corrected_palc_p128.json'
julia --startup-file=no --threads=1 --project=. benchmark/pnjl/benchmark_magnetic_palc.jl magnetic magnetic.json 3

# main 对应目录，使用同一个 benchmark 脚本。
julia --startup-file=no --threads=1 --project=. benchmark/pnjl/benchmark_magnetic_palc.jl nonmagnetic nonmagnetic.json 3
```

### 完整路径的实际成本

磁场统一为 `T=5 MeV, μB=500–1700 MeV, p128/段, n_max=64, zeta256, pz_max=40`；`ds=.002, dsmax=.006, Newton tol=1e-10`。目标点与上轮相同。无磁场为 `p24×t8, ξ=0, ds=.01, dsmax=.05, Newton tol=1e-8`。

| 场景 | 目标数 / 目标根数 | PALC 顶点 | 路径求解时间 | 其中 PALC | 另加全顶点诊断重建 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 无磁场，T120，μq=0:15:450 MeV | 31 / 31 | 59 | 0.05883 s | 0.04149 s | — |
| 无磁场，T120，μq=296:.25:318 MeV | 89 / 111 | 23 | 0.10410 s | 0.01569 s | — |
| 无磁场，T200，μq=0:15:450 MeV | 31 / 31 | 47 | 0.04864 s | 0.03092 s | — |
| 磁场，eB=.05 GeV² | 8 / 12 | 448 | 10.44276 s | 10.06660 s | 1.22933 s |
| 磁场，eB=.155 GeV² | 8 / 22 | 530 | 5.56814 s | 5.30523 s | .60921 s |
| 磁场，eB=.4 GeV² | 8 / 18 | 628 | 4.11208 s | 3.94716 s | .45935 s |

磁场三组的三次时间范围分别为 `10.38479–10.55740 s`、`5.55891–5.57243 s`、`4.10312–4.21701 s`。PALC 占求解时间约95–96%；锚点为 `.19892/.08318/.04551 s`，目标校正为 `.15243/.15902/.10320 s`。每组锚点只有一次尝试，目标校正分别12/22/18次，失败/恢复均为零。

| eB (GeV²) | PALC 普通残差调用 | 嵌套 AD Jacobian 调用 | PALC 平均每保存顶点成本 |
| --- | ---: | ---: | ---: |
| .05 | 2590 | 816 | 22.47 ms |
| .155 | 3236 | 998 | 10.01 ms |
| .4 | 3524 | 1079 | 6.29 ms |

无磁场每顶点约 `.66–.70 ms`；磁场每顶点调用 Jacobian 约1.7–1.9次，没有随总耗时成百倍增加。强场路径的顶点和调用数更多，总时间却更低，主要原因是 Landau 热积分需要处理的有效层数更少。

三条 p128 路径与上轮根状态逐项相同，16/14/8次回弯及全部52个预选根保留，最大目标残差 `1.09e-9`。无磁场细窗口111根保留，全部案例每条路径只用一个锚点，目标校正无恢复。磁场求解分配量约29–40 MiB，三次 warm 样本内部 GC 均为零；当前优先级由算术成本决定。

### 热积分、磁真空项与求导

对 `eB=.155,p128` 的完整低温 PALC 做独立采样剖析：2447个有效栈样本中，1808个经过热积分（73.9%）、623个经过 ζ′（25.5%）、1501个经过 Polyakov 分布（61.3%，属于热积分的子集）；线性代数为1个样本。比例是调用栈的包含关系，不能相加；Windows 采样具有统计波动。

以下微基准在固定自洽态测量，时间单位为 ms。Jacobian 是五场 gap 的导数，即势的二阶导数；磁真空和热项通过原公式分项求导，分项梯度和完整 gap 的最大差为 `5.50e-16`。

| 状态 | gap 残差 | AD Jacobian | 中心 FD Jacobian | 热项 Hessian | 磁真空 Hessian |
| --- | ---: | ---: | ---: | ---: | ---: |
| T5，μB1500，eB=.05 | 2.404 | 14.906 | 26.501 | 13.481 | 1.058 |
| T5，μB1500，eB=.155 | .982 | 6.186 | 10.783 | 4.947 | 1.076 |
| T5，μB1500，eB=.4 | .513 | 3.441 | 5.528 | 2.291 | 1.090 |
| T5，μB1000，eB=.155，七根中最大 Mu 态 | .201 | 1.546 | 2.243 | .463 | 1.015 |
| T120，μB930，eB=.155 | 8.803 | 47.589 | 108.097 | 47.054 | 1.022 |
| T200，μB930，eB=.155 | 22.265 | 120.362 | 275.961 | 120.730 | 1.016 |

各微基准独立计时，分项中位数不严格相加。分项 AD 的类型和编译路径也与整体不同，表格用于识别数量级。AD 比 FD 快1.45–2.29倍，二者 Jacobian 相对范数差为 `7.93e-11–2.12e-9`。目前磁场继续采用 AD 有明确的单机依据。

磁场的 `n_max=64` 只是循环上限：T5、μB1500 时三味实际进入积分的层数之和为35/12/5，分段节点总数为15616/6016/2816。低密大质量态只有640个热积分节点，磁真空求导因而成为主要部分。

为与无磁场作同 T、μ 参考，T120/T200 取 `μq=310 MeV`。磁场按现有自动层数策略固定为148/386，热积分节点总数67328/172544；无磁场按已有测试预算24×8，三味共576个节点。

| 温度 | 无磁场 gap / AD Jacobian | 磁场 gap / AD Jacobian | 当前预算下耗时比（gap / AD） |
| --- | ---: | ---: | ---: |
| T120 | .02577 / .23863 ms | 8.80262 / 47.58923 ms | 342 / 199 |
| T200 | .02672 / .23108 ms | 22.26539 / 120.36165 ms | 833 / 521 |

这解释了高温未必是磁场路径最便宜的起点：热尾覆盖的 Landau 层数随温度增加，现有每段p128会在大量高层继续使用完整阶数。上述倍率描述现有预算与入口，等误差的最低成本需在后续积分预算优化时比较。

### 已量到的预算收益与优化顺序

| 已有精度支持的低温路径 | p128 | p64 | 提速 | 根数 / 回弯数 |
| --- | ---: | ---: | ---: | ---: |
| eB=.155 | 5.56814 s | 3.53871 s | 1.57× | 22 / 14，保持一致 |
| eB=.4 | 4.11208 s | 2.87547 s | 1.43× | 18 / 8，保持一致 |

p64 与原p128目标根最大距离 `4.37e-14`，两种预算的延拓调用数也相同。此收益适用于已经比较的两条T5路径；默认p128保留到更广消费者的导数与根集合核验完成。

后续建议按以下顺序推进：

1. **先优化热积分内的代数与导数。** `_scaled_polyakov_terms` 目前每个粒子/反粒子方向调用四次指数，节点变换又计算 `exp(t)` 和 `exp(-t)`，嵌套Dual算术放大了这些成本。优先用稳定的单指数及其幂复用、复用每层常数；进一步利用每味只依赖 `(Mf,Φ,Φbar)` 的结构，在小维度上求导后通过质量映射组装五场残差和Jacobian。保留移动节点的导数及低温尾部稳定性。
2. **处理磁真空导数。** 当前为三味、每味256个 ζ 节点反复求导，还包含 arctan 的18项展开。先缓存与质量无关的变换节点/权重；重点评估沿用上轮独立核已经交叉验证的 loggamma/digamma 导数关系，直接计算磁真空一阶/二阶导数。完整势仍可沿用当前表达式。这一项对低密路径尤其有价值。
3. **按目标误差减少积分工作量。** 两条T5路径的p64已有1.43–1.57倍实测收益；高温则应按Landau层/能量段分配阶数，减少指数衰减尾部的过度采样。验收同时看gap、Jacobian和目标根，求导过程中固定离散层数。
4. **再调PALC步长。** 当前磁场 `dsmax=.006`，明显小于无磁场 `.05`。可在平滑段扩大步长，在转折邻域保留已验证分辨率；以全部预选根、回弯及独立残差检验收益。不能用两组现有顶点数直接预言能提速多少。
5. **整理重复重建。** `calculate_magnetic_rho` 先重算Ω组件再计算密度，诊断 `_state_record` 又单独计算Ω与gap。复用质量、配置及热积分pair有明确代码依据；全顶点输出目前额外增加约.46–1.23秒。该部分属于输出效率改善，核心延拓的95–96%成本仍由前几项处理。

以上除p64预算对照与AD/FD微基准外，其余均为基于采样和代码的优化候选，尚无实现后的提速承诺。下一轮可以先做保持公式及当前预算的热核/磁真空求导优化，再复用本次三条路径和同温微基准验证性能与数值稳定性。本轮没有修改核心算法、默认扫描策略或内部数值基线，也没有重算完整相图。
