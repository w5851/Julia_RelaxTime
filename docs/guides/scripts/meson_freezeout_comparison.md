# 化学冻结线的 q=0 外推比较

入口 `scripts/analysis/relaxtime/run_charged_gbu_freezeout_comparison.jl` 默认比较
`cleymans2006` 与 `andronic2018` 的文献中心参数，在每个碰撞能量对应的 `(T, muB)`
直接求解 quark-only BQS 背景，再计算 `K_plus` 与 `pi_plus` 的 GBU 密度。
冻结线共用 [FreezeoutProfiles API](../../api/models/scans/FreezeoutScan.md)，
可用独立 TOML 加入同 family 的新参数组。

可显式选择 `ihrg2019`、`frg2024_mudep`、`frg2024_muindep` 和
`functional_qcd2026_mean`。新参数的来源和适用范围见
[补充文献比较](../../analysis/relaxtime/meson_freezeout_defaults_20261017_analysis/literature_candidates_20261017.md)。
FRG 原文未拟合 7.7 GeV，也没有施加奇异数中性；不能直接对这些 profile
沿用 CLI 的整个默认低能网格并称为文献覆盖。建议显式指定 11.5–200 GeV 的
共同实验能量，7.7 GeV 作为单列诊断。2026 方案本轮使用原文 mean 行。

该入口复用 PR #325 的 `q0_lambda_reference`；`lambda=omega+(mu1-mu2)`、
`lambda0=sqrt(lambda^2-q^2)` 仅用于 `lambda>=q`，其余相位为零。Bose 权重使用
外部 `omega`，有限窗口下限为 `1e-5 fm^-1`。`nS=0`、`nQ/nB=0.4`、`xi=0`，
不加入介子反馈。结果是局域平衡的指定介子密度比，没有实验接受度或共振衰变链。

```powershell
julia --project=. scripts/analysis/relaxtime/run_charged_gbu_freezeout_comparison.jl `
  --profiles cleymans2006,andronic2018 --sqrts-list 7.7,11.5,19.6,200 `
  --resolution standard --output data/outputs/results/relaxtime/meson_density/example_defaults
```

不传 `--sqrts-list` 使用脚本中固定的能量网格，在 horn 邻域加密。能量单位为 GeV；
点按能量降序延续背景 seed，每个参数化重新开始。`--profile-dir` 可指定冻结线
配置目录。`--config` 控制现有 charged-GBU 背景合同；其有限-q生产数值设置不会
被本分析的谱积分分辨率修改。

| resolution | q nodes | mesh | cut / tail nodes | omega nodes | qmax (fm^-1) |
| --- | ---: | ---: | ---: | ---: | ---: |
| screening | 8 | 64 | 32 / 32 | 64 | 8 |
| standard | 16 | 128 | 64 / 64 | 96 | 8 |
| fine | 32 | 256 | 96 / 96 | 128 | 8 |

这些是本 q0 比较的分辨率档位，不能仅凭档位名称宣称数值收敛。需用单独批次的
相同能量比较检查排序稳健性。有限窗口结果也不等于零下限极限或完整有限-q生产验收。
该工具不会改变 `Models.run_charged_gbu_freezeout_scan` 的现有生产默认。

输出包括 `freezeout_ratios.csv`、逐点 `points/<profile>/energy_NNN.json`、
`run.json`、`manifest.json` 及实际源码/配置快照。CSV 记录能量、T、muB、两种密度
（fm^-3）、无量纲比值、背景残差、状态、端点警告和失败原因。
原子 checkpoint 用 SHA-256 校验；`--resume` 仅允许相同源码、配置、能量、profile
和分辨率，保留成功与失败记录。失败比值为空值/NaN，不能当零或跨缺口连线。

`--parameter-policy` 默认为 `published_central`。跨文献组合或参数带扫描必须显式
使用 `literature_envelope_exploration`；该标识写入 run/manifest，并参与恢复身份，
不能把一次默认参数扫描以不同参数政策恢复。本参数只记录科学来源口径，不改变
积分或背景求解。旧批次的原文字符串元数据保留，数值读取器仍兼容。

跨文献研究的范围与预先固定的评分合同见
[参数包络研究](../../analysis/relaxtime/meson_freezeout_defaults_20261017_analysis/parameter_envelope_20261018.md)。
用以下入口生成新的配置目录，再将 `profiles/` 传给本 Julia runner：

```powershell
python scripts/analysis/relaxtime/freezeout_parameter_envelope.py prepare --output NEW_DESIGN_DIRECTORY
```

设计文件保存各参数单位、原文来源、采样种子、全部参数组及拒绝原因。
`assess --designs ... --runs ... --output ...` 读取同分辨率的冻结批次并对固定观测点
评分；`refine --designs ... --assessment ... --output ...` 生成一轮轴向邻域配置。
它们拒绝覆盖已有目录，且没有统计置信度或全局优化保证。
`horn-map` 使用已计算的热力学峰点调整能量映射，并显式检查原参数包络；其输出
仍必须经 Julia 直接计算核验。`run_freezeout_envelope_refinement.jl` 接收选择文件、
数值计划和分片编号 `0/1`，以同一既有 runner 完成 standard/fine 批次；分片按
整条参数曲线划分，保持每条曲线内部按能量降序延续求解。

绘图和评分只读取已冻结数值与实验输入。实验按体系、中心度、快度和电荷口径标记；
Au+Au 是主比较集合，Pb+Pb 的 horn 作为辅助集合。默认选择指标为 Au+Au 发布点上的
等权 RMSE，同时报告 horn 能区、较高能区和峰位/峰高，不用稠密理论网格的点数加权。
具有公共系统误差的实验不能按独立点构造正式拟合优度或参数置信区间。

历史默认比较使用 `data/inputs/meson_freezeout_experiment_20261017/`，可从已冻结的
原始表重建；生成器会核对来源 hash、电荷与中心度，不需要联网：

```powershell
python scripts/analysis/relaxtime/prepare_freezeout_experiment_data.py `
  --output data/outputs/results/relaxtime/meson_density/experimental_ratios_rebuilt.csv
```

参数包络研究显式传入 `data/inputs/meson_freezeout_experiment_20261018/experimental_ratios.csv`。
新版只将 `NA49_17.27` 的中心度修正为 `0-5%`，数值和全部评分字段不变；
唯一修订及父版本 hash 在同目录 `sources.json`。旧生成器仍复现 17 日输入，
上述目录的区分不代表已完成全部实验元数据的逐项审定。

完整能量扫描完成后生成对照图；若使用自选能量列表，必须包含全部主比较及 NA49
实验能量，缺点会报错，不会用插值代替直接计算：

```powershell
python scripts/analysis/relaxtime/plot_freezeout_comparison.py `
  --run-dir data/outputs/results/relaxtime/meson_density/example_defaults_fine `
  --reference-run-dir data/outputs/results/relaxtime/meson_density/example_defaults_standard `
  --output-dir data/outputs/figures/relaxtime/meson_density/example_defaults_png_review
```

`--reference-run-dir` 可省略。提供时核对两批的源码、物理配置和积分支撑相同，
逐点比较共同坐标的比值差，保留失败状态变化；新的细网格点不会增加实验评分权重。
这检验分辨率变化下的排序稳定性，不代替尾部和端点极限验收。

图像目录只保存彩色/灰度 PNG 与 `plot_manifest.json`。同名 `results` case 保存
比较指标、逐点残差、可选的分辨率差异、图注和绘图源码快照；可用 `--result-dir`
指定该目录。已有任一 case 都拒绝覆盖。交付为 `png_review`，图例位于同一子图内，
不会自动更新论文的 current 图层。

默认比值图固定 `cleymans2006` 为蓝色实线、`andronic2018` 为橙色虚线，子集／重排
保持身份；未登记参数组须先显式添加样式，不能循环两套编码。比值与解析 T–muB
渲染器使用[局部样式及显式字体策略](../../api/plotting/PlotStyle.md)，退出后恢复调用者
样式，并在 manifest 核对曲线编码及实际普通／数学字体文件。

参数包络另用两个只读分析入口汇总与绘图：

```powershell
python scripts/analysis/relaxtime/summarize_freezeout_envelope.py `
  --screen-assessment ALL_SCREEN_ASSESSMENT/assessment.json `
  --design SELECTED_DESIGN/design.json `
  --standard-runs STANDARD_RUN `
  --fine-runs FINE_RUN_0 FINE_RUN_1 `
  --baseline-runs DEFAULTS_FINE_RUN LITERATURE_DEFAULTS_FINE_RUN `
  --experiment data/inputs/meson_freezeout_experiment_20261018/experimental_ratios.csv `
  --output NEW_SUMMARY_DIRECTORY

python scripts/analysis/relaxtime/plot_freezeout_envelope.py `
  --summary-dir NEW_SUMMARY_DIRECTORY `
  --output-dir data/outputs/figures/relaxtime/meson_density/NEW_ENVELOPE_PNG_CASE
```

汇总核对冻结批次及源码快照 hash、物理合同、选定参数与坐标，并将旧默认方案
重新按同一九个 Au+Au 观测点评分。`summary.json` 保存最优样本、完整排名、
同坐标 standard/fine 差异、失败和适用限制；`fine_curves.csv` 与
`selected_parameters.csv` 分别保存直接计算曲线和五参组合。
排序稳定性与差异如实报告，档位名称不替代收敛证明。

灰色带来自固定的完整初筛曲线集合；蓝色带来自主 RMSE 不超过初筛最小值
1.10 倍的完整曲线集合。最优曲线使用加密后的直接计算。初筛带不会因选中三条
曲线完成 fine 核验而自动获得同等精度；两种带均没有统计置信度。
任何已有汇总或图像 case 都拒绝覆盖，修订图应使用新的同级目录。

相同参数组的 T–muB 图使用独立的解析准备阶段；绘图阶段只读准备好的表：

```powershell
python scripts/analysis/relaxtime/plot_freezeout_thermodynamic_envelope.py prepare `
  --summary-dir NEW_SUMMARY_DIRECTORY --output NEW_TMU_TABLE_DIRECTORY

python scripts/analysis/relaxtime/plot_freezeout_thermodynamic_envelope.py render `
  --prepared-dir NEW_TMU_TABLE_DIRECTORY `
  --output-dir data/outputs/figures/relaxtime/meson_density/NEW_TMU_PNG_CASE
```

`prepare` 沿用冻结的完整/近最优成员，在共同 muB 处计算各条解析冻结线温度，
默认覆盖 0–500 MeV。`thermodynamic_band.csv` 保存温度边界与取极值的参数组；
`energy_anchors.csv` 保存 4.87、7.7、19.6、200 GeV 的能量位置。
`preparation.json` 记录全部成员、原始能区、与已有计算坐标的核验和解析延伸口径。
多项式 T(muB) 只依赖 a、b、c，d、e 仍参与能量标记。
解析延伸不产生密度结果，也不改变原 ratio 失败状态。这里的 Andronic 仍为
原文默认参数，不能将此图当作 Andronic 函数族也已调优的比较。

## Andronic 函数族的同构探索

Andronic 候选准备入口为 `scripts/analysis/relaxtime/andronic_parameter_envelope.py`。
默认设计配置 `config/physics/freezeout/parameter_domains/andronic_exploration_v1.toml`
通过 SHA-256 绑定 `andronic_literature_envelope_v2.toml`，不自行扩展文献域。

```powershell
python scripts/analysis/relaxtime/andronic_parameter_envelope.py prepare `
  --output path/to/new-andronic-design
```

生成的 `design.json`、`parameters.csv`、`profiles/*.toml` 与来源快照分别记录固定
五参数、角色、输入 hash 和数值分片。首轮含四个未计算的完整文献中心组、
32 个分层内部样本与 10 个轴向边界探针。2018 默认组复用既有基准；Elfner 与其重复，
JHEP 只有温度三参数，均不生成新的完整密度候选。

使用本页相同的 Julia runner，选择 `--resolution screening`、生成的 profile 目录及
`--parameter-policy literature_envelope_exploration`，按 `design.json` 内的 17 个
`screening_energies_GeV` 直接求解。`parts` 给出两份不重叠的 profile 名单，可作为两个独立进程运行。

```powershell
python scripts/analysis/relaxtime/andronic_parameter_envelope.py assess `
  --designs path/to/new-andronic-design/design.json `
  --runs path/to/screen-0 path/to/screen-1 `
  --experiment data/inputs/meson_freezeout_experiment_20261018/experimental_ratios.csv `
  --output path/to/new-assessment
```

`assess` 复用多项式研究的九个固定 Au+Au 评分点、等权 RMSE、1.10 倍近最优阈值与
完整曲线成带规则，逐点校验本批 Andronic 参数产生的 T、muB；失败不能缩小评分集合。
`refine`、`horn-map`、`select` 均接受 `--assessment`、`--designs`、`--output`，
另可用 `--config` 指定设计配置。细化只做预定的一轮五轴正负 5% 步长；
`horn-map` 尝试最多三组域内的 `E'=k E, A'=A-ln(k)/w`，其中 `k=x_peak/7.7`，
保持完整 T(muB) 几何。没有域内可行映射时报告零候选，不扩大参数域。

`select` 保留 RMSE 领先组和近最优 horn 形状备选，再按 RMSE 补足三组。
使用 `run_freezeout_envelope_refinement.jl` 对相同参数组进行 standard/fine 核验。
`summarize_freezeout_envelope.py` 从设计识别函数族和参数名，复用既有默认参数 fine 结果。
数值背景、密度核和有限 Bose 窗口沿用本页合同。

多阶段本地批次可使用 `run_freezeout_envelope_worker.jl EXECUTION_DIR PART` 保持数值 runner 已加载，
`PART` 取 0 或 1。每个进程依序读取原子发布的 `worker_PART_NNN.json` 请求，
每批仍独立执行既有 `run_scan` 并保留源码快照、run、点文件和 manifest。
请求包括 `schema=freezeout_envelope_worker_request_v1`、`action=scan`、`output`、`profiles`、
`profile_dir`、`energies_GeV` 和 `resolution`；`action=stop` 结束该进程。
完成记录绑定请求 hash 与 scan identity。此方式减少重复编译，不改变候选或密度公式。
`run_andronic_envelope_pipeline.py --prefix andronic_envelope_20261018` 是这次研究批次的协调入口：
它要求同前缀首轮设计/计算和两个 worker 已启动，按固定顺序创建后续同级结果目录，
在数值汇总与解析表完成后停止；图像渲染仍单独执行。
汇总时可用 `--other-family-summary path/to/polynomial-summary/summary.json` 对照既有多项式结果：
脚本校验冻结文件、实验输入、九点集合与 fine 物理合同，并从旧 fine 曲线重算原评分。
`other_family_comparison` 保存源 summary hash、参数、评分和比较范围，不改写旧产物。

`plot_freezeout_envelope.py` 与 `plot_freezeout_thermodynamic_envelope.py` 继续采用上节入口，
从 summary 自动识别函数族。T-muB 图对与 ratio 图相同的固定成员，在共同 muB 网格上逐项取界；
Andronic 通过 `x=(D/muB-1)/E` 解析投影，muB=0 取 T_lim。
与多项式不同，D、E 参与 Andronic 的 T(muB) 几何，只有特定联动变换才保持几何不变。
解析延伸不生成新密度点，旧 2018 ratio 的失败缺口仍然保留。

## 现有研究证据的恢复

默认比较、文献筛选和两族参数探索的同一研究链保存在
[`freezeout_review_20261022`](../../../data/outputs/results/relaxtime/meson_density/freezeout_review_20261022/README.md)。
CSV、summary、参数表和 PNG 可直接读取；需要逐点记录、源码快照或批次复核时，
先在仓库根目录恢复数值和检索证据 ZIP：

```powershell
python scripts/analysis/relaxtime/restore_freezeout_review_evidence.py verify
python scripts/analysis/relaxtime/restore_freezeout_review_evidence.py restore
```

恢复入口核对 SHA-256、CRC、仓库相对路径和已有文件冲突；只恢复缺失文件，
不同字节明确报错。原始绝对路径只表示当时运行位置，当前复核使用本 checkout
对应的同名批次目录。多项式文献 catalog 的 `source_root` 相对仓库根解析。
