# 热平均优化与 mode A 原始异常验证

作者于 2026-10-22 授权尝试性能优化，并以原始生产数据验证非一阶异常。当前工作树为 `D:/w/meson-freezeout-defaults`，分支 `codex/meson-freezeout-defaults`，基准 HEAD `698feaf8fada670fde81c93bab6093b50093935f`。实际运行源码由各批次的不可变快照记录。

## 结论

- `thermal_performance_v1` 已完成 36 组候选计时；`native_performance_v1` 完成 18 组实际接口计时。预计算和严格等质量情形的相对速度简化已接入当前工作分支的实际热平均调用路径，最大同点相对差异 `8.66e-15`。这不表示代码已经合并或部署。
- 多表共享求积仍在诊断层。对各插值区间取并集后，在原 5D 节点累加非负重心权重，可等价积分线性 σ 和线性 xσ；不使用能量箱中心近似。代数测试及实际表的小网格逐查询核验通过。
- 阈值采样升级尚未更改默认缓存。`anomaly_selection_v1/selection.json` 选出 11 个窗口、46 个背景、249 个背景/过程组合；覆盖 mode A 的历史平滑窗口、同相分支残余窗口和快速上升对照。固定原背景，每个窗口均限制在同一相分支内，未跨越一阶边界。
- 两个原批次保留 19 项成功结果，分组投影继续计算余下 230 项，全部成功重放。初轮有 7 项未满足自适应或参考精度要求，经两轮针对性加密后，249 项全部通过；11 个窗口没有缺失任务。
- 四处孤立凹陷得到明显纠正，连续快速变化保留。自适应与独立细参考重建的六种 τ 最大相对差异为 `9.8247e-5`。这里检验的是缓存采样层，不能据此认定所有非一阶变化同源，也不构成整套生产数据的替换。

最终 [summary.json](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/summary.json)、[逐道结果](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/rates.csv)、[六种 τ 重建](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/tau_reconstruction.csv) 和 [来源审计](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/audit.json) 均已生成。原始生产来源文件的哈希保持不变。

## 实际调用路径的性能

`_omega_integral_5d` 在每次调用内预计算单粒子能量、占据数及相对夹角。仅在 `mi == mj` 时把相对速度写成 `2sqrt((s-4m²)/s)`；近似相等质量继续走原表达式。节点、支持区间、缓存/扣除策略、band 输出和原累加顺序保持不变。

以下为三次交错计时的中位数，单位秒；缓存固定为原生产缓存。

| μB [MeV], αT, ξ, 通道 | 原墙钟 | 新墙钟 | 原 CPU | 新 CPU |
| --- | ---: | ---: | ---: | ---: |
| 0, 1, 0.37, uū→uū | 11.8482 | 6.1742 | 10.9688 | 6.0625 |
| 450, 1, 0.26, uū→uū | 11.0472 | 10.4967 | 10.8438 | 8.8281 |
| 0, 1, 0.37, us→us | 10.1306 | 7.4662 | 10.0625 | 7.3438 |

CPU 中位数分别降低约 45%、19%、27%；墙钟受机器负载影响，不能保证所有背景固定比例加速。此处只计一次热平均，不包括缓存构建、Julia 启动/JIT、背景求解、其他通道及最终输运动量积分，因此不能直接与此前“建缓存＋热平均”的 13–20 秒比较。

对应检查：原热平均单元测试 97 项、新生产调用同点/边界/band 检查 74 项、共享求积检查 65 项通过。每个异常验证任务也通过了完整历史节点的原结果重放。旧 `tau_xi_probe` 文件哈希检查和固定 τ 的输运积分测试不能单独证明此次热平均优化正确，未用它们替代数值重放。

## 诊断求积的角度合并

`GroupedThermalProjection.jl` 把相对夹角向量完全相同的 `(cosθ_i,cosθ_j)` 节点归为一组，然后合并组内完全相同的 `cosΘ` 所对应的 φ 权重。占据数仍逐个极角节点计入，不假设其对称性，也不移动或近似任何 s 查询自变量。默认 `angle20/phi36` 从 14400 个角度组合减少为 2855 个不同几何组合。

这利用的是当前截面缓存只以 s 为查询自变量的实际实现。若以后保留更多外部方向依赖，需要重新推导可合并条件。该实现目前只用于分析求积，未改变公开热平均 API 的默认执行方式。

`grouped_performance_v2` 的三次交错计时中位数如下。两种实现同时计算原缓存、自适应表及两级参考，计时范围不同于上一表的单次生产热平均。

| μB [MeV], ξ, 通道 | 原投影墙钟 [s] | 分组墙钟 [s] | 原 CPU [s] | 分组 CPU [s] |
| --- | ---: | ---: | ---: | ---: |
| 0, 0.37, uū→uū | 12.4301 | 2.5316 | 12.2969 | 2.5156 |
| 450, 0.26, uū→uū | 12.2372 | 2.0719 | 12.1094 | 1.9844 |
| 0, 0.37, us→us | 12.9362 | 2.3949 | 12.8125 | 2.3750 |

首次严格对照的 `1.25e-11` 差异定位于原长序列普通累加。独立全节点 Kahan 求和显示三个情形原主积分的相对舍入偏差分别约 `−1.43e-11`、`−1.80e-11`、`−4.63e-12`。分组主积分改用补偿求和后，与该参考以及原投影的其他三项积分的最大相对差异为 `5.11e-15`。验收仍用 `1e-11`，原数据重放仍用 `1e-8`；没有放宽误差要求。40 项补充检查覆盖非对称自定义网格、不同质量、阈值直算以及 raw/xσ 表示。

## 输入与验证设计

用户给定的图清单位于 `D:/Desktop/Julia_RelaxTime/data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2/plot_manifest.json`。其原始 `mode_a_scan.csv`（910 行）及来源绘图清单的 SHA-256 已匹配。对应的 `mode_a_diag.csv`（38220 行）记录各道密度、平均截面与逆弛豫时间贡献。

窗口内选择对物种总逆弛豫时间的局部曲率贡献至少 `5e-4` 的过程，并至少包含排名前两道。未选过程完整保留在 `channel_plan.csv`；后续 τ 重建是明确标注的部分通道替换。此筛选不能证明未选过程的绝对误差为零。

46 个背景均记录为 `converged=true`，最大背景残差范数 `1.2631e-9`。逐道 `multiplicity × density × rate` 与其记录贡献完全一致，贡献之和也与诊断表 total 一致。诊断表 total 与主表 inverse τ 最大相对差异 `2.83e-12`，来自两份历史表本身；本轮以主表原 τ 为锚点，增量替换选道贡献，保留该微小历史差异。

热平均维持 `p128/angle20/phi36`，截面 `t24`。比较原网格重放、局部 xσ 加尾部 σ 自适应表、独立 log-x/√x 固定加密的 192/384 级参考。阈值和整个截面支持区间分别检查；自适应目标 `1e-4` 只控制缓存误差，不代表外层积分、环函数或完整输运已经收敛。失败与未收敛结果必须保留，并单独加密复核。

本轮保持原数值口径：各向同性 A 用于有效耦合、full-RS A 用于极化、B₀ 分布修正仍为旧一阶 ξ 形式。改进 B₀ 端点是独立方向，本批不把其变化混入采样对照。

接受标记分别为：原平均截面重放相对差异 `≤1e-8`，相邻两级热平均参考变化 `≤2e-4`，自适应表与细参考差异 `≤3e-4`。初轮采用 192/384 级参考；未通过项提高至 384/768、768/1536 或 1536/3072 级，具体层级逐行记录。这些是本轮诊断门槛，不是完整输运系数的总误差界。

## 各窗口的结果

以下四处原来的下降凹陷在保持背景和外层节点不变的情况下，恢复为单调变化；没有对 τ(ξ) 做平滑。中心点数值用于展示数值误差的量级，单位 fm。

| μB [MeV] | αT | ξ | 物种 | 原 τ | 重算 τ | 观察 |
| --- | --- | --- | --- | ---: | ---: | --- |
| 0 | 1.0 | +0.37 | u | 2.145957 | 2.285158 | 孤立凹陷基本消失 |
| 450 | 1.0 | +0.26 | u | 1.674240 | 3.798967 | 大幅孤立凹陷基本消失 |
| 450 | 1.1 | +0.35 | u | 1.408881 | 1.472461 | 下降凹陷消失，底层快速上升和弯曲保留 |
| 0 | 1.2 | −0.47 | s | 0.566420 | 0.683790 | 孤立凹陷基本消失 |

其余窗口的快速变化没有被此层采样修正消除：μB=0、αT=1 的 ξ≈−0.22 对照；μB=0、αT=1.1 的 ξ≈+0.48；μB=450、αT=1 的 ξ≈−0.20；以及 μB=900 的四个同相分支窗口。详细数值见 [window_summary.csv](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/window_summary.csv)。这些是保留快速变化的数值证据，本轮没有重新求介子极点，不能单凭它们给出新的 Mott 归因。

图像按项目目录习惯放在 `data/outputs/figures/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/`。每个窗口分别展示原曲率更大的夸克或反夸克物种；选择依据为原始数据，全部六种 τ 均保存在 CSV 中。图例位于子图内，连线仅连接计算节点。

- [夸克对照 PNG](../../../../data/outputs/figures/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/anomaly_windows.png) · [PDF](../../../../data/outputs/figures/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/anomaly_windows.pdf)
- [反夸克对照 PNG](../../../../data/outputs/figures/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/anomaly_antiquark_windows.png) · [PDF](../../../../data/outputs/figures/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/anomaly_antiquark_windows.pdf)

## 加密核验及生产替换判断

初轮的 2 项参考未收敛、5 项自适应误差超标都位于 μB=900、αT=1 的窗口。`anomaly_refinement_v1` 复用保存的直接截面值，把自适应目标收紧至 `1e-5`，并提高参考层级。7 项自适应结果全部达标，但其中 5 项的固定参考仍需加密。`anomaly_refinement_v2` 使用 1536/3072 级参考后，这 5 项的参考变化最大 `2.1522e-5`，自适应偏差最大 `5.6706e-5`，均达标。

这里暴露了两项真实限制：`x` 加权的局部误差目标不等于完整热权重下的误差目标；192/384 两层接近，也可能暂时漏掉窄结构。因此必须同时比较独立设计的自适应表和固定参考。本轮没有放宽验收容差，也没有以曲线更光滑作为收敛证据。每次被更新的结果保存在 [superseded_rates.csv](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/anomaly_analysis_v1/superseded_rates.csv)，原批次文件保持不变。

当前实现选择是：保留已经接入实际调用路径的预计算和等质量速度简化；角度合并与自适应缓存继续作为可复现诊断工具。采样策略设为生产默认前，应补齐完整 17 道、真实热积分误差控制及缓存/API 兼容检查；B₀ 端点、t 积分和外层热积分的收敛仍需分开处理。当前未改动这些数值默认，未重新生成整套生产数据或触发 Actions。

## 保留的试算与执行记录

- `anomaly_pilot_v1/v2` 保留早期性能试算的停止状态。v2 的建表与投影完成，但异质 callable tuple 的高节点联合核验异常缓慢；随后改为小网格实际表等价核验，并将 tuple 查询递归展开。未以这些失败试算给出物理结论。
- `anomaly_pilot_v3` 的 ξ=0.37 两道通过全部检查。完整支持区间的 uū→uū 细参考 `0.5067039024 fm²` 与上一轮只改善局部阈值得到的约 `0.5036461 fm²` 有差别，说明远离阈值的缓存也有小幅误差。
- `anomaly_full_s0_v1/s1_v1` 共保留 19 项完整结果；其 `complete=false` 表示原全批计划被暂停。`rates.csv`、输入及源码快照保留，继续批次只补缺失任务。
- `grouped_performance_v1` 保留严格浮点对照未通过的状态。`grouped_performance_v2` 已完成 18 组计时，并以独立全节点补偿求和确认误差来源；相关结果见上表。
- `anomaly_grouped_s0_v1/s1_v1` 分别完成 112、118 项，跳过此前已完成的 9、10 项；没有重复计入旧结果。两轮精度补算分别处理 7、5 项。
- 四组 Julia 检查合计 276 项断言通过；脚本入口检查覆盖 313 个 Julia 文件，依赖策略检查通过。重建公式、缺失任务和不同参考层级的汇总行为另做了针对性检查。图像已实际查看，历史输入及运行源码按各自快照核验。

所有本轮产物位于 [结果目录](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/)。历史数据、后处理显示层和 Actions 均未替换。
