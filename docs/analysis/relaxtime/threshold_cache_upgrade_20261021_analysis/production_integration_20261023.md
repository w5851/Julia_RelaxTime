# 输运新默认的生产接入与外部对齐

工作树为 `D:/w/meson-freezeout-defaults`，分支 `codex/meson-freezeout-defaults`，
HEAD 基底为 `dffe4070eb5c0788e9c8a1e227cf40978b51ce3e`。
数值实现尚未提交，各批次 `source/` 与 manifest 保存实际代码。
当前进度与待办由[生产接入任务](../../../dev/active/2026-10-23_输运采样生产接入与重跑合同.md)维护。
此前代表集的 1% 收敛资格和完整性能比较见[采样验收报告](production_sampling_20261022.md)。

本报告中的数值验收和 Julia--Fortran 对齐均记录完整 \(A^\xi\) 切换前的生产候选；
2026-10-24 已将有效耦合从各向同性 \(A^{(0)}\) 切换为完整 \(A^\xi\)，
因此下文结果只作为历史验收与重跑基线，不能替代切换后的独立收敛、联合验证和
同点对齐。

## 配置和实现

两个生产扫描 CLI 默认使用 `thermal_weighted_v1`，统一背景、A、热平均、缓存预算及末端积分。
旧固定网格、扣除/加回、w0cdf/PCHIP 与无身份 CSV 缓存加载已从生产实现移除。
历史实验和生产数据仍按各自冻结源码重放。
相同相对夹角的运算合并与热权重自适应采样共同使用；两者原理及局部/完整加速口径见采样报告。

`TransportSamplingProfiles.solve_profile_transport` 按已验收的 17 道顺序重置极化缓存、建表，
再由真实 workflow 消费缓存。采样预算耗尽或无有效背景时保留失败，不能写入有效 τ。
续算绑定实际数值源码、模型、profile 和物理配置。

## 新背景的相边界与温度锚点

`phase_revalidation_v5` 完成九条 mode-A 路径的双向续接，910 个计划点全部通过，
另有 18 次冷种子审计全部通过。比较 25 fm⁻¹、p96/cos24 与 30 fm⁻¹、p128/cos32。
质量最大相对变化为 `2.80351545e-6`，熵为 `2.97450138e-5`。
冷种子/续接解的最大质量相对变化为 `5.56849e-9`，最大 Ω 绝对差为 `3.19744e-14`。

| μB [MeV] | 锚点定义 | 新 T [MeV] | 相对旧参考移动 [MeV] |
| --- | --- | ---: | ---: |
| 0 | 轻凝聚温度导数绝对峰 | 200.08981 | +0.00091 |
| 450 | 同上 | 182.72820 | +0.08808 |
| 900 | 两分支等 Ω | 125.73600 | +0.00194 |

三点均通过独立 0.1 MeV 锚点检查。crossover 搜索分辨率为 0.002 MeV，
细化差为零表示在该搜索分辨率内未分辨移动，不能视作物理误差严格为零。
μB=900、αT=1 的固定温度一阶边界在两档背景下均为 ξ≈0，
二分结果 `−1.78814e-10`，括区间为 `[−0.003,+0.003]`。
共存两侧点继续分别保留。本检查覆盖生成的生产路径，不替代全局相图核验。

v3/v4 的锚点或部分背景记录、中断及源码快照保留。
v5 修复了根求解迭代预算未透传的问题；每次 Newton 或 trust-region 尝试使用同一 40 步上限，
残差和物理性判据保持原标准。`phase_v5_independent_check.json` 复核实际点集、误差和当前 322 个相关源码/配置文件。

## Julia–Fortran 四个 ξ=0 点

`fortran_alignment_v2` 完成 8 组背景、68 个通道和 4 组最终量比较，零失败。
采用原始四个学位论文读图点，未通过调节积分上限拟合目标，外部标准仍为 5%。

| μB [MeV] | T [MeV] | 论文读图 η/s | 旧 Fortran η/s | 新 Julia η/s | 新 Julia 对读图偏差 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 190.32399 | 0.32200669 | 0.32661787 | 0.32625212 | +1.3184% |
| 0 | 250.61110 | 0.09975537 | 0.09910711 | 0.09954296 | −0.2129% |
| 800 | 130.05398 | 4.00685901 | 4.17319201 | 4.16792134 | +4.0197% |
| 800 | 317.62728 | 0.30046150 | 0.29773692 | 0.30895846 | +2.8280% |

四点均通过。新 Julia 对旧 Fortran 的 η/s 相对差异依次为
`−0.11198%、+0.43977%、−0.12630%、+3.76894%`。
匹配 Fortran 常数（包括其舍入的 ħc）及 15 fm⁻¹ 背景时，质量、Polyakov 环、熵的最大相对差异为 `5.69243e-9`。
冻结 Fortran 背景与 τ 后，密度、17 道组合得到的 τ、末端 η 重放最大相对差异为 `6.13177e-10`。
新生产配置下的 68 道平均率对旧 Fortran 的最大相对变化为 1.02949%；这里同时含新背景与下游求积变化，不能全部归为缓存误差。

最高温点的差异以尾部为主：冻结背景与 τ，仅将末端 15 fm⁻¹ 改为半无穷，η 增加 5.20299%；
新背景的熵增加 1.18275%，新 τ 也有小幅变化，最终 η/s 对 Fortran 增加 3.76894%。
`η/s` 的分子、分母应分别比较，不能把最终百分比全部解释为弛豫时间或散射率误差。
旧 Fortran 四点的 5% 读图符合性不意味着其所有网格已收敛；内部 1% 收敛仍由独立及联合比较验收。

## 当前生产入口复算

`production_replay_all_v2` 已在当前源码上调用真实 `run_phase_guided_scan`，
完成十二点、108 项最终量及 204 道采样诊断，零失败，进程退出码为 0。
相对原已验收候选的最大变化为 `1.45729785e-6`（0.00014573%），满足 0.01% 接入要求；
相对独立参考的最大变化为 **0.6863883%**，满足 1%。
产物检查逐点核对实际 CSV、有效配置、收敛状态及通道预算；其 322 个数值源码/配置文件与相界核验及当前工作树一致。

### 背景输入与严格回归基线

首次对局部分支种子候选 v1 运行回归，结果为 464 项通过、九项失败。
输运定点 117 项、截面定点 15 项及一阶稳定分支 27 项全部通过；九项失败属于生产入口的两个点：

| 点 `(μB [MeV], αT, ξ)` | 超过 v1 内部容差的量 | 最大相对变化 |
| --- | --- | ---: |
| `(0, 1, +0.37)` | τ_u/d/ubar/dbar、η/s、ζ/s、σ/T | `1.45730e-6` |
| `(900, 1.1, 0)` | τ_ubar/dbar | `1.21935e-6` |

原候选保留局部分支种子，生产入口重新选择稳定解。首点质量最大差约 `7.33e-10 fm⁻¹`，
原候选的 u/d 质量差约 `1.38e-10 fm⁻¹`，生产解精确相等；第二点质量最大差约 `3.73e-11 fm⁻¹`。
逐道比较定位到首点的 `uubar_to_ddbar` 和第二点的 `udbar_to_udbar`、`dubar_to_dubar`。

`run_transport_background_replay_diagnostic.jl` 对三个通道分别固定原候选背景、生产背景、重复生产背景，
并将原候选的 u/d 质量取均值作为受控扰动，共完成十二次建表。
冻结背景后各自复现原率，最大相对差为 `2.22e-16`；三个通道的同输入重复结果均逐值相同。
仅消除原背景的微小 u/d 不等后，率相对生产值的差分别为 `4.30e-9`、`4.46e-10`、`4.46e-10`。
这些证据把漂移关联到背景求解入口及近等质量输入的数值敏感性；质量取均值只用于诊断，未加入生产算法。

据此建立 `baseline_thermal_weighted_transport_v2.csv`，明确保护稳定背景生产入口。
十二点坐标及独立参考不变，另以 `qualified_candidate` 保留 v1 原候选；v1 文件和首次失败日志继续保留。
内部回归仍为 `rtol=1e-6, atol=1e-8`，独立参考仍为 1%，未以放宽容差接受漂移。
完整 v2 回归 314 项通过。随后在一个新进程中经实际 CLI `main()` 重新计算两个受影响点，
省略 `--sampling-profile`，每次启动独立扫描，不使用上一扫描的续接解。
两点各 22 项默认入口检查及现有单点回归 28 项全部通过；十八项最终量与 v2 **逐值相同**。
另行核对两份实际产物、配置、源码身份及 34 道预算，均通过。全过程未再修改数值源码。

## Runner 输入

[runner 输入](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/runner_inputs_v1/runner_inputs.json)
使用新 case `thermal_weighted_v1_reanchored_20261023`：六点先导分四个 jobs，暂定整批为 910 点、54 片，
每片最多 20 点。输入坐标及温度已逐点核对新相界计划，包含高温负 ξ 端点及一阶共存侧。

当前数值实现尚未提交，未指定计算 SHA，因此 `pilot_input_ready=false`、`full_dispatch_ready=false`。
指定 SHA 时，输入生成器会核对提交的实际数值文件集合及内容，不能用只满足 40 位格式的历史 SHA 替代当前实现。
固定实际计算提交后，先做六点 runner 先导及产物检查，再依据完整时间、启动和峰值内存确定整批分片。

本 runner 输入仅含 mode-A。mode-B 历史 T={120,160,200}、μB={0,450,900} MeV 的 909 点范围已核对，
新背景下的分支/共存侧及批次尚未准备。[生产参数留档](production_plan_20261023.md)
保存九条具体温度、各层积分、缓存预算、W 资源使用和 runner 剩余项，结构化 JSON 关联实际输入与配置。
长数值任务优先使用已经连通的 W；GitHub runner 的时间和内存资格仍在 GitHub 检查。

本次 `production_replay_all_v2/timings.csv` 是验证运行的时间记录，非受控性能基准。
其中 `(μB,αT,ξ)=(0,1.2,-0.5)` 耗时约 1566 秒，已列入先导；不能直接把暂定 20 点分片作为最终用时预算。
专门的完整路径计时及角度合并加速口径仍以[采样验收报告](production_sampling_20261022.md)为准。

## 直接证据

结果根为 `data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/`：

- [相界与背景状态](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/phase_revalidation_v5/status.json)、同目录 `anchors.csv`、`background.csv`、`cold_seed_audits.csv`、`traces/`。
- [Fortran 对齐状态](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/fortran_alignment_v2/status.json)、[四点最终量](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/fortran_alignment_v2/observables.csv)、同目录 `background.csv`、`channels.csv`、`model_parameters.json`、`inputs/fortran/`。
- [生产复算状态](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/production_replay_all_v2/status.json)、同目录 `comparisons.csv` 和各点 `points/`。
- [默认入口重复校验](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/production_repeat_defaults_v1/verification.json)、同目录 `points/`、`inputs/`、`source/`；`production_integration_checks_v1/production_default_repeat_artifacts_v1_check.json` 核对两份真实产物。
- `production_integration_checks_v1/production_artifacts_all_v2_check.json` 保留十二点真实产物核验，`background_drift_review_v1.json` 汇总三个通道的固定背景诊断。
- 固定背景的完整状态、率、表和源码分别在 `production_background_replay_v1`、`production_background_replay_muB900_udbar_v1`、`production_background_replay_muB900_dubar_v1`。
- `production_integration_checks_v1/fortran_alignment_v2_summary.json` 保留对齐汇总及输入文件哈希。
- `production_integration_checks_v1/python_production_checks_v9.log` 的 34 项检查通过；包括结果行未收敛、缺失收敛标记而汇总仍称成功时的拒绝检查。`actionlint_v2.log` 验证两个修改的工作流通过；`actionlint_v3.log` 另验证已去除过时参数的 manual workflow，共三个 workflow 通过。
- `unit_core_selected_v1.log` 的 284 项及 `integration_core_selected_v1.log` 的 37 项通过；原 v1 回归失败保存在 `regression_production_selected_v1.log`，各阶段退出码在 `verification_sequence_status_v1.json`。
- `production_integration_checks_v1/production_default_repeat_and_regression_v2.log` 保留 v2 的 314 项、两点默认入口的 44 项及现有单点回归的 28 项通过记录。
- [W 手动入口检查](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/remote_manual_preflight_v1/reviewed_result.json)
  核对两点小型自定义 profile、23/23 日志、退出码 0、34 道采样诊断及运行前后 879 个源文件身份一致。
  原包装的单独标记存在 soft-scope 记录错误，原文件与更正依据均保留；生产收敛与 runner 时间资格不从该小型检查推定。

当前未启动 GitHub Actions，未覆盖历史生产 CSV 或论文图。
