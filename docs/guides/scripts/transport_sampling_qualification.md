# 输运采样的联合验证、续算与完整计时

这些入口用于算法验收，不启动 Actions、不覆盖原生产 CSV。
共享热平均实现及两个生产扫描 CLI 已切换到 `thermal_weighted_v1`，旧固定网格、
扣除/加回及无身份 CSV 缓存加载实现已清理；本地验证状态由
[生产接入任务](../../dev/active/2026-10-23_输运采样生产接入与重跑合同.md)记录。
验收问题、独立分层结果和当前结论见
[生产采样分析](../../analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/production_sampling_20261022.md)。
底层接口见 [AverageScatteringRate](../../api/relaxtime/transport/AverageScatteringRate.md)。

切换完整 A^xi 耦合前的历史汇总为 `threshold_cache_upgrade_20261021/joint_complete_summary_v1`：
十二点、108 项通过 1%，两点各两次完整计时的输出逐值复现候选。
本机单点中位数为 346.53 和 424.53 秒；p48 已接入默认生产配置。
切换前源码的十二点生产复现、910 点新背景及四点 Julia–Fortran 对齐已通过；
v2 回归和两点全新默认入口重复检查完成，runner 输入及后续实测由当前接入任务维护。
这些数值验收不覆盖 2026-10-24 的完整 A^xi 耦合切换；新口径验证状态以当前接入任务为准。
完整配置、适用范围与成本预算以分析报告为准。

完整 A^xi 口径的十二点联合验收已在 quote-remote 通过，108 项最终量的最大相对变化为
0.137533%，408 条通道率完整；[返回复核](../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_joint_representatives_p48_v1/returned/reviewed_result.json)
保留实际源码包与输入身份。独立 cheap/kernel/outer 分层也已通过，A/t 最大相对变化为
1.49490e-6，外层最大为 0.00386317；
[分层返回复核](../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_layer_validation_setup_v2/returned/reviewed_result.json)
与联合结果分别保留。真实生产入口的十二点、108 项最终量与 204 道诊断也已通过，
四个 ξ=0 Fortran 点通过外部 5%；结果及适用域见
[产物复核](../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_production_replay_v3/returned/reviewed_result.json)。
v3 回归 314/314、跨目录回放身份 13/13 以及五次真实默认入口的 45 项量和 85 道诊断均通过；
[默认入口复核](../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_baseline_acceptance_setup_v1/returned/reviewed_result.json)
记录 quote-remote 的两个稳态中位时间 143.087 秒、147.331 秒。该时间不作为 GitHub runner 成本。
mode-A 的三个新锚点、910 个背景点及 18 次冷种子审计也已通过，见
[背景复核](../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_phase_validation_v1/returned/reviewed_result.json)。

在计算工作树根目录使用 Julia 1.12.5 的根 `Project.toml`。输入目录须包含已选定并保留来源的
`review_points.csv` 和 `channel_plan.csv`。所有输出目录必须是新的空目录。
输出保留实际输入、源码/配置快照和 manifest；未提交代码以快照为准，HEAD 只标识基底提交。

在无 `.git` 的冻结源码包中，设置 `THRESHOLD_SOURCE_MANIFEST` 为包的 `package_manifest.json`。
同目录须保留 `source.zip`；入口核对 ZIP 的 `zip_sha256`、`source_files` 中的逐文件 SHA256，
并要求本次源码快照中的文件全部列入包清单。缺失或不匹配时拒绝启动。
输出保存包清单副本及其 SHA256，`source_archive` 记录归档身份与 `head_base`；
此时顶层 `head`/`branch` 留空，未提交包的 `calculation_commit` 也必须为空。
未设置该变量时继续使用实际 Git 工作树的分支和 HEAD。
真实 phase-guided 入口还需要默认 `data/reference/pnjl/issue130_phase_reference_v2/` 运行时输入。
仅打包源码不足以启动该入口；可在新的执行副本中附加逐文件校验的参考输入包，
同时保留源包与参考输入包的独立身份，不修改已经冻结的源码目录。
后续收敛、数值回归、复算与计时在 `quote-remote` 执行。

## 联合验证

```text
run_transport_joint_validation.jl INPUT_DIR NEW_OUTPUT_DIR
    [pilot|representatives SHARD SHARDS REUSE_DIR|- CANDIDATE_P REFERENCE_DIR|-]
```

- `pilot` 检查 μB=0、αT=1、ξ=−0.22；`representatives` 检查分析报告列出的十二个点。
- 分片编号从 0 开始。`REUSE_DIR` 只能复用已通过、输入与数值配置匹配的完整点。
- `REFERENCE_DIR` 只复用已通过运行的独立参考；核对输入、参考配置、数值源码和实际参考函数。
- `CANDIDATE_P` 显式覆盖候选外层动量节点。省略时仍是原先导的 96，不能假定已改为 48。
- 每点消费全部 17 道 native 缓存，检查六种 τ、η/s、ζ/s、σ/T；最终相对变化目标是 1%。
  缓存插值预算为 0.1%，独立参考 σ 表的加权 L1 预算为 0.02%。
- `status.json:joint_target_passed` 才是该批的最终结果。一个完成的通道或一个输出文件不能替代它。

例如，显式运行 48/36/36 的完整先导：

```powershell
$base = 'data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021'
julia --startup-file=no --project=. --heap-size-hint=1000M `
  scripts/analysis/relaxtime/run_transport_joint_validation.jl `
  "$base/convergence_selection_v1" "$base/joint_example_p48" pilot 0 1 - 48
```

τ 外层和其归一化密度沿用半无穷映射 scale=10；本轮 A 和末端输运的映射 scale=5。
这些是不同积分层的设置。背景热上限是 25 fm⁻¹，物理截面 Λ 截断保持不变。

## 完整 A^xi 的分层验证

完整 A^xi 切换后的分层入口为：

```text
run_transport_full_rs_layer_validation.jl INPUT_DIR NEW_OUTPUT_DIR cheap|kernel|outer
```

三个阶段使用新的空目录，在 quote-remote 串行执行。`cheap` 重建十二点的候选与参考背景，
以 RS 变量代换检查内部及 workflow 密度、粒子对归一化误差，并独立加密二维 A 和末端输运矩。
`kernel` 在八个敏感背景/通道组合分别比较 A32/64/128、t24/48/96 及 A/t 联合加密。
`outer` 在三个敏感通道独立增加动量、极角、方位角、密度节点，改变半无穷映射尺度并做联合加密。
后两阶段为每档重新生成独立截面表，流式消费实际热节点；不读取历史各向同性耦合的率或表。
检查 `status.json:all_cases_complete` 和 `layer_target_passed`，并保留失败与逐通道变化。
分层通道的 1% 筛查、0.1% 密度归一化预算和最终九项量的 1% 联合验收分别记录。
末端 Gram 矩阵比较适用于该背景上任意正 τ 和 bulk 括号；参考积分仍须加密检查。
重建背景的 A 要求逐位一致。有效耦合 helper 的 `@fastmath` 内联可能产生机器舍入差，
`identity.csv` 单独记录其最大相对差并限制为 `8eps(Float64)`；这不是物理收敛容差。

## mode-B 背景网格复核

```text
run_transport_mode_b_background_validation.jl INPUT_DIR NEW_OUTPUT_DIR
```

该入口已在 quote-remote 通过 909 点及 94 次冷种子审计，见
[mode-B 背景复核](../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_mode_b_background_setup_v1/returned/reviewed_result.json)。
固定 T={120,160,200} MeV、
μB={0,450,900} MeV、ξ=−0.5:0.01:0.5，共 909 点，不消费 mode-A 的温度锚点。
使用与 mode-A 相同的两档背景参数，分别做双向分支续接；每条曲线在两端、四分位点、
中点以及发现的相变符号括区端点做独立冷种子审计。质量和熵变化上限仍为 1%，
冷种子质量相对差为 1e-6、热力学势绝对差为 1e-8。失败及未解析的精确共存点保留，不自动删点。
`phase_plan_passed` 只表示该离散网格的背景验收；`transport_qualified=false` 明确保留
输运代表点收敛和真实入口复算的独立要求，也不将边界括区当作全局相图。

通过背景复核后，`prepare_transport_mode_b_qualification.py BACKGROUND_DIR NEW_INPUT_DIR`
从实际稳定分支的轨迹提取代表点，覆盖九条曲线的 ξ=−0.5、三个反向端点、三个各向同性控制点，
并加入每个一阶共存括区的两个网格端点。输入保留背景、轨迹、源码和选择理由。
`run_transport_mode_b_qualification.jl INPUT_DIR NEW_OUTPUT_DIR STAGE [JOINT_SUMMARY]`
复用现有 cheap/kernel/outer、联合及真实生产复算方法，`STAGE` 为
`cheap|kernel|outer|pilot|joint|production`；各阶段仍须独立通过原有门槛。
代表点准备器的 6 项测试、Julia 选择器的 7 项断言及 cheap/kernel/outer 三层已通过；
17 点联合运行的 153 项最终量及 578 条通道率全部通过，最大相对变化为 0.0720757%；
真实入口的 17 点、153 项量和 289 道诊断也已通过；相对候选最大变化为 6.29522e-9。
[完整返回复核](../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_mode_b_qualification_setup_v1/returned/reviewed_result.json)
已在 quote-remote 完成，归档及 1584 份源码快照、15 份输入快照已取回核验。
此资格覆盖该代表点集合，整个 909 点输运网格尚未计算；背景、代表点和整批产物的状态分别保留。

## 中断后的检查点恢复

```text
resume_transport_joint_validation.jl SOURCE_RUN NEW_OUTPUT [SOURCE_RUN NEW_OUTPUT ...]
```

该入口用于已经开始全部目标点、并完成其 native candidate 的运行。它核验输入和数值源码，
检查已存散射率的误差预算、完整通道集合、密度及九项比较，保留已完成点。
未完成的 reference 先重建并核对原背景，再只计算缺失通道。可以一次串行恢复多个分片。
原目录和中断状态不修改；新 manifest 记录恢复来源和 CSV 哈希。

以下情形明确拒绝恢复：存在未处理的数值失败、native candidate 尚未完整消费缓存、
部分目标点尚未开始、输入/公式/数值配置发生变化，或检查点数据不完整。
这些情形需要针对实际缺项重新运行对应点；不能补写一个成功状态来跳过验证。

## 生产基线与只读回放

当前 full-RS v3 基线候选由已验收的真实入口产物生成，v1/v2 及其旧数值承诺保留。
旧基线审查先区分实际数值差异与身份核验失败；不能通过放宽容差或直接刷新 CSV 消除失败。
只读回归回放允许源码根目录搬迁，但配置文件的仓库相对路径须相同且列入源码指纹。
核验时仅恢复原配置绝对位置用于重建身份哈希，所有数值参数、采样选项和源码哈希仍精确匹配，
并核对保存的身份内容。当前生产续算接口保持原有合同。

## 完整单点计时

当前默认入口使用：

```text
benchmark_transport_production_entrypoint.jl INPUT_DIR NEW_OUTPUT_DIR BASELINE_CSV [REPEATS]
```

此入口位于 `benchmark/relaxtime/`，使用根项目环境，在 quote-remote 执行。
它对 μB=0/αT=1/ξ=−0.22 和 μB=900/αT=1/ξ=−0.04 调用真实 `run_phase_guided_scan`，
不传 `--sampling-profile`，从而检查实际默认值。首次完整调用作为预热单独保留；随后每点至少重复两次，
记录完整 wall/CPU 时间、GC 和累计分配量。累计分配量不等于峰值内存，W 的计时不代替 runner 先导。
每次调用核验 17 道诊断，九项最终量须复现提供的生产基线（rtol=1e-6、atol=1e-8），
并满足其独立参考的 1% 门槛。检查 `status.json:default_entrypoint_passed`；原基线和历史结果不覆盖。

切换前的新旧算法成本对照使用以下历史入口：

```text
benchmark_transport_threshold_complete.jl INPUT_DIR NEW_OUTPUT_DIR [REPEATS CANDIDATE_P]
```

计时点是 μB=0/αT=1/ξ=−0.22 与 μB=900/αT=1/ξ=−0.04，默认两次交错重复。
比较旧缓存的原角度路径、旧缓存的 grouped 路径和 thermal-weighted 候选。
含旧算法的该次计时应从 `complete_timing_p48_v1/source` 冻结源码快照重放，
不能在当前已移除旧 API 的源码树直接执行旧路径。原始数据、计时和源码快照保持不变；
历史消融及旧策略对照脚本遵循同一规则，当前生产只维护热权重实现。
计入背景、17 道建表与热平均、bulk 导数和末端输运，显式 GC 计入总用时。
候选按 `REQUIRED_RATE_KEYS_FOR_TAU` 顺序逐道重置极化缓存并构建 native 表，
随后把 `cs_caches` 传入实际 workflow 消费；生产适配应保留这一已验证的构建顺序和缓存生命周期。
旧路径的建表由两种角度方法共享，各自总成本由共同组件与独立投影时间组成。

`warmup.json` 单列预热及其 JIT 时间；稳态计时不含这一预热和 Julia 进程/模块加载。
`timings.csv` 同时给出 wall/CPU 用时，`channels.csv` 给出建表分项，
`observables.csv` 保留完整输出。旧路径是成本基准，不能作为新算法的精度参考。
必须与其他 Julia 数值任务分开运行；分片并发期间的诊断耗时不用于生产预算。

## 汇总

```powershell
python scripts/analysis/relaxtime/summarize_transport_joint_validation.py `
  "$base/joint_native_pilot_p48_v1" `
  "$base/joint_native_representatives_p48_s0_recovered_v1" `
  "$base/joint_native_representatives_p48_s1_recovered_v1" `
  --expected-points 12 --benchmark "$base/complete_timing_p48_v1" `
  --outdir "$base/joint_complete_summary_v1"
```

汇总核对完整状态、数值配置/源码、输入、点位、九项量、相对变化和通过标记；
拒绝重复点、缺失计时点或重复次数、非有限用时。输出 `summary.json` 和 `comparisons.csv`。
计时的九项实际输出也与联合验收比较，记录候选复现差异并要求仍满足同一 1% 参考。
退出码 1 表示联合目标未达到，应查看缺项或失败，不能改低预期点数来取得通过。

## 共用密度节点的独立核验

候选和联合参考的热平均内部密度都使用 p32。仅比较这两档完整结果不能排除密度中的共同误差。
`verify_transport_joint_density.py` 先复现既有 Julia 密度记录，再用 RS 变量代换
`n(ξ)=n(0)/sqrt(1+ξ)` 和独立半无穷高阶径向积分检查新背景的内部归一化。
结果中的粒子对误差上界只针对固定正分子下的密度分母，不能替代完整散射率收敛。
本轮 24 组背景、144 项通过，证据保存在 `$base/joint_density_normalization_v1/summary.json`。

```powershell
python scripts/analysis/relaxtime/verify_transport_joint_density.py `
  "$base" "$base/joint_density_example"
```

此入口读取本轮已命名的联合分片和既有密度对照，输出目录同样要求全新；不启动 Julia 补算。

十二点资格只覆盖所述代表集；两点本机计时不等于 GitHub runner 实测，也不是全扫描用时上界。
生产入口、相边界/温度锚点及 Julia–Fortran 对齐结果见
[生产接入报告](../../analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/production_integration_20261023.md)。
正式重跑的分片依据 mode-A 六点及 mode-B 五点的 runner 先导校正，实测记录见下文生产参数留档。

## 生产入口的复现与新背景核验

版本化配置的 API、积分层次及失败语义见
[TransportSamplingProfiles](../../api/relaxtime/transport/TransportSamplingProfiles.md)。
以下检查采用新的空输出目录，保留源码快照；不启动 GitHub Actions。

```text
run_transport_production_validation.jl INPUT_DIR NEW_OUTPUT_DIR [pilot|remaining|all SHARD SHARDS]
run_transport_phase_revalidation.jl INPUT_DIR NEW_OUTPUT_DIR [ANCHOR_SOURCE]
run_transport_fortran_alignment.jl INPUT_DIR FORTRAN_EVIDENCE_DIR NEW_OUTPUT_DIR
```

生产复现调用真实 `run_phase_guided_scan`，经 mode-B 固定温度单点计划进入 `run_scan`，
复用 `joint_complete_summary_v1` 的冻结独立参考。最早的单点先导直接调用 `run_scan`，其源码快照独立保留。
它检查 17 道采样诊断及九项最终量，要求对参考仍在 1% 内；相对原已验收候选的复现差异上限是 0.01%。
`pilot` 为 μB=0、αT=1、ξ=−0.22；`remaining` 覆盖其余 11 点。这里的时间包含入口加载/JIT，
不能替代已有专门计时的性能结论。内存受限的本机应串行运行 Julia 数值批次；并发受干扰的记录单独标注。

相边界入口先使用新模型重新求 μB=0、450、900 MeV 的 ξ=0 温度锚点，再验证整个 mode-A 计划。
crossover 锚点采用轻夸克凝聚对温度导数的绝对峰；一阶锚点采用两分支等 Ω。
独立参考同时提高背景上限和节点。每条路径从相反端点分别续接两分支，
中点另用冷种子及邻近种子池审计稳定解；所有新坐标都检查稳定分支、质量和熵。
每次根求解预算为 40，Newton 与 trust-region 共用此预算，残差和物理性门槛不变；
`traces/` 保留每个方向的状态、质量、Ω 和残差。
可选 `ANCHOR_SOURCE` 只复用已通过锚点，并要求输入、配置及除当前分析脚本外的源码快照字节一致；
算法或相关源码改变后应重新核验锚点。
存在一阶边界时另保存固定温度下的 ξ 根和参考根。此检查只覆盖所生成的生产路径，不认证全局相图。

Julia–Fortran 入口读取旧 Fortran 的四个 `point.csv` 和 17 道 `w_ij*.dat`，
从归一化前后的率恢复密度，验证 τ 的通道组合及固定背景/τ 的末端 η。
背景比较单列匹配 Fortran 常数和有限 15 fm⁻¹ 的层次，以及新生产配置；
最后比较新配置的各道率、τ 和 η/s。外部读图 5% 与内部 1% 收敛目标分开报告。
旧 Fortran 的阈值上移和均匀网格不会被当作精确参考，原记录保持不变。

## 默认入口和回归迁移

`scan_relaxation_times_vs_T.jl`、`scan_transport_vs_xi_T150_muB800.jl` 现在是共享
`run_gap_transport_scan.jl` 的参数包装，输出使用统一的 gap CSV schema。
前者保留 `--mub-list`、`--out`、`--diagnostics-out`；后者保留固定 T/μB 和 ξ 区间/列表别名。
`run_manual_relaxation_scan_workflow.jl` 的数值设置统一通过 `--sampling-profile <name|path>`，
不再单独维护节点和有限末端上限选项。自定义配置不自动继承默认配置的收敛资格。

旧 mode-B v1 基线 CSV 保留为历史结果。完整 A^xi 切换前，
`baseline_thermal_weighted_transport_v1.csv` 保存联合验收时保留局部分支种子的候选；
`baseline_thermal_weighted_transport_v2.csv` 保护真实生产入口重新选择稳定背景后的行为。
两个点的九项量在首次 v1 回归中略超容差。固定两组背景分别重放三个相关通道后，各自复现原率，
同输入重复结果一致；因此单独建立生产入口基线，保留 v1、失败记录及诊断，未修改数值实现或容差。
完整 A^xi 切换后的当前回归使用 `baseline_thermal_weighted_transport_v3.csv`，
314/314 项及五次真实默认入口已通过；v2 的数值差异与原始失败记录继续保留。

```text
export_thermal_weighted_transport_baseline.py QUALIFIED_DIR INPUT_DIR NEW_BASELINE.csv
    [--production-replay REPLAY_DIR [--source-archive SOURCE.zip]]
```

导出器位于 `scripts/dev/`。当前生产基线使用 `--production-replay`，从真实的十二点产物读取数值，
核对源码/配置、原候选和独立参考；每项相对候选须在 0.01% 内、相对参考在 1% 内。
输出以 `candidate` 保存生产值，以 `qualified_candidate` 保留原候选，`reference` 保持不变；
已有 CSV/JSON 不会被覆盖。省略该选项时仍导出局部分支种子候选，不能当作生产入口的数值承诺。
导出成功后还需审查差异，并对受影响点做全新的默认入口复算。

无 `.git` 的生产复现使用已核验的 `source_archive`：导出器重新检查包清单副本、完整 ZIP 的
SHA256、归档逐文件哈希和结果源码快照，再要求实际扫描产物的 `git_commit` 为 `unknown`。
当前数值源码、配置、点集和通道预算仍须匹配。基线元数据保留归档身份及独立的 `head_base`，
不将工作树基底填作计算提交。归档搬迁后可用 `--source-archive` 指定同一份 ZIP；字节必须一致。

内部回归仍采用 `rtol=1e-6, atol=1e-8`，另检查最终量相对独立参考在 1% 内。
core 检查 ξ=+0.37 的代表点，full 检查十二点。

如果已运行当前源码的完整生产复现，可设置 `TRANSPORT_REGRESSION_REPLAY_DIR` 指向其结果目录，
让回归测试读取真实产物、核对当前数值源码/配置身份和通道预算，避免重复数值计算。
未设置时回归执行真实 phase-guided 入口；身份不一致时拒绝复用。

## Runner 输入与实际产物检查

```text
prepare_transport_production_inputs.py PHASE_DIR NEW_OUTPUT_DIR
    --replay-dir REPLAY_DIR [--replay-dir OTHER_REPLAY_DIR ...]
    --case-name NEW_CASE [--calculation-sha FULL_SHA] [--shard-points 20]
check_transport_production_artifact.py RESULT_DIR
    --profile PROFILE_TOML --calculation-sha FULL_SHA --tier dry-run|production
```

输入生成器要求整个新背景计划、锚点及十二个复现点通过，重新计算九项量的误差，
核对 profile 一致性、当前数值源码及冻结快照、已验收的点集和参考值。输出六点 runner 先导及暂定整批分片。
指定 `--calculation-sha` 时还验证该对象确为本地可读的提交，其实际数值文件集合及内容与验收代码一致。
Git 常规换行转换和无关文档修改不影响此检查；未知提交、旧实现或未提交的数值文件会被拒绝。
未选择含已审阅实现的完整提交时 `pilot_input_ready=false`；
该输入生成器写出的 `full_dispatch_ready` 始终为 false；成本复核通过后，另生成可执行的批次计划，保留原始输入。
生成输入不会触发工作流。

workflow 通过 `benchmark/relaxtime/transport_runner_timing.jl` 原样调用扫描入口；
`TRANSPORT_RUNNER_TIMING_PATH` 指定新的 `runtime_timing.toml`。该文件分别记录加载和扫描的
耗时、GC、JIT/重新编译秒数及累计分配字节，并记录实际 CPU、可见内存和 Julia/BLAS 线程数。
JIT 已包含在各阶段耗时中，不再次加到总时间；Julia 1.10 缺少编译计时字段时明确标为不可用，
当前 workflow 固定的 Julia 1.12.5 提供这些字段。扫描入口抛错仍以失败退出，并尽可能保留计时诊断。
`entrypoint_returned=true` 只表示入口正常返回，数值成功仍由下面的产物检查决定。

外层 GNU time 的 `process_time.txt` 包含完整扫描进程的 wall/CPU 时间和峰值 RSS（KiB），
包括 Julia 启动；依赖安装和其他 job 开销从 Actions step 记录读取。
累计分配不是峰值内存，单点/小批扫描进程的峰值也不是整项 job 所有步骤的峰值。

产物检查读取真实 CSV、诊断、有效配置和源码合同，要求结果行明确 `converged=true`，检查每点九项量与 17 道误差预算，
并将源码哈希与当前检出的计算源比较。续算时仅消费与成功结果行 `run_id` 相同的通道诊断；
旧失败尝试继续保留。dry-run 只证明计划生成，不能替代数值完成。

下载先导后，在 quote-remote 运行
`python -B scripts/analysis/relaxtime/review_transport_runner_pilot.py PILOT_CASE`。
case 须包含绑定输入、派发清单、固定计算源码包，以及每个 run 的 API 记录、原始 ZIP 和解压产物。
复核器按生产合同递归核对 `src/`、`scripts/relaxtime/`、PNJL 配置与物理配置下的 Julia/TOML 文件，
另核对根 Project/Manifest；这与资格输入生成器的较窄脚本范围不同。
`reviewed_pilot.json` 中的 `pilot_artifacts_passed` 和 `runner_cost_qualified` 均为 true 后才准备剩余批次。
每点预算为同模式先导整段 scan 最大耗时的两倍，每个 job 另预留两倍已观察到的非 scan 开销，
以 180 分钟目标确定分片。`steady_point_seconds=null` 表示当前包装没有独立的逐点稳态实测。

## 生产留档与远端检查

当前[生产参数留档](../../analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/production_plan_20261023.md)
记录 profile、各层积分、九条温度、两种模式的输入与归档要求。
mode-A 的 910 点和 mode-B 的 909 点输入已绑定计算提交
`2e5c411c716d5cbc524c3fc6cfb335fe05b39444` 和标签 `transport-full-rs-20261026-v1`。
七个先导 jobs 的 11 点、187 道诊断与成本复核已通过；剩余各 904 点分别按最多 14/15 点分为 72/63 片。
[实际复核与批次](../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/returned/reviewed_pilot.json)
已取回核验，剩余批次尚未派发。

长数值任务优先在 W 运行。使用现有 SSH 配置的 `quote-remote`、BatchMode 公钥和主机校验，
源码包解压到任务隔离目录，采用便携 Julia 1.12.5；保留源码包、有效参数和输出身份。
远端初期单个 Julia 数值进程，避免把逻辑线程数量直接换算成并发进程数。
Windows PowerShell 通过短的 UTF-16LE `EncodedCommand` 加载审核过的脚本内容，不修改系统执行策略；
较长的脚本先上传，避免命令行长度限制。

`relaxtime-manual-scans.yml` 的数值入口同样使用 `--sampling-profile`，不再传递旧节点或积分模式选项。
run-id 隔离目录中的结果、图、执行命令、失败日志和 `/usr/bin/time -v` 随 artifacts 保存。
W 上的小型自定义 profile 集成检查只验证入口与 sidecars；生产收敛沿用独立/联合验收，
runner 分片用时与峰值内存使用上述 GitHub runner 实测，后续批次继续保留相同计时和产物合同。
