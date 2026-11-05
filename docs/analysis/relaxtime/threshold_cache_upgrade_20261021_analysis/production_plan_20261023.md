# 新采样生产参数与整批执行（2026-10-23）

## 当前状态与交付范围（2026-11-05）

本阶段的输运重计算、正式数值验收及计划内有限窗口分析已完成，
[PR #331](https://github.com/w5851/Julia_RelaxTime/pull/331)已合并到 `bad5f98ce99ffd81a9ff0f60336a6c3c32a3d44a`。
适用CI已通过，包括修复后的benchmark script smoke；SSH公钥权限已通过实际推送核验。
工程交付范围与恢复方式见[交付入口](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/engineering_delivery_20261104_v1/README.md)。

后续交付在 `D:/w/transport-v15-paper`、`codex/transport-v15-paper` 推进：v15的76张PDF已完成，
原PNG逐字节复现，PDF字体、矢量、尺寸及三张论文组合图视觉检查通过；作者已授权正式晋升、current及论文替换。
详见[v15正式交付](../../../../data/outputs/results/relaxtime/transport/phase_guided/publication_clean_v15_pdf/README.md)。
论文已接入所需三张输运PDF和最小数据依据，15页预览、41项引用和逐页视觉核验通过；新正文为 `agent_review`，独立相图沿用原来源。
消费者迁移后已退役557张旧PNG及论文279份旧附件，均保留已有Git恢复来源；旧树39787个冗余未跟踪文件已按外部归档清单清理，8处tracked残余已核对保留去向并恢复，旧树Git状态干净。当前会话仍占用 `D:/w/meson-freezeout-defaults`，该worktree待切换会话目录后回收。
既有外部归档继续复用，不重复保存整套数值证据。
历史记录中的 pending 只描述生成时的阶段，不推翻后续晋升记录。

| 对象 | 当前状态与入口 |
| --- | --- |
| 计算实现 | 完整 A^xi、thermal-weighted 采样、cm 相对速度和 v3 回归已验收；计算源码固定为 `2e5c411c716d5cbc524c3fc6cfb335fe05b39444` |
| 原网格 | 142 个 jobs、1819 点、30923 道诊断完成并正式晋升；8 项准入证据见[晋升记录](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_full_grid_v1/formal_promotion_v1/promotion.json) |
| 最新正式数值 | [η 加密版本](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_eta_augmented_20261029_v1/README.md)：1867 点、31739 道记录；1810 个原值、9 个原坐标重算值及48 个新增点 |
| 当前正式图层 | [v15 PDF及原PNG引用](../../../../data/outputs/results/relaxtime/transport/phase_guided/publication_clean_v15_pdf/README.md)，对应1867点；[current](../phase_guided_transport/publication_clean_current.json)已切换 |
| 论文交付 | 三张v15输运组合图及独立相图组成四图审阅稿；新科学解释与文字状态见论文 `writing_package/00_status.md`，不继承旧正文批准 |
| Mott 分析 | π、K、η、σ 的计划内检查完成；13 个窗口及175 个参考求值、过程/物种分解见[Mott 分析](mott_local_bumps_20261028.md) |
| σ 最新补检 | 连续近阈混合标量平方积分的左侧 h^(-1/2)、右侧对数段及误差对照完成；完整 τ 和末端输运量的高阶导数未据此认证 |

仍保留的研究问题按下表限定，不作为既定数值交付的未完成计算。

| 分类 | 未决问题及再次启动条件 |
| --- | --- |
| research | 定量使用两处 π 最窄尖峰时，补更细 ξ 及独立外层/t/B₀ 积分检查；现有 s 表加严不能认证连续导数极限 |
| research | 要宣称完整 τ、η、ζ、σ_el 均为 C¹ 且非 C²，须逐项证明背景、其他通道、边界和非解析系数不抵消；当前解析结论有条件，σ 数值对象为局部积分 |
| research | 混合谱支未决、严格退耦和采样间双穿越/切触零尚未全面排除；“未见 Mott”只在已查扫描分辨率下成立 |
| research | 认证极细窗口高阶导数时，处理生产行列式保护的人工曲率、离散热权重/硬边界及右侧最内层精度；此次未改变生产核 |
| independent | [A 解析降维](full_rs_a_reduction_20261025.md)原型已验，生产接入与端到端收益仍独立待办 |
| independent | v15 PDF及current已交付；论文15页审阅稿编译通过，旧PNG、论文旧附件及旧树冗余副本已退役。保留v15渲染依赖和唯一恢复链，原会话worktree待切换cwd后回收 |

以下保留各批次的实施与验收依据；涉及历史图层或准备阶段的语句不代替上表的当前状态。

## 实施与验收依据

此前的本地数值验收针对各向同性有效耦合口径；本轮已将生产候选切换为完整
\(A^\xi\) 有效耦合；2026-10-25 已按同一套门槛通过独立分层与十二点联合验证，
十二点真实生产入口及四个 ξ=0 Fortran 点也已通过。旧 v2 差异已审查，v3 回归 314/314、
跨目录回放 13/13、五次真实默认入口及 mode-A 的 910 点新背景验收已通过。
mode-B 的 909 点及 94 次冷种子审计、cheap/kernel/outer 三层通过；17 点联合运行的 153 项最终量、
578 条通道率全部通过，最大相对变化为 0.0720757%。真实入口的 17 点、153 项量及 289 道诊断通过，
相对候选最大变化为 6.29522e-9；完整返回包及 1584 份源码快照已取回并核验。
Actions 的七个先导 jobs、十一点及后续 135 个新分片均已完成并验收。
2026-10-26 的[完整网格验收](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_full_grid_v1/final_review_v1/full_grid_acceptance.json)
通过 mode-A 的 910 点、mode-B 的 909 点和全部 30923 道诊断；完整背景与局部变化审计通过，
合并表、逐点来源及原始归档已取回并核对哈希。[生产接入任务](../../../dev/archived/2026-10-23_输运采样生产接入与重跑合同.md)保留完整实施记录，当前状态以本页顶部为准。

本批于 2026-10-28 按作者明确授权完成[数值晋升](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_full_grid_v1/formal_promotion_v1/promotion.json)
及 [publication_clean_v14 图层正式化](../../../../data/outputs/results/relaxtime/transport/phase_guided/publication_clean_v14_pdf/formalization.json)。
交付 76 对彩色／灰度 PNG 和 76 份矢量 PDF，补齐固定 T／FIG4；原 1819 点无平滑、删点或局部替换。
μB=900 的 π 与 η Mott 窗口补查在 quote-remote 完成，现有积分加密证据支持保留其快速响应。
旧版本从 current 入口退役，原数据、图像和冻结审计保持字节；物理清理见[清单](../../../../data/outputs/results/relaxtime/transport/phase_guided/publication_clean_v14_pdf/legacy_cleanup/README.md)。

2026-10-29按作者同意的方案完成[η局部加密正式数值版本](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_eta_augmented_20261029_v1/README.md)：
1810个原值不变、9个原坐标重算、48个新ξ点，共1867点及31739道记录。
回放、接缝与逐点选源审计均在quote-remote通过，既有源码和回归资格保持适用。
v15采用v14形式的PNG已由作者审阅，PDF当时按作者指示暂缓，已于2026-11-05完成交付和晋升。2026-10-31完成的其他Mott有限窗口诊断涵盖全部主网格及175个参考点、同端点精度、ξ分割和同宽度过程／物种分解；K也存在斜率响应。[当时审阅记录](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/mott_other_slopes_20261029_v1/scientific_review_v1.json)保留原始范围。11月4日收束的σ连续近阈积分补检、条件性光滑性分析及仍未认证的完整输运量结论见[Mott分析](mott_local_bumps_20261028.md)。这些诊断不替换正式数据或current。

2026-10-25 的 [A 解析降维原型评估](full_rs_a_reduction_20261025.md)在 quote-remote 通过；
现有 A 单次仅需约 25–33 微秒，暂不把其生产接入列为本轮验收前置。
本轮候选继续使用二维完整 A^xi；收敛、数值回归、同点复算和计时均在 quote-remote 运行，
其中 W 的计时只用于 W 的任务规划。本机不执行这些数值任务；GitHub runner 的成本按后述先导单独验收。

切换前的参数记录保留于[production_parameters.json](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/production_plan_v1/production_parameters.json)，
保存当时配置、相界、Fortran 对齐参数和 runner 输入的路径及 SHA256。
完整 A^xi 口径的实际配置和坐标以新验收 manifest 及
`full_rs_phase_validation_v1/returned/phase_revalidation/sampling_plan.csv` 为准。

## 本轮验收与交付状态

切换前口径的十二点真实生产入口、108 项最终量及 204 道采样诊断曾通过，最大独立参考差异为
0.6863883%；这些结果作为历史验收保留。切换为完整 \(A^\xi\) 后，必须重新确认
最终量 1% 目标、采样诊断、背景锚点和 Julia–Fortran 同点对齐，不能沿用旧基线作为新口径的收敛证明。
完整结果见[接入与外部对齐报告](production_integration_20261023.md)。

| 准备项 | 当前状态 |
| --- | --- |
| 新采样默认、旧实现清理、生产入口 | 完整 \(A^\xi\) 的独立分层、十二点联合、真实入口、四点外部对齐、v3 基线及默认入口验收通过 |
| phase-guided、manual scans、ci-integration workflow 的源码合同 | 已修正，三个 workflow 通过 actionlint；手动 argv 的实际数组构建检查通过 |
| W 隔离环境的手动入口集成检查 | 两点、23/23 断言及 34 道诊断通过；运行前后源码身份一致，结果已取回 |
| 实际计算提交及 GitHub 可取到的源码版本 | `2e5c411c716d5cbc524c3fc6cfb335fe05b39444` 与固定标签 `transport-full-rs-20261026-v1` 已发布并核验；322 份数值文件及 159 份冻结文件符合最终提交预览 |
| 两种模式的 runner 输入 | 输入已绑定计算提交；mode-A 六点/四个 jobs、mode-B 五点/三个 jobs 全部完成 |
| GitHub runner 先导时间、启动/JIT、峰值内存及真实产物核验 | W 复核通过：11 点、187 道诊断、实际参数、来源及计时均通过；返回包已取回核验 |
| mode-A 整批分片 | 六个先导点加 72 片，完整 910 点已成功取回、验收并转正 |
| mode-B 新背景、共存分支及代表点覆盖，再生成批次 | 五个先导点加 63 片，完整 909 点已成功取回、验收并转正 |
| 持久归档、整批完整性核验及新图层 | 1819 点来源及背景审计通过；publication_clean_v14 原始 PNG/PDF 与 FIG4 完整交付，current 已更新 |

两种模式的 135 个新分片已全部成功并取回，实际用量见下文整批统计。
本批沿用作者直接执行、取回并验收整批的授权，没有追加正式长度先导。

## 物理参数与积分配置

采用 PNJL `default` 模型和 `default` 物理配置：Λ=602.3 MeV，mud0=5.5 MeV，ms0=140.7 MeV，
GΛ²=1.835，KΛ⁵=12.36；log Polyakov 势的 T0=210 MeV。
μu=μd=μs=μB/3，传播子 ξ 跟随热力学 ξ，计算 bulk。
ħc=197.3269804 MeV·fm。末端 η/s、ζ/s、σ/T 和六种 τ 使用同一实际背景。

采样配置为[thermal_weighted_v1.toml](../../../../config/models/pnjl/transport_sampling/thermal_weighted_v1.toml)，
SHA256 为 `fd7e9a2ef1c3a313e7621be0d3213ce0cf03ce3f48a3f64505066166ffb2afca`。

| 层次 | 动量积分 | 节点与尺度 |
| --- | --- | --- |
| 背景热积分 | 有限 [0,25] fm⁻¹ | p96/cos24；残差≤1e−8 |
| A | 半无穷 [0,∞) | p32/cos16，映射尺度 5 fm⁻¹；极化 full-RS、耦合 full-\(A^\xi\) |
| τ 热平均 | 半无穷 [0,∞) | p48/cos36/phi36，映射尺度 10 fm⁻¹；截面 t24 |
| 热平均内部数密度 | 半无穷 [0,∞) | p32/cos36，映射尺度 10 fm⁻¹ |
| workflow 数密度 | 有限 [0,25] fm⁻¹ | p64/cos8 |
| η、ζ、σ 末端积分 | 半无穷 [0,∞) | p32/cos16，映射尺度 5 fm⁻¹ |

25 fm⁻¹ 是背景热积分和 workflow 数密度的上限；散射支持 cutoff 为 Λ≈3.0522942113 fm⁻¹。
对每个入射/出射道，截面只在由阈值 \(s_{\rm bo}\) 到 cutoff 动量确定的
\(s_{\rm up}\) 内有支持，区间外贡献为零。半无穷映射的尺度不代表积分上限，
也不移除该散射支持限制；它与 NJL 真空正则化分属不同层次。

热权重缓存采用 rtol=1e−3、atol=1e−12、seed_order=48、weighted_window=0.6，
每道预算为 max_evaluations=40000、max_nodes=12000、max_rounds=40、max_events=30000000。
按已验收的 17 道顺序逐道重置极化缓存；预算耗尽保留失败点，不输出有效 τ。
每次根求解尝试最多 40 步。缓存局部容差、最终量 1% 收敛、外部读图 5% 和内部回归容差各有独立用途。

## 生产范围

mode-A 固定 μB={0,450,900} MeV，αT={1,1.1,1.2}，T=αT·Tphase(μB,ξ=0)。
ξ=−0.5:0.01:0.5；μB=900 MeV、αT=1 用 ξ=−0.003、+0.003 替换 ξ=0，保留共存两侧，
该曲线 102 点，其余八条各 101 点，总计 910 点。

| μB [MeV] | 复核 Tphase [MeV] | αT=1 的 T | αT=1.1 的 T | αT=1.2 的 T |
| --- | ---: | ---: | ---: | ---: |
| 0 | 200.0898076 | 200.0898076 | 220.0987884 | 240.1077691 |
| 450 | 182.7282025 | 182.7282025 | 201.0010227 | 219.2738430 |
| 900 | 125.7360003 | 125.7360003 | 138.3096003 | 150.8832003 |

mode-B 新背景及输入覆盖 T={120,160,200} MeV、μB={0,450,900} MeV、ξ=−0.5:0.01:0.5，
九条各 101 点，共 909 点。T=120 MeV、μB=900 MeV 的一阶边界在 [−0.13,−0.12] 内，
离散网格没有精确共存点；边界两侧保留并完成代表点输运验收。
两种模式的完整 1819 点输运已完成并通过本批产物验收；背景与既有合格背景的
最大质量/熵相对差为 `2.17057e-8`/`1.91574e-8`，可区分稳定分支不一致为零。
表中的 Tphase 表示当前配置下重新核验的温度锚点，本批背景证据覆盖这些生产路径。

## 计算资源与 runner 先导

长数值任务优先使用 W。现有连接是 `quote-remote`，配置文件为
`D:/Desktop/CodexBountyLab/executions/cross-chain-quote-observer-2026-10-15/tunnel-setup-2026-10-16/ssh_config`。
沿用配置中的公钥、BatchMode 和主机校验，经堡垒机连接 W；远端身份为 `W\82306`。
无需修改现有 SSH 配置或另用失效的 `remote-julia` 别名。

W 为 i5-14400F、10 核/16 逻辑线程、约 16 GiB 内存，探测时可用约 8.4 GiB。
初期一次只运行一个 Julia 数值进程。使用隔离目录
`D:/CodexCompute/Julia_RelaxTime/transport-preflight-20261023-v1`，便携 Julia 固定为 1.12.5，
系统默认 1.11.6 保持原状。源码包及逐文件身份见
[package_manifest.json](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/remote_preflight_setup_v1/package_manifest.json)。
W 用于独立收敛补算、背景审计、机制诊断及非 Actions 集成检查。
本次手动扫描集成检查使用两点小型自定义 profile，验证实际 CSV/sidecars，不能作为生产配置的收敛检查。
检查已经通过，实际两个点为 `(T,μB,ξ)=(150,0,0)、(190,0,−0.5)`，bulk 按测试要求关闭。
[远端结果复核](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/remote_manual_preflight_v1/reviewed_result.json)
核对 23/23 测试日志、退出码 0、实际成功行和 34 道采样诊断，源包中的 879 个文件在运行前后身份一致。
完整测试进程约 412 秒，包含首次包预编译和入口/JIT，不能按两个点均摊后作为生产用时。

原包装脚本受 Julia 顶层 soft scope 影响，把单独摘要的通过标记写为 false；原记录保留，
独立复核明确更正依据。后续包装使用 Ref，并已单独检查标记写出，无须重复两点数值计算。
远端本次采用冻结源码包，manifest 的 `git_commit=unknown` 如实表示没有计算提交；不填造假的 SHA。

GitHub runner 固定 Julia 1.12.5、Julia/OpenBLAS 各单线程、heap hint 1200M、job 上限 360 分钟。
mode-A 先导使用相同生产配置，以四个 jobs 检查六点：
`(μB,αT,ξ)=(0,1,−0.22)、(0,1,+0.37)、(0,1.2,−0.5)、(450,1.1,+0.35)、(900,1,−0.04)、(900,1,−0.003)`。
mode-B 先导以三个 jobs 检查五点：`(T,μB,ξ)=(120,0,+0.5)`、
`(120,900,−0.5/−0.13/−0.12)`、`(200,900,−0.5)`；其中低温高密度 job 实际跨过已验证的一阶边界。
每片计算负载暂按约 180 分钟规划，为启动和波动留空间；20 点/片不是最终预算。

先导必须测完整路径时间、启动/JIT、峰值内存并检查实际产物。
W 的时间可用于 W 的任务规划，GitHub runner 的整批预算仍从 runner 先导确定。
workflow 已接入原样调用 CLI 的计时包装，记录 `runtime_timing.toml` 中的加载/扫描及其 JIT，
以及外层 `process_time.txt` 中的完整扫描进程时间和峰值 RSS；job 安装等开销另从 Actions 记录获取。
2026-10-26 的 workflow 静态检查通过；包装在 W 上的 29 项单测和两点真实 dry-run 也已通过，
1639 份源码文件在运行前后不变，返回包已取回并核验；实际 runner 验收见下表。

[runner 复核器](../../../../scripts/analysis/relaxtime/review_transport_runner_pilot.py)修订版的七项行为检查已在 W 通过。
它将实际点集、17 道预算、计算/工作流来源、完整进程时间及 RSS 与具体派发请求核对。
当前包装没有逐点稳态计时；分片计算预算采用先导整段扫描最大用时的两倍作为每点预算，
不除以先导点数，同时为每个 job 留出两倍已观察到的非扫描开销。
以 180 分钟规划目标校正每片点数，并要求先导峰值 RSS 不超过主机内存的 80%。
这是一组先导样本支持的保守规划，不是全网格每一点的成本上界；JIT 已包含在阶段耗时中，不重复相加。

实际计算和 workflow 均固定到 `2e5c411c716d5cbc524c3fc6cfb335fe05b39444`，
执行 ref 为 `transport-full-rs-20261026-v1`。七份原始 API 记录、ZIP、解压产物和派发请求保留于
[先导目录](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/dispatch_manifest.json)。

| 模式与先导曲线 | run | 点数 | 完整扫描进程 [s] | Actions job [s] |
| --- | --- | ---: | ---: | ---: |
| A，μB=0、αT=1 | [40224161650](https://github.com/w5851/Julia_RelaxTime/actions/runs/40224161650) | 2 | 373.16 | 809 |
| A，μB=0、αT=1.2 | [40224178198](https://github.com/w5851/Julia_RelaxTime/actions/runs/40224178198) | 1 | 272.83 | 718 |
| A，μB=450、αT=1.1 | [40224194586](https://github.com/w5851/Julia_RelaxTime/actions/runs/40224194586) | 1 | 164.78 | 452 |
| A，μB=900、αT=1 | [40224212631](https://github.com/w5851/Julia_RelaxTime/actions/runs/40224212631) | 2 | 358.26 | 732 |
| B，T=120、μB=0 | [40224228254](https://github.com/w5851/Julia_RelaxTime/actions/runs/40224228254) | 1 | 210.12 | 649 |
| B，T=120、μB=900 | [40224250574](https://github.com/w5851/Julia_RelaxTime/actions/runs/40224250574) | 3 | 347.60 | 641 |
| B，T=200、μB=900 | [40224267466](https://github.com/w5851/Julia_RelaxTime/actions/runs/40224267466) | 1 | 228.98 | 630 |

[完整复核](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/returned/reviewed_pilot.json)
通过 11 点、187 道诊断与 327 份递归生产源码身份，最大进程 RSS 为 1,808,244,736 字节，
约占可见主机内存的 10.8%。runner CPU 包含 znver3、znver5 和 sapphirerapids，预算采用同模式最慢的先导。
第一版离线复核误用了较窄的 322 文件资格输入范围；新增的五份合同文件均在原冻结包中且哈希一致。
v2 校正枚举范围后通过，原复核器及失败证据保留；数值源码和七份原始计算产物没有改变。

| 模式 | 最大先导 scan [s] | 每点规划预算 [s] | 每 job 预留开销 [s] | 剩余点/片数 | 最多点/片 |
| --- | ---: | ---: | ---: | --- | ---: |
| A | 341.929 | 683.859 | 957.964 | 904 / 72 | 14 |
| B | 327.099 | 654.198 | 940.242 | 904 / 63 | 15 |

[批次清单](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/returned/remaining_batches/batches_manifest.json)
绑定复核报告和每份请求的 SHA256；mode-A 的六点、mode-B 的五点以原 artifact 身份复用。
复用点与剩余点无交集，并集分别等于完整的 910/909 点网格。
按上述余量相加，[剩余批次规划](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_runner_pilot_v1/cost_summary.json)
合计 371.615 runner 小时，单片最大规划时间分别为 175.533/179.220 分钟。
这是一份保守资源计划；实际运行用量见下节，费用由账户计费记录确定。

## 完整批次实际结果

135 个新分片全部成功，连同七个先导的完整网格在 quote-remote 通过逐点、17 道预算、
源码/配置、完整性、背景和相邻变化审计。[最终结果与合并表](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_full_grid_v1/README.md)
保留逐 run 原始证据及逐点来源。10 项超过相邻变化筛查线的观测量对均位于已验证的一阶边界，
非边界筛查项为零；470 条有限 τ 比值标记核对后原样保留。

| 运行集合 | jobs | 点数 | 累计 job 小时 | 最长 job [s] |
| --- | ---: | ---: | ---: | ---: |
| A 新分片 | 72 | 904 | 34.0192 | 2558 |
| B 新分片 | 63 | 904 | 32.6736 | 2523 |
| 复用先导 | 7 | 11 | 1.2864 | 809 |
| 全部 | 142 | 1819 | 67.9792 | 2558 |

[实际计时汇总](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_full_grid_v1/final_review_v1/runner_accounting.json)
记录 135 个新 jobs 合计 66.6928 runner 小时，中位 job 为 1912 秒，最长为 42.63 分钟，
均低于 180 分钟规划目标和 360 分钟上限。扫描进程峰值 RSS 为 1,945,735,168 字节
（约 1.81 GiB），占对应主机可见内存约 11.61%。RSS 不覆盖所有 job 步骤，JIT 已包含在阶段耗时中；
累计 job 时间不等于账户账单，也不替代逐点稳态实测。

两次 POST 网络中断在完整运行列表核对缺失后使用原冻结请求单次补提交，
每个请求最终恰好对应一个 run。首个新 job 至最后一个结束历时 4 小时 58 分 10 秒，
包含派发恢复间隙。冻结请求、已成功数值、原始失败记录和返回归档均保留。

## 原批次留档与阶段完成条件

每片保留 calculation SHA、workflow SHA、run ID、profile、实际生效配置、计划 CSV、结果 CSV、
17 道诊断、失败证据、执行命令和完整进程时间/峰值内存。
在 Actions artifacts 的 30 天保留期内取回到新的持久 case，核对源码/配置、完整点集、重复点和失败点后合并。
新数值已用于生成与 v13 同构的原始 PNG 审阅图包：
[full_rs_raw_20261027_png_review](../../../../data/outputs/results/relaxtime/transport/phase_guided/full_rs_raw_20261027_png_review/README.md)。
75 对彩色/灰度 PNG、原始顶点和全部图像质量检查通过，17 项定向绘图测试通过；
仅两处已审计一阶不连续断线，1819 点和全部质量提示原样保留。
作者于 2026-10-28 大致审阅后要求检查 μB=900 的两处局部小凸起。
[Mott 局部诊断](mott_local_bumps_20261028.md)已在 quote-remote 完成：21 点的 π 反向跨阈分别位于
ξ∈(−0.05,−0.04) 和 (0.19,0.20)，K 跨阈稍早，混合模态未见新跨阈；A 升阶后区间不变。
随后按作者要求重点检验平均散射率热积分：六点的 24 项 s 插值和 12 项主导通道外层积分
均在 quote-remote 通过本次判据，返回包 406 个文件已核对。两处局部最大值和原网格峰位保留，
凸起高度分别变化 −0.132%/+0.587%，电导率最大修正 0.0169%；所查积分误差不足以解释凸起。
通道分解与局部收敛支持阈值相关结构。π/K 分离假说的已有数据初查显示：K 跨阈时总电导率
仍下降，π 跨阈时明显上升；峰后六种 τ 增加而电导率回落，严格归因须同时区分末端权重变化，
当前按作者讨论方向收束局部数值排查，以 π 阈值相关响应承接论文中的公式链条与物理解释，
π/K 分离保留为次要线索；Δξ 网格相对位置影响表观折角仍为合理但未验证的解释。
更细 ξ、末端权重分解和分振幅对照仅为更强论断的可选补证。
A 解析降维仍是独立待办，不改变本批二维 A 口径。原始数值与 PNG 保持冻结，
不将作者的初步视觉反馈解释为 PDF 接受记录或论文资格。
上述原始 PNG 批次保留生成时的审阅状态；后续 v14、v15 PDF 与图层晋升已完成，具体以本页顶部当前状态为准。

本轮准备完成的判据是：实际可取到的计算版本、commit 绑定输入、runner 成本与产物先导通过，
mode-A/mode-B 批次及归档方式全部具体可执行。这些准备条件已满足：绑定先导输入可用，
经成本复核后的独立批次计划为 `full_dispatch_ready=true`。原始输入草案保持冻结。
135 片已由[冻结队列](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_full_grid_v1/queue_plan.json)
执行完毕，调度时最多保持 20 个未结束任务；这是本次控制器的并发上限。
[收集状态](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_full_grid_v1/collection_status.json)
记录全部成功和已返回状态。完整性、数值来源、逐点诊断与最终背景审计已在 quote-remote 通过；
全网格数值验收、正式晋升和论文使用保持分别记账。
