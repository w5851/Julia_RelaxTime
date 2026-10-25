# 新采样生产参数与重跑准备（2026-10-23）

此前的本地数值验收针对各向同性有效耦合口径；本轮已将生产候选切换为完整
\(A^\xi\) 有效耦合；2026-10-25 已按同一套门槛通过独立分层与十二点联合验证，
十二点真实生产入口及四个 ξ=0 Fortran 点也已通过。旧 v2 差异已审查，v3 回归 314/314、
跨目录回放 13/13、五次真实默认入口及 mode-A 的 910 点新背景验收已通过。
mode-B 的 909 点及 94 次冷种子审计、cheap/kernel/outer 三层通过；17 点联合运行的 153 项最终量、
578 条通道率全部通过，最大相对变化为 0.0720757%。真实入口的 17 点、153 项量及 289 道诊断通过，
相对候选最大变化为 6.29522e-9；完整返回包及 1584 份源码快照已取回并核验。
Actions 重跑准备尚未全部完成。当前已准备 mode-A 的 910 点及 mode-B 的 909 点输入，
共七个先导 jobs、十一点；尚未提交、推送或触发 Actions。
当前进度由[生产接入任务](../../../dev/active/2026-10-23_输运采样生产接入与重跑合同.md)维护。

2026-10-25 的 [A 解析降维原型评估](full_rs_a_reduction_20261025.md)在 quote-remote 通过；
现有 A 单次仅需约 25–33 微秒，暂不把其生产接入列为本轮验收前置。
本轮候选继续使用二维完整 A^xi；收敛、数值回归、同点复算和计时均在 quote-remote 运行，
其中 W 的计时只用于 W 的任务规划。本机不执行这些数值任务；GitHub runner 的成本按后述先导单独验收。

切换前的参数记录保留于[production_parameters.json](../../../../data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/production_plan_v1/production_parameters.json)，
保存当时配置、相界、Fortran 对齐参数和 runner 输入的路径及 SHA256。
完整 A^xi 口径的实际配置和坐标以新验收 manifest 及
`full_rs_phase_validation_v1/returned/phase_revalidation/sampling_plan.csv` 为准。

## 已完成的验收与仍需准备的事项

切换前口径的十二点真实生产入口、108 项最终量及 204 道采样诊断曾通过，最大独立参考差异为
0.6863883%；这些结果作为历史验收保留。切换为完整 \(A^\xi\) 后，必须重新确认
最终量 1% 目标、采样诊断、背景锚点和 Julia–Fortran 同点对齐，不能沿用旧基线作为新口径的收敛证明。
完整结果见[接入与外部对齐报告](production_integration_20261023.md)。

| 准备项 | 当前状态 |
| --- | --- |
| 新采样默认、旧实现清理、生产入口 | 完整 \(A^\xi\) 的独立分层、十二点联合、真实入口、四点外部对齐、v3 基线及默认入口验收通过 |
| phase-guided、manual scans、ci-integration workflow 的源码合同 | 已修正，三个 workflow 通过 actionlint；手动 argv 的实际数组构建检查通过 |
| W 隔离环境的手动入口集成检查 | 两点、23/23 断言及 34 道诊断通过；运行前后源码身份一致，结果已取回 |
| 实际计算提交及 GitHub 可取到的源码版本 | 已准备独立暂存区的提交范围预览，322 份数值文件及 121 份冻结文件检查通过；计时补充后仍须刷新预览，尚未形成实际提交 |
| 两种模式的 runner 输入 | mode-A 六点/四个 jobs、mode-B 五点/三个 jobs；`calculation_sha=null`，`pilot_input_ready=false` |
| GitHub runner 先导时间、启动/JIT、峰值内存及真实产物核验 | 待执行 |
| mode-A 整批分片 | 暂定最多 20 点/片、54 片；待 runner 实测后校正 |
| mode-B 新背景、共存分支及代表点覆盖，再生成批次 | 909 点背景、94 次审计、三层独立、17 点联合及真实入口通过；909 点输入已准备，暂定 54 片 |
| 持久归档、整批完整性核验及新图层 | 新 case 和归档合同已记录，待实际结果 |

这些是准备状态，不能以本地验收或 W 的运行时间代替 GitHub runner 的成本验收。

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

| μB [MeV] | 新 Tphase [MeV] | αT=1 的 T | αT=1.1 的 T | αT=1.2 的 T |
| --- | ---: | ---: | ---: | ---: |
| 0 | 200.0898076 | 200.0898076 | 220.0987884 | 240.1077691 |
| 450 | 182.7282025 | 182.7282025 | 201.0010227 | 219.2738430 |
| 900 | 125.7360003 | 125.7360003 | 138.3096003 | 150.8832003 |

mode-B 新背景及输入覆盖 T={120,160,200} MeV、μB={0,450,900} MeV、ξ=−0.5:0.01:0.5，
九条各 101 点，共 909 点。T=120 MeV、μB=900 MeV 的一阶边界在 [−0.13,−0.12] 内，
离散网格没有精确共存点；边界两侧保留并完成代表点输运验收。
两种模式的 1819 点输入已具体准备，当前输运验收仍是代表点覆盖，不是 1819 点输运已完成。

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
1639 份源码文件在运行前后不变，返回包已取回并核验。此时尚未获得 GitHub runner 成本资格。

[runner 复核器](../../../../scripts/analysis/relaxtime/review_transport_runner_pilot.py)的六项行为检查已在 W 通过。
它将实际点集、17 道预算、计算/工作流来源、完整进程时间及 RSS 与具体派发请求核对。
当前包装没有逐点稳态计时；分片计算预算采用先导整段扫描最大用时的两倍作为每点预算，
不除以先导点数，同时为每个 job 留出两倍已观察到的非扫描开销。
以 180 分钟规划目标校正每片点数，并要求先导峰值 RSS 不超过主机内存的 80%。
这是一组先导样本支持的保守规划，不是全网格每一点的成本上界；JIT 已包含在阶段耗时中，不重复相加。

## 留档与完成条件

每片保留 calculation SHA、workflow SHA、run ID、profile、实际生效配置、计划 CSV、结果 CSV、
17 道诊断、失败证据、执行命令和完整进程时间/峰值内存。
在 Actions artifacts 的 30 天保留期内取回到新的持久 case，核对源码/配置、完整点集、重复点和失败点后合并。
新数值通过审阅后再产生新的论文图层；历史生产与当前图层继续可追溯。

本轮准备完成的判据是：实际可取到的计算版本、commit 绑定输入、runner 成本与产物先导通过，
mode-A/mode-B 批次及归档方式全部具体可执行。当前 `pilot_input_ready=false`、`full_dispatch_ready=false`。
