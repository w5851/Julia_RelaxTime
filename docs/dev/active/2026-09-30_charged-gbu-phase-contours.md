# charged GBU 相图介子数密度等高图任务单

## 1. 背景与目标

组会要求用当前最新 charged infinite-thermal GBU 路线，在固定 quark-only
`FixedMuBConservedCharges` BQS 背景上构造 `(T, μB)` 相图中的介子数密度与
`K+/pi+` 等高图，为后续评估化学冻结线参数化提供诊断依据。图中应叠加参考
的一阶相变线、crossover 线和当前默认化学冻结线；Bose 支撑不安全或背景求解
失败的网格点必须显式记录并 mask，不得裁零或插值掩盖。

本任务在 task ledger 中是 `independent` 后续；不替换已 accepted 的公式闭合主线。

## 2. 范围与非目标

- 固定 `rho_Q/rho_B=0.4`、`rho_S=0`，不改 PNJLCore。
- 使用 `ChargedGBUResearchWorkflow.background` 与当前无限热 kernel；筛选阶段
  使用低分辨率设置，不把筛选结果称为完整 production。
- 至少生成 `pi_plus`、`K_plus` 和 `K+/pi+`；失败点保存失败原因。
- 网格坐标使用 MeV：横轴 `muB_MeV`，纵轴 `T_MeV`。
- 不包含介子反馈、强子化、共振衰变或实验最终产额；不更新 production baseline。
- 不在本任务中通过调参制造 horn，不裁掉负密度或失败点。

## 3. Benchmark 口径

上一版 benchmark 在每个独立 Julia 进程内调用旧的有限热
`CausalGBUResearch.density`，且首次 `Models.solve` 的 JIT 包含在背景计时中。
因此旧的约 40--55 s/点和数小时外推不能用于当前路线，保留为 historical
diagnostic。用户提出的“背景应为 ms 级”也不能直接假定：当前是 8 维
BQS 联合求解，并会尝试多个内部 seed；同进程实验显示第一次调用约几十秒，
warm/seed continuation 通常约 0.2--2 s，仍需远程重复测量。

新的 `benchmark/relaxtime/bench_charged_gbu_contour_point.jl` 已改为当前
`ChargedGBUResearchWorkflow` 无限热 kernel，使用 v4 schema 分别记录进程首调用、
同进程 no-seed 重复、seed 首调用与重复、低分辨率 density 首调用与重复。
默认背景重复 5 次、密度重复 3 次，报告均值/中位数/范围、编译/GC/分配字段和失败。
进程启动和 include/load 不计入内部样本；workflow 另外保存 process wall time。
每点都测密度。`screening` 值有限不代表完整 local/Mott/topology/eta/q-order gates；
另一个 workload 直接调用原 `channel_density`，读取原 production 配置并运行全部
smoke 门禁。只有门禁通过的重复样本可用于成功产点成本估算。
shell 的单位是 fm^-2，外层 q 权重积分后密度才是 fm^-3。
首次 seed 调用独立预热，重复时固定相同前一点 seed。背景状态、Omega、残差和 seed
随样本保留；密度固定使用该点无 seed 首调用的背景，不把分支变化混进积分性能。
compile 字段仅覆盖 `@timed` 内可见的工作，不能声称包含全部调用者 inference，也不
从总耗时中机械扣除；报告 warm 样本中观察到编译的次数。

本机仅做了一次短实现探针（Julia 1.12.5、同进程两个 continuation 点、低分辨率）：
冷首背景约 `30.28 s`，同点 warm 无 seed/seed 约 `1.31/1.43 s`，一个 continuation
点约 `0.40 s`，四个通道的 warm screening shell 约 `0.47--0.64 s/channel`。
这组数只证明计时分层和当前无限热入口可运行，不用于决定远程扫描步长；runner
硬件、重复样本和完整 production gates 尚未测量。

必须在 GitHub Actions 上运行重复样本后，才决定二维扫描步长；本机不运行完整
二维网格。

## 4. 建议初始范围

第一阶段保留已讨论的矩形：

```text
T      = 40:10:220 MeV
muB    = 0:50:800 MeV
```

这 19×17 个点只是工作范围候选，步长未由新 benchmark 决定前不启动全图。此前
文档把相变参考 CSV 中的 `mu_transition_MeV` 误读成 `muB`；该字段实际是 `muq`
（图中另有 `muB=3*muq` 转换），而且历史参考线是 equal-flavor 背景，不能直接
称作当前 BQS 的相变线。若要完整覆盖该历史线，应另行扩大到约 `muB=1200 MeV`
并明确背景差异。当前下游默认的
`data/reference/pnjl/issue130_phase_reference_v2/accepted/` 中，xi=0 Maxwell
参考覆盖 T=5--131 MeV、muB=873.586--1083.118 MeV（muB=3*muq）；若字面要求
包含整条参考一阶线，候选范围应为 T=5--220、muB=0--1200 MeV，而不是只扩大
muB 上限。该 equal-flavor 参考不认证 BQS 相变位置；BQS 附近的分支选择仍须复核。
本轮只进行代表点 benchmark，不暗中扩大或启动网格，步长待定。

## 5. 本轮交付入口

- `benchmark/relaxtime/bench_charged_gbu_contour_point.jl`：源码快照、前后 hash、
  冷/warm/seed 与两种密度 workload；每点原子 checkpoint，NaN 序列化为 null。
- `tests/unit/relaxtime/test_charged_gbu_contour_benchmark_contract.jl`：纯合成
  统计、失败保留、参数验证、取消与计时字段测试，不含速度阈值。
- `.github/workflows/relaxtime-charged-gbu-contour-benchmark.yml`：当前无限热
  kernel 的冷/warm/continuation benchmark artifact。screening 与 production
  两个 job，各用 Julia 1.12.5 和单线程，在包 instantiate/precompile 后运行。
  immutable checkout、run id/attempt、源码/环境快照、stdout、process wall time
  和失败时 artifact 均保留。初次只对当前 topic 分支窄路径 push 触发，手动入口
  可重复；不为执行新 workflow 自动合并 main。

本地 scanner、scan workflow 和 scan contract 仍是未提交草稿，未发布或触发。
恢复时 source/config 身份、逐点 checksum、失败重跑和聚合完整性还需补齐。

## 6. 后续任务

### M0：可靠 benchmark

- [x] 修正无限热入口、重复采样、JIT/GC 字段、源码前后 hash 和纯合成统计测试。
- [ ] 在 GitHub Actions 上完成至少一组 warm 重复与 continuation 样本。
- [ ] 报告冷首 JIT、warm 背景、seed continuation、density probe 的均值/中位数/
      范围，并据此决定 `T`、`muB` 步长和 shard 数。

### M1：远程筛选

- [ ] 用确定的步长启动 diagnostic shards。
- [ ] 合并点文件，验证每个网格 key 唯一、source hash 一致、失败点未被当作零。

### M2：图形与参考线

- [ ] 生成 `n_pi+`、`n_K+`、`K+/pi+` 热图/等高图。
- [ ] 叠加明确标注单位、BQS、diagnostic 状态的参考线；另画 mask/failed 图。

### M3：局部加密与完整验收

- [ ] 在冻结线、相变参考线、低温高 μB 区域和 ratio 梯度大的单元加密。
- [ ] 对冻结线代表点、等高线 extrema 和 mask 两侧点运行完整 production gates。
- [ ] 比较 screening 与完整 gates 点，记录是否改变 horn/单调性判断。

### M4：证据治理

- [ ] 生成 manifest、输入/源码 hash、行数、NaN/Inf/失败点摘要和图形清单。
- [ ] 只有作者审核后才提升为 accepted research artifact；不自动更新 baseline。

## 7. 风险与回退

- BQS solver 在一阶线附近可能有多支；保存 seed、残差和失败点，不拼接共存相。
- Bose 支撑不安全时 mask 并保留 gate 信息，不裁零。
- 低分辨率伪结构必须由 M3 完整 gates 点复核。
- GitHub Actions benchmark 未运行前，不给出可靠墙钟或步长结论。
