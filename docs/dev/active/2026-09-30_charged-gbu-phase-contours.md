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

本机不承担数值扫描。GitHub Actions run `36832410627`（Julia 1.12.5、单线程、
screening 六点）给出：进程首背景 `15.98 s`，同点 warm no-seed/seed 中位数
`0.84/0.87 s`；其余点 warm 背景约 `0.15--1.73 s`。低分辨率四通道 screening
density 每点约 `1.02--1.06 s`，benchmark 命令的 process wall time 为 `142.46 s`；
包 instantiate/precompile 是 workflow 的独立步骤，不计入这些点成本。

同一 run 的 production job 在 60 分钟上限前只完成首点：四通道完整 smoke gate
通过，density 首调用 `842.82 s`、一次重复 `837.58 s`，第二点尚未完成，任务被取消。
因此 screening 秒级成本不能外推 full production。

随后 GitHub Actions run `36840120193` 对代表点 `(T,muB)=(165.923,23.525) MeV`
进行两次完整 production density 重复：背景 cold `16.26 s`、warm no-seed/seed
`0.848/0.844 s`；density 首调用 `835.05 s`，两次重复中位数 `830.11 s`
（`829.67--830.56 s`），三者及四通道均通过 `full_smoke_production_gates`，
首调用与重复差约 `0.59%`。该结果把首调用初始化与 production 稳态成本分开，
但仍只是一个代表点，不构成完整网格收敛。

上述结果说明背景 warm 求解不是历史意义上的毫秒级单一 gap solve：当前入口是
8 维 `FixedMuBConservedCharges` 联合求解，含多候选 seed、BQS 热力学积分和 AD
Jacobian；production density 还包括 q-order、local/Mott/topology、eta 与 tail
等完整门禁。`@timed` 的 compile 字段只覆盖表达式内部可观察的编译，不能用来从
总耗时机械扣除 JIT。

screening benchmark 已足以支持先按 `T=40:10:220 MeV`、`muB=0:50:800 MeV`
启动 diagnostic 网格；正式 production 仍必须按单点/单通道 shard 运行。二维网格
和所有数值任务均通过 GitHub Actions，不在本机运行。

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
  和失败时 artifact 均保留。workflow 只保留 `workflow_dispatch`，因此每次运行
  都可指定网格、分片和 tag；不为执行新 workflow 自动合并 main。

scanner、scan workflow 和 scan contract 已提交到 topic branch，并通过 GitHub Actions
run `36849353973` 完成一次 4-shard diagnostic scan。scan workflow 现在只保留
`workflow_dispatch`；本轮用于触发扫描的 topic-branch `push` 已移除，避免后续提交
自动重跑整张网格。aggregate job 还会用 solver-free Python 绘图器生成
`n_pi+`、`n_K+`、`K+/pi+` 热图和三联失败/mask 图；绘图器只接受精确 screened
rows，失败行保留为 mask，不插值、不零填充。远程 artifact 保留 point JSON、CSV、
manifest、聚合索引和 PNG/plot manifest；本地只下载到临时审计目录，不写入仓库。

当前绘图器还默认叠加化学冻结线、crossover 和 Maxwell 参考线，并在 plot
manifest 中记录 `rho_Q/rho_B=0.4`、`rho_S=0`、无介子反馈、单位及 `muB=3 muq`
坐标转换。phase reference 是 equal-flavor 历史参考，只作导向线，不能宣称为当前
BQS 相变边界；超出 screening 域的线段保持空白。候选选择器
`select_charged_gbu_contour_candidates.py` 读取合并 screening CSV，按冻结线、
`K+/pi+` 梯度极值、mask 边界和相变参考线选取稀疏点，失败点不派发。新的
`relaxtime-charged-gbu-contour-full-gate.yml` 按候选点并行运行四通道完整
local/Mott/topology/eta/tail/q-order gate，每点独立保留 provenance；其结果仍为
diagnostic，不更新 production baseline。GitHub 的 `workflow_dispatch` 只会在默认
分支已注册 workflow 后提供该入口，因此本轮 topic branch 的实际 12 点 dispatch
使用了仓库已有的 benchmark workflow `production` workload；新 workflow 保留为
合并后可复用的稀疏入口。

## 6. 后续任务

### M0：可靠 benchmark

- [x] 修正无限热入口、重复采样、JIT/GC 字段、源码前后 hash 和纯合成统计测试。
- [x] 在 GitHub Actions 上完成 warm 重复、seed continuation 与代表点完整
      production 重复（runs `36832410627`、`36840120193`）。
- [x] 报告冷首 JIT、warm 背景、seed continuation、screening density 与 full-gate
      density 的中位数/范围；初始 diagnostic 网格采用 `10/50 MeV` 步长，使用
      远程 T-row shards。production 步长不由 screening 成本推断。

### M1：远程筛选

- [x] 用确定的步长启动 diagnostic shards（`T=40:10:220 MeV`、
      `muB=0:50:800 MeV`，4 个 T-row shards）。
- [x] 合并点文件，验证 323/323 个网格 key 唯一、每个 shard 的源码 hash
      一致，且失败点未被当作零；分片点数为 `85/85/85/68`，321 点为
      `screened`，2 点（`T=160 MeV, muB=200/250 MeV`）因 `K±` 的
      `static instability or unresolved Bose endpoint` 为 `gate_failed`。
- [x] 保留 run `36849353973` 的远程 artifacts 作为 diagnostic evidence；该
      结果不是 full production、不是 baseline，也不授权自动绘图或物理晋升。

### M2：图形与参考线

- [x] 通过 scan aggregate job 生成 `n_pi+`、`n_K+`、`K+/pi+` 热图以及
      三联失败/mask 图；run `36854992530` 的 plots artifact manifest 确认
      `323` 行、`321 screened`、`2 gate_failed`，四张 PNG 均非空且带 SHA-256。
      绘图器为 solver-free，失败点不插值、不零填充。
- [x] 叠加明确标注单位、BQS、diagnostic 状态的参考线；参考线的 equal-flavor
      与 `mu_q -> mu_B` 转换在 plot manifest 中显式记录。
- [x] mask/failed 图已随上述 aggregate artifact 生成。

### M3：局部加密与完整验收

- [x] 由 solver-free 选择器在冻结线、相变参考线覆盖域、mask 边界和 ratio 梯度
      大的单元抽取候选点；未覆盖的相变线段显式计数，不伪造网格点。
- [x] 通过现有 benchmark workflow 的 `production` workload 在 GitHub Actions 运行
      12 个稀疏代表点（每点四个 charged channel）；run IDs 为
      `36863447044`、`36863456109`、`36863466376`、`36863474205`、
      `36863484426`、`36863494390`、`36863503078`、`36863513408`、
      `36863524728`、`36863533978`、`36863547938`、`36863558667`。
      这些点覆盖冻结线、梯度极值、mask 两侧和 crossover 参考线；候选源为
      screening run `36862279055` 的 manifest SHA-256
      `1562a6fc32d8aa48531749db874a36e3316eb5ab0fa391362481de6687f68ffe`。
      12/12 points、48/48 channels 均 `accepted`，`final_order=16`，没有把失败
      点转成零或成功值；production 总耗时范围为 `667.4--851.7 s/point`。
- [x] 对同一 12 点比较 screening 与 full-gate 的 `K+/pi+`：最大绝对相对差约
      `2.91%`，平均绝对相对差约 `0.885%`。因此 screening 可用于候选筛选，不能
      取代 full gate；该对照未改变当前局部的 ratio 单调性/高梯度区域判断，也不
      构成全网格收敛或 production baseline 更新授权。

### M4：证据治理

- [x] 绘图器和候选选择器生成 manifest、输入/源码 hash、行数、NaN/Inf/失败点摘要
      和图形/候选清单。
- [ ] 只有作者审核后才提升为 accepted research artifact；不自动更新 baseline。

### M5：授权扩域与加密 screening（2026-10-02）

用户明确授权在 GitHub Actions 运行 `T=40:5:220 MeV`、
`muB=0:25:1200 MeV`，共 `37×49=1813` 点。四通道、quark-only BQS 背景和
单点 screening 设置保持不变；这只是外部扫描网格加密，不提高每点内部积分
精度，不运行额外 benchmark，也不全量运行 production gates。

用户指定沿用无 BQS 约束的 equal-flavor 相变/CEP 参考，不另行定位 BQS CEP。
图和 manifest 仍区分“密度计算的 BQS 背景”与“相变导向线的 equal-flavor 背景”。
当前范围覆盖历史 CEP 附近，但 T 下限仍为 40 MeV，不声称覆盖到 T=5 MeV 的
完整一阶参考线。

本轮绘图新增三个带数值标签的等高线版本，同时保留三个原值热图和失败/mask
图。热图无插值；等高线仅在四角均为 `screened` 的单元内进行 display-only
取线，`corner_mask=false`，不跨失败区、不平滑、不外推、不裁负值。正密度若跨
两数量级以上，选取七个对数间隔的等高值（原值和热图颜色尺度不变）；ratio
以及存在非正密度的情况使用线性等高值。具体等高值、有效/阻断单元数和
生成器 hash 写入 manifest。图例移到数据 axes 外，PNG 输出 600 dpi。

该旧分支尚无当前 `figure_production_v2` 公共层；本轮图仍是独立的 diagnostic
PNG 审查产物，不声称通过论文级字体/布局测量、公共 validator 或矢量交付验收。
没有为作图迁移其他主线，也不修改生产数值配置。

- [x] 新增等高线、mask 保持和原值保留的合成 profile/PNG 契约测试。
- [x] 在 topic branch 推送并 dispatch 四个 T-row shards；run `36957893472`，
      source SHA `eb77d0616be4da6597df01e6159ad69b6eb19b89`；CLI 网格参数为
      `--t-grid 40:220:5 --muB-grid 0:1200:25`（start:stop:step）。
      `workflow_dispatch` 实际成功，无需添加 push 触发或合并 main。虽然当前
      默认分支的该文件 GET 返回 404，但历史已注册的 scan workflow 本轮可派发；
      不将此前新 full-gate workflow 的 404 限制泛化为所有 topic workflows。
- [x] 核对远程网格唯一性、点/失败数和首版 PNG/CSV 输出 hashes：
      run `36957893472` 全部 jobs success，`1813` 唯一网格键，`1805 screened`、
      `8 gate_failed`；每个分片 point count 为 `490/441/441/441`。原始数据只下载
      到 `D:/Temp/charged-gbu-contour-fine-36957893472/`，不提交数值 bulk。
- [ ] 下载并视觉审核热图、等高线和 mask；不自动晋升正式研究/论文产物。

本轮新增 mask 点为 `(T,muB)=(85,800),(95,550),(95,575),(105,675),`
`(145,375),(145,400),(150,375),(150,400) MeV`。背景和 pi 通道成功，
K± 均报 `static instability or unresolved Bose endpoint`；这不是凝聚的唯一诊断。
screening 的 BQS continuation seed 会随 muB 步长改变，不能假定加密前后的 mask
位置恒定；本轮没有为解释 mask 变化启动额外背景/端点复测。

另有 `15` 个有限输出的 K+ 负值，出现在右上角 T=200--220 MeV、
muB=1100--1200 MeV 的部分网格：最小 `n_K+=-8.6194608535e-4 fm^-3`，
最小 `K+/pi+=-0.0281890510`。这些仍是未通过全门禁的 signed diagnostic，
不裁零、不强行归因；绘图以叉号标记负值，而不把 `screened` 解释为 positivity
或 production acceptance。pi± 与 K- 的本轮 screened 总密度没有负值。

首版视觉检查发现 mask 标题重叠及密集等高标签。修正只涉及显示：最长有效路径
上每级至多一个标签，按弧长错开位置，短路径不贴标签；mask 图留出标题空间。
新增 `replot_run_id` workflow input，从旧 run 的冻结 shards 重绘，scan job 跳过。
manifest 分开数值 source SHA/run 与 postprocess SHA/run，不能把重绘 commit
写成数值 source commit。该入口用于本轮排版修复，不启动新的数值采样。

Python 绘图/候选契约 `15/15` 通过，包括负值保留、短路径标签、replot 跳过扫描
和 source/postprocess 分离；未改物理核、生产配置或容差。全仓数值 regression
未运行，因为本轮没有这些数值实现改动。

首个冻结重绘 run `36960286457` 已成功，scan job 实际 `skipped`，CSV hash
与数值 run 完全相同。进一步标签复核用最终画布上的 text bbox 排除相交标签；
只隐藏重叠的文字，不隐藏等高线、不删数值。标签不使用 inline 路径裁切，避免
逐次贴标签改变路径后把下一个标签吸附到邻近等高级；这是显示修正，不是物理
phase branch 或积分方式的修正。

## 7. 风险与回退

- BQS solver 在一阶线附近可能有多支；保存 seed、残差和失败点，不拼接共存相。
- Bose 支撑不安全时 mask 并保留 gate 信息，不裁零。
- 低分辨率伪结构必须由 M3 完整 gates 点复核。
- GitHub Actions 已给出 screening 的远程墙钟与 M1 网格完整性证据；full
  production 仍按单点/单通道 shard 运行，不能由 screening 成本外推。
