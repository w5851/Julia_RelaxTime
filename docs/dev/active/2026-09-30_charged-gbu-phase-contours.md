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
- [x] 下载并完成 agent PNG 视觉审核、最终输出 hashes 与冻结 CSV 比对；
      等待作者审核，不自动晋升正式研究/论文产物。

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

最终 PNG-review 由冻结重绘 run `36960810064` 生成（全部执行 jobs success，
scan 为 `skipped`），postprocess SHA
`104d5c311e44a99883926c404a574498807a76a2`；数值 source 仍为
run `36957893472` / `eb77d0616be4da6597df01e6159ad69b6eb19b89`。
输出为三张原值热图、三张等高线图和一张 mask 图，7/7 PNG SHA-256 复核通过。
最终目录为 `D:/Temp/charged-gbu-contour-fine-36957893472/final-plots/`。

- 冻结数值 CSV SHA-256：
  `3a77588502ccc9728f08ab49cf4bb59f1840d709cad229502c561b3210cd50ea`，
  三次绘图输入/派生 CSV 字节一致，没有重跑或修改数值。
- 最终 plot manifest SHA-256：
  `571ae0b35bf8ec4f1a6971a1235d0d9b2ca1d4f253888536ac88884241c87622`。
- 原始 Actions 从创建到完成约 20 分钟；四个 scan job 总墙钟为
  `15m30s/13m42s/18m55s/13m09s`，包含环境准备，不是单点 density benchmark。
- 候选选择器生成 44 个 solver-free 候选，只生成清单，未派发 full gates。
- 本轮的 source/config/point/figure 都保留 provenance，数据没有进入 Git；
  Actions artifacts 保留 30 天，不等于永久研究档案，作者接受后再治理保留。

### M6：冻结线附近的二维 ratio 调整地图（2026-10-02）

用户授权仅使用既有数据快速作图，不重新扫描、不拟合实验、不更改默认冻结线。
新增 `scripts/analysis/relaxtime/plot_charged_gbu_freezeout_neighborhood.py`，从 M5
冻结 CSV 和 manifest 验证 hash、唯一键、行数、status、ratio 与两通道密度的一致性。
绘图选取 `0<=muB<=750 MeV`、`|T-Tfo(muB)|<=20 MeV` 的探针带：
`242 screened`、`6 failed`，有效 ratio 范围 `0.279063--0.744012`。
这只是已采样邻域，不是冻结线置信带或全精度验收。

局部颜色使用线性尺度 `0.25--0.75`，等高值步长 `0.05`；原始值不改变、不裁切。
热图无插值；等高线仍只在四个有效角点的单元内 display-only 取线，
`corner_mask=false`，不跨失败或未选区域。图中默认线的数字为碰撞能量 GeV；
`a±10 MeV` 示意线只改变 `a_GeV`，保持 `b,c,d,e` 和能量到 muB 的映射不变。
示意线是坐标探针，不是新产额或实验拟合。

此旧分支缺少 v2 公共层，生成命令通过 `--plotting-support-root` 只读复用
`D:/Desktop/Julia_RelaxTime` 的当前 v2 层，不迁移或修改该工作树。该支持工作树
存在其他任务改动，因此 manifest 同时保存其 HEAD、dirty 标记和所消费模块/profile
的精确 hashes；不能仅用 HEAD 代表其内容。输入和支持文件生成前后 hash 相同。

最终新 sibling：
`D:/Temp/charged-gbu-contour-fine-36957893472/freezeout-guidance-v3/`。
`freezeout_neighborhood_ratio.png` 已经 agent 视觉审核；公共 v2 validator 通过，
600 dpi、最小测量字高约 `2.096 mm`、无裁切/文字重叠、图例位于 axes 外。
保持 `audit/png_review`、`manuscript_eligible=false`，等待作者审核。
新合成/输入契约测试 `6/6`；本轮绘图与候选测试合计 `21/21`；所复用的公共
plotting contract/quality 测试 `31/31`。不运行数值 regression，因为未修改计算核。

复现入口：上述脚本以 `--csv <M5 final-plots/contour_points_merged.csv>`、
`--source-manifest <M5 final-plots/plot_manifest.json>`、
`--plotting-support-root <含当前v2层的只读工作树>`、`--output-dir <新目录>` 运行。
已有输出目录被拒绝覆盖；本轮未 push、dispatch Actions、创建 PR 或修改 production。

#### M6 组会显示修订：±50 MeV（2026-10-02）

用户明确要求扩大到冻结线 ±50 MeV、移除参数探针和能量数字，并强调本图用于
组会，不执行论文级正式化。新增同一脚本的 `--presentation --band-MeV 50`
分支；原 v2 审查模式保留。组会分支不加载论文级公共层、不申请 APS/矢量资格，
但继续验证输入 hash、状态、ratio、mask 和不覆盖历史目录。

显示 `muB=0--750 MeV`、已有 `T=40--220 MeV` 网格中的选定邻域：
`608 screened`、`7 failed`，ratio 范围 `0.0475296--1.12918`。带的低温端
低于 40 MeV 时截在原输入域，不补算或外推。默认冻结线改用红虚线以区别黑色
等高线；没有 a±10MeV 曲线、没有能量标签。全部等高值使用 0.05 步长，色标
也按 0.05 标注。标签在一次 inline-clabel 调用中直接嵌在对应路径上，生成时
逐一核对其文本与选取的原等级一致；长曲线补重复标签，不跨 mask 重连。

组会最终 PNG/manifest 位于
`D:/Temp/charged-gbu-contour-fine-36957893472/freezeout-guidance-band50-meeting-v2/`。
画布 `12.2×8.8 in`、220 dpi，仅作组会显示审核；已 agent 视觉检查，不声称
v2 最终尺寸/投稿合格。相关新契约 `9/9` 与既有绘图/候选契约 `15/15` 通过，
原始数值 CSV hash 不变。未改变默认配置、计算核、全精度 gate 或数值 baseline。

组会标签简化版保留在同根目录 `freezeout-guidance-band50-meeting-v3/`：
去掉标签白色底框，只保留 0.7 pt 细描边；标签和色标数字每 0.10 一次，
等高线仍为 0.05。与 meeting-v2 的输入 CSV hash、全部 contour levels、
color limits、608 screened 和 7 failed 均相同。11/11 相关契约测试通过，
agent 已审核 PNG；这是显示修订，不是新数值计算。

### M7：同网格 q=0 外推远程重算（2026-10-02）

用户授权在 GitHub Actions **新计算** M5 的同一完整网格，生成与 M6
`meeting-v3` 同构的冻结线邻域组会图；不是从 finite-q 数据改标签。

- 网格 `T=40:5:220 MeV`、`muB=0:25:1200 MeV`，1813 点、4 个 T-row shards；
  CLI 为 start:stop:step，即 `40:220:5`、`0:1200:25`。
- 同一 quark-only BQS continuation、四通道、无限 PNJL 热目标、双线真空截断；
  同一 `mesh/cut/tail/omega/q=64/32/32/64/8`、`qmax=8 fm^-1`。
- 新的 `charged_gbu_q0_reference.jl` 仅复用当前无限热 profile，不调用历史
  有限热 `bubble_at`。沿用此前对照的 `q0_lambda_reference`：
  `lambda0=sqrt((omega+shift)^2-q^2)`，仅 `omega+shift>=q`；其他区域相位零
  是该外推近似的定义，不是删除直接 finite-q 谱。
- 根在 q=0 解析 gap 中独立计数，再映射为 `hypot(lambda_root,q)-shift`。
  Bose 用外部 omega，GBU 分部积分保留 `-n*g(omega_threshold)`；单位 fm^-3。
  保留同一 shell 支撑、阈值/Mott 不确定性和 UV 门槛，不 fold/翻符号/裁零。
- 显示同一 `muB=0--750`、冻结线 ±50 MeV、已采样 `T=40--220`；同一
  12.2×8.8 in、220 dpi、0.05 等高线、0.10 无底框标签、红虚线，无参数探针
  或能量数字。固定色标 `[0,1.15]`；若超界，显示 overflow 色及色标端箭头，
  不改数值。mask 由新计算的失败结果决定，不能复制旧 mask 强求一致。

- [x] 实现显式 route selector、source hash 和 resume 隔离；默认 finite-q 不变。
- [x] workflow 接入新数值 route 和 solver-free 同构组会图。
- [x] 完成合成根/坐标/GBU 边界、route 混用拒绝和绘图不变约束验证。
- [x] 推送最小源码变更并新 dispatch；核对数值 SHA、参数和实际 scan jobs。
- [x] 审核1813个唯一网格点、四分片、背景一致性、失败/负值和输入输出 hashes。
- [x] 下载并视觉审核组会 PNG；报告与 finite-q 的诊断差异，不晋升 production。

新增 Julia 测试仅用合成 profile/纯代数；真实数值运行全部留在 Actions。
生产配置、PNJLCore、上游求解器与 baseline 不修改；不在本轮合并分支。

本轮 focused Julia `68/68`、Python `30/30` 通过。首轮 Julia 的唯一失败是测试
将 Float64 相位与精确 Irrational π 严格比较；改为同类型断言后通过，未改物理
核或容差。docs consistency、formula-route、script entrypoint/governance、
data-output guard、task-ledger preflight 与 `git diff --check` 均通过。
全仓 active-doc 检查报告五份既有 2026-07 文档超龄（>60d）；这些文档未修改，
不在本轮擅自归档。完整数值 regression/full gates 未运行：本轮新增的是独立
screening 近似，生产实现/config/baseline 与 finite-q 分支的求值调用不变。

#### M7 远程结果与边界

GitHub Actions run `36971423991` 实际执行四个新数值分片和 aggregate，全部
jobs success；数值 source SHA 为 `e6b3de025d59213e480927ebffa5f188e3e71ebd`，
未使用 `replot_run_id`。1813 个唯一 key、分片 `490/441/441/441` 和同一
背景/config/infinite kernel 源码 hashes 已核对。新的 route、坐标约定和 source
hash 进入 scan/plot manifest，不把旧 finite-q CSV 伪装成外推结果。

- 全域 `733 screened / 1080 gate_failed`，无 background_failed。
- K+ 的 1073 点报 `static instability or unresolved Bose endpoint`；另 7 点
  K± 都报 `q0 reference onset is not positive real`。不强行将这些错误唯一归因
  于凝聚、额外根或未收敛，也不放宽门槛。
- 同构邻域选定615点：`301 screened / 314 failed`，可显示 ratio 范围
  `0.0221157040--0.4313014934`。固定 `[0,1.15]` 色标无超界；0.05 等高线
  和 0.10 标签保持。agent PNG 视觉审核通过，但大片失败域使本图不能作为
  已闭合的完整冻结线调参地图。
- 本轮 screened 全域四通道均无负密度；全域 ratio 仍可达 `2.33762697`。
  “screened”不等于完整门禁或模型可靠性，不能只展示正的成功值宣称近似通过。

背景解不是逐字完全相同：4点的解向量最大绝对差大于 `1e-8`，为
`(40,1100),(80,775),(85,800),(105,675) MeV`，全部在新扫描中失败；后两点
旧图也失败。其余共同 screened 点的解向量最大绝对差 `4.32166525e-11`，
残差最大差 `1.76973495e-13`。不能把那4点的差异算作 q0/finite-q 近似误差，
也不能因源码相同就假定 solver 在所有多解处选择同支。

新增一组纯合成 Landau profile 测试（6项），显示原 q=0 在
`lambda=shift` 的外部静态占据零点，在 `0<q<shift` 的内部 invariant 映射后
一般会落到 `sqrt(shift^2-q^2)`，不再自动是零点。它证明该近似不保证端点
条件，不证明每一个真实失败点都是这一机制。单文件含新测试 `40/40` 通过，
没有修改 gate、物理核或真实扫描数值。

本轮沿用历史 q0 的**坐标近似**，但不是逐字沿用历史有限热数值密度算法：
`causal_gbu_dense_execution.jl` 用有限热核和非零 `s.lower` 的导数积分；当前
适配器用无限热核和当前零端点条件下的 GBU 分部积分。历史参考算法只检查
lambda0=0 的 onset，不等价于检查外部 omega=0。故历史十点成功记录不能直接
证明本轮全部 BQS 背景的外推闭合，不为复现旧值引入有限 Bose 截断。

输出只保留于 `D:/Temp/charged-gbu-contour-q0-36971423991/`：

- `inputs/`：四个原始分片（points JSON、CSV、manifest）。
- `plots/`：全域合并 CSV、热图、等高图和 mask 图。
- `meeting/freezeout_neighborhood_ratio.png`：本轮同构组会图，未覆盖旧图。
- `audit_report.json`：solver-free 网格/source/背景/失败/图像 hash 对照。

合并 CSV SHA-256：
`6b98c1b4c0d26a10c73ae6347e00448a08882ebe6bfa559f4895e96c1c69d401`；
组会 PNG SHA-256：
`69b33eca5ef0168afafc17d158311147a3c5a223f2db699fa24c701134e1bf2c`。
数值 bulk 不提交，未开/合并 PR、未修改 main 或 production 默认。
后续若需要填满外推图，必须另行评审其低频/化学势坐标与详细平衡，而不是
移除检查、静默删除 Landau 段或复制早期带积分下限的结果。

最终展示图在 `meeting-review-v2/`：只将原组会图脚注的“no rescan”改成
“no fit”，避免误解为没有新运行数值。该 solver-free 重绘的输入 CSV hash、
selection/mask/ratio/color/contour/label 字段与远程 meeting 图一致，未重扫；
最终 PNG SHA-256 为
`0d86c71cf2c446e4b59e39f0777b723548eb7ad18d4b8d16f0e142c24f39fce9`。
Python `30/30` 复测通过，agent 已完成最终 PNG 审核。Actions 从 dispatch 到完成
约17分11秒，包含环境准备；不是新的单点 benchmark 结论。

### M8：授权的小型外推端点复测（2026-10-02）

用户授权在 Actions 对约三个**已有背景**进行小型只读复测；不求新背景、
不重扫整图、不改容差。此项为本独立任务的 `research` follow-up，不替换
accepted 主线，也不授权修复、删谱、改变外推约定或填补 mask。

新增 `audit_charged_gbu_q0_endpoint.jl`，固定输入为历史十点对照中的 7.7 GeV、
3 GeV 背景，以及 M7 run `36971423991` 的 `(145,375) MeV` onset 失败背景。
前者读取已 hash-bound 的 `backgrounds.csv`；后者只通过保存的八维 seed
代数恢复质量、化学势和 KMT coupling，不调用上游求解。选定输入字节及
输入/源码哈希一同进入新 artifact；拒绝覆盖已有诊断目录。

- K+ 主探针、pi+ 和 K- 对照；检查原 q0 的 `lambda=shift`、`lambda=0`，
  以及各 q 外推的 `sqrt(shift^2-q^2)`。分别保存 Re/Im F、相位、GBU 权重和
  原门槛的两个条件；不把合并错误文本直接译成凝聚。
- 三档 `mesh/cut/tail=64/32/32,128/64/64,256/96/96`，并在关键点用原谱核
  PV 128/256 阶与解析 q0 on-shell cut 核对；区分插值误差和真实非零 cut。
- 同一坐标比较旧 `Lth=24,mesh=512` 核；在相同冻结背景上区分核更新和
  端点门禁/积分算法变化。原有 `1e-8` 静态虚部门槛不变。
- 对出现非零低频外推相位的壳层，记录 `omega_min=1e-3...1e-7 fm^-1`
  的有限窗口导数、bulk 与上下边界项；仅作端点探针，不运行完整 density。
- 保留直接有限 q 的 K+ 静态对照。此小样本不能外推全部1080个mask的原因，
  也不能仅凭 Re F 为负认证新的平衡相或物理凝聚。

scan workflow 增加默认关闭的 `endpoint_audit`。开启后 scan/aggregate 实际应
跳过，仅下载旧分片并执行独立小型 job；普通扫描和重绘默认语义保持不变。
本机仅运行纯代数/合成 profile 单元与工作流契约，真实核探针全部留在 Actions。

- [x] 聚焦单元/工作流契约和治理验证：Julia `110/110`、Python `30/30`；
      docs/formula-route/script/data-output/ledger 与 diff 检查通过。
      首次测试捕获 CSV 不接受 `nothing`，改为仅在输出层明确表示 `missing`，
      未改核或容差。全仓 active-doc 检查仍为五份既有七月文档超龄，未擅自归档。
- [x] 选择性提交诊断源码，dispatch 小型 audit；核对 scan/aggregate skipped。
- [x] 审核实虚部、原 cut、节点对照、旧核/新核与 IR 边界证据并报告结论。

#### M8 远程数值证据与诊断结论

Actions run `36976406317` / source SHA
`3f6c8e2cdc3df4182ae9b74d86bef7ee4ccca471` 成功；只有 prepare 和小型
endpoint-audit 执行，scan/aggregate 均 skipped。artifact 明确记录三背景、零新
背景求解、未计算完整密度、容差未改、无 production 晋升。含351行三档端点、
117行历史核对照、9行直接有限q静态对照、60行有限IR窗口；17个输出哈希逐一
核对通过。数值只保留于
`D:/Temp/charged-gbu-q0-endpoint-audit-36976406317/endpoint-output/`，未进入Git。
manifest SHA-256 为
`e3afd753ed9c3fab9d177984d9900c9b55fe49847df28da0381646767982de9a`。
远程两文件合成契约 `88/88` 通过；核探针步骤约34秒，不是完整密度benchmark。

下表的逆传播子为无量纲 `F=1-4K*Pi`；取当前8节点网格首壳层
`q=0.1588405740098553 fm^-1`。原q0静态位置是 `lambda=shift`，外推后外部
`omega=0` 对应 `lambda0=sqrt(shift^2-q^2)`。两个历史背景中均位于Landau cut。

| 保存背景 | 原q0静态Re F（原PV256） | 外推静态lambda0 (fm^-1) | 外推静态Im F（原cut） | mesh256插值Im误差 |
| --- | --- | --- | --- | --- |
| 7.7 GeV：T=139.6121、muB=421.6499 MeV | +0.15663527，Im=0 | 0.40093476 | +1.58803995e-6 | -2.09431e-10 |
| 3 GeV：T=79.95694、muB=719.0764 MeV | +0.14766760，Im=0 | 0.73204038 | +1.09509399e-8 | -2.42181e-14 |

三档虚部均超过原 `1e-8` 门槛。原q0 on-shell cut不含径向求积误差；原PV
128/256在这些位置的实部差约机器精度。用旧有限热Lth24核在**同一映射位置**
得到Im F为 `1.58798752e-6 / 1.09509300e-8`，同样失败；直接有限q原核静态
Im F分别约 `4.23e-20 / -6.75e-23`，实部正。故这两个样本的mask不是换成
无限热核造成，也不是粗插值误差或静态Re F变负；是保留的内部lambda外推
没有保持外部Bose原点的占据零点。该结论针对已复测样本，不把全1073个合并
错误逐点标成同一种原因。

历史24节点首壳层 `q=0.01925112` 的7.7 GeV原cut Im F也为 `3.12258e-8`，
并非新8节点碰巧选中异常点。3 GeV在此较小q处虚部低于数值门槛，但不严格为零；
不能将“某q通过容差”提升为所有q的严格端点证明。

**历史十点为何成功：**它检查原q0静态点与内部origin，不检查每个boosted
静态点，且使用 `omega_min=1e-5 fm^-1` 的导数积分。此次在相同历史背景、
相同旧核上已验证端点条件仍缺失。相位端点很小（7.7/3 GeV在首8节点壳层
分别约 `-1.0103e-5 / -7.3890e-8`），GBU权重三次抑制，不能据此宣称旧ratio
数值上已有巨大误差；但也不能把有限窗口成功当作零端点极限闭合。

数学上，令 `W=delta-sin(2delta)/2`，在有限光滑窗口有
`integral(-g')W/pi = integral g*W'/pi - [gW]_a^b/pi`。
若真实映射端点 `W(0)=W0!=0`，则bulk积分的主导项为 `T*W0/(pi*a)`；
不能省去下边界并沿用零端点的部分产额式。7.7 GeV首8节点壳层在
`omega_min=1e-5 -> 1e-7 fm^-1` 时，bulk为
`-1.96818e-14 -> -1.97886e-12 fm^-2`，下边界同时增长，导数仅约
`9.44e-17 -> 1.58e-16 fm^-2`。这是一个局部低频窗口，**不是完整密度误差界**。
IR表刻意复用历史相位权重差分算法；极小权重的导数平台还受浮点抵消影响，
不把它解释为严格收敛。端点非零的判断以原cut与稳定GBU权重为依据。

**另一类失败：**保存的 `(145,375) MeV` 背景，K±在内部origin的原PV256
Re F约 `-0.04038549`；K+原q0外部静态Re F约 `-0.04427347`，直接有限q首壳层
也约 `-0.00829499`。旧Lth24核和三档插值同样负，故不是仅外推引入或插值翻号。
原q0静态虚部为零不代表静态实部正。这支持当前正常charged支撑门禁应拒绝
该**保存背景**，不证明实验中的凝聚或新的稳定平衡相。该seed的奇异质量排序
`Mu≈3.3021 > Ms≈1.8979 fm^-1`、`phi_s>0`也已保留；残差小不等于平衡支
已认证，后续背景选择/稳定性需要另行授权评审，本轮没有重求或改选背景。

下一步为 `research` 评审：外推近似在有限化学势/Landau低频处是否有与当前
产额定义兼容的推导；另独立复核负静态背景的支选择与稳定性。不能通过取消
端点门禁、裁掉Landau、改变phase分支或复制旧窗口值来填图。直接有限q默认、
原核、PNJLCore、生产配置、baseline及已有PNG/CSV均保持不变；未开/合并PR。

### M9. q0 有限窗口、数值警告与同背景重计算（2026-10-04）

用户明确授权按前述建议修改、重计算和绘图，要求警告保留在数值中、不要出现在
新等高图中。此授权接续 M8 的方法评审，覆盖 q0 参考的有限窗口处方；本任务仍为
独立 research/screening 工作，不替换 accepted 直接有限 q 或 production baseline。

- 默认 `endpoint_policy=finite_window`、`omega_lower_inv_fm=1e-5`。将外推端点
  非零虚部降为诊断；仍拒绝原 q0 onset、负静态实部、非正根、异常拓扑及阈值/UV。
- 每个 cut 计算有限窗口的导数积分，采用减去常数 `W(a)` 后的分部积分恒等式，
  不显式差分相位，不丢弃上下边界；使用对数分段和稳定权重差。保留
  `strict_zero_limit` 供历史对照，不能将有限窗口结果标成零下限收敛结果。
- 同网格 `T=40:5:220, muB=0:25:1200 MeV`。使用 run `36971423991` 保存的
  1813 个背景 seed，代数恢复质量/耦合并验证背景源码，避免重新选支。
- Actions 先在 `(140,425),(80,725),(160,0),(145,375) MeV` 四背景比较
  下限 `1e-4/1e-5/1e-6`、omega 64/128、profile 64/128、q 8/16 和严格旧口径。
  新窗口/omega 的局部比较事先采用 `1e-3` 相对加 `1e-12 fm^-3` 绝对判据；
  谱/q 分辨率独立报告，不据此宣称全图 production 收敛。
- CSV/点 JSON 保存端点虚部、相位、权重、警告计数、边界项和处方。
  `--clean-contours` 隐去警告符号/失败图例，失败点仍为空值且禁止跨空缺等高线。
  新 PNG 保留组会图层定位及当前冻结线，历史 CSV/PNG/manifest 不覆盖。
- 本地验证层为公式/合成 profile、背景恢复及 CLI/绘图契约；真实核、固定背景密度
  比较和全网格重新积分均在 Actions。稳定 Models/API 与背景方程未变，无需改 API。

- [x] 聚焦 Julia 测试分两组 `212/212`、`139/139`（含重复公式测试），Python `32/32`；
      docs/formula-route/script-entrypoints/data-output/ledger 与 diff 检查通过。
      初次夹具使用了非接口要求的 Vector 输入，改回 meanfield-state phi 类型；
      工作流条件测试同步新增独立 audit 分支，未改任何数值容差。
- [ ] Actions 四背景敏感性审计，记录失败保留与误差。
- [ ] Actions 同背景 1813 点密度重算，核验哈希和新增可用点。
- [ ] 生成无警告标记的全域/冻结线邻域等高图，核验并交付 PNG。

## 7. 风险与回退

- BQS solver 在一阶线附近可能有多支；保存 seed、残差和失败点，不拼接共存相。
- Bose 支撑不安全时 mask 并保留 gate 信息，不裁零。
- 低分辨率伪结构必须由 M3 完整 gates 点复核。
- GitHub Actions 已给出 screening 的远程墙钟与 M1 网格完整性证据；full
  production 仍按单点/单通道 shard 运行，不能由 screening 成本外推。
