# PNJL S 形临界判据可行性审计 v1

## 状态与范围

- overall verdict：`diagnostic_candidate`
- rho-support cascade verdict：`synthetic_candidate`
- full-window multi-evidence verdict：`performance_rejected`
- production authority：`false`
- reference promotion：`not_applicable`
- 目标：验证密度支撑区间能否把当前离散 secant 检测器漏掉的窄 S 形路由到少量固定 `rho` 补点，并建立准确度、CEP 外推与性能审计合同。
- 非目标：不修改 phase production 默认路径，不重新计算全温区相图，不晋升 reference，不把 BifurcationKit 加入根依赖。

本分析只使用解析 cusp 曲线族

\[
\mu(x)=x^3+a x,
\]

其中 `a < 0` 有 S 形，`a = 0` 为临界 cubic，`a > 0` 单调。解析真值只用于 feasibility 和测试，不作为 PNJL 物理结论。

## 判定合同

原型输出四种状态：

- `resolved_s_shape`：全区间/中心窗口、3/5 阶和奇偶留出拟合都给出两侧正斜率、中间负斜率的
  `+ -> - -> +` 拓扑，且独立定向验证至少连续 3 个 secant 确认同一拓扑。
- `near_critical`：稳定拟合提示隐藏负斜率或最小斜率接近零，但缺少独立认证。
- `resolved_monotone`：稳定正导数拟合与覆盖完整区间的独立、更密验证同时成立。
- `unresolved`：离散 topology、拟合窗口/阶次、留出误差或独立验证之间不一致。

`g=min(dmu/drho)` 只是一项预警证据。正的离散 secant 或正的拟合最小斜率均不能单独生成 `resolved_monotone`。

完整的全区间/中心窗口、3/5 阶拟合路径保留为保守对照，不进入建议的 hot path。新级联合同为：

1. 现有 coarse 点已经给出持续 `+ -> - -> +` 时直接返回，新增函数评估和 polynomial fit 都为 0；
2. 否则只在当前低斜率密度支撑与上一温度 spinodal 密度区间的交集内路由补点；
3. 合并窗口内已有 coarse 点与最多 12 个新增固定 `rho` 点，只做一次 local cubic；
4. `resolved_s_shape` 同时要求采样 secant 的持续 `+ -> - -> +` 与 cubic 拓扑一致；
5. `resolved_monotone` 只允许在上一温度 prior 已限定窗口、实际 secant 与 cubic 均保有正 margin 时产生；
6. 其他情况保持 `near_critical` 或 `unresolved`，不得靠拟合符号单独升级。

这里利用 cusp 正规形的不同消失尺度：接近 CEP 时，spinodal 密度宽度满足
`Delta rho ~ (T_CEP-T)^(1/2)`，而化学势高度满足更快消失的
`Delta mu ~ (T_CEP-T)^(3/2)`。因此密度区间可用于路由，不把微弱的化学势高度当作唯一信号。

## 解析用例结果

完整数值见 [case_results.csv](tables/case_results.csv)。当前 7/7 用例符合预先声明的状态：

| 用例 | 当前 secant baseline | 多证据结果 | 说明 |
|---|---|---|---|
| 明显 S 形 | 检出 | `resolved_s_shape` | 离散 fold、跨窗口拟合与独立验证一致 |
| 欠采样窄 S 形 | 漏检 | `near_critical` | 只预警，不冒充真值 |
| 窄 S + 定向验证 | 漏检 | `resolved_s_shape` | 独立密集验证确认持续的 `+ -> - -> +` 拓扑 |
| 临界 cubic | 未检出 S | `near_critical` | 最小斜率在临界容差内 |
| 未验证单调曲线 | 未检出 S | `unresolved` | 正 secant 不足以认证单调 |
| 已验证单调曲线 | 未检出 S | `resolved_monotone` | 独立覆盖与正 margin 同时成立 |
| 带噪单调曲线 | 误检 S | `unresolved` | 留出误差和拟合稳定性阻止误升级 |

该结果支持“用跨窗口拟合做临界预警、用独立定向求解做认证”的算法结构，但完整多证据实现的分类成本过高，已标记为 `performance_rejected`。
全局单调下降曲线和只有一个负 secant 的验证曲线另有确定性回归，均不能进入 `resolved_s_shape`。

## rho 支撑级联与 CEP 外推

[cascade_case_results.csv](tables/cascade_case_results.csv) 中 7/7 用例符合预声明状态。明显 S 形走 fast path，新增 0 点；两个 coarse 漏检的窄 S 形各新增 8 点、1 次 cubic 后得到采样与拟合一致的持续拓扑；临界 cubic 保持 `near_critical`；未认证单调、带噪单调和全局下降反例均未误升级。

[temperature_sequence_results.csv](tables/temperature_sequence_results.csv) 使用 13 个 coarse `rho` 点跟踪同一解析 cusp 的 10 个递增温度。明显 S 形阶段直接返回；进入欠采样区后复用上一温度 spinodal 密度中心/宽度，每个温度新增 8 点。7 个低温 resolved 宽度单调收缩，临界点为 `near_critical`，临界点上方只有在 prior 窗口内采样 secant 与 cubic 都为正时才得到 `resolved_monotone`。

对 7 个 resolved 宽度拟合 `Delta rho^2 = m T + c`，
[cep_extrapolation_summary.csv](tables/cep_extrapolation_summary.csv) 给出
`T_CEP=1.0268768`，合成真值为 `1.0`，绝对误差 `0.0268768`，通过预声明的 `0.05` feasibility gate。这只证明密度宽度外推在受控正规形上可用，不是 PNJL CEP 结果。

## 性能边界

[accuracy_performance_frontier.csv](tables/accuracy_performance_frontier.csv) 只比较同进程解析曲线分类的工作计数与耗时。表中显式记录解析函数评估、polynomial fit 和 fit holdout 评估；`pnjl_residual/jacobian/newton/anchor/branch` 保持缺失，因为本 PR 没有运行 PNJL solver。多证据方法的价值在于减少错误 topology 决策，而不是声称曲线分类本身比 secant 更快。

本次 500 次重复记录中，当前 secant baseline 约 `4.10 us/case`，rho-support fast path 约 `33.54 us/case`，完整 rho-support cascade 约 `59.48 us/case`，旧多证据路径约 `742.74 us/case`。因此级联在本次解析分类中仍约为 baseline 的 `14.5x`，不能据此宣称 hot-path 性能已通过；但相对旧多证据路径约快 `12.5x`。

更关键的 solver 工作量代理是新增采样数：500 次重复的级联只使用 `12000` 个 targeted 点和 `1500` 次 cubic，而旧路径使用 `1002500` 个独立 validation 点、`25500` 次 polynomial fit 和 `183000` 个 holdout 评估。折算每个默认解析用例，级联平均新增 `24/7` 个点，单一 ambiguous 用例实际为 8 点且硬上限为 12；fast path 为 0。真实 PNJL 中每个 fixed-rho 点包含昂贵 equilibrium solve，故必须由后续同物理 pilot 测量，而不能把这里的解析 wall-time 或函数计数直接外推为 runner-minutes。

[process_startup_timings.csv](tables/process_startup_timings.csv) 分离记录 analysis module include、首次 suite 调用与同进程 warm 调用。首次 suite 可能包含完整 `Models` 加载；该成本不是 PNJL equilibrium 或 PALC continuation 成本。

启动探针单独记录 module include、首次 suite 调用、warm full-window 路径与 warm rho-support cascade。首次调用包含 `Models` 冷加载、编译及 suite 自检，数值受本机进程状态影响，不作为固定性能基线。因此后续 Actions pilot 必须让同一 shard 复用一个 Julia 进程，并分别报告启动与数值工作量。

PALC 历史记录 `125.4 s cold / 1.37 s warm` 被保留为 `historical_noncomparable`：两个数字来自不同场景，前者包含 BifurcationKit 编译，均不能与本解析 benchmark 或完整 phase production 直接比较。现有正式化结论仍是 `ready_for_production_replacement=false`。

真实 PNJL 性能 gate 必须记录：

- equilibrium/residual 调用数；
- Jacobian 构造与 Newton 迭代数；
- anchor solve、分支数和 continuation 点数；
- 单 shard 临界路径与总 runner-minutes；
- 冷启动和同进程 warm 成本。

缺少这些计数时只能保留 diagnostic verdict。

## 后续 pilot

下一阶段应在 GitHub Actions 上运行窄窗口、同物理口径 pilot，而不是启动 C3/O1 或 formal production：

1. 固定旧 canonical CEP 邻域的代表性 `xi=-0.5,0,0.5`。
2. 每个 xi 在同一进程共享 thermo/solver 配置与初始曲线。
3. 比较当前 dense secant/Maxwell、rho-support cascade，以及仅作 oracle 的 full-window 多证据或 isolated fold 求解；PALC 只在性能证据支持时保留。
4. 对所有方法执行相同 ground-state pressure governance。
5. 只有准确度不低于高分辨率 oracle 且同精度总成本不高于当前路径时，才可成为 production candidate；更准确但明显更慢则保持 diagnostic oracle。

## 复现

生成分析包：

```powershell
julia --project=. scripts/analysis/pnjl_criticality_feasibility.jl `
  --output-dir=docs/analysis/pnjl/criticality_feasibility_v1 `
  --repetitions=500
```

随后运行启动/热路径探针；该命令写入同一 evidence package 并刷新 manifest：

```powershell
julia --project=. scripts/perf/pnjl_criticality_feasibility_probe.jl `
  --output=docs/analysis/pnjl/criticality_feasibility_v1/tables/process_startup_timings.csv `
  --repetitions=100
```

第一条命令会移除旧的启动表，避免把不同代码版本的 timing 混入新 case/performance 表；第二条命令重新生成启动表并把两个 producer script、README 和全部表格的 SHA256 写回 manifest。

聚焦单测：

```powershell
$env:UNIT_FILES='models/test_pnjl_criticality_feasibility.jl'
julia --project=. tests/unit/runtests.jl
```

## 证据边界

claim 级证据见 [claim_ledger.csv](tables/claim_ledger.csv)。当前没有可直接用于论文的 PNJL CEP 或性能结论；解析曲线只证明算法合同在受控反例上可行。

manifest 的 `repo_head` 是生成时的基线 HEAD；本分析包与 producer 同属一个尚未提交的 PR，因此同时记录
`generation_worktree_dirty=true`，并以 `producer_scripts[*].sha256` 绑定实际执行代码。不能只用基线 HEAD 替代 producer hash。
