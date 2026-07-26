# PNJL S 形临界判据可行性审计 v1

## 状态与范围

- verdict：`diagnostic_candidate`
- production authority：`false`
- reference promotion：`not_applicable`
- 目标：验证多证据临界判据能否在不把拟合最小斜率当作真值的前提下，发现当前离散 secant 检测器漏掉的窄 S 形候选，并建立准确度与性能审计合同。
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

该结果支持“用跨窗口拟合做临界预警、用独立定向求解做认证”的算法结构，但不支持把解析曲线速度或准确率外推到 PNJL。
全局单调下降曲线和只有一个负 secant 的验证曲线另有确定性回归，均不能进入 `resolved_s_shape`。

## 性能边界

[accuracy_performance_frontier.csv](tables/accuracy_performance_frontier.csv) 只比较同进程解析曲线分类的工作计数与耗时。表中显式记录解析函数评估、polynomial fit 和 fit holdout 评估；`pnjl_residual/jacobian/newton/anchor/branch` 保持缺失，因为本 PR 没有运行 PNJL solver。多证据方法的价值在于减少错误 topology 决策，而不是声称曲线分类本身比 secant 更快。

本次 500 次重复记录中，当前 secant baseline 约 `8.75 us/case`，多证据路径约 `1.09 ms/case`，曲线判别层约慢 `125x`；后者同时计入解析验证曲线、全区间/中心窗口 3/5 阶拟合和留出检查。该比值只适用于内存中的解析曲线，不能外推为 PNJL phase solve 比值。

[process_startup_timings.csv](tables/process_startup_timings.csv) 分离记录 analysis module include、首次 suite 调用与同进程 warm 调用。首次 suite 可能包含完整 `Models` 加载；该成本不是 PNJL equilibrium 或 PALC continuation 成本。

本次探针中 module include 约 `6.46 s`，首次 suite 调用约 `161.70 s`，同进程 warm 路径约 `4.77 ms/case`。首次调用包含 `Models` 冷加载、编译及 suite 自检，数值受本机进程状态影响，不作为固定性能基线。因此后续 Actions pilot 必须让同一 shard 复用一个 Julia 进程，并分别报告启动与数值工作量。

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
3. 比较当前 dense secant/Maxwell、局部多证据判据、isolated multi-branch PALC 或直接 fold 求解。
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
