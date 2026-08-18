# C2 三维相图 v4 视觉闭合审计

## 范围与输入

本包只做 C2 结果的 solver-free 后处理，固定输入为 run `31862752226`、calculation SHA `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48` 和 postprocess SHA `fd359e792a89beb5ab12349bba761dc58ee16761`。没有调用 equilibrium solver，没有写入 reference，也没有改变 Maxwell、CEP 或三态判据。

坐标约定为 `x=mu_q`、`y=xi`、`z=T`。三维图只消费 `boundary_*.csv`、`crossover_*.csv`、`cep_*.csv` 和 `phase_grid_convergence_*.csv`；原始输入 hash 见 `figures/plot_manifest.json`。

## 物理筛选合同

偏导响应峰是数值候选，不自动拥有 crossover 物理含义。同一 `(xi, mu_q)` 不得同时标记为 crossover 和 Maxwell；在 `mu_q > mu_CEP` 的一阶区，响应峰保留在 `tables/crossover_physical_filter.csv`，但不绘制为 crossover。v4 因此绘制 1157 个保留响应峰，排除 331 个高于 CEP 代理的响应峰。

CEP 当前仍是温度 bracket 投影；中点仅为显示 landmark，不是 confirmed CEP。原生采样跨越 CEP 筛选边界时不跨 gap 补成物理面。

## 大 xi 处橙色面缺口

`tables/crossover_cep_sampling_gap.csv` 对每个 xi 比较最后一个 `mu_q <= mu_CEP` 的响应峰和第一个 `mu_q > mu_CEP` 的原生响应峰。大 xi 的直接证据如下：

| xi | CEP 化学势代理 (MeV) | 最后保留 mu_q (MeV) | 首个排除 mu_q (MeV) | 原生间隔 (MeV) |
| ---: | ---: | ---: | ---: | ---: |
| 0.3500 | 332.8530 | 324.3781 | 349.3303 | 24.9522 |
| 0.4000 | 338.1579 | 325.7758 | 350.8354 | 25.0597 |
| 0.4500 | 343.3177 | 327.1276 | 352.2912 | 25.1637 |
| 0.5000 | 348.3991 | 328.4366 | 353.7009 | 25.2644 |

因此，靠近 CEP 的下一个响应峰已经落在一阶区并被物理规则排除，而不是被错误画成 crossover。缺口的第一层归因是离散 `mu_q` 采样间隔约 25 MeV；这项归因由逐 xi 的原生点表直接支持。它不证明 CEP 已闭合，也不证明 solver 在缺口内没有响应峰。

此外，C2 boundary 在部分切片的相邻温度间隔为 10--16 MeV，超过 v3 的 8 MeV 三角化上限。这会造成 Maxwell 面的显示空洞，但属于 support/网格显示问题，不等同于物理面消失。

## v4 显示策略与边界

v4 是 `visualization-only`：所有 finite/converged Maxwell boundary 行（包括自动 geometry/interpolation 未闭合行）统一使用蓝色；不显示 unresolved 状态颜色，也不显示高于 CEP 的灰色叉号。为改善视觉连续性，v4 只把 boundary 三角化的 display-only `max_dT` 明确设为 16 MeV；不生成新数值行，不修改 CSV，不改变任何 gate。

`tables/maxwell_surface_point_status.csv` 仍保留 5424 条 unresolved grid 状态，manifest 标明 `closure_visualization_only=true`。因此 v4 不能用于 phase-reference promotion、strict 论文图或数值收敛结论。正式/strict 图仍必须断开 support gap，并保留 unresolved 语义。

## Maxwell 面空缺的分解

这里要区分三种网格：

1. `rho` 网格（Stage B `0.00625`、Stage C targeted `0.003125`）决定单个外部点的 Maxwell 几何精度；
2. 外部 `T/xi` 自适应网格决定 boundary surface 的 support；
3. crossover 的 16 个 `mu_q` 响应峰采样是另一套网格，不能用来推断 Maxwell 的外部 support。

C2 boundary CSV 共有 6886 行，且这些行都是最终 `converged=true` 的 Maxwell 结果。相比之下，`phase_grid_convergence` 有 22791 行，其中 5424 行 unresolved。按唯一 `(xi,T)` 坐标比较，约 3572 个 grid 坐标没有对应 boundary 行；其中约 2928 个坐标的主要原因是 `stable_no_s_shape`，表示该处本来就是 monotone 区，不应绘制 Maxwell。其余缺口主要来自 `rho_geometry_not_converged`、`geometry_tolerance_exceeded` 和 `hybrid_stage_c_not_converged`，这些是证书未闭合，不能直接画成 Maxwell。

一个直接例子是 `xi=0.5`：`T=17` 的 rho 细化记录在 level 3/4 报告 `rho_geometry_not_converged`/`hybrid_stage_c_not_converged`，所以 boundary CSV 没有该点；而 `T=107.0625` 及更高附近的记录出现 `stable_no_s_shape`，这是进入 crossover/monotone 侧后的物理缺失，不是丢失的 Maxwell 点。也就是说，C2 的细 `rho` 采样并不保证外部 Maxwell 面在所有 `(xi,T)` 上都有可绘制的 certified row。

v4 只把已有 finite/converged boundary 行统一为蓝色，并用显式 16 MeV display-only 三角化上限减少绘图空洞；它不会为 monotone 或 unresolved 坐标生成 Maxwell 行。

## 结论

当前图中不存在“同一 `(xi, mu_q)` 同时为 crossover 和 Maxwell”的物理双标记；橙色端点附近的空缺主要是 CEP 物理筛选叠加原生 `mu_q` 采样稀疏造成的可审计 gap。v4 适合作者检查整体拓扑和面的大致相对位置，不能替代逐点证书或 CEP 重验。
