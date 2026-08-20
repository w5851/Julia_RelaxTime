# C2 targeted manual review v1

这是一个 solver-free 的作者审核包，不是 production/reference 结果。它固定读取
source run `31941614867` 的 9 个 regression target，calculation SHA 为
`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`，workflow/postprocess SHA 为 `aa9f3a70a821df863dd9f4fbbc8a89b053cdbf8c`。

## 9 个 rho-mu 曲线

`curves/` 保存 source artifact 中的 hybrid 和 independent full-fine oracle 曲线副本；
`figures/` 为每个点的完整范围加局部 S 形放大图，`nine_point_local_overview.png`
提供 3x3 总览。`tables/point_summary.csv` 保留分类、候选数、交点数、geometry 和
局部坐标范围。所有曲线点都经过 finite/converged 和重复 rho 检查。

当前观测是：9 个点的 oracle 均为唯一三交点一阶候选，而 hybrid 均为
`ambiguous_near_critical`，原因仍是 Stage-B 到 Stage-C 的 geometry certificate 未闭合；
本包不把 oracle 标签写回 hybrid。

## 三个 CEP 超限

`tables/cep_current_failures.csv` 保留 C2 blocking audit 的原始三行：
`xi=0.125, 0.39375, 0.5`，当前 bracket 宽度均为 `0.125 MeV`，超过硬门禁
`0.1 MeV`。`tables/cep_refinement_plan.csv` 给出每个 bracket 的 midpoint：
`126.25`、`113.5` 和 `107.0 MeV`。这些 midpoint 需要在相同 calculation SHA
下进行新的数值切片；本包没有调用 solver，因此不声称 CEP 已闭合。

## 作者审核边界

作者可逐点确认曲线中的 S 形、三交点、Maxwell 水平线及 hybrid/oracle 差异，并在
后续 overlay 决策中记录接受、保持 ambiguous 或要求进一步数值细化。当前包不修改
容差、support ranking、production label、reference 或 transport。
