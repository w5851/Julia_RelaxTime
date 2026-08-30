# Issue #130 RS accepted-primary matched convergence v1

## 结论

本包记录 accepted-primary runtime 下的 p104 -> p128 收敛证据。它是独立的
solver numerical evidence，不是对旧 strict-era production result/figure 的重写。

作者已接受以下语义：

- 旧 strict workflow 生成的 `prod_v2` raw result 和既有 figure 保持不变，继续作为已审核的
  正式历史产物；
- accepted-primary p104/p128 只作为当前 runtime 合同下的网格收敛证明和 reference 对照；
- accepted 与 strict 的差异首先归因于 phase-reference 输入层变化，不能把它重新套入旧 strict
  workflow 后再宣称为 accepted-primary 收敛；
- 当前差异在作者认可的物理范围内，因此不需要根据 accepted-primary 结果替换已有 raw 或 figure；
- 若未来需要发布完全基于 accepted-primary 的正式 raw/figure，应创建新的 versioned
  promotion/import case，不能覆盖旧 strict-era case。

## 固定输入与范围

- calculation SHA：`3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- workflow/contract SHA：`0448a2df5574ab55dd74b93ba5e359afd696b035`
- phase reference：`phase_reference_layer=accepted`、`phase_reference_mode=runtime`
- topology：`direct_coexistence`；保留 mode-A 的 `xi=+/-0.003` 近零两侧语义
- mode-A：`fixed-muB-phase-scaled`
- mode-B：`fixed-T-sparse-muB`
- p104：`p_num=104, tau_angle_nodes=18, tau_phi_nodes=30, tau_n_sigma=20, sigma_grid_n=440`
- p128：`p_num=128, tau_angle_nodes=20, tau_phi_nodes=36, tau_n_sigma=24, sigma_grid_n=560`
- 外部 artifact 根目录：
  `D:/Desktop/Julia_RelaxTime_issue130_artifacts/rs_accepted_primary_convergence_20260830/`

## 完整性

每档 resolution 有 30 个 shard：15 个 mode-A、15 个 mode-B；p104 与 p128 各有
`1819` 个 scan rows，合计 `60` 个 `run_manifest.json` 和 `3638` 个 scan rows。
每个 resolution 的 `channel_diagnostics.csv` 合计 `76398` 行，两档合计 `152796` 行。
所有 run manifest 的 `summary.error_count=0`、`skipped_count=0`，对应
`failed_points.csv` 均为零行；本包不复制这些全量数值曲线。

## 数值解释边界

accepted-primary p104 -> p128 的总体变化符合当前应用级收敛判断。此前 comparator 标出的
局部大相对差异保留为诊断：三个点落在已知传播子分母敏感窗口，另一个约 `1.35%` 的点位于
一阶/非一阶相结构切换附近，不能把它解释为传播子分母近零。它们不构成 solver failure、
非有限输出或 phase-reference 合同回归。

这份结论不等同于“strict 与 accepted 在逐行数值上相同”。两者使用不同 reference view，
在 mode-A 中 phase anchor 还会通过 `T=alpha_T*T_phase_base` 进入输运输入；因此存在小幅漂移
是合同预期的一部分。当前任务只判断该漂移是否在作者认可范围内，不以旧 strict 输出覆盖新
accepted evidence，也不以 accepted evidence追溯改写旧 provenance。

## 正式产物决策

| asset | 当前处理 | 说明 |
| --- | --- | --- |
| strict-era `prod_v2` raw | 保留 | 不覆盖、不改名，作为既有正式历史结果 |
| strict-era figure | 保留 | 不因 accepted-primary 收敛证据重绘 |
| accepted-primary p104 | 保留在外部 artifact | 低阶端收敛证据，不导入正式 result tree |
| accepted-primary p128 | 保留在外部 artifact | 当前 accepted runtime 的高阶端证据，不导入正式 result tree |
| accepted-primary formal replacement | 不执行 | 只有未来明确需要 accepted 正式 case 时另立 versioned PR |

## 验证与限制

本包只做 manifest/结果完整性和语义审计，不在本机调用 PNJL equilibrium solver，不改变
Maxwell、equilibrium solver、输运核、容差或 direct-coexistence route。旧 strict workflow
不重新执行；重新执行会重新施加 strict 过滤/fallback 语义，不能验证本包的 accepted-primary
合同。

详细字段、输入路径和 claim 绑定见同目录的 `manifest.json`、`tables/source_inventory.csv`、
`tables/comparison_scope.csv` 和 `tables/claim_ledger.csv`。
