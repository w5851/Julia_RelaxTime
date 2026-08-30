# C1 `(xi,T,mu)` 诊断相面

本目录基于 C1 artifact `31762201725`（calculation SHA `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`）生成，不调用 solver。

图中蓝色半透明曲面是 `boundary_*.csv` 的 Maxwell 共存化学势，橙色半透明曲面是 `crossover_*.csv` 的 peak crossover。坐标顺序为 `(xi, T, mu)`。

本图只是诊断可视化，不是 phase-reference 晋升结果。C1 的 `phase_grid_convergence` 共 11375 条记录，其中 1629 条未闭合；这些区域没有被强行插值成物理结论。C1 CEP 仍为 ambiguous，因此未绘制 CEP 曲线。

输入审计、哈希、三角剖分数量、派生点表和 claim ledger 见 `manifest.json`、`figures/plot_manifest.json`、`tables/surface_point_summary.csv` 和 `tables/claim_ledger.csv`。
