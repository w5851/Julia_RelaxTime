# C1 `(mu_xi_T)` 诊断相面

本目录基于 C1 artifact `31762201725`（calculation SHA `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`）生成，不调用 solver。

图中蓝色半透明曲面是 `boundary_*.csv` 的 Maxwell 共存化学势，橙色半透明曲面是 `crossover_*.csv` 的 peak crossover。坐标顺序为 `(mu, xi, T)`。

本图只是诊断可视化，不是 phase-reference 晋升结果。C1 的 `phase_grid_convergence` 共 11375 条记录，其中 1629 条未闭合；这些区域没有被强行插值成物理结论。C1 CEP 仍为 ambiguous，因此没有绘制单值 CEP 曲线；`figures/c1_cep_temperature_brackets_vs_xi.png` 仅以半透明带显示每个 ξ 的 `[T_bracket_low, T_bracket_high]`，中线明确标为 diagnostic only。

`crossover` 是对每个采样 `(xi,mu)` 独立在温区内寻找 `peak`，原始 CSV 没有用 Maxwell 面做物理遮罩。因此 `tables/crossover_vs_maxwell_audit.csv` 只做几何关系审计，不把一阶区域内的响应峰自动解释成第二张物理 crossover 面。`tables/crossover_maxwell_endpoint_separation.csv` 给出每个 ξ 的 Maxwell 最高温度与 crossover 最低温度差；它用于检查两张面是否在同一温度域重叠。

输入审计、哈希、三角剖分数量、CEP bracket 图元数据、派生点表、Maxwell 支持上边界、crossover 关系审计、端点温度分离和 claim ledger 见 `manifest.json`、`figures/plot_manifest.json`、`figures/c1_cep_temperature_brackets_vs_xi.json`、`tables/surface_point_summary.csv`、`tables/maxwell_surface_support_audit.csv`、`tables/crossover_vs_maxwell_audit.csv`、`tables/crossover_maxwell_endpoint_separation.csv` 和 `tables/claim_ledger.csv`。
