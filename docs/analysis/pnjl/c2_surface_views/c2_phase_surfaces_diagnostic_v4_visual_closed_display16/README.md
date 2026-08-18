# C2 diagnostic phase surfaces v4 (visual closed)

本包的 x 轴为夸克化学势 `mu_q`，y 轴为 `xi`，z 轴为 `T`。Maxwell 与 crossover 的物理筛选遵循同一条规则：同一 `(xi, mu_q)` 不同时赋予两种物理含义；在一阶区 `mu_q > mu_CEP` 的偏导响应峰只作为诊断数据，不绘制为 crossover。

本 v4 版本为 `visualization-only`：所有有限/converged Maxwell boundary 行（包括自动 geometry/interpolation 未闭合行）统一以蓝色绘制，隐藏状态颜色和高于 CEP 的灰色响应峰。它不生成合成数值行、不改变源曲线；仅使用 manifest 明确记录的 16 MeV display-only 三角化上限改善视觉连续性，不把未闭合证书升级为闭合。
`tables/crossover_physical_filter.csv` 记录响应峰筛选，`tables/crossover_cep_sampling_gap.csv` 记录每个 xi 切片最后保留点到首个被排除点的原生 mu 间隔；大 xi 处的视觉缺口主要由该离散间隔造成，不能跨 gap 隐式补线。`tables/maxwell_surface_point_status.csv` 保留原始 grid 状态，即使 v4 图中不再用颜色区分。CEP 仍保留 bracket 和 midpoint 语义，midpoint 不是 confirmed CEP。

主图见 `figures/c2_phase_surfaces_mu_xi_T_visual_closed.png`，证据表见 `tables/claim_ledger.csv`。verdict 为 `visualization_only_closed_surface_ready`，不产生 reference promotion。
