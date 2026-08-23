# C2 diagnostic phase curves v5 (no triangulation)

本包的 x 轴为夸克化学势 `mu_q`，y 轴为 `xi`，z 轴为 `T`。Maxwell 与 crossover 的物理筛选遵循同一条规则：同一 `(xi, mu_q)` 不同时赋予两种物理含义；在一阶区 `mu_q > mu_CEP` 的偏导响应峰只作为诊断数据，不绘制为 crossover。

本 v5 模式完全禁用三角化：Maxwell 与 crossover 只绘制原生有序 support 的相邻线段，超过显式采样间隔门限的 gap 保持为空，不填面、不跨 gap、不生成合成点。该图用于区分数据 support 缺失与三角化造成的视觉空洞，不替代逐点 Maxwell/geometry 证书。
`tables/crossover_physical_filter.csv` 记录响应峰筛选，`tables/crossover_cep_sampling_gap.csv` 记录每个 xi 切片最后保留点到首个被排除点的原生 mu 间隔；大 xi 处的视觉缺口主要由该离散间隔造成，不能跨 gap 隐式补线。`tables/maxwell_surface_point_status.csv` 保留原始 grid 状态；`tables/grid_unresolved_diagnostics.csv` 记录 unresolved 行是否保留细化层级/几何诊断字段，以及是否存在同坐标的 converged boundary 行。这里的 `left/right/midpoint` 是细化层级元数据，不是 Maxwell 化学势 bracket。CEP 仍保留 bracket 和 midpoint 语义，midpoint 不是 confirmed CEP。

主图见 `figures/c2_phase_surfaces_mu_xi_T_no_triangulation.png`，证据表见 `tables/claim_ledger.csv`。verdict 为 `diagnostic_no_triangulation_ready`，不产生 reference promotion。
