# C2 diagnostic phase surfaces v3

本包用于作者从全局相图检查相结构拓扑，不替代逐点 Maxwell/geometry 证书。x 轴为源数据的夸克化学势 `mu_q`（图上记作 `mu_q`），y 轴为 `xi`，z 轴为 `T`。Maxwell 面使用 C2 finite/converged boundary 输出。crossover 面只保留每个 `xi` 切片中 `mu_q <= mu_CEP` 代理的响应峰；`mu_q > mu_CEP` 的响应峰仍保留为灰色审计点，但不再赋予 crossover 物理含义。CEP 保留温度 bracket 和中点语义。

`tables/crossover_physical_filter.csv` 记录每个响应峰的物理筛选，`tables/crossover_maxwell_endpoint_separation.csv` 显示各 xi 切片的两面端点温度差，避免把离散端点间隔误读为双侧 crossover。C2 grid 未闭合状态叠加在 Maxwell 面上，颜色分别对应 rho geometry、temperature interpolation 和 xi interpolation；图中没有用插值填补数据缺口，也没有把 CEP 中点写成物理 CEP。当前 verdict 为 `global_diagnostic_surface_ready_for_author_review`，不产生 reference promotion。

主图见 `figures/c2_phase_surfaces_mu_xi_T_diagnostic.png`，状态表见 `tables/maxwell_surface_point_status.csv`，证据边界见 `tables/claim_ledger.csv`。
