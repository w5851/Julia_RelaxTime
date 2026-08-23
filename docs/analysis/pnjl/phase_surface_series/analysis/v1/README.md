# C2 diagnostic phase surfaces v1

本包用于作者从全局相图检查相结构拓扑，不替代逐点 Maxwell/geometry 证书。x 轴为 `mu_B`，y 轴为 `xi`，z 轴为 `T`。Maxwell 面使用 C2 finite/converged boundary 输出；crossover 面把 `mu_q` 转换为 `mu_B=3*mu_q`；CEP 使用最后一阶点的化学势投影，并保留温度 bracket 和中点语义。

C2 grid 未闭合状态叠加在 Maxwell 面上，颜色分别对应 rho geometry、temperature interpolation 和 xi interpolation；图中没有用插值填补数据缺口，也没有把 CEP 中点写成物理 CEP。当前 verdict 为 `global_diagnostic_surface_ready_for_author_review`，不产生 reference promotion。

主图见 `figures/c2_phase_surfaces_muB_xi_T_diagnostic.png`，状态表见 `tables/maxwell_surface_point_status.csv`，证据边界见 `tables/claim_ledger.csv`。
