# PNJL CEP endpoint-local production shadow v4

本目录是 full 24-anchor endpoint-local shadow 的诊断 evidence 包，不是
reference promotion。`full_hybrid_candidate` 只表示当前 gate 通过，仍需作者审核
代表性 rho–mu 曲线、三个 endpoint 证书和成本，再决定是否进入 C0→C1→C2。

- calculation SHA: `8fb86255c7894004cb1c52bbb03f9f53a4828411`
- postprocess SHA: `6ceaa1ea73ca684f7333556c5d9e0492bc567899`
- numerical source run: `31067984700`
- aggregate replay run: `31068641378`
- approved deep run: `31065343978`
- verdict: `full_hybrid_candidate`
- evidence state: `final`
- figure policy: `independent_rho_mu_zoom_with_phase_markers_v2`
- no-support controls use the display-only smooth-window policy
  `smooth_window_rho_mu_zoom_v1`

仓库只导入聚合表和代表性 PNG；完整 `curve_points.csv` 保留在 Actions/local artifact，
通过 `tables/curve_index.csv`、source manifest 和 SHA 追溯。三态物理语义、Maxwell
容差、geometry gate、equilibrium solver、旧 reference 和 transport 均未被本包改写。
