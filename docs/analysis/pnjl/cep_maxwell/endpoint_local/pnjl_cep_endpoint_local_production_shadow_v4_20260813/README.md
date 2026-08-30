# PNJL CEP endpoint-local production shadow v4

本目录是 full 24-anchor endpoint-local shadow 的诊断 evidence 包，不是
reference promotion。`full_hybrid_candidate` 只表示当前 gate 通过，仍需作者审核
代表性 rho–mu 曲线、三个 endpoint 证书和成本，再决定是否进入 C0→C1→C2。

- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- postprocess SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- numerical source run: `31713534102`
- aggregate replay run: `31714535418`
- approved deep run: `31710995191`
- verdict: `full_hybrid_candidate`
- evidence state: `final`
- figure policy: `independent_rho_mu_zoom_with_phase_markers_v3_tight_envelope`
- first-order local rho padding: `max(4% of phase/support envelope, 0.01 rho)`
- local mu-axis display padding: `max(8% of local mu span, 0.0002 MeV)`
- no-support controls use the display-only smooth-window policy
  `smooth_window_rho_mu_zoom_v1`

仓库只导入聚合表和代表性 PNG；完整 `curve_points.csv` 保留在 Actions/local artifact，
通过 `tables/curve_index.csv`、source manifest 和 SHA 追溯。三态物理语义、Maxwell
容差、geometry gate、equilibrium solver、旧 reference 和 transport 均未被本包改写。
