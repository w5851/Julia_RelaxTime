# Phase-reference limited evidence audit v1

本包对 C2 中 452 个 `ok`、272 个 `hybrid_stage_c_not_converged` 和 39 个 classification-transition rho 记录进行限定范围的 solver-free raw-oracle 审计；763 行去重后为 761 个坐标，其中 2 个坐标含多个 unresolved level。每个坐标从 Zenodo 归档读取完整 `rho=0:0.003125:4` 曲线，检查 1,281 个点、finite/converged、重复 rho 和 `mu_B(rho)` 的离散 `+→−→+` 拓扑。

当前 verdict 为 `raw_curve_coverage_complete_diagnostic_only`。`s_shape_present` 只说明独立 raw 曲线存在 S 形拓扑；它不等于 Maxwell 面积/geometry 证书，也不会覆盖 hybrid/C2 标签。因此本包不触发 reference promotion、不重跑 C0/C1/C2，也不调用 equilibrium solver。

详细结果见 `tables/raw_curve_audit.csv`、`tables/reason_shape_summary.csv`、`tables/shape_threshold_sensitivity.csv`、`tables/claim_ledger.csv`、`decision.json` 和代表性 `figures/`。
