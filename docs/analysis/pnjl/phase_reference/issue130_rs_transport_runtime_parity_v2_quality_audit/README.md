# Issue #130 RS numerical pilot quality audit v1

这是对 numerical pilot source artifact 的 solver-free 后处理审计；不调用 equilibrium solver，不改写数值 CSV。

- source run: `32684074876`
- workflow head: `27e9642d431ba7afd845f2b187f77c0fbbe3be4d`
- calculation SHA: `3c5f6b3c9bd535cff7657364dadb2efc31f2ea48`
- original in-run verdict: `pilot_solver_or_curve_failure`
- replay classification: `pilot_pair_complete_with_common_quality_warnings_diagnostic_only`
- common `tau_u_ubar_ratio_high` warnings: `5`

## 结论

5 个质量警告点在 candidate runtime 与 legacy 两套 reference 中完全成对出现，且原因相同；两套运行均 19/19 点完成，无 failed point、NaN/Inf、重复键或 equilibrium nonconvergence。
这些点的 `tau_u/tau_ubar` 比值超过既有 `scan_quality` 阈值 6，是 transport 质量审查标记，不等价于 solver failure。高 `muB=900 MeV` 下 `n_u/n_ubar` 强烈不对称，并由 u/ubar 通道贡献共同体现；该归因只支持到“共同质量警告”，不升级为传播子机制结论。

## Reference 差异为何会传播到输运

本 pilot 的 `mode_a_fixed_muB_phase_scaled` 并非固定同一温度后只切换 phase label。plan 先从 reference 插值得到 `T_phase_base_MeV`，再使用 `T=alpha_T*T_phase_base_MeV` 调用 equilibrium 与 transport。candidate/legacy 的 anchor 温度因此可以不同；即使最终 phase kind/structure 相同，连续的 T、平衡态、质量、Polyakov loop、密度和散射率仍会略有差异。

## 证据边界

本包支持：workflow 完成、共同质量警告、reference-specific transport regression 未见明显信号。它不支持：删除 legacy reference、宣称全域 RS numerical convergence、或把质量警告自动视为可忽略的物理异常。

详见 `collector_replay/`、`tables/quality_warning_summary.csv`、`tables/channel_attribution.csv`、`tables/anchor_drift_summary.csv` 和 `claim_ledger.csv`。
