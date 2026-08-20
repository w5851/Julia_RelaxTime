# Phase-guided transport p128 xi001 tau-first 突变分析

本分析包只消费仓库内已入库的 `first_canonical_v1_p128_xi001_validated_anchored_prod_v1` 正式产物，不重跑 production，不修改主结果 CSV。分析主线是弛豫时间 `tau_*`；`eta_over_s`、`zeta_over_s` 等输运系数只作为 `tau` 与热力学量组合后的下游响应来讨论。

## 输入核对

| mode_key | scan_rows | channel_rows | failed_rows | xi_count | xi_step | plot_count | has_antiquark_tau_fields | has_zeta_over_s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mode_a | 909 | 38178 | 0 | 101 | 0.01 | 36 | True | True |
| mode_b | 909 | 38178 | 0 | 101 | 0.01 | 36 | True | True |


正式图像已补齐反夸克弛豫时间：`tau_ubar_vs_xi.png`、`tau_dbar_vs_xi.png`、`tau_sbar_vs_xi.png`。本轮也补出 `zeta_over_s_vs_xi.png`，便于把 `eta_over_s/zeta_over_s` 作为下游趋势响应一起审阅。

## Tau 突变总体分类

自动检测覆盖 mode A/B 全部 `tau_u,tau_d,tau_s,tau_ubar,tau_dbar,tau_sbar` 曲线。相邻 xi 步长为 0.01，候选突变由 log-step 的 robust outlier 规则筛出，再按同一 panel/series 的相邻 species 聚合成窗口。

| cause_verdict | count |
| --- | --- |
| channel_rate_spike_supported | 2 |
| upstream_first_order_branch_jump_supported | 2 |
| weak_or_broad_tau_variation_candidate | 6 |


最强窗口如下：

| window_id | plot_panel | plot_series | target_xi | affected_tau_fields | max_tau_step_factor | cause_verdict | top_channel_evidence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| mode_b_T200p0_muB450p0_xip0p31_1 | T200.0 | muB450.0 | 0.31 | tau_d;tau_dbar;tau_sbar;tau_u;tau_ubar | 36.44 | channel_rate_spike_supported | dbar:dubar_to_dubar(0.66);ubar:dubar_to_dubar(0.66);dbar:uubar_to_ddbar(0.17);ubar:uubar_to_ddbar(0.17);dbar:uubar_to_uubar(0.17) |
| mode_b_T120p0_muB900p0_xim0p09_1 | T120.0 | muB900.0 | -0.09 | tau_d;tau_dbar;tau_s;tau_sbar;tau_u;tau_ubar | 7.181 | upstream_first_order_branch_jump_supported | not_primary;upstream_driver=m_u;max_background_rel_step=0.5927 |
| mode_a_muB450p0_alpha1p0_xip0p26_2 | muB450.0 | alpha1.0 | 0.26 | tau_d;tau_dbar;tau_u;tau_ubar | 4.24 | channel_rate_spike_supported | ubar:uubar_to_ddbar(0.43);dbar:uubar_to_ddbar(0.43);ubar:uubar_to_uubar(0.47);dbar:uubar_to_uubar(0.47);u:uubar_to_ddbar(0.32) |
| mode_a_muB900p0_alpha1p0_xip0p00_1 | muB900.0 | alpha1.0 | 0 | tau_d;tau_dbar;tau_s;tau_sbar;tau_u;tau_ubar | 3.856 | upstream_first_order_branch_jump_supported | not_primary;upstream_driver=m_u;max_background_rel_step=0.4785 |


## 机制解释边界

- 一阶或上游分支窗口：若 `phase_reference_kind/phase_structure` 指向一阶邻域，并且 `m_u`、`Phi`、`s_fm3inv` 等背景量在同一 xi 邻域快速变化，则把 tau 突变归因到上游平衡态/分支变化。用户指出的 `mode A, muB=900, alpha_T=1.0, xi≈0` 属于这一类；虽然 `phase_curr` 字符串没有跳变，但 `m_u` 和熵密度有阶跃，因此不能只用 `phase_curr` 标签判断。
- 背景平滑但 tau 单点下探的窗口：直接近因是 channel diagnostics 中少数通道的 rate/contribution 局部尖峰。例如 `mode A, muB=450, alpha_T=1.0, xi=0.26` 由 `uubar_to_uubar` 与 `uubar_to_ddbar` 放大主导；`mode B, T=200, muB=450, xi=0.31` 由 `udbar_to_udbar/dubar_to_dubar` 放大主导。
- 对这些 channel-rate 尖峰，是否能写成“传播子分母近零”取决于下方逐窗口 denominator-chain verdict；没有通过 rate 复现或局部收敛补证的窗口仍保留为机制候选。


## 非一阶 channel-rate spike 的 denominator-chain 补证

本轮对 `channel_rate_spike_supported` 的两个非一阶窗口做定点深拆，并把 6 个 `weak_or_broad_tau_variation_candidate` 窗口全部纳入 denominator-chain 检查；`eta_over_s/zeta_over_s` 仍作为 tau 的下游响应，不进入根因判定。

| window_id | plot_panel | plot_series | observable | primary_species | mechanism_verdict | dominant_channels | dominant_denominator_branch | max_rate_reproduction_rel_error | denominator_sigma_alignment | upstream_branch_flag |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mode_b_T200p0_muB450p0_xip0p31_1 | T200.0 | muB450.0 | tau_dbar | dbar | small_denominator_supported | dubar_to_dubar;uubar_to_uubar;uubar_to_ddbar | simple_1m4KPi | 0.0 | true | false |
| mode_a_muB450p0_alpha1p0_xip0p26_2 | muB450.0 | alpha1.0 | tau_ubar | ubar | small_denominator_supported | uubar_to_uubar;uubar_to_ddbar;dubar_to_dubar | mixed_detM | 2.979989057790282e-14 | true | false |
| mode_a_muB0p0_alpha1p1_xip0p48_1 | muB0.0 | alpha1.1 | tau_dbar | dbar | small_denominator_supported | dubar_to_dubar;uubar_to_uubar;uubar_to_ddbar | simple_1m4KPi | 5.626880713523278e-13 | true | false |
| mode_a_muB900p0_alpha1p1_xip0p18_1 | muB900.0 | alpha1.1 | tau_dbar | dbar | small_denominator_supported | dubar_to_dubar;uubar_to_uubar;uubar_to_ddbar | simple_1m4KPi | 1.137760542769262e-14 | true | false |
| mode_a_muB450p0_alpha1p1_xip0p35_1 | muB450.0 | alpha1.1 | tau_dbar | dbar | small_denominator_supported | dubar_to_dubar;uubar_to_uubar;uubar_to_ddbar | simple_1m4KPi | 0.0 | true | false |
| mode_a_muB900p0_alpha1p2_xip0p49_1 | muB900.0 | alpha1.2 | tau_dbar | dbar | small_denominator_supported | dubar_to_dubar;uubar_to_uubar;uubar_to_ddbar | simple_1m4KPi | 4.401128194660395e-15 | true | false |
| mode_a_muB450p0_alpha1p0_xim0p20_1 | muB450.0 | alpha1.0 | tau_dbar | dbar | small_denominator_supported | dubar_to_dubar;uubar_to_uubar;uubar_to_ddbar | simple_1m4KPi | 8.077694440113023e-16 | true | false |
| mode_b_T200p0_muB0p0_xim0p21_1 | T200.0 | muB0.0 | tau_dbar | dbar | small_denominator_supported | dubar_to_dubar;uubar_to_uubar;uubar_to_ddbar | simple_1m4KPi | 0.0 | true | false |


- `mode_a_muB450p0_alpha1p0_xip0p26_2`：在 denominator-chain 证据层支持小分母机制。`uubar_to_uubar/uubar_to_ddbar` 贡献覆盖主导份额，`sigma(s)` 峰和 mixed `detM` 峰都落在近阈值 band；修正机制脚本口径后，直调 rate 与 channel diagnostics 的复现误差为机制表中的机器精度量级。该窗口的上游背景量平滑，因此近因不是一阶相变或上游分支突跳；但本轮没有完成额外 high-rate gate，论文表述应保留 scope limit。
- `mode_b_T200p0_muB450p0_xip0p31_1`：在 denominator-chain 证据层也支持小分母机制。`dubar_to_dubar/uubar_to_uubar/uubar_to_ddbar` 的 `sigma(s)` 峰和 simple `1-4KΠ` 峰在近阈值 band 对齐，且上游背景量平滑；修正机制脚本的 A-builder 口径后，直调 `average_scattering_rate` 与生产 `channel_diagnostics.csv` 的 rate 复现误差为 `0`。
- 此前 mode B 约 `0.47` 的 rate 复现误差来自诊断脚本 bug：机制脚本手工重建传播子 A 场时没有使用 production workflow 的 `a_builder` 配置。当前机制脚本已显式记录并使用 `p_nodes=16,p_max=20.0,cos_nodes=4,use_aniso=true`；这不是 production 数据 bug。
- 6 个 weak/broad 窗口全部完成深拆，均显示 simple `1-4KΠ` 小分母支持；新增的 4 个窗口是 `mode_a_muB450p0_alpha1p1_xip0p35_1`、`mode_a_muB900p0_alpha1p2_xip0p49_1`、`mode_a_muB450p0_alpha1p0_xim0p20_1` 和 `mode_b_T200p0_muB0p0_xim0p21_1`。
- 本机尝试的局部 high-rate convergence gate 未在可控时间内完成；`local_rate_convergence_gate.csv` 因此只保留表头。当前结论是 denominator-chain 补证，不是新的 production-grade 收敛证明。


## 下游输运系数响应

`eta_over_s` 和 `zeta_over_s` 不是本轮突变根因入口。定性上，`eta`、`zeta` 是带有 `tau` 权重的动量积分，`eta_over_s`、`zeta_over_s` 还要除以熵密度。因此：

- 在 channel-rate spike 窗口，背景熵密度通常平滑，`tau` 单点下探会直接拖低 `eta_over_s` 与 `zeta_over_s`，表现为同 xi 的下游凹陷。
- 在一阶/上游分支窗口，`tau` 常上升，同时 `s_fm3inv` 下降；除以熵密度后，`eta_over_s/zeta_over_s` 的跳变会比 `eta/zeta` 本身更显著。
- 若某个下游 ratio 与 tau 方向不完全一致，应优先检查熵密度和对应 transport numerator，而不是把 ratio 曲线直接当成散射机制证据。

## 关键图

- `docs/analysis/relaxtime/phase_guided_transport_p128_xi001_analysis/figures/tau_jump_windows_overview.png`
- `docs/analysis/relaxtime/phase_guided_transport_p128_xi001_analysis/figures/downstream_response_windows.png`

## 产物表

- `tables/input_inventory.csv`：输入行数、hash、图像数量和字段核对。
- `tables/tau_jump_step_candidates.csv`：逐 species 的相邻 xi 突变候选。
- `tables/tau_jump_window_summary.csv`：聚合窗口、机制 verdict、上游背景量和主导通道摘要。
- `tables/tau_jump_channel_attribution.csv`：窗口内 top channel 的 contribution/rate 局部变化。
- `tables/downstream_transport_response_summary.csv`：`eta_over_s/zeta_over_s/sigma_over_T` 对 tau 突变的下游响应关系。
- `tables/mechanism_window_candidates.csv`：非一阶 channel-rate spike 的 denominator-chain 深拆候选。
- `tables/mechanism_window_summary.csv`：已深拆窗口的机制 verdict。
- `tables/denominator_chain_summary.csv`、`tables/denominator_chain_band_table.csv`、`tables/denominator_ds_samples.csv`：传播子分母、`sigma(s)` 与 rate band 的局部证据。
- `tables/local_rate_reproduction_mismatch_root_cause.csv`：mode B 早期 rate 复现偏差的诊断脚本口径根因。
- `tables/local_rate_convergence_gate.csv`：本机局部 high-rate gate 输出；当前仅有表头，表示未在本轮可控预算内完成。
- `tables/claim_ledger.csv`：可写入论文或需要作者确认的 claim 账本。

## 作者确认项

- xi001 新中间点继承 p128 积分参数与旧 xi=0.05 锚点比较证据，但非锚点局部尖峰还没有逐点 p104->p128 gate。
- 目前的机制支持来自 production 口径 rate 复现与 denominator-chain 对齐；若要把这些局部结构写成更强的 production-grade 收敛结论，仍应补可完成的局部高精度 convergence gate。
