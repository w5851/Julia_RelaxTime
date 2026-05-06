# charged/freeze-out literature workflow reproduction 最小结果说明

更新日期：2026-05-04

## 1. 运行口径

- target: `data\outputs\results\relaxtime\literature\meson_density_targets\blaschke2019col_kminus_piminus_mu_pi_100_fig4_right.csv`
- freezeout profile: `default`
- path profile: `baseline_freezeout`
- flavor chemical profile: `default`
- meson chemical profile: `blaschke2019_mu_pi_100`
- regime: `phase_shift_gbu`
- workflow output: `data\outputs\results\relaxtime\meson_density\freezeout_validation\blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default\workflow_scan.csv`
- comparison output: `data\outputs\results\relaxtime\meson_density\freezeout_validation\blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default\comparison_vs_target.csv`

## 2. 当前结果定位

这是一条 **workflow reproduction / staged comparison**，不是 strict literature reconstruction，也不是正式 validation gate。

当前固定的是：

- `sqrt(s_NN)` 直接取自 literature target 的横轴采样点
- `sqrt(s_NN) -> (T, mu_B)` 使用本仓库 freeze-out baseline profile
- path strategy 通过显式 path profile 作用在 baseline 之上
- flavor-level `mu_u, mu_d, mu_s` 使用显式 flavor chemical profile
- `mu_pi` 使用显式 meson chemical profile

## 3. 数值摘要

- points: `48`
- missing rows: `0`
- converged rows: `48`
- target range: `[ 0.0, 0.1692842263342661 ]`
- model range: `[ 0.000974, 0.021497 ]`
- mean abs diff: `0.10069692442550236`
- max abs diff: `0.1477872263342661`
- mean rel diff: `0.8088425593717752`
- max rel diff: `0.8730123859410721`

## 4. 解读约束

- 该比较当前主要用于确认 charged/freeze-out workflow 与生产输出链已闭环可运行。
- 当前不对相对误差设置 pass/fail 门槛，也不把该结果直接纳入 regression gate。
- 若与文献量级仍有系统差异，应优先把它解释为 path/regime 语义差异，而不是直接回退 workflow 入口契约。
