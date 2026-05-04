# charged/freeze-out literature target 最小结果说明

更新日期：2026-05-03

## 1. 运行口径

- target: `data\outputs\results\relaxtime\literature\meson_density_targets\blaschke2019col_kminus_piminus_mu_pi_100_fig4_right.csv`
- freezeout profile: `default`
- meson chemical profile: `blaschke2019_mu_pi_100`
- regime: `phase_shift_current`
- workflow output: `data\outputs\results\relaxtime\meson_density\freezeout_validation\blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_current_default\workflow_scan.csv`
- comparison output: `data\outputs\results\relaxtime\meson_density\freezeout_validation\blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_current_default\comparison_vs_target.csv`

## 2. 当前结果定位

这是一条 **结果层 staged comparison**，不是 stitched path 的严格回归真值。

当前固定的是：

- `sqrt(s_NN)` 直接取自 literature target 的横轴采样点
- `sqrt(s_NN) -> (T, mu_B)` 使用本仓库 freeze-out baseline profile
- `mu_pi` 使用显式 meson chemical profile
- 尚未引入 flavor-level `mu_s` 现象学链

## 3. 数值摘要

- points: `48`
- missing rows: `0`
- converged rows: `48`
- target range: `[ 0.0, 0.1692842263342661 ]`
- model range: `[ 0.007191, 0.014988 ]`
- mean abs diff: `0.10555874333243385`
- max abs diff: `0.1542962263342661`
- mean rel diff: `0.8638957885518251`
- max rel diff: `1.6462184929385186`

## 4. 解读约束

- 该比较当前主要用于确认 charged/freeze-out workflow 已闭环可运行。
- 若与文献量级仍有系统差异，下一阶段优先检查 `mu_s` / path strategy，而不是先改当前 baseline workflow 入口。
