# charged/freeze-out literature target 最小结果说明

更新日期：2026-05-04

## 1. 运行口径

- target: `data\outputs\results\relaxtime\literature\meson_density_targets\blaschke2019col_kminus_piminus_mu_pi_100_fig4_right.csv`
- freezeout profile: `default`
- path profile: `freezeout_plus_constT_proxy`
- flavor chemical profile: `default`
- meson chemical profile: `blaschke2019_mu_pi_100`
- regime: `stable`
- workflow output: `data\outputs\results\relaxtime\meson_density\freezeout_validation\blaschke2019col_kminus_piminus_mu_pi_100_stable_constTproxy\workflow_scan.csv`
- comparison output: `data\outputs\results\relaxtime\meson_density\freezeout_validation\blaschke2019col_kminus_piminus_mu_pi_100_stable_constTproxy\comparison_vs_target.csv`

## 2. 当前结果定位

这是一条 **结果层 staged comparison**，不是 stitched path 的严格回归真值。

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
- model range: `[ 0.0, 0.140068 ]`
- mean abs diff: `0.015277246937880863`
- max abs diff: `0.030110226334266105`
- mean rel diff: `0.1278990769416856`
- max rel diff: `0.22853257892728876`

## 4. 解读约束

- 该比较当前主要用于确认 charged/freeze-out workflow 已闭环可运行。
- 若与文献量级仍有系统差异，下一阶段优先检查 `mu_s` / path strategy，而不是先改当前 baseline workflow 入口。
