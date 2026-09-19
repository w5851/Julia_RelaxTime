# Issue #130 RS `publication_clean_v2_residual_smoothed`

## 目的与边界

本包是现有 `publication_clean_v2` 的第二层 display-only 派生。它只消费 v2 的
`publication_clean_points.csv`，不调用 equilibrium/transport solver，不修改 raw CSV、
production registry、v2 原目录或 `main`。

本轮处理的是 v2 仍可见的局部显示结构：小分母敏感窗口的三点肩部、mode-A
`muB=900` 的 tau 局部斜率鼓包，以及 `alpha_T=1.1` tau 曲线的两段快速上升。
这不是传播子有限宽度、分母截断或新的物理计算。

## 显示算法

对每个显式 recipe 窗口，保留左右 anchor 的 v2 display 值；在正值曲线的
`log(y)` 空间构造单调 Hermite 趋势，再按 recipe 中的 `strength` 与原 v2 值作
几何混合。`strength=1` 用于孤立小突变，较小强度用于连续 branch response。
候选曲线还必须通过逐曲线 roughness gate：若窗口内最大的相邻 log-y 斜率跳变
没有下降，则该 observable 保留 v2 值，避免对原本平滑的曲线制造新拐点。

实际渲染的一阶 gap 是硬保护区：窗口不得进入 gap 内部或跨越 gap。branch-local
窗口只允许在同一个 `phase_curr` 分支内处理，因此不会把两个相态拼成一条曲线。

## 输入 provenance

- source artifact: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v2`
- source manifest SHA256: `c39bd3dd34869bd4c1b9c76f53f86fbfe677b2e27deaebb5ff2c9a7aa1558fea`
- source v2 point rows: `21828`
- recipe: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_v2_residual_smoothing/tables/residual_smoothing_windows.csv`
- recipe SHA256: `71464033d1fe282fddd4c075f8c136d047ce15cf5ef5746a64277b0631cae499`
- derived solver called: `false`
- canonical/raw data modified: `false`
- manuscript eligible: `false`

## 结果摘要

- recipe windows: 14
- applied windows: 14
- retained v2 windows: 0
- applied observable curves: 103
- roughness-guard retained curves: 7
- residual-smoothed point rows: 329
- mean absolute local log residual: `0.0114411` -> `0.00544369`
- point residual comparison: 257 improved / 72 worsened
- publication figures: 72
- audit figures: 2

`tables/publication_clean_points.csv` 同时保留 `raw_value`、`v2_clean_value` 和最终
`clean_value`。改变点还记录 `smoothing_window`、锚点、拟合方法、phase gate 状态及
`local_residual_before/after`；所有改变点都在 `tables/residual_point_audit.csv` 中逐点记录，
逐曲线 roughness 决策见 `tables/residual_curve_audit.csv`。

## 当前处理决策与作者审阅项

- roughness gate 保留 v2 的曲线数：7；这些曲线不应继续套用同一窗口。当前保留项为：`mode_b_T200p0_muB0p0_xim0p21_window:tau_s; mode_b_T200p0_muB0p0_xim0p21_window:tau_sbar; mode_b_T200p0_muB900p0_xim0p10_window:zeta_over_s; mode_b_T200p0_muB900p0_xim0p10_window:tau_u; mode_b_T200p0_muB900p0_xim0p10_window:tau_d; mode_a_muB450p0_alpha1p0_xip0p26_window:tau_s; mode_a_muB450p0_alpha1p0_xip0p26_window:tau_sbar`。
- applied 曲线中相对改动达到 5% 的项目列在下表；它们只适合作者确认，不应自动视为投稿最终值：

| window_id | observable | max relative change | roughness ratio |
| --- | --- | ---: | ---: |
| mode_a_muB900p0_alpha1p0_tau_quark_branch_ubar_dbar | tau_ubar | 0.09589 | 0.6425 |
| mode_a_muB900p0_alpha1p0_tau_quark_branch_ubar_dbar | tau_dbar | 0.09589 | 0.6425 |
| mode_a_muB900p0_alpha1p1_tau_rise_one_ubar_dbar | tau_ubar | 0.07412 | 0.4634 |
| mode_a_muB900p0_alpha1p1_tau_rise_one_ubar_dbar | tau_dbar | 0.07412 | 0.4634 |
| mode_b_T200p0_muB900p0_xim0p10_window | tau_sbar | 0.05799 | 0.005285 |
| mode_b_T200p0_muB0p0_xim0p21_window | tau_ubar | 0.05469 | 0.2088 |
| mode_b_T200p0_muB0p0_xim0p21_window | tau_dbar | 0.05469 | 0.2088 |
| mode_b_T200p0_muB0p0_xim0p21_window | tau_u | 0.05469 | 0.2088 |
| mode_b_T200p0_muB0p0_xim0p21_window | tau_d | 0.05469 | 0.2088 |

- `muB=900, alpha_T=1.0/1.1` 的反夸克 tau 是连续 branch response，处理目标只是降低局部斜率尖锐度；若作者不接受约 5--10% 的显示改动，应降低对应 recipe strength 或直接保留 v2。
- `mode B / T=200, muB=900 / tau_sbar` 的窗口虽可降低图面肩部，但 high-rate convergence gate 仍为空；在补 gate 前只能保留为 channel-rate candidate，不能以平滑图面替代机制证据。
- `tables/review_adjustment_map.csv` 中的 author-review 插值是输入侧审阅记录，不代表本层已自动应用；任何采用都应重新生成并核对 manifest。
- 一阶 gap、端点和跨分支连接继续禁止填补；连续宽响应不应通过扩大窗口被抹平。

## 解释边界

1. 这些值仍是 display-only 派生值，不是 solver 重算、收敛证明或物理正则化。
2. 小分母窗口的机制证据仍应引用原有 mechanism audit；本包只改变论文图的显示形状。
   特别是 `mode B / T=200 MeV / muB=900 MeV / tau_sbar` 仍是尚未完成
   high-rate convergence 的 channel-rate candidate，不能因图形变平而升级为已证实机制。
3. `muB=900, alpha_T=1` 的 `xi=-0.003/+0.003` 和 mode-B `T=120, muB=900`
   的 `xi=-0.13/-0.12` 一阶端点继续使用 v2 raw endpoint，gap 不被填充。
4. 任何定量峰值、临界行为或输运系数精确拟合仍必须使用 raw/v2 evidence，不能用
   本层的平滑值替代。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v2_residual_smoothed.py
python -m pytest tests/unit/python/test_phase_guided_publication_clean_v2_residual_smoothed.py
```
