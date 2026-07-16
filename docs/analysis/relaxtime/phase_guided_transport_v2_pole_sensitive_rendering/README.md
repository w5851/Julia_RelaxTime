# Phase-guided transport v2 极点敏感派生显示审计

## 范围与结论

本分析包消费 `first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1`，不修改 production CSV、canonical figure 或 `production_registry.json`。它把旧 xi001 denominator-chain 机制表迁移到 v2，并生成**仅供作者审阅**的非破坏性派生显示候选。

当前结论：8 个窗口通过小分母机制迁移门槛；2 个一阶/上游分支窗口被硬保护。虚线桥接是视觉指南，不是新计算值，也不能静默替换正式图。

## 输入与迁移门槛

| mode | scan rows | diagnostic rows | max tau v1/v2 drift | max rate v1/v2 drift |
| --- | ---: | ---: | ---: | ---: |
| mode_a | 909 | 38178 | 0.00016650389 | 0.00053765629 |
| mode_b | 909 | 38178 | 8.2279794e-06 | 0.00037216686 |

迁移门槛为 `0.001`。v2 的能量语义修改不改变弛豫时间定义；实际 tau 和 channel-rate 漂移均显著低于门槛，因此允许继承已有 denominator-chain 定点证据。该门槛不替代尚未完成的局部 high-rate convergence gate。

## Pole-sensitive display mask 候选

| window | panel | series | xi | denominator branch | eligible observables |
| --- | --- | --- | ---: | --- | --- |
| mode_b_T200p0_muB450p0_xip0p31_1 | T200.0 | muB450.0 | 0.31 | simple_1m4KPi | tau_u;tau_d;tau_ubar;tau_dbar;tau_sbar;eta;eta_over_s;zeta;zeta_over_s;sigma;sigma_over_T |
| mode_a_muB450p0_alpha1p0_xip0p26_2 | muB450.0 | alpha1.0 | 0.26 | mixed_detM | tau_u;tau_d;tau_ubar;tau_dbar;eta;eta_over_s;zeta;zeta_over_s;sigma;sigma_over_T |
| mode_a_muB0p0_alpha1p1_xip0p48_1 | muB0.0 | alpha1.1 | 0.48 | simple_1m4KPi | tau_u;tau_d;tau_ubar;tau_dbar;eta;eta_over_s;zeta;zeta_over_s;sigma;sigma_over_T |
| mode_a_muB900p0_alpha1p1_xip0p18_1 | muB900.0 | alpha1.1 | 0.18 | simple_1m4KPi | tau_ubar;tau_dbar |
| mode_a_muB450p0_alpha1p1_xip0p35_1 | muB450.0 | alpha1.1 | 0.35 | simple_1m4KPi | tau_u;tau_d;tau_ubar;tau_dbar;eta;eta_over_s;zeta;zeta_over_s;sigma;sigma_over_T |
| mode_a_muB900p0_alpha1p2_xip0p49_1 | muB900.0 | alpha1.2 | 0.49 | simple_1m4KPi | tau_ubar;tau_dbar |
| mode_a_muB450p0_alpha1p0_xim0p20_1 | muB450.0 | alpha1.0 | -0.20 | simple_1m4KPi | tau_ubar;tau_dbar |
| mode_b_T200p0_muB0p0_xim0p21_1 | T200.0 | muB0.0 | -0.21 | simple_1m4KPi | tau_u;tau_d;tau_ubar;tau_dbar;sigma;sigma_over_T |

规则：

1. tau 只有在旧分析的 `affected_tau_fields` 中才进入 mask；
2. transport 只有在目标点相对相邻两点均值的绝对 log 偏离不小于 `0.03` 时才进入 mask；
3. raw 点始终以橙色叉号保留；solid guide 在 eligible 点断开，并以虚线连接两侧真实邻点；
4. `linear_neighbor_guide_value` 只用于画图，不写回 production，也不声称是物理值。

## 一阶窗口硬保护

| window | panel | series | xi | phase structure |
| --- | --- | --- | ---: | --- |
| mode_b_T120p0_muB900p0_xim0p09_1 | T120.0 | muB900.0 | -0.09 | first_order_possible |
| mode_a_muB900p0_alpha1p0_xip0p00_1 | muB900.0 | alpha1.0 | 0.00 | first_order_possible |

这两个窗口保留 raw 曲线和跳变，不应用桥接或 mask。保护依据是背景质量、Polyakov loop 和熵密度的同步快速变化，而不是单纯依赖 phase 字符串标签。

## 派生图

- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_b_T200p0_muB450p0_xip0p31_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB450p0_alpha1p0_xip0p26_2.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB0p0_alpha1p1_xip0p48_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB900p0_alpha1p1_xip0p18_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB450p0_alpha1p1_xip0p35_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB900p0_alpha1p2_xip0p49_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_a_muB450p0_alpha1p0_xim0p20_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/pole_sensitive_mode_b_T200p0_muB0p0_xim0p21_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/first_order_protected_mode_b_T120p0_muB900p0_xim0p09_1.png`
- `docs/analysis/relaxtime/phase_guided_transport_v2_pole_sensitive_rendering/figures/first_order_protected_mode_a_muB900p0_alpha1p0_xip0p00_1.png`

每个 pole-sensitive 图包含 primary tau 与 `eta_over_s/zeta_over_s/sigma_over_T` 四个局部面板。每个一阶图包含相同的下游比值与 `tau_u`，用于核对跳变未被删除。

## 证据边界与作者判断

- `supported`：v2 输入完整、无 failed/NaN/负 rate，v1→v2 tau/rate 迁移门槛通过；一阶窗口保护规则明确。
- `supported_with_scope_limit`：8 个窗口已有 denominator-chain 与生产 rate 复现证据，但没有逐窗口新的 high-rate convergence gate。
- `author_check`：是否允许把带显式标注的派生图用于论文展示；是否需要计算层有限宽度/极点正则化和新 production slug。
- 本包不把小分母结构直接定性为随机数值噪声，也不把虚线桥接升级为物理预测。

## 复现

```powershell
python scripts/analysis/relaxtime/build_phase_guided_pole_sensitive_rendering.py
```

关键表：

- `tables/input_inventory.csv`
- `tables/window_classification.csv`
- `tables/rendered_point_audit.csv`
- `tables/pole_sensitive_mask.csv`
- `tables/first_order_protection.csv`
- `tables/claim_ledger.csv`
- `figures/plot_manifest.json`
