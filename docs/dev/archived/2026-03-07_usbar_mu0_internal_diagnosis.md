---
title: usbar→usbar 在 μ_B=0 下的内部口径排查记录
archived: true
original: docs/dev/active/2026-03-07_usbar_mu0_internal_diagnosis.md
archived_date: 2026-03-08
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# usbar→usbar 在 μ_B=0 下的内部口径排查记录

**记录日期**：2026-03-07  
**范围**：`usbar_to_usbar`，优先聚焦 `μ_B=0`、`T=210/250 MeV`  
**目的**：避免重复研究；沉淀当前诊断脚本、结果路径与已验证结论。

## 1. 当前问题定义

- 文献比较显示 `usbar_to_usbar` 是 `σ(√s)` 文献验证中的主要系统性异常过程。
- `μ_B=0` 下存在两类不同失配：
  - `T=210 MeV`：模型峰位偏右且峰高过大。
  - `T=250 MeV`：模型峰位基本对齐，但峰高整体偏大。

## 2. 已完成诊断脚本

- [scripts/debug/diag_usbar_mu0_meson_poles.jl](scripts/debug/diag_usbar_mu0_meson_poles.jl)
  - 输出 `T=210/250 MeV, μ_B=0` 的 `K`、`sigma_K` 极点质量、宽度、阈值和 gap。
- [scripts/debug/diag_usbar_mu0_channel_decomp.jl](scripts/debug/diag_usbar_mu0_channel_decomp.jl)
  - 在总截面同口径 `t` 积分下拆分 `sP/sS/tS/tP/interference`。
- [scripts/debug/diag_usbar_T210_width_coupling.jl](scripts/debug/diag_usbar_T210_width_coupling.jl)
  - 比较零宽度与有限宽度 `K` 分母零点位置。
- [scripts/debug/diag_usbar_T210_polarization_breakdown.jl](scripts/debug/diag_usbar_T210_polarization_breakdown.jl)
  - 比较 `B0(+k0)`、`B0(-k0)` 与对称平均后的 `Π_us^P` 和分母。
- [scripts/debug/diag_usbar_T210_b0_terms_sensitivity.jl](scripts/debug/diag_usbar_T210_b0_terms_sensitivity.jl)
  - 展开 `B0 = term1 - term2 + term3 - term4`，并做 `A_sum`、`K4567_plus`、`Re B0` 的单参数灵敏度反推。
- [scripts/debug/diag_usbar_term1_internal_compare.jl](scripts/debug/diag_usbar_term1_internal_compare.jl)
  - 对 `term1 = tilde_B0(:plus, -λ, k=0, m_u, m_s, μ_u)` 的内部量做细分诊断。
  - 同时输出当前 hybrid 实现与历史 quadgk 口径的对比，用来排除积分器策略本身造成的系统性偏差。
- [scripts/debug/diag_usbar_term1_counterfactuals.jl](scripts/debug/diag_usbar_term1_counterfactuals.jl)
  - 对 `term1` 做反事实拆分：固定几何输入只换 thermo，或固定 thermo 只换几何输入。
  - 用来区分 `term1` 差异到底主要来自 `E0/p(E0)/λ` 这类运动学量，还是来自 `distribution_value_b0(:plus, ...)` 承载的热分布权重。
- [scripts/debug/diag_usbar_term1_thermo_components.jl](scripts/debug/diag_usbar_term1_thermo_components.jl)
  - 把 thermo 进一步拆成 `T`、`μ_B/3`、`Φ/Φbar` 三块，量化谁在主导 `term1` 变化。
- [scripts/debug/diag_usbar_mu0_T210_T250_terms_compare.jl](scripts/debug/diag_usbar_mu0_T210_T250_terms_compare.jl)
  - 在同一 `μ_B=0` 家族内直接对比 `T=210` 与 `T=250` 的 `term1/term2/term4/B0/denom`，避免 `μ_B=800` 的化学势效应污染主判断。
- [scripts/debug/diag_usbar_term1_pnjl_vs_fd.jl](scripts/debug/diag_usbar_term1_pnjl_vs_fd.jl)
  - 在固定几何和 thermo 点位上，把 `distribution_value_b0(:plus, ...)` 的 PNJL 权重与简化 Fermi-Dirac 权重并排比较。
- [scripts/debug/diag_usbar_term4_internal_compare.jl](scripts/debug/diag_usbar_term4_internal_compare.jl)
  - 对 `term4 = tilde_B0(:minus, -λ, k=0, m_s, m_u, μ_s)` 做和 `term1` 同口径的内部展开，并用 quadgk 交叉检查数值器影响。
- [scripts/debug/diag_usbar_b0_pnjl_vs_fd.jl](scripts/debug/diag_usbar_b0_pnjl_vs_fd.jl)
  - 把 PNJL vs FD 的比较从单个 `term1` 扩展到完整 `B0 = term1 - term2 + term3 - term4` 四子项。
- [scripts/debug/diag_usbar_minus_branch_compare.jl](scripts/debug/diag_usbar_minus_branch_compare.jl)
  - 直接比较 `term2` 和 `term4` 的 `:minus` 分支内部量，确认哪一支真的携带主导差异。
- [scripts/debug/diag_usbar_minus_distribution_values.jl](scripts/debug/diag_usbar_minus_distribution_values.jl)
  - 在 `E0` 处并排输出 `:minus` 分支的 PNJL 与 FD 权重，用来量化 `1 - \bar f(E)` 本身的差异。
- [scripts/debug/diag_usbar_term4_weight_replacement.jl](scripts/debug/diag_usbar_term4_weight_replacement.jl)
  - 保持 strange `term4` 的几何结构不变，只替换 `:minus` 权重为 `PNJL / blend / FD`，直接观察 `term4`、`Re B0` 和 `Re denom` 如何回落。
- [scripts/debug/diag_antiquark_distribution_components.jl](scripts/debug/diag_antiquark_distribution_components.jl)
  - 在 `term4` 的 `E0` 处展开 PNJL `antiquark_distribution` 的分子/分母结构，并与简化 FD 的 `\bar f(E)` 并排比较。
- [scripts/debug/diag_usbar_term4_counterfactuals.jl](scripts/debug/diag_usbar_term4_counterfactuals.jl)
  - 在 `T210/T250 literature_peak` 间做 `geom × thermo` 交叉替换，定量拆分 strange `term4` 的几何放大与热权重放大。
- [scripts/debug/diag_antiquark_distribution_parameter_attribution.jl](scripts/debug/diag_antiquark_distribution_parameter_attribution.jl)
  - 在固定 `E0` 处把 `antiquark_distribution` 拆成 `x=(E+μ)/T` 与 `Φ/\bar Φ` 两类来源，量化谁在主导 `\bar f(E)` 的变化。
- [scripts/debug/diag_usbar_curve_counterfactuals.jl](scripts/debug/diag_usbar_curve_counterfactuals.jl)
  - 在完整 `σ(usbar→usbar)` 曲线上复用原始 `t` 积分与 blocking 口径，只局部替换 strange `sP/K` 通道中的 `term4` 结构。
  - 同时输出基线曲线、`term4` 权重替换曲线和 `term4` 几何替换曲线，用于直接比较峰位与峰高如何响应。

## 3. 结果文件

- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/T210_muB0_usbar_decomp.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/T210_muB0_usbar_decomp.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/T250_muB0_usbar_decomp.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/T250_muB0_usbar_decomp.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/T210_muB0_b0_terms_sensitivity.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/T210_muB0_b0_terms_sensitivity.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_chain_compare_T150_T210.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_chain_compare_T150_T210.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term1_internal_compare.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term1_internal_compare.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term1_counterfactuals.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term1_counterfactuals.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term1_thermo_components.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term1_thermo_components.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_mu0_T210_T250_terms_compare.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_mu0_T210_T250_terms_compare.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term1_pnjl_vs_fd.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term1_pnjl_vs_fd.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_internal_compare.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_internal_compare.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_b0_pnjl_vs_fd.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_b0_pnjl_vs_fd.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_minus_branch_compare.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_minus_branch_compare.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_minus_distribution_values.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_minus_distribution_values.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_weight_replacement.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_weight_replacement.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/antiquark_distribution_components.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/antiquark_distribution_components.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_counterfactuals.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_counterfactuals.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/antiquark_distribution_parameter_attribution.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/antiquark_distribution_parameter_attribution.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_curve_counterfactuals.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_curve_counterfactuals.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_curve_counterfactuals_peak_summary.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_curve_counterfactuals_peak_summary.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_geometry_components.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_geometry_components.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_geometry_component_ratios.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_geometry_component_ratios.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_regular_bin_contributions.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_regular_bin_contributions.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_regular_factor_decomposition.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_regular_factor_decomposition.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_kinetic_geometry_breakdown.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_kinetic_geometry_breakdown.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_pplusp0_breakdown.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_pplusp0_breakdown.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_E0_gap_components.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_E0_gap_components.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_E0_gap_ratios.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_term4_E0_gap_ratios.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_threshold_origin_components.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_threshold_origin_components.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_threshold_origin_swaps.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_threshold_origin_swaps.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_threshold_origin_comparisons.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_threshold_origin_comparisons.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_k0_source_audit.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_k0_source_audit.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_ms_branch_origin.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_ms_branch_origin.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_ms_branch_origin_summary.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_ms_branch_origin_summary.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_baseline_peak_alignment.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_baseline_peak_alignment.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_baseline_peak_alignment_summary.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_baseline_peak_alignment_summary.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_gap_stationarity_gradients.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_gap_stationarity_gradients.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_gap_stationarity_counterfactuals.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_gap_stationarity_counterfactuals.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_peak_formation_chain.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_peak_formation_chain.csv)
- [data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_peak_formation_chain_summary.csv](data/outputs/results/relaxtime/cross_section/usbar_mu0_channel_decomp/usbar_peak_formation_chain_summary.csv)

## 4. 已验证结论

### 4.1 通道分解

- `T=210 MeV`：峰几乎完全由 strange `s` 道 `P` 通道主导。
  - 文献峰附近 `~470 MeV`：`sigma_s_P ≈ 18.75 mb`，`sS`、`t` 道与干涉几乎可忽略。
  - 模型峰附近 `~498-500 MeV`：`sigma_s_P ≈ 49.0 mb`，其余项仍很小。
- `T=250 MeV`：峰由 strange `s` 道 `P+S` 共同抬高。
  - `~582 MeV`：`sigma_s_P ≈ 5.29 mb`，`sigma_s_S ≈ 4.46 mb`。
  - `t` 道和干涉项几乎不贡献峰形。

### 4.2 极点与宽度耦合

- `T=210 MeV`：
  - 阈值 `≈ 467.876 MeV`
  - 有宽度 `K` 极点质量 `≈ 508.038 MeV`
  - 零宽度 `Re[1-4KΠ]=0` 的估计位置 `≈ 515.581 MeV`
- 结论：宽度耦合不是把峰推高的根因；有限宽度反而把零点从 `~515.6 MeV` 拉低到 `~508.0 MeV`。

### 4.3 对称化排查

- 在 `μ_B=0`、`T=210 MeV`、`k=0` 的关键点上，`B0(+k0)` 与 `B0(-k0)` 基本一致。
- `Re[denom(+k0)] = 0`、`Re[denom(-k0)] = 0`、`Re[denom(avg)] = 0` 都落在 `≈ 515.581 MeV`。
- 结论：`num_s_quark=1` 的 `±k0` 对称化不是当前 `μ_B=0` 峰位偏高的主因。

### 4.4 B0 四子项与灵敏度

- `T=210, literature peak ~469.614 MeV`：
  - `term1..4 re ≈ [-0.675, 1.485, 0.0224, -4.316]`
  - `Re B0 ≈ 2.1786` 主要由 `-term4` 和 `-term2` 结构控制，`term3` 很小。
  - `A_sum = A_u + A_s ≈ -23.5771`
  - `prefactor ≈ -2.8315`
  - `Re Π ≈ 1.1302`
  - `Re[1 - 4KΠ] ≈ 0.0742`
- 若仅允许单参数变化并要求文献峰位处 `Re[1 - 4KΠ] = 0`：
  - `A_sum` 需再降低约 `10.1%`
  - `K4567_plus` 需增大约 `8.0%`
  - `Re B0` 需增大约 `38.6%`
- 解释：
  - 要把共振从 `~508-516 MeV` 拉回到文献峰位 `~470 MeV`，最“经济”的单参数改动不是 `Re B0`，而是 `K4567_plus` 或 `A_sum` 量级上的约 `8%-10%` 调整。
  - 但 `Re B0` 仍是控制 `Π_real(k0)` 随能量推进的关键项；只是以单参数口径看，它不是最小改变量方向。


### 4.20 曲线级反事实：把单点结论推进到完整 `σ(usbar→usbar)`

- 为了把前面的 `B0/denom` 结论直接推到可观测峰形，新增了 `diag_usbar_curve_counterfactuals.jl`。
- 这支脚本复用了原始总截面计算口径：
  - `t` 道、`sS` 道、干涉项和 blocking 因子全部保持原实现；
  - 只在 strange `sP/K` 通道局部替换 `Π_us^P` 里的 `term4`。
- 当前比较了三类反事实：
  - `baseline`：原始实现；
  - `weight_blend50 / weight_fd`：保持几何不变，只把 strange `term4` 的 `:minus` 权重从 PNJL 向 FD 混合或完全替换；
  - `geometry_peer`：保持 base thermo 不变，只把 strange `term4` 的几何/奇点结构换成另一条温度曲线对应的几何输入。

- `T210, μ_B=0`：
  - `baseline`：峰位 `≈ 500 MeV`，峰高 `≈ 49.29 mb`；文献峰附近 `470 MeV` 处总截面 `≈ 18.85 mb`。
  - `weight_blend50`：整条曲线显著塌陷，扫描区间内最高只剩 `≈ 19.44 mb`，且文献峰附近只剩 `≈ 1.73 mb`。
  - `weight_fd`：抑制更强，扫描区间内最高只剩 `≈ 7.40 mb`，文献峰附近只剩 `≈ 0.61 mb`。
  - `geometry_peer`：把 `T250` 的几何结构嫁接到 `T210` 后，曲线也整体塌陷，扫描区间内最高 `≈ 8.51 mb`，文献峰附近 `≈ 3.70 mb`。

- `T250, μ_B=0`：
  - `baseline`：峰位 `≈ 572 MeV`，峰高 `≈ 9.81 mb`；文献峰 `582 MeV` 处 `≈ 9.75 mb`。
  - `weight_blend50`：峰高降到 `≈ 8.16 mb`，峰位被推到更高能量端 `≈ 662 MeV`。
  - `weight_fd`：峰高进一步降到 `≈ 5.61 mb`，峰位移到 `≈ 598 MeV`。
  - `geometry_peer`：把 `T210` 的几何结构嫁接到 `T250` 后，峰位左移到 `≈ 554 MeV`，峰高暴涨到 `≈ 28.44 mb`；在文献峰 `582 MeV` 处也被抬到 `≈ 23.02 mb`。

- 这组曲线级结果说明：
  - strange `term4` 的 `:minus` 权重确实不是只在单点 `B0/denom` 上有效，它会直接控制完整 `σ` 曲线的整体强度；
  - 但在 `T210 ↔ T250` 的家族比较里，真正把峰形从“一条较弱曲线”搬运成“一条更强且更靠左的曲线”的，仍然是几何/奇点结构，而不是 thermo 权重本身；
  - 换句话说，`weight lever` 决定“这条曲线能被压低多少”，`geometry lever` 决定“这条曲线更像 T210 还是更像 T250”。
### 4.5 A → G → K 链条与 T150/T210 对比

- `calculate_G_from_A` 是线性的：
  - `G^f = -(N_c / 4π²) * m_f * A_f`
- `K4567_plus` 也是简单线性映射：
  - `K4567_plus = G_fm2 + 0.5 * K_fm5 * G_u`
- 因此 `A_u` 的变化只通过 `G_u` 线性传入 `K4567_plus`，没有隐藏的非线性放大。

- `T=210, μ_B=0`：
  - `A_sum ≈ -23.5771`
  - `G_u ≈ 0.3356`
  - `K4567_plus ≈ 0.20479`
  - anomaly shift `0.5 * K * G_u ≈ 0.00783`，只占 `G_fm2` 的约 `3.97%`
- `T=150, μ_B=800`：
  - `A_sum ≈ -23.9318`
  - `G_u ≈ 0.4009`
  - `K4567_plus ≈ 0.20631`
  - anomaly shift `≈ 0.00935`，只占 `G_fm2` 的约 `4.75%`

- 两点直接对比：
  - `A_sum` 只差约 `1.5%`
  - `G_u` 相差约 `19.4%`
  - 但 `K4567_plus` 只差约 `0.74%`
- 结论：
  - `K4567_plus` 由常数项 `G_fm2` 主导，`A/G` 链对它的调制相对有限。
  - 在 `T=210` 文献峰位处，前述灵敏度分析要求约 `+8%` 的 `K4567_plus` 变化才足以把 `Re[1-4KΠ]` 拉到零；而当前 `T150` 与 `T210` 间实际只观察到 `0.74%` 的变化。
  - 这使“effective coupling 映射本身偏强”作为主因的可能性显著下降。

### 4.6 T150, μ_B=800 与 T210, μ_B=0 的 term1..term4 结构对比

- `T=210, μ_B=0` threshold：
  - `term re ≈ [-0.680, 1.496, 0.0225, -4.330]`
  - `Re B0 ≈ 2.176`
- `T=150, μ_B=800` threshold：
  - `term re ≈ [-1.504, 1.467, 0.0228, -4.222]`
  - `Re B0 ≈ 1.274`

- 共同点：
  - 两点的 `term3` 都很小，可忽略。
  - 两点都由 `term4` 主导主要正贡献（通过 `-term4` 进入 `Re B0`），`term2` 提供主要抵消。
  - 就“哪一支主导结构”而言，它们属于同一类 us P 通道结构。

- 差异点：
  - `T150, μ_B=800` 的 `term1` 负贡献显著更大，导致 `Re B0` 仅约 `1.27`；
  - `T210, μ_B=0` 的 `term1` 负贡献较小，因此 `Re B0` 抬升到约 `2.15-2.18`。
- 结论：
  - 两点属于“同一结构类型、不同数值平衡”的情况，而不是两种完全不同的 B0 机制。
  - `μ_B=0, T=210` 的异常更像是 us P 通道里 `term1/term4/term2` 的平衡把 `Re B0` 整体抬高了约 `40%`，而不是 `A` 或 `K4567_plus` 单独失真。

### 4.7 term1 内部结构与积分器排查

- 新增 `diag_usbar_term1_internal_compare.jl` 后，对 `term1` 直接展开了：
  - `λ_term1`
  - `denominator_term`
  - `E0`
  - `E0 - Emin`
  - `f0 = p(E0) * f^+(E0)`
  - 当前 hybrid 路径里的 `regular`、`logterm`、`pv_integral`
  - 与旧 quadgk 路径的并排结果

- hybrid 与 quadgk 的对比结果：
  - `T=210, μ_B=0, literature peak ~469.614 MeV`：
    - `term1_hybrid ≈ -0.6749195`
    - `term1_quadgk ≈ -0.6749187`
    - 相对差仅 `~1.15e-6`
  - `T=210, μ_B=0, threshold`：
    - `term1_hybrid ≈ -0.6802743`
    - `term1_quadgk ≈ -0.6804300`
    - 相对差 `~2.29e-4`
    - 这是极点正贴端点时的 gap fallback 分支，仍远小于解释峰位失配所需的量级
  - `T=150, μ_B=800, pole_mass ~525.700 MeV`：
    - `term1_hybrid ≈ -1.446066`
    - `term1_quadgk ≈ -1.446064`
    - 相对差 `~1.10e-6`

- 结论一：
  - `term1` 的异常不是由当前 hybrid 求积策略造成的。
  - 至少在当前关注点位上，hybrid 与 quadgk 几乎完全一致；端点 fallback 也只有 `1e-4` 级相对差。

- 结论二：`T210` 的 `term1` 之所以明显没有 `T150` 那么负，主要来自内部量本身，而不是数值器误差。
  - `T=210, μ_B=0, literature peak`：
    - `E0 - Emin ≈ 7.52e-3 fm^-1`
    - `f0 ≈ 2.889e-2`
    - `regular ≈ 6.329e-1`
    - `logterm ≈ 1.702e-1`
    - `regular + logterm ≈ 8.031e-1`
  - `T=150, μ_B=800, pole_mass`：
    - `E0 - Emin ≈ 4.09e-2 fm^-1`
    - `f0 ≈ 1.589e-1`
    - `regular ≈ 1.265`
    - `logterm ≈ 6.616e-1`
    - `regular + logterm ≈ 1.926`

- 对比解释：
  - 两点的 `|λ_term1|` 只相差约 `12%`，不足以单独解释 `term1` 从 `-1.446` 抬到 `-0.675`。
  - 真正显著的差异在于 `T210` 的奇点 `E0` 更贴近下端点 `Emin`，导致：
    - `p(E0)` 更小，`f0` 更小；
    - 主值积分里的 `regular + logterm` 组合也只有 `T150` 的约 `42%`。
  - 因此 `T210` 的 `term1` 偏小，当前更像是 `E0` 的区间位置和对应分布权重带来的物理输入差异，而不是“积分器把 term1 算浅了”。

### 4.8 term1 反事实拆分：几何输入 vs 热分布输入

- 为了继续区分 `term1` 的来源，新增了 `diag_usbar_term1_counterfactuals.jl`：
  - `geom210_lit`：固定 `T=210, μ_B=0, literature peak` 的几何输入（`m_u, m_s, λ, E0` 等）
  - `geom150_pole`：固定 `T=150, μ_B=800, pole_mass` 的几何输入
  - 然后在两套几何上分别代入：
    - 各自原始 thermo
    - 对方的 thermo

- 四组结果如下：
  - `geom210_lit + thermo210_lit`：
    - `term1 ≈ -0.6749`
    - `p0 ≈ 7.23e-2`
    - `dist(E0) ≈ 0.3996`
    - `regular + logterm ≈ 0.8031`
  - `geom210_lit + thermo150_pole`：
    - `term1 ≈ -1.6253`
    - `p0` 不变
    - `dist(E0) ≈ 0.8653`
    - `regular + logterm ≈ 1.9340`
  - `geom150_pole + thermo210_lit`：
    - `term1 ≈ -0.5850`
    - `p0 ≈ 1.884e-1`
    - `dist(E0) ≈ 0.3715`
    - `regular + logterm ≈ 0.7792`
  - `geom150_pole + thermo150_pole`：
    - `term1 ≈ -1.4461`
    - `p0 ≈ 1.884e-1`
    - `dist(E0) ≈ 0.8437`
    - `regular + logterm ≈ 1.9262`

- 关键量化结论：
  - 固定 `geom210`，只把 thermo 从 `210/0` 换成 `150/800`：
    - `term1` 从 `-0.6749` 变成 `-1.6253`
    - 变化量 `≈ -0.9504`
  - 固定 `geom150`，只把 thermo 从 `150/800` 换成 `210/0`：
    - `term1` 从 `-1.4461` 变成 `-0.5850`
    - 变化量 `≈ +0.8611`
  - 相比之下，固定 thermo 只换几何：
    - 在 `thermo210` 下，几何变化只带来 `≈ +0.0899`
    - 在 `thermo150` 下，几何变化只带来 `≈ +0.1792`

- 解释：
  - `term1` 的主导差异来自 thermo，而不是 geometry。
  - 换句话说，真正把 `term1` 从 `~ -1.45` 压到 `~ -0.67` 的，不是 `E0` 更贴近端点这件事本身，而是这套点位上的 `distribution_value_b0(:plus, ...)` 权重显著更小。
  - 几何量并非完全不重要，但它更像次级修正：在当前两组对比点上，量级只有热分布效应的大约 `10%` 到 `20%`。

- 这也修正了上一轮的表述：
  - 上一轮说“更像 `E0` 贴端点后的权重/几何位置变化”。
  - 现在更准确的说法应是：`E0` 贴端点只是背景条件，真正主导 `term1` 强弱的是 thermo 输入下的分布权重大小；几何位置变化只是在这个基础上的次级调制。

### 4.9 thermo 再拆分：T、μ、Φ/Φbar 谁主导 term1

- 在确认“thermo 比 geometry 更主导”后，又新增了 `diag_usbar_term1_thermo_components.jl`，把 thermo 拆成：
  - 仅换 `T`
  - 仅换 `μ_B/3`
  - 仅换 `Φ/Φbar`
  - 以及它们的组合

- 对 `geom210_lit` 的结果：
  - `base210`: `term1 ≈ -0.6749`
  - `T_only_150`: `term1 ≈ -0.4961`，变化 `+0.1789`
  - `mu_only_800`: `term1 ≈ -1.6233`，变化 `-0.9484`
  - `Phi_only_150`: `term1 ≈ -0.5811`，变化 `+0.0939`
  - `full150`: `term1 ≈ -1.6253`，总变化 `-0.9504`

- 对 `geom150_pole` 的结果：
  - `base210`: `term1 ≈ -0.5850`
  - `T_only_150`: `term1 ≈ -0.4132`，变化 `+0.1718`
  - `mu_only_800`: `term1 ≈ -1.4562`，变化 `-0.8712`
  - `Phi_only_150`: `term1 ≈ -0.4944`，变化 `+0.0906`
  - `full150`: `term1 ≈ -1.4461`，总变化 `-0.8611`

- 这组结果说明：
  - `μ_B/3` 的化学势位移几乎单独重现了全部 thermo shift。
  - 单独降低 `T` 并不会把 `term1` 拉得更负，反而会让它更不负。
  - 单独把 `Φ/Φbar` 换成 `T150, μ_B=800` 的值，也只带来较小且同方向的“变浅”效应。
  - 因而 `full150` 比 `base210` 更负，几乎完全是由 `μ` 项主导；`T` 和 `Φ/Φbar` 只是在其上做小幅反向修正。

- 当前阶段最重要的解释更新：
  - 之前用 `T150, μ_B=800` 作为“结构正常点”去对比 `T210, μ_B=0`，在 `term1` 层面会被很强的化学势效应主导。
  - 因此，这组对比虽然能证明 `K4567_plus` 不是主因、也能证明积分器不是主因，但它已经不适合作为继续定位 `μ_B=0` 异常根因的主对照组。
  - 若要继续定位 `μ_B=0` 的“异常”，下一步更应该在同一 `μ_B=0` 家族内比较，比如：
    - `T=210` vs `T=250`；
    - 或固定 `μ_B=0` 只扫 `T`；
    - 而不是继续把 `μ_B=800` 的点当作主参考。

### 4.10 同 μ_B=0 家族对照：T210 vs T250 的 term1/term2/term4

- 新增 `diag_usbar_mu0_T210_T250_terms_compare.jl` 后，直接在 `μ_B=0` 条件下对比了 `T=210` 和 `T=250` 的 `term` 结构。

- 结果摘要：
  - `literature_peak` 点：
    - `T210 @ 469.614 MeV`：
      - `term1 ≈ -0.6749`
      - `term2 ≈ 1.4850`
      - `term4 ≈ -4.3161`
      - `Re B0 ≈ 2.1786`
      - `Re denom ≈ 0.0742`
    - `T250 @ 582.275 MeV`：
      - `term1 ≈ -0.2669`
      - `term2 ≈ 0.8883`
      - `term4 ≈ -2.8302`
      - `Re B0 ≈ 1.7240`
      - `Re denom ≈ 0.0946`
  - `pole_mass` 点：
    - `T210 @ 508.204 MeV`：`Re B0 ≈ 2.1496`
    - `T250 @ 566.498 MeV`：`Re B0 ≈ 1.7613`

- 对 `T250 - T210` 的差分解释：
  - 在 `literature/model/pole` 这些相关点位上：
    - `term1` 确实在 `T250` 更“不负”，即数值更接近 0；
    - `term2` 也显著减小，意味着抵消项变弱；
    - 但最强的变化来自 `term4` 的绝对值显著下降。
  - 以 `literature_peak` 为例：
    - `Δ term1 ≈ +0.408`
    - `Δ term2 ≈ -0.597`
    - `Δ term4 ≈ +1.486`
    - 最终 `Δ Re B0 ≈ -0.455`

- 结论：
  - 在同一 `μ_B=0` 家族内，`T210` 的 `Re B0` 较高，不是因为它的 `term1` 比 `T250` 更异常地偏浅。
  - 相反，若只看 `term1`，`T250` 反而更浅。
  - 真正把 `T210` 的 `Re B0` 抬得更高的主因，是 `term4` 在 `T210` 的正贡献更强；`term2` 和 `term1` 的变化只是在其上做部分抵消。
  - 这说明当前“只盯 `term1`”已经不够，下一轮需要把注意力扩展到 `term4 = tilde_B0(:minus, -λ, m_s, m_u, μ_s)` 的 thermo/几何来源。

### 4.11 PNJL 权重 vs 简化 Fermi-Dirac 权重

- 新增 `diag_usbar_term1_pnjl_vs_fd.jl` 后，在固定点位上比较了：
  - `distribution_value_b0(:plus, ...)` 当前使用的 PNJL 权重
  - 去掉 `Φ/Φbar` 修正后的简化 FD 权重

- 结果摘要：
  - `T210, literature_peak`：
    - `PNJL dist(E0) ≈ 0.3996`
    - `FD   dist(E0) ≈ 0.4182`
    - `term1_PNJL ≈ -0.6749`
    - `term1_FD   ≈ -0.7780`
    - `Δ term1 (FD - PNJL) ≈ -0.1030`
  - `T210, model_peak`：
    - `Δ term1 ≈ -0.1060`
  - `T250, literature_peak`：
    - `term1_PNJL ≈ -0.2669`
    - `term1_FD   ≈ -0.3338`
    - `Δ term1 ≈ -0.0668`

- 结论：
  - PNJL 修正确实会把 `term1` 再压浅一截，即比简化 FD 更不负。
  - 但这个效应的量级只有约 `0.07` 到 `0.11`，显著小于：
    - `T210` 与 `T250` 在 `Re B0` 上的差异量级；
    - 更远小于此前跨 `μ_B` 对照里由化学势带来的大幅变化。
  - 因此，问题不太像“只要把 PNJL 改成 FD 就能解释 `μ_B=0` 的 usbar 失配”。
  - 更准确的说法是：PNJL 修正提供了一个次级但方向一致的浅化效应；它值得保留在可疑列表里，但目前不足以单独解释主异常。

### 4.12 term4 内部结构与积分器排查

- 新增 `diag_usbar_term4_internal_compare.jl` 后，对 `term4` 做了与 `term1` 相同的内部展开与 `hybrid/quadgk` 对照。

- 数值器结论先行：
  - `term4` 的 `hybrid` 与 `quadgk` 也高度一致。
  - 代表性点位：
    - `T210, literature_peak`：相对差约 `2.67e-6`
    - `T250, literature_peak`：相对差约 `1.94e-5`
    - 端点 fallback 的 threshold 点也仅约 `1.4e-4` 到 `1.5e-4`
  - 这说明 `term4` 的大幅差异同样不是数值积分器造成的。

- 内部量对比：
  - `T210, literature_peak`：
    - `term4 ≈ -4.3161`
    - `E0 - Emin ≈ 1.29e-3 fm^-1`
    - `dist(E0) ≈ 0.9035`
    - `regular ≈ 4.6690`
    - `logterm ≈ 0.4668`
    - `regular + logterm ≈ 5.1359`
  - `T250, literature_peak`：
    - `term4 ≈ -2.8302`
    - `E0 - Emin ≈ 3.16e-1 fm^-1`
    - `dist(E0) ≈ 0.8426`
    - `regular ≈ 2.7765`
    - `logterm ≈ 1.3991`
    - `regular + logterm ≈ 4.1757`

- 解释：
  - `T210` 的 `term4` 更强，既不是数值器假象，也不只是 `dist(E0)` 的小幅差异。
  - 真正更关键的是 strange 支这边：
    - `E0` 更贴近下端点；
    - 同时 `dist(E0)` 仍维持在很高水平；
    - 两者一起把 `regular + logterm` 推到了比 `T250` 更大的量级。
  - 这与上一节的同家族结论一致：当前 `μ_B=0` 家族里，主导 `Re B0` 差异的确实是 `term4`。

### 4.13 完整 B0 的 PNJL vs FD 对照

- 新增 `diag_usbar_b0_pnjl_vs_fd.jl` 后，把 PNJL 与简化 FD 的比较扩展到了四子项整体。

- 结果非常关键：
  - `T210, literature_peak`：
    - `PNJL`: `term1 ≈ -0.6749`, `term2 ≈ 1.4850`, `term4 ≈ -4.3161`, `Re B0 ≈ 2.1786`
    - `FD`:   `term1 ≈ -0.7780`, `term2 ≈ 0.2945`, `term4 ≈ -0.3971`, `Re B0 ≈ -0.6434`
    - `Δ(FD-PNJL)`: `term1 ≈ -0.1030`, `term2 ≈ -1.1905`, `term4 ≈ +3.9189`, `Re B0 ≈ -2.8220`
  - `T250, literature_peak`：
    - `PNJL`: `Re B0 ≈ 1.7240`
    - `FD`:   `Re B0 ≈ -0.1581`
    - `Δ(FD-PNJL) ≈ -1.8821`

- 这组结果比单独看 `term1` 更有信息量：
  - PNJL 修正对 `term1` 只是次级影响；
  - 但对 `term2` 和尤其 `term4`，影响是主导量级；
  - 一旦把 PNJL 换成简化 FD，整个 `Re B0` 甚至会翻到负值。

- 结论：
  - 如果只看 `term1`，会低估 PNJL 修正在完整 `B0` 里的作用。
  - 对当前 `usbar_to_usbar` 的 strange `s` 道来说，真正高度依赖 PNJL 口径的是 `:minus` 分支主导的项，尤其 `term4`，其次是 `term2`。
  - 因此，现在更值得怀疑的不是“PNJL 修正有没有影响”，而是：
    - 当前 `distribution_value_b0(:minus, ...) = 1 - antiquark_distribution(...)` 这一路的物理口径是否与目标文献/旧实现一致；
    - 或者 strange `:minus` 分支在当前 PNJL 参数下确实把 `B0` 推得过强。

### 4.14 term2 / term4 的 :minus 分支内部对照

- 新增 `diag_usbar_minus_branch_compare.jl` 后，直接比较了 `term2` 和 `term4` 在 `T=210/T=250, μ_B=0` 下的内部量。

- 首先，一个很关键的结构事实出现了：
  - 在当前关注的 `literature/model/pole` 点位上，`term2` 没有落到 `k=0` 奇点分支，因此没有对应的 `E0` 残数入口。
  - `term4` 则明确落在奇点分支上，且内部量可稳定展开。
  - 这说明即便 `term2` 对 `Re B0` 有数值贡献，它和 `term4` 并不是同一类“共振增强入口”。

- `term4` 的同家族对比：
  - `literature_peak`：
    - `T210`: `dist(E0) ≈ 0.9035`, `regular + logterm ≈ 5.1359`, `term4 ≈ -4.3161`
    - `T250`: `dist(E0) ≈ 0.8426`, `regular + logterm ≈ 4.1757`, `term4 ≈ -2.8302`
    - 差分：
      - `Δ dist(E0) ≈ -0.0609`
      - `Δ (regular + logterm) ≈ -0.9602`
      - `Δ term4 ≈ +1.4859`
  - `pole_mass`：
    - `T210`: `regular + logterm ≈ 5.0695`
    - `T250`: `regular + logterm ≈ 4.2233`
    - `Δ term4 ≈ +0.9946`

- 解释：
  - `T210` 的 `term4` 更强，不只是因为 `dist(E0)` 稍大；
  - 更关键的是它在 strange `:minus` 支上的主值积分主体 `regular + logterm` 本身也更大。
  - 所以 `term4` 的强化来自两部分共同作用：
    - `1 - \bar f(E0)` 权重更大；
    - strange 支对应的奇点/几何结构也更能放大主值积分。

### 4.15 :minus 分布口径核对

- 代码与参考口径的一致性：
  - [src/relaxtime/OneLoopIntegrals.jl](src/relaxtime/OneLoopIntegrals.jl) 已明确写出：
    - `:plus => f^+(+E, μ)`
    - `:minus => f^+(-E, μ) = 1 - f^-(+E, μ)`
  - 历史 quadgk 对照脚本 [docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md](docs/dev/archived/2026-01-22_RelaxTime_Fortran_Comparison_Scripts.md) 也沿用同一 `distribution_value_b0(sign_flag, ...)`。
  - Mathematica 参考推导 [docs/reference/mathematica/B0.md](docs/reference/mathematica/B0.md) 中，B0 的第二、四项确实对应 `f^-` 型分布；因此把 `f^+(-E, μ)` 重写为 `1 - f^-(E, μ)` 在形式上是自洽的。

- 但形式自洽不等于量级无问题。新增 `diag_usbar_minus_distribution_values.jl` 后，在 `E0` 处直接比较了 `:minus` 权重：
  - `T210, literature_peak, term4`：
    - `PNJL minus ≈ 0.9035`
    - `FD minus ≈ 0.1294`
    - 差值 `≈ -0.7740`
    - 比值 `FD/PNJL ≈ 0.143`
  - `T250, literature_peak, term4`：
    - `PNJL minus ≈ 0.8426`
    - `FD minus ≈ 0.1826`
    - 差值 `≈ -0.6599`
    - 比值 `≈ 0.217`

- 这意味着：
  - 对 `term4` 而言，当前 PNJL 的 `1 - \bar f(E)` 权重比简化 FD 大了约 `4.6` 到 `7.0` 倍。
  - 因此，当前问题已经不再是“`:minus` 公式写法是否正确”这么粗的层面；
  - 真正的问题更像是：在当前 `T/μ/Φ/\bar Φ` 参数下，这个自洽公式给出的 strange `:minus` 权重是否物理上过强，或者是否与文献作图时采用的分布口径不同。

### 4.16 strange term4 权重替换实验

- 为了把“可疑”进一步推进到“因果性证据”，新增了 `diag_usbar_term4_weight_replacement.jl`：
  - strange `term4` 的几何输入、奇点位置、`regular + logterm` 都保持不变；
  - 只把 `:minus` 权重从当前 PNJL 权重连续替换成 `blend_25 / blend_50 / blend_75 / FD`；
  - 然后重算 `term4`、完整 `Re B0` 与 `Re denom`。

- `T210, literature_peak ~469.614 MeV`：
  - `PNJL`：`term4 ≈ -4.3161`，`Re B0 ≈ 2.1786`，`Re denom ≈ 0.0742`
  - `blend_25`：`term4 ≈ -3.3363`，`Re B0 ≈ 1.1989`，`Re denom ≈ 0.1605`
  - `blend_50`：`term4 ≈ -2.3566`，`Re B0 ≈ 0.2191`，`Re denom ≈ 0.2469`
  - `blend_75`：`term4 ≈ -1.3769`，`Re B0 ≈ -0.7606`，`Re denom ≈ 0.3332`
  - `FD`：`term4 ≈ -0.3971`，`Re B0 ≈ -1.7403`，`Re denom ≈ 0.4195`

- `T210, model_peak ~497.448 MeV`：
  - `PNJL`：`term4 ≈ -4.0384`，`Re B0 ≈ 2.1606`，`Re denom ≈ 0.0293`
  - `FD`：`term4 ≈ -0.3599`，`Re B0 ≈ -1.5179`，`Re denom ≈ 0.4326`

- `T250, literature_peak ~582.275 MeV`：
  - `PNJL`：`term4 ≈ -2.8302`，`Re B0 ≈ 1.7240`，`Re denom ≈ 0.0946`
  - `blend_50`：`term4 ≈ -1.5866`，`Re B0 ≈ 0.4804`，`Re denom ≈ 0.3366`
  - `FD`：`term4 ≈ -0.3430`，`Re B0 ≈ -0.7632`，`Re denom ≈ 0.5787`

- 结论：
  - 只改 strange `term4` 的 `:minus` 权重，不改几何结构，就足以让 `Re B0` 发生符号级别翻转，并把 `Re denom` 从近共振小正值推回到显著更大的正值。
  - 这比“PNJL vs FD 的相关性”更强，已经是直接的因果性数值证据：当前异常的主导入口就是 strange `term4` 的 `:minus` 权重。

### 4.17 antiquark_distribution 分子分母展开

- 为了继续解释“为什么 `1 - \bar f(E)` 会这么大”，新增了 `diag_antiquark_distribution_components.jl`。
- 它在 `term4` 的 `E0` 处直接输出：
  - `x = (E + μ) / T`
  - `exp(-x), exp(-2x), exp(-3x)`
  - PNJL 分子 `\bar Φ e^{-x} + 2Φ e^{-2x} + e^{-3x}`
  - PNJL 分母 `1 + 3\bar Φ e^{-x} + 3Φ e^{-2x} + e^{-3x}`
  - 最终 `antiquark_distribution` 与 `1 - antiquark_distribution`
  - 对照的简化 FD `\bar f(E)` 与 `1 - \bar f(E)`

- `T210, literature_peak`：
  - `x ≈ 1.9059`
  - `exp(-x) ≈ 0.1487`
  - `numerator ≈ 0.1291`
  - `denominator ≈ 1.3375`
  - `anti_PNJL ≈ 0.0965`
  - `anti_FD ≈ 0.8706`
  - `minus_PNJL ≈ 0.9035`
  - `minus_FD ≈ 0.1294`

- `T250, literature_peak`：
  - `x ≈ 1.4985`
  - `exp(-x) ≈ 0.2235`
  - `numerator ≈ 0.2576`
  - `denominator ≈ 1.6364`
  - `anti_PNJL ≈ 0.1574`
  - `anti_FD ≈ 0.8174`
  - `minus_PNJL ≈ 0.8426`
  - `minus_FD ≈ 0.1826`

- 解释：
  - 当前问题不是 PNJL 让 `1 - \bar f(E)` “略微偏大”，而是 PNJL 让 `\bar f(E)` 本身在 relevant `E0` 处变得极小。
  - 在这些点位上，FD 反夸克分布约为 `0.82` 到 `0.87`，但 PNJL 只有 `0.10` 到 `0.16`；因此 `1 - \bar f(E)` 被直接抬到 `0.84` 到 `0.90`。
  - 就当前数值看，主要现象可以概括为：PNJL 分布把 strange `term4` 的 `:minus` 权重压向接近 1，而这正是 `Re B0` 被显著抬高的直接来源。

### 4.18 strange term4：几何放大 vs 热权重放大

- 为了回答“当前 `term4` 的增强到底主要来自几何还是热权重”，新增了 `diag_usbar_term4_counterfactuals.jl`。
- 它固定 `T210/T250 literature_peak` 两个代表点，做四组 `geom × thermo` 组合：
  - `geom210 + thermo210`
  - `geom210 + thermo250`
  - `geom250 + thermo210`
  - `geom250 + thermo250`
- 其中：
  - `geometry` 包含 `k0, λ, m_u, m_s, E0` 以及奇点位置相关结构；
  - `thermo` 包含 `T, μ, Φ, \bar Φ`，即进入 `distribution_value_b0(:minus, ...)` 的权重口径。

- 结果如下：
  - `geom210 + thermo210`：
    - `term4 ≈ -4.3161`
    - `Re B0 ≈ 2.1786`
    - `Re denom ≈ 0.0742`
  - `geom210 + thermo250`：
    - `term4 ≈ -4.1442`
    - `Re B0 ≈ 2.0067`
    - `Re denom ≈ 0.0893`
    - 相对 `geom210 + thermo210`，`Δ term4 ≈ +0.1719`
  - `geom250 + thermo210`：
    - `term4 ≈ -2.9443`
    - `Re B0 ≈ 1.8382`
    - `Re denom ≈ 0.0724`
    - 相对 `geom210 + thermo210`，`Δ term4 ≈ +1.3717`
  - `geom250 + thermo250`：
    - `term4 ≈ -2.8302`
    - `Re B0 ≈ 1.7240`
    - `Re denom ≈ 0.0946`
    - 相对 `geom210 + thermo210`，`Δ term4 ≈ +1.4859`

- 这组数值非常关键：
  - `T250 - T210` 的完整 `term4` gap 约为 `+1.4859`；
  - 若只换 thermo，不换 geometry，只能解释其中约 `+0.1719`；
  - 若只换 geometry，不换 thermo，已经能解释其中约 `+1.3717`。

- 结论：
  - 在 `T210 ↔ T250` 这条 `μ_B=0` 家族比较线上，strange `term4` 的差异主要由几何/奇点结构决定，热权重只是次级修正。
  - 按当前数值量级看，几何项解释了约九成以上的 `term4` gap，而 thermo 权重大约只解释一成左右。
  - 这和前一轮“只换 strange `term4` 的 :minus 权重就能显著改变 `B0/denom`”并不矛盾：
    - PNJL 权重对“绝对峰强度/接近共振程度”仍是直接因果杠杆；
    - 但在 `T210` 与 `T250` 的相互比较里，哪一个更强，主要还是由几何放大决定。

### 4.19 antiquark_distribution：x vs Φ/Φbar 归因

- 为了回答“PNJL 反夸克分布为什么这么小，以及到底是 `x=(E+μ)/T` 还是 `Φ/\bar Φ` 在主导”，新增了 `diag_antiquark_distribution_parameter_attribution.jl`。
- 在固定 `E0` 的前提下，它分别替换：
  - 只换 `x`
  - 只换 `Φ/\bar Φ`
  - `x` 与 `Φ/\bar Φ` 一起换
  - 同时保留 `FD same x` 作为参照

- `geom210_lit`（`E0` 固定为 `T210 literature_peak` 对应点）结果：
  - `base210`：`anti ≈ 0.0965`，`minus ≈ 0.9035`
  - `x <- 250`：`anti ≈ 0.1301`，变化 `+0.0335`
  - `Φ <- 250`：`anti ≈ 0.1079`，变化 `+0.0113`
  - `full <- 250`：`anti ≈ 0.1434`，变化 `+0.0468`
  - `FD same x`：`anti ≈ 0.1294`，`minus ≈ 0.8706`

- `geom250_lit`（`E0` 固定为 `T250 literature_peak` 对应点）结果：
  - `base250`：`anti ≈ 0.1574`，`minus ≈ 0.8426`
  - `x <- 210`：`anti ≈ 0.1210`，变化 `-0.0365`
  - `Φ <- 210`：`anti ≈ 0.1436`，变化 `-0.0138`
  - `full <- 210`：`anti ≈ 0.1088`，变化 `-0.0486`
  - `FD same x`：`anti ≈ 0.1826`，`minus ≈ 0.8174`

- 归因解释：
  - 在 `T210 ↔ T250` 之间，`\bar f(E)` 的变化主要由 `x=(E+μ)/T` 驱动；
  - `Φ/\bar Φ` 的变化方向一致，但量级明显更小；
  - 就 `geom210` 这组数值看，`x` 单独已经解释了大约七成的 `anti` 回升，`Φ/\bar Φ` 只贡献约四分之一量级，其余是非线性交叉项。

- 进一步看“为何 PNJL 比 FD 更小”：
  - `geom210`：
    - `FD same x anti ≈ 0.1294`
    - `PNJL anti ≈ 0.0965`
    - PNJL 分子 `≈ 0.1291`，小于 `exp(-x) ≈ 0.1487`
    - PNJL 分母 `≈ 1.3375`，大于 `1 + exp(-x) ≈ 1.1487`
    - 即：`T210` 处既有分子压低，也有分母抬高，两者都在压制 `\bar f(E)`
  - `geom250`：
    - `FD same x anti ≈ 0.1826`
    - `PNJL anti ≈ 0.1574`
    - PNJL 分子 `≈ 0.2576`，反而大于 `exp(-x) ≈ 0.2235`
    - 但 PNJL 分母 `≈ 1.6364`，远大于 `1 + exp(-x) ≈ 1.2235`
    - 即：`T250` 处主要是分母抬高在压制 `\bar f(E)`

- 结论：
  - 若问“PNJL 反夸克分布为何比 FD 小”，当前更准确的答案不是单一句“因为 `Φ/\bar Φ` 小”，而是：
    - `T210/T250` 之间的相对变化主要由 `x` 驱动；
    - PNJL 相对 FD 的额外压低，则更多体现为 Polyakov 修正把分母系统性抬高，`T210` 时分子压低也同时参与，`T250` 时则主要是分母效应占优。

### 4.21 geometry lever 再拆：不是 logterm 主导，而是 regular 主值积分主导

- 为了继续回答“geometry 到底通过什么组成块放大 `term4`”，新增了 `diag_usbar_term4_geometry_component_attribution.jl`。
- 这支脚本在前一轮 `geom × thermo` 的基础上继续拆开：
  - `1 / |λ_term4|`
  - `E0 - Emin`
  - `p(E0)`
  - `dist(E0)`
  - `f0 = p(E0) * dist(E0)`
  - `log_kernel`
  - `regular`
  - `logterm`
  - `regular_share = regular / (regular + logterm)`

- 在固定 thermo 的几何对比下，结果几乎相同：
  - `geometry_effect_under_thermo210`：
    - `|term4_210| / |term4_250| ≈ 1.4659`
    - `inv|λ| ratio ≈ 1.2399`
    - `(regular + logterm) ratio ≈ 1.1823`
  - `geometry_effect_under_thermo250`：
    - `|term4_210| / |term4_250| ≈ 1.4643`
    - `inv|λ| ratio ≈ 1.2399`
    - `(regular + logterm) ratio ≈ 1.1810`

- 更关键的是各组成块的方向：
  - `regular ratio ≈ 1.62 - 1.63`
  - `logterm ratio ≈ 0.315 - 0.316`
  - `p0 ratio ≈ 0.0689`
  - `dist(E0) ratio ≈ 1.014 - 1.017`
  - `f0 ratio ≈ 0.070`
  - `log_kernel ratio ≈ 4.51`
  - `E0 - Emin ratio ≈ 4.08e-3`

- 解释：
  - `T210` 的奇点确实比 `T250` 更贴近下端点，`E0 - Emin` 只有后者的大约 `0.4%`；
  - 这会让 `log_kernel` 变大约 `4.5` 倍，但由于 `p(E0)` 极小，`f0` 反而只剩后者的约 `7%`，因此 `logterm` 并没有被放大，反而只有 `T250` 的约三分之一；
  - 真正被放大的是 `regular` 主值积分本体，它在 `T210` 下比 `T250` 大约 `1.6` 倍；
  - 再叠加 `1 / |λ|` 的约 `1.24` 倍增益，才得到最终 `|term4|` 的约 `1.46` 倍增强。

- 因而当前可以把 geometry lever 说得更精确：
  - 它不是“奇点更贴边，所以 logterm 更大”；
  - 而是“奇点/阈值结构改变了 principal-value regular 部分的整体形状，同时较小的 `|λ|` 又进一步放大了它”。

### 4.22 regular 能区分箱：T210 的增强来自宽阈值区，不是单个极窄尖点

- 为了继续回答“`regular` 到底是哪一段能区在抬高 `T210`”，新增了 `diag_usbar_term4_regular_bin_contributions.jl`。
- 这支脚本把 `regular_integrand(E)` 在 `ΔE = E - Emin` 上分箱积分：
  - `[0, 0.005]`
  - `[0.005, 0.02]`
  - `[0.02, 0.05]`
  - `[0.05, 0.1]`
  - `[0.1, 0.2]`
  - `[0.2, 0.4]`
  - `[0.4, 0.8]`
  - `[0.8, 1.6]`
  - 余下尾部

- `T210, literature_peak`：
  - `regular ≈ 4.6690`
  - `logterm ≈ 0.4668`
  - `ΔE ≤ 0.2` 已累计约 `28.1%` 的 regular
  - `ΔE ≤ 0.4` 已累计约 `42.5%`
  - `ΔE ≤ 0.8` 已累计约 `64.4%`
  - 单个最大分箱是 `[0.8, 1.6]`，贡献约 `34.2%`

- `T250, literature_peak`：
  - `regular ≈ 2.7765`
  - `logterm ≈ 1.3991`
  - `ΔE ≤ 0.2` 只累计约 `14.3%`
  - `ΔE ≤ 0.4` 只累计约 `26.3%`
  - `ΔE ≤ 0.8` 累计约 `48.0%`
  - 单个最大分箱同样是 `[0.8, 1.6]`，贡献约 `39.8%`

- 更直接地比较各分箱的绝对 regular 贡献：
  - `[0, 0.005]`：`T210 / T250 ≈ 8.85`
  - `[0.005, 0.02]`：`≈ 5.34`
  - `[0.02, 0.05]`：`≈ 3.83`
  - `[0.05, 0.1]`：`≈ 3.02`
  - `[0.1, 0.2]`：`≈ 2.46`
  - `[0.2, 0.4]`：`≈ 2.03`
  - `[0.4, 0.8]`：`≈ 1.69`
  - `[0.8, 1.6]`：`≈ 1.44`

- 解释：
  - `T210` 的 regular 放大确实从极近阈值区开始最强，但它不是只靠 `E0` 邻域的一个极窄尖点；
  - 从 `ΔE \lesssim 0.8 fm^-1` 的整段阈值上方能区，`T210` 都系统性地比 `T250` 更大；
  - 到最远尾部时，这个优势才明显消失。

- 这使当前结论进一步收敛为：
  - strange `term4` 的 geometry 放大，本质上是 `T210` 把更大的一整段阈值邻近能区都推到了 higher-regular-contribution 分支，而不只是把 singular log piece 局部拉尖。

### 4.23 regular 因子分解：T210 与 T250 的 family gap 几乎全是动量项

- 为了继续回答“在这些 regular 分箱里，究竟是 `p(E)` 还是 `1-\bar f(E)` 在主导差异”，新增了 `diag_usbar_term4_regular_factor_decomposition.jl`。
- 这支脚本使用了 regular integrand 的精确分解：

$$
\frac{p(E)w(E)-p_0 w_0}{E-E_0}
=
\frac{(p(E)-p_0)w(E)}{E-E_0}
+
\frac{p_0\,(w(E)-w_0)}{E-E_0},
$$

其中：
- `kinetic term = ((p(E)-p0) * w(E)) / (E-E0)`
- `weight term = p0 * (w(E)-w0) / (E-E0)`
- `w(E) = distribution_value_b0(:minus, E, ...)`

- 总体分解结果：
  - `T210, literature_peak`：
    - `regular ≈ 4.6690`
    - `kinetic ≈ 4.6614`
    - `weight ≈ 7.65e-3`
    - 即：weight 只占 regular 的约 `0.16%`
  - `T250, literature_peak`：
    - `regular ≈ 2.7765`
    - `kinetic ≈ 2.5923`
    - `weight ≈ 1.842e-1`
    - 即：weight 大约占 regular 的 `6.6%`

- 分箱层面看，`T210` 更极端：
  - 在 `[0, 0.005]` 到 `[0.8, 1.6]` 这些分箱里，`kinetic/full` 基本都在 `99.8%` 左右；
  - 对应的 `weight/full` 只在 `2.8e-4` 到 `2.0e-3` 量级。

- `T250` 下，weight 项虽然比 `T210` 更可见，但仍是次级：
  - 近阈值分箱 `weight/full` 大约在 `5%` 到 `7%`；
  - 中高能区也大致维持在 `5.8%` 到 `7.2%`；
  - `kinetic/full` 仍然稳定在 `92.8%` 到 `95%` 左右。

- 解释：
  - 这一步把前面的结论再压缩了一层：
    - PNJL `:minus` 权重当然仍然是决定绝对峰强能否被压低的重要杠杆；
    - 但在 `T210` 和 `T250` 的 family gap 里，regular 本身为什么会差这么多，答案几乎完全是运动学项 `p(E)-p0` 在主导，而不是 `w(E)-w0` 在主导。
  - 换句话说：
    - `weight lever` 更像“绝对强度旋钮”；
    - `regular family lever` 则几乎就是“动量/阈值运动学旋钮”。

### 4.24 kinetic 几何再拆：主导不是 `E+E0`，而是 `1 / (p(E)+p0)`

- 为了直接回答“为什么 `T210` 的动量项 `(p(E)-p0)` 会在整段阈值上方能区系统性更大，而 `T250` 没有”，新增了 `diag_usbar_term4_kinetic_geometry_breakdown.jl`。
- 关键恒等式是：

$$
\frac{p(E)-p_0}{E-E_0} = \frac{E+E_0}{p(E)+p_0}.
$$

- 也就是说，kinetic 项其实是：

$$
\frac{(p(E)-p_0)w(E)}{E-E_0} = w(E)\,\frac{E+E_0}{p(E)+p_0},
$$

- 因而现在可以把它再拆成两部分：
  - 分子项：`E + E0`
  - 分母项：`1 / (p(E) + p0)`

- 在共同的 `ΔE = E - Emin` 分箱上，结果很清楚：
  - bin `[0, 0.005]`：
    - `kinetic_ratio ≈ 9.318`
    - `slope_ratio ≈ 8.286`
    - `num_ratio ≈ 1.165`
    - `invden_ratio ≈ 7.114`
  - bin `[0.005, 0.02]`：
    - `kinetic_ratio ≈ 5.647`
    - `num_ratio ≈ 1.164`
    - `invden_ratio ≈ 4.316`
  - bin `[0.02, 0.05]`：
    - `kinetic_ratio ≈ 4.071`
    - `num_ratio ≈ 1.163`
    - `invden_ratio ≈ 3.119`
  - bin `[0.05, 0.1]`：
    - `kinetic_ratio ≈ 3.217`
    - `num_ratio ≈ 1.162`
    - `invden_ratio ≈ 2.476`
  - bin `[0.1, 0.2]`：
    - `kinetic_ratio ≈ 2.634`
    - `num_ratio ≈ 1.158`
    - `invden_ratio ≈ 2.044`
  - bin `[0.2, 0.4]`：
    - `kinetic_ratio ≈ 2.180`
    - `num_ratio ≈ 1.152`
    - `invden_ratio ≈ 1.718`
  - bin `[0.4, 0.8]`：
    - `kinetic_ratio ≈ 1.823`
    - `num_ratio ≈ 1.141`
    - `invden_ratio ≈ 1.476`
  - bin `[0.8, 1.6]`：
    - `kinetic_ratio ≈ 1.541`
    - `num_ratio ≈ 1.123`
    - `invden_ratio ≈ 1.302`

- 这组比值给出的答案非常直接：
  - `E + E0` 这一路只提供了大约 `12%` 到 `16%` 的放大；
  - 真正把 secant slope 拉大的，是 `1 / (p(E)+p0)` 这一路，它在近阈值区给出约 `2x` 到 `7x` 的增强。

- 进一步做反事实替换后，这个判断更稳：
  - 在 `T210` 权重固定时，若只把分母项换成 `T250` 的 `1/(p+p0)`，kinetic 会大幅塌陷；
  - 若只把分子项换成 `T250` 的 `E+E0`，kinetic 只会小幅下降。
  - 例如：
    - bin `[0, 0.005]`：`actual T210 ≈ 1.152e-1`，`T210 with D250 ≈ 1.619e-2`，`T210 with N250 ≈ 9.886e-2`
    - bin `[0.2, 0.4]`：`actual T210 ≈ 6.723e-1`，`T210 with D250 ≈ 3.917e-1`，`T210 with N250 ≈ 5.835e-1`

- 解释：
  - 现在已经可以把这句“为什么”说得很具体：
    - `T210` 的 kinetic 项在整段阈值上方都更大，主因不是 `E+E0` 更大，而是 `p(E)+p0` 在同样 `ΔE` 下系统性更小；
    - 这又与前面的结论完全一致，因为 `T210` 的 `p0` 极小，意味着它始终待在更接近阈值的 low-momentum branch 上。
  - 换句话说，broad-threshold enhancement 的核心几何特征已经进一步收缩为：
    - `T210` 的 secant-slope denominator `p(E)+p0` 在整段相关能区都更小，因此 `(p(E)-p0)/(E-E0)` 系统性更大。

### 4.25 `p(E)+p0` 再拆：主导不是 `p(E)` 分支，而是异常小的 `p0`

- 为了继续回答“`p(E)+p0` 为什么在 `T210` 整段阈值上方都更小”，新增了 `diag_usbar_term4_pplusp0_breakdown.jl`。
- 这次不再看整个 secant slope，而是直接对逆分母

$$
\frac{1}{p(E)+p_0}
$$

做两类反事实：
  - 只换 `p0`：保留 `T210` 的 `p(E)` 分支，只把 `p0` 替成 `T250` 的值；
  - 只换 `p(E)`：保留 `T210` 的 `p0`，只把 `p(E)` 分支替成 `T250` 的值。

- 在共同 `ΔE = E - Emin` 分箱上，原始比值已经给出一个很强的先验：
  - `p_ratio = <p(E)>_210 / <p(E)>_250 ≈ 1.097 - 1.132`
  - `p0_ratio = p0_210 / p0_250 ≈ 0.0689`
- 也就是说：
  - `T210` 的阈值上方 `p(E)` 本身并不更小，反而系统性更大约 `10%` 到 `13%`；
  - 真正异常的是 `p0`，它只有 `T250` 的约 `6.9%`。

- 直接看逆分母反事实，结论更明确：
  - bin `[0, 0.005]`：
    - `invden_ratio ≈ 7.114`
    - 只换 `p0` 后，`p0swap_recovery ≈ 0.139`
    - 只换 `p(E)` 后，`pswap_recovery ≈ 1.064`
  - bin `[0.005, 0.02]`：
    - `invden_ratio ≈ 4.316`
    - `p0swap_recovery ≈ 0.227`
    - `pswap_recovery ≈ 1.095`
  - bin `[0.02, 0.05]`：
    - `invden_ratio ≈ 3.119`
    - `p0swap_recovery ≈ 0.311`
    - `pswap_recovery ≈ 1.107`
  - bin `[0.05, 0.1]`：
    - `invden_ratio ≈ 2.476`
    - `p0swap_recovery ≈ 0.388`
    - `pswap_recovery ≈ 1.112`
  - bin `[0.1, 0.2]`：
    - `invden_ratio ≈ 2.044`
    - `p0swap_recovery ≈ 0.466`
    - `pswap_recovery ≈ 1.114`
  - bin `[0.2, 0.4]`：
    - `invden_ratio ≈ 1.718`
    - `p0swap_recovery ≈ 0.550`
    - `pswap_recovery ≈ 1.113`
  - bin `[0.4, 0.8]`：
    - `invden_ratio ≈ 1.476`
    - `p0swap_recovery ≈ 0.636`
    - `pswap_recovery ≈ 1.107`
  - bin `[0.8, 1.6]`：
    - `invden_ratio ≈ 1.302`
    - `p0swap_recovery ≈ 0.720`
    - `pswap_recovery ≈ 1.095`

- 这组结果的物理含义很直接：
  - 若把 `T210` 的 `p0` 换成 `T250` 的值，`1/(p(E)+p0)` 会立刻塌到原来的 `14%` 到 `72%`；
  - 若只把 `p(E)` 分支换成 `T250` 的增长分支，`1/(p(E)+p0)` 几乎不塌，甚至会微幅上升到原来的 `1.06x` 到 `1.11x`。

- 因而现在可以把分母来源说得更精确：
  - `T210` 的逆分母优势并不是来自“阈值上方的 `p(E)` 分支本身更慢”；
  - 恰好相反，在同一 `ΔE` 下，`T210` 的 `p(E)` 还略大一些；
  - 真正把 `p(E)+p0` 压小的，是那个异常小的 `p0`，也就是 strange `term4` 奇点几乎贴在阈值上，使整个相关能区都被锁在 low-offset denominator 分支上。

- 换句话说，当前这条几何证据链已经进一步闭合为：
  - `E0 - Emin` 极小
  - `=> p0` 极小
  - `=> p(E)+p0` 在整段阈值上方能区都偏小
  - `=> 1/(p(E)+p0)` 系统性更大
  - `=> (p(E)-p0)/(E-E0)` 系统性更大
  - `=> regular` 的 kinetic 项系统性更大
  - `=> strange term4` 在 `T210` 更强。

### 4.26 `E0-Emin` 归因闭环：主导是 `k0-(m_s+m_u)` 的阈值超额被压扁

- 为了继续回答“为什么 `E0-Emin` 会小到把 `p0` 压成 `0.07` 量级”，新增了 `diag_usbar_term4_E0_gap_attribution.jl`。
- 对 `term4` 的 `k=0` 奇点，源码里的闭式关系可以直接写成：

$$
E_0 = -\frac{\lambda^2 + m^2 - m'^2}{2\lambda}, \qquad \lambda = -k_0,
$$

因此对于 strange `term4`（`m=m_s, m'=m_u`），有精确分解：

$$
E_0 - E_{\min}
= E_0 - m_s
= \frac{\bigl(k_0-(m_s+m_u)\bigr)\bigl(k_0-(m_s-m_u)\bigr)}{2k_0}.
$$

- 也就是说，`gap = E0-Emin` 可以被拆成两块：
  - 阈值超额：`Δ_thr = k0 - (m_s + m_u)`
  - 不对称因子：`A = (k0 - (m_s - m_u)) / (2k0)`

- 实际结果如下：
  - `k0=geom210_lit | mass=geom210_lit`：
    - `gap ≈ 1.289e-3`
    - `Δ_thr ≈ 8.805e-3`
    - `A ≈ 1.464e-1`
    - `p0 ≈ 7.231e-2`
    - `threshold ≈ 2.37107`
  - `k0=geom250_lit | mass=geom250_lit`：
    - `gap ≈ 3.162e-1`
    - `Δ_thr ≈ 1.2868`
    - `A ≈ 2.457e-1`
    - `p0 ≈ 1.049`
    - `threshold ≈ 1.66399`

- 直接比较 `T210/T250` 实际比值：
  - `gap_ratio ≈ 4.077e-3`
  - `Δ_thr_ratio ≈ 6.843e-3`
  - `A_ratio ≈ 5.959e-1`
  - `p0_ratio ≈ 6.892e-2`

- 这组数值说明：
  - `E0-Emin` 之所以缩到 `T250` 的 `0.4%`，主导不是不对称因子 `A`，它只额外给出一个约 `0.6x` 的压缩；
  - 真正决定性的是 `Δ_thr = k0-(m_s+m_u)`，它已经先缩到了 `0.7%` 量级。

- 再做单独反事实替换，结论更稳：
  - 只换质量阈值，保留 `k0_210`：
    - `gap` 会从 `1.289e-3` 跳到 `1.322e-1`
    - 相当于放大 `≈ 102.6x`
  - 只换 `k0`，保留 `m_s/m_u`：
    - `gap` 会从 `1.289e-3` 跳到 `1.245e-1`
    - 相当于放大 `≈ 96.6x`

- 这说明：
  - 在 `T210 literature_peak` 这组几何里，`E0-Emin` 的异常小并不是某一个单独变量的偶然效应；
  - 它来自“更低的 `k0`”和“更高的 strange 阈值 `m_s+m_u`”两边同时把 `Δ_thr` 压到几乎为零；
  - 而 `m_s-m_u` 进入的那条不对称因子只是在这个基础上再做次级压缩。

- `p0` 也随之闭环了。因为：

$$
p_0^2 = E_0^2 - m_s^2 = (E_0-m_s)(E_0+m_s),
$$

实际比值满足：
  - `p0_ratio ≈ 6.892e-2`
  - `sqrt(gap_ratio) ≈ 6.385e-2`
  - `sqrt((E0+m)_ratio) ≈ 1.079`

- 也就是说：
  - `p0` 极小几乎完全是 `gap = E0-Emin` 极小的直接后果；
  - `E0+m_s` 这一边只给了约 `8%` 的修正，远不是主因。

- 到这一步，几何链条已经进一步收口为：
  - `T210 literature_peak` 的 `k0` 本身更低；
  - 同时 `T210` 的 strange 阈值 `m_s+m_u` 更高；
  - 两者共同把 `Δ_thr = k0-(m_s+m_u)` 压到接近零；
  - 于是 `E0-Emin` 极小，继而 `p0` 极小；
  - 再由 `p(E)+p0` 小、`1/(p(E)+p0)` 大，把 strange `term4` 的 kinetic regular 整段放大出来。

### 4.27 strange threshold 本身的来源：`T210` 的高阈值主要还是 strange 质量没掉下来

- 为了继续回答“为什么 `T210` 的 `m_s+m_u` 会这么高”，新增了 `diag_usbar_threshold_origin.jl`。
- 这支脚本不再停留在 `E0-gap` 几何层，而是直接回到平衡解与 `K` 介子工作流，输出：
  - `m_u, m_s`
  - strange threshold `m_s + m_u`
  - `K` 极点质量与 `K - (m_s+m_u)` gap
  - 以及 `u/s` 单独替换后的 threshold counterfactual

- 实际结果如下：
  - `T210, μ_B=0`：
    - `m_u ≈ 67.893 MeV`
    - `m_s ≈ 399.983 MeV`
    - `m_s + m_u ≈ 467.876 MeV`
    - `K mass ≈ 508.204 MeV`
    - `K gap ≈ 40.327 MeV`
  - `T250, μ_B=0`：
    - `m_u ≈ 16.120 MeV`
    - `m_s ≈ 312.230 MeV`
    - `m_s + m_u ≈ 328.350 MeV`
    - `K mass ≈ 566.498 MeV`
    - `K gap ≈ 238.148 MeV`

- 因而 `T210 - T250` 的 strange threshold gap 为：
  - `Δ(m_s+m_u) ≈ 139.526 MeV`
  - 其中
    - `Δm_u ≈ 51.773 MeV`，占约 `37.1%`
    - `Δm_s ≈ 87.753 MeV`，占约 `62.9%`

- 只做单独替换时：
  - 只把 `m_u` 从 `T210` 换成 `T250`，`T210` threshold 会下降 `≈ 51.773 MeV`
  - 只把 `m_s` 从 `T210` 换成 `T250`，`T210` threshold 会下降 `≈ 87.753 MeV`

- 这说明：
  - `T210` 的高 strange threshold 不是由 light 端单独抬起来的；
  - 两个组分质量都还偏大，但更主导的是 strange 端本身，即 `m_s` 仍然显著维持在较高 constituent-mass 分支上；
  - 因而“为什么 threshold 高”这条链现在可以更具体地说成：
    - `T210` 的 light sector 还没像 `T250` 那样充分掉下去；
    - 但更关键的是 strange sector 掉得更慢，它贡献了约六成的 threshold 抬升；
    - 这也是为什么单靠 light mass 下降并不能把 `Δ_thr` 恢复到 `T250` 的量级。

### 4.28 `k0` 来源审计：当前 debug 锚点与 validation 锚点不是同一口径

- 为了继续回答“为什么所选 `k0` 会把 `Δ_thr` 压得这么小”，新增了 `diag_usbar_k0_source_audit.jl`。
- 这支脚本把四类能量锚点并排放到同一张表里：
  - 当前 debug 脚本里手写的 `manual_lit` / `manual_model`
  - `usbar_curve_counterfactuals_peak_summary.csv` 里的 baseline 扫描峰
  - `solve_gap_and_meson_point` 给出的 `K` pole mass
  - validation CSV 里的 literature `sqrt(s)` 锚点

- `T210, μ_B=0`：
  - `threshold ≈ 467.876 MeV`
  - `pole mass ≈ 508.204 MeV`
  - `baseline scan peak ≈ 500.000 MeV`
  - `manual_model ≈ 497.448 MeV`
  - `manual_lit ≈ 469.614 MeV`
  - `validation literature sqrt(s) ≈ 950.744 MeV`
  - 相对 threshold 的超额分别是：
    - baseline peak：`≈ +32.124 MeV`
    - manual model：`≈ +29.572 MeV`
    - manual lit：`≈ +1.738 MeV`
    - validation：`≈ +482.868 MeV`

- `T250, μ_B=0`：
  - `threshold ≈ 328.350 MeV`
  - `pole mass ≈ 566.498 MeV`
  - `baseline scan peak ≈ 572.000 MeV`
  - `manual_model/manual_lit ≈ 582.275 MeV`
  - `validation literature sqrt(s) ≈ 639.268 MeV`
  - 相对 threshold 的超额分别是：
    - baseline peak：`≈ +243.650 MeV`
    - manual model/manual lit：`≈ +253.925 MeV`
    - validation：`≈ +310.918 MeV`

- 这组审计结果给出两个重要更新：
  - 第一，当前 debug 文档里反复使用的 `469.614 / 582.275 / 497.448 MeV` 并不是直接从 validation CSV 继承的锚点，而是独立手写进诊断脚本的比较点；
  - 第二，`T210` 的 near-threshold 几何放大并不只是 `manual_lit ≈ 469.614 MeV` 这个极端点的人为产物，因为 baseline 扫描峰本身也只有 `threshold + 32 MeV`；
  - 但 `manual_lit` 确实把这种 near-threshold 敏感性进一步推到了最强，它只比 threshold 高 `1.74 MeV`，因此会把 `Δ_thr`、`E0-Emin` 与 `p0` 压到最极端；
  - 相比之下，当前 validation CSV 里的 literature 锚点是另一套口径，离 threshold 远得多，因此不应直接拿来和本文档里的 `manual_lit` 混成同一个“文献峰位”。

- 因而现在可以把 `k0` 这一侧的判断说得更精确：
  - 若问“为什么 `T210` 模型侧会靠近 strange threshold”，答案是：即便不看手工 literature 点，只看 baseline 扫描峰，它也确实比 threshold 只高约 `32 MeV`，这属于模型内生的 near-threshold 峰位；
  - 若问“为什么本轮 `E0-gap` 证据链会显得特别极端”，答案则是：当前使用的 `manual_lit` 锚点又把这个 near-threshold 结构额外往 threshold 边上推进了一步；
  - 因此，后续如果要把 root-cause 结论和 validation 口径完全打通，就必须显式区分：
    - internal debug anchor
    - model scan peak
    - validation literature anchor
    这三类点位，不能再把它们混叫成同一个 `literature_peak`。

### 4.29 `m_s` 高分支来源审计：当前 tested seed 族里并没有第二个也收敛的低 `m_s` 分支

- 为了继续回答“`T210` 的 strange sector 为什么会把 `m_s` 留在这么高的分支”，新增了 `diag_usbar_ms_branch_origin.jl`。
- 这支脚本不再只看当前 workflow 默认种子，而是对同一物理点并排尝试：
  - `hadron`
  - `quark`
  - `weak_chiral_conf`
  - `ht_guess_0p8 / 0p9 / 0p95`
- 它记录每个 seed 的：
  - 是否收敛
  - `omega`
  - `phi_s`
  - `Φ/Φbar`
  - `M_s`
  - 以及同一配置下是否真的出现多条已收敛 branch

- `T210, μ_B=0` 的结果是：
  - 收敛的 seed 只有 `quark` 与三组 high-temperature guess
  - 它们全部落到同一个解：
    - `M_s ≈ 399.983 MeV`
    - `phi_s ≈ -1.65447`
    - `Φ ≈ Φbar ≈ 0.65229`
    - `omega ≈ -23.3118`
  - 这些收敛结果之间的 spread 几乎为零：
    - `ms_span ≈ 0`
    - `phi_s_span ≈ 2.1e-11`
    - `Phi_span ≈ 1.9e-11`
  - 相反，`hadron` 与 `weak_chiral_conf` 都没有收敛到物理解；它们停在：
    - 更高的 `M_s ≈ 633 MeV`
    - 更低的 `Φ ≈ 0.074`
    - 但残差仍高达 `~7.18e-1`

- `T250, μ_B=0` 也给出同样结构：
  - 收敛的 seed 仍然只有 `quark` 与 high-temperature guess
  - 全部落到同一个解：
    - `M_s ≈ 312.230 MeV`
    - `phi_s ≈ -1.10291`
    - `Φ ≈ Φbar ≈ 0.76227`
    - `omega ≈ -27.8823`
  - 收敛解的 spread 也几乎为零：
    - `ms_span ≈ 0`
    - `phi_s_span ≈ 4.1e-12`
    - `Phi_span ≈ 3.3e-12`
  - `hadron` 与 `weak_chiral_conf` 同样只给出未收敛的高 `M_s`、低 `Φ` 停点。

- 这组结果说明：
  - 当前 tested seed 族里，并没有出现“另一个也收敛的低 `m_s` 分支”可供 selector 误选；
  - 换句话说，`T210` 的高 `m_s` 不是 `default_omega_selector` 或 workflow 默认 seed 偶然挑错 branch 的结果；
  - 当前在 `T210` 真正稳定收敛下来的，就是一条 `Φ≈0.65`、`M_s≈400 MeV` 的单一 branch；
  - 因而若还要继续追“为什么 `m_s` 高”，优先级已经不在 seed/selector，而应转到：
    - gap 方程本身
    - `phi_s` 的自洽平衡位置
    - 以及这一组 `Φ/Φbar` 与 strange condensate 的绑定结构。

### 4.30 baseline 峰相对 pole / `Re denom=0` / blocking 的位置：blocking 不是主导，`T210` 峰更像 near-pole 但 pre-zero 的峰

- 为了继续回答“为什么 baseline peak 本身会贴 strange threshold”，新增了 `diag_usbar_baseline_peak_alignment.jl`。
- 这支脚本对 baseline 曲线逐点输出：
  - `sigma_total`
  - `blocking`
  - `Re[1-4KΠ]`
  - `Im[1-4KΠ]`
  - 并和 `K` pole mass 一起汇总到同一张 summary 表里。

- `T210, μ_B=0`：
  - `threshold ≈ 467.876 MeV`
  - `pole mass ≈ 508.204 MeV`
  - `baseline peak ≈ 500.000 MeV`
  - `Re denom=0 ≈ 515.581 MeV`
  - 相对位置为：
    - peak 相对 threshold：`≈ +32.124 MeV`
    - peak 相对 pole：`≈ -8.204 MeV`
    - peak 相对 `Re denom=0`：`≈ -15.581 MeV`
  - blocking 变化很平：
    - peak：`≈ 0.5753`
    - pole：`≈ 0.5835`
    - `Re denom=0`：`≈ 0.5916`

- `T250, μ_B=0`：
  - `threshold ≈ 328.350 MeV`
  - `pole mass ≈ 566.498 MeV`
  - `baseline peak ≈ 572.000 MeV`
  - `Re denom=0 ≈ 688.133 MeV`
  - 相对位置为：
    - peak 相对 threshold：`≈ +243.650 MeV`
    - peak 相对 pole：`≈ +5.502 MeV`
    - peak 相对 `Re denom=0`：`≈ -116.133 MeV`
  - blocking 同样只缓慢变化：
    - peak：`≈ 0.5998`
    - pole：`≈ 0.5960`
    - `Re denom=0`：`≈ 0.6694`

- 这组结果带来三个直接更新：
  - 第一，blocking 在 peak / pole / `Re denom=0` 附近都只是在 `~0.58-0.67` 的平缓变化，不像决定峰位的主导量；
  - 第二，`T210` baseline peak 明显处在 “threshold 和 pole 之间”，并且已经靠近 pole，但仍然落在 `Re denom=0` 之前；
  - 第三，`T250` baseline peak 则几乎贴着 pole，却离 `Re denom=0` 很远，这说明在 `T250` 家族里峰位更像是由宽峰/虚部结构控制，而不是由 `Re denom` 零点直接钉住。

- 因而现在可以把“模型峰为什么贴阈值”再说得更精确：
  - `T210` 的模型峰确实不是 blocking 人为压出来的；
  - 它更像是 near-threshold strange geometry 先把 s-channel strength 提前抬起来，于是峰在到达 pole/零点之前就已经形成；
  - 但这条“提前抬起”的强度链依然要回到前面已经锁定的 strange `term4` 几何/奇点放大，而不是回到 `t` 道或 blocking。

### 4.31 gap 驻点分量归因：`T210` 的高 `m_s` 是 vacuum 与 Polyakov 势共同压住的自洽平衡，而不是 branch 误选

- 为了把“为什么当前 gap 方程会把 `T210` 固定在高 `m_s` branch 上”直接写成驻点条件，新增了 `diag_usbar_gap_stationarity_attribution.jl`。
- 这支脚本做了两件事：
  - 对平衡解处的 `Ω = χ + U_poly + Ω_vac + Ω_therm` 分别求梯度，直接看 `phi_s / Φ / Φbar` 驻点是哪些项在互相平衡；
  - 再做 counterfactual，把 `phi_s` 和 `Φ/Φbar` 单独替换成另一温度点的值，比较 `Δχ / Δpoly / Δvac / Δtherm / ΔΩ`。

- 在 `T210, μ_B=0` 的平衡点：
  - `∂Ω/∂phi_s` 的分量平衡为：
    - `chi_grad ≈ -1.3245`
    - `vac_grad ≈ +1.5493`
    - `therm_grad ≈ -0.2248`
    - 总和 `≈ 0`
  - 也就是说：
    - chiral 项与 thermal 项都在把 `phi_s` 往“更不负、更低 `m_s`”的方向推；
    - 真正把它拉回高 `m_s` branch 的主导项是 vacuum 项。
  - `∂Ω/∂Φ` 与 `∂Ω/∂Φbar` 则是：
    - `poly_grad ≈ +1.7510`
    - `therm_grad ≈ -1.7510`
    - 总和 `≈ 0`
  - 即 `Φ` 驻点主要是 Polyakov 势和 thermal 项两者直接对顶。

- 对 `T210` 再做 counterfactual：
  - 若只把 `phi_s` 换成 `T250` 的值：
    - `M_s` 会从 `399.98 MeV` 掉到 `314.24 MeV`
    - `Δchi ≈ -0.6107`
    - `Δvac ≈ +0.8027`
    - `Δtherm ≈ -0.1296`
    - 最终 `ΔΩ ≈ +0.06235`
  - 若只把 `Φ/Φbar` 换成 `T250` 的值：
    - `Δpoly ≈ +0.7027`
    - `Δtherm ≈ -0.3809`
    - 最终 `ΔΩ ≈ +0.3218`
  - 若两者一起换：
    - `ΔΩ ≈ +0.3653`

- 这说明在 `T210`：
  - 把 strange condensate 放松到 `T250` 那样，虽然 `chi` 和 `therm` 都会更有利，但会被更大的 vacuum 代价反超；
  - 把 `Φ` 抬到 `T250` 那样，thermal 的确更有利，但 Polyakov 势惩罚更大，净效应仍然升高 `Ω`；
  - 因而当前 `T210` 的高 `m_s` 与较低 `Φ≈0.65` 不是 branch/seed 误选，而是同一个自洽最低 `Ω` 平衡点上的绑定结构。

- `T250, μ_B=0` 给出对偶结论：
  - `∂Ω/∂phi_s` 的平衡为：
    - `chi_grad ≈ -0.8696`
    - `vac_grad ≈ +1.3200`
    - `therm_grad ≈ -0.4504`
    - 总和 `≈ 0`
  - 相比 `T210`，`T250` 的 thermal 拉低 `phi_s` 的作用更强，已经接近 vacuum 的抵消量级。
  - 对 `Φ` 来说：
    - `poly_grad ≈ +3.7011`
    - `therm_grad ≈ -3.7011`
    - 说明 `T250` 的 `Φ≈0.76` 是 thermal 项更强之后才稳定下来的更去禁闭平衡。
  - 若把 `T250` 的 `Φ` 压回 `T210` 的值：
    - `Δpoly ≈ -0.3286`
    - `Δtherm ≈ +0.8237`
    - 最终 `ΔΩ ≈ +0.4951`
  - 这说明在 `T250`，较高 `Φ` 之所以稳定，不是因为 Polyakov 势喜欢它，而是因为 thermal 项强烈偏好它。

- 因而现在可以把“`T210` 的高 `m_s` 为什么稳定”说得更精确：
  - `T210` 不是还有另一条低 `m_s` 稳定 branch 被错过；
  - 当前这条 branch 之所以稳定，是因为在 `phi_s` 方向上，vacuum 项仍然压过 thermal 去恢复手征的倾向；
  - 同时在 `Φ` 方向上，Polyakov 势仍把系统压在比 `T250` 更低的 loop 值上；
  - 两者共同把 `T210` 固定在“较低 `Φ` + 更负 `phi_s` + 更高 `m_s`”的同一套自洽平衡里。

### 4.32 峰形成链再收口：`T210` 的 pre-pole 峰来自 propagator 先见顶，而不是 blocking 或 numerator 先掉头

- 为了继续回答“为什么 `T210` 的 baseline 峰会在 pole 和 `Re denom=0` 之前就形成”，新增了 `diag_usbar_peak_formation_chain.jl`。
- 这支脚本把 strange `sP` 的积分链拆成：
  - `|D_s^P|^2`
  - `numerator_sP = s_12^- s_34^-`
  - `phase_kinematic_factor`，其中包含：
    - `t` 积分区间长度
    - blocking
    - `1 / (16π s_12^+)`
    - 也就是把 propagator 拿掉后的整体系数
  - 最终 `sigma_s_P`

- `T210, μ_B=0` 的结果是：
  - `sigma_s_P` 峰位：`≈ 500 MeV`
  - `|D_s^P|^2` 峰位：`≈ 494 MeV`
  - `K` pole：`≈ 508.204 MeV`
  - `Re denom=0`：`≈ 515.581 MeV`
  - `phase_kinematic_factor` 与 `numerator_sP` 的最大值都出现在当前扫描上边界 `≈ 620 MeV`

- 这组位置关系非常关键：
  - `|D_s^P|^2` 本身已经先于 pole 见顶，且比 pole 早约 `14.2 MeV`；
  - `sigma_s_P` 则落在它之后、但仍早于 pole，大约只比 pole 低 `8.2 MeV`；
  - 与此同时，phase/kinematic 因子在当前扫描窗内还在继续单调上升，根本没有提前掉头。

- 因而 `T210` 的 pre-pole 峰不是因为：
  - blocking 提前变小；
  - 或 numerator/phase-space 提前塌掉；
  - 恰好相反，phase 因子还在往高能端推峰。

- 更准确的机制是：
  - propagator 强度 `|D_s^P|^2` 自己先在 pole 之前滚到最大；
  - phase/kinematic 因子随后继续上升，把最终 `sigma_s_P` 的峰从 `494 MeV` 再往右推到 `500 MeV`；
  - 但这个右推还不足以把峰推到 pole，更不可能推到 `Re denom=0`。

- `T250, μ_B=0` 也支持同一解释框架：
  - `sigma_s_P` 峰位：`≈ 562 MeV`
  - `|D_s^P|^2` 峰位：`≈ 520 MeV`
  - `K` pole：`≈ 566.498 MeV`
  - `Re denom=0`：`≈ 688.133 MeV`
  - `phase_kinematic_factor` 与 `numerator_sP` 仍在扫描上边界 `≈ 720 MeV` 达到最大

- 这说明在 `T250`：
  - propagator 峰和 pole 的间距更大；
  - 但由于 phase 因子在更高能区继续推升，最终 `sigma_s_P` 峰会被拖到接近 pole 的位置；
  - 只是即便如此，它仍远早于 `Re denom=0`，因此也不是“实部零点直接锁峰”的机制。

- 因而现在可以把“为什么 `T210` 会在 pole 之前成峰”再收口一层：
  - 主因已经不是 blocking；
  - 也不是 numerator/phase-space 提前塌陷；
  - 更直接地说，是 strange `sP` propagator 的强度 `|D_s^P|^2` 本身就在 pole 之前先滚到最大，而 phase 因子只是在其后把截面峰轻微往右推。
- 这也给出下一步最直接的新问题：
  - 为什么 `|D_s^P|^2` 自身会在 `T210` 家族里于 pole 之前先见顶；
  - 这一步已经不再是截面层问题，而更像是 `Re/Im denom` 竞争与宽峰结构本身的传播子层问题。

## 5. 当前判断

- 当前最像的问题不在 `t` 道、不在 `σ_K` 主导、不在 `±k0` 对称化、也不在宽度耦合本身。
- 真正需要继续排查的是 strange `s` 道 `K` 通道中：
  - `A_u + A_s` 是否整体偏大；
  - `B0` 的四子项组合是否和旧实现/目标口径有系统差异；
  - 其中在同 `μ_B=0` 家族内，当前更值得优先追的已经从 `term1` 转向 `term4`，因为它才是拉高 `T210` 的 `Re B0` 的主驱动项。
  - `term1` 这边目前新增证据表明：
    - hybrid 数值积分误差基本可排除；
    - PNJL 修正确实会让 `term1` 更浅，但量级偏次级；
    - 跨 `μ_B` 对照里最强的 shift 来自 `μ_B/3`，因此不宜再把 `μ_B=800` 当主控制组。
  - 完整 `B0` 视角下，PNJL 口径对 `term2/term4` 的影响却是主导量级，尤其 `term4`；这使 `distribution_value_b0(:minus, ...)` 成为新的高优先级排查入口。
  - 进一步细看后，当前高优先级入口可以再收窄：`term2` 在相关点位上并不走奇点增强，而 `term4` 明确走奇点增强且 `:minus` 权重大得异常，因此最优先目标是 strange `term4`。
  - 最新两支脚本进一步把证据链闭合了：
    - strange `term4` 的 `:minus` 权重替换本身就足以显著拉低 `Re B0` 并抬高 `Re denom`；
    - 而这个权重之所以异常强，是因为 PNJL `antiquark_distribution` 在 relevant `E0` 处只有 `0.10` 到 `0.16`，远低于 FD 的 `0.82` 到 `0.87`。
  - 再往下拆后，目前可以更精确地区分两种“主因”：
    - 若比较 `T210` 和 `T250` 哪一个更强，主因是 strange `term4` 的几何/奇点放大；
    - 若比较 PNJL 权重与更弱权重（blend/FD）会不会显著拉低峰强，则答案仍然是会，且这是直接因果杠杆。
  - 最新的曲线级反事实把这两个判断都推进到了完整 `σ` 层：
    - 只换 `term4` 权重时，`T210` 曲线会从 `~49 mb` 直接塌到 `~7-19 mb`；
    - 只换几何时，`T250` 曲线会被抬到 `~28 mb` 并左移到 `~554 MeV`，而 `T210` 曲线则会被明显压低。
  - 最新的 geometry 细拆又进一步说明：
    - `T210` 的优势不是 logterm 更大，而是 regular 主值积分本身更大；
    - 并且这个 regular 优势不是单点尖峰，而是覆盖了 `ΔE \lesssim 0.8 fm^-1` 的整段阈值上方能区。
  - 最新的 factor decomposition 又进一步说明：
    - 这部分 regular gap 几乎完全来自动量项 `(p(E)-p0)w(E)/(E-E0)`；
    - 分布权重项 `p0(w(E)-w0)/(E-E0)` 对 `T210` 基本可以忽略，对 `T250` 虽可见但仍是次级修正。
  - 最新的 kinetic geometry breakdown 又进一步说明：
    - 动量项之所以大，不是因为 `E+E0` 明显更大；
    - 真正主导的是 secant slope 分母 `p(E)+p0` 更小，而这一路在近阈值区提供了远大于分子项的增强。
  - 最新的 `p(E)+p0` 再拆又进一步说明：
    - 这个更小的分母也不是由 `p(E)` 阈值增长分支主导；
    - 直接主导项是异常小的 `p0`，即 strange 奇点与阈值的极小间隙本身。
  - 最新的 `E0-Emin` 解析归因又进一步说明：
    - 这个极小间隙的直接来源，是 `k0-(m_s+m_u)` 的阈值超额几乎被压扁；
    - `m_s-m_u` 对应的不对称因子只是在其上做次级压缩。
  - 最新的 threshold-origin 脚本又进一步说明：
    - `T210` 的高 strange threshold 不是单纯由 light 端造成；
    - `m_s` 本身仍维持得更高，贡献了约 `63%` 的 `T210-T250` threshold gap，`m_u` 贡献约 `37%`。
  - 最新的 `m_s` branch audit 又进一步说明：
    - 在当前 tested seed 族里，并没有第二个也收敛的低 `m_s` 分支；
    - `T210` 的高 `m_s` 因此不是 selector 偶然挑错 branch，而是当前求解器稳定收敛到的单一 branch 本身就维持在高 strange constituent-mass 上。
  - 最新的 gap stationarity attribution 又进一步说明：
    - 在 `T210` 的 `phi_s` 方向上，真正把系统压回高 `m_s` branch 的主导平衡项是 vacuum；
    - 在 `Φ` 方向上，则是 Polyakov 势与 thermal 项的对顶，而且当前 `T210` 仍被 Polyakov 势压在较低 loop 值上；
    - 因而“高 `m_s` + 低 `Φ`”是同一个自洽最低 `Ω` 平衡点的绑定结构，而不是独立偶然量。
  - 最新的 `k0` source audit 又进一步说明：
    - 当前 debug 链条里的 `manual_lit/manual_model` 不是 validation CSV 的同一套锚点；
    - 其中 `T210 manual_lit` 只比 threshold 高约 `1.74 MeV`，会把 near-threshold 几何敏感性推到最强；
    - 但就算退回 baseline 扫描峰，`T210` 仍只比 threshold 高约 `32 MeV`，因此“模型峰本身贴 strange threshold”这件事仍然成立。
  - 最新的 baseline peak alignment 又进一步说明：
    - `T210` baseline peak 落在 threshold 与 pole 之间，且先于 `Re denom=0` 出现；
    - `T250` baseline peak 则几乎贴 pole 但远早于 `Re denom=0`；
    - 两个家族里 blocking 都只缓慢变化，因此它不是决定峰位位置的主导杠杆。
  - 最新的 peak formation chain 又进一步说明：
    - `T210` 的 `sigma_s_P` pre-pole 峰不是因为 phase-space 或 blocking 提前掉头；
    - 相反，phase/kinematic 因子还在继续上升，真正先见顶的是 propagator 强度 `|D_s^P|^2` 本身；
    - 这把问题继续收口到了传播子层的 `Re/Im denom` 竞争，而不再只是截面层乘子分解。
  - `antiquark_distribution` 本身的参数归因也已基本明确：
    - `T210/T250` 间的相对变化主要由 `x=(E+μ)/T` 驱动；
    - PNJL 相对 FD 的额外压低则更多体现为分母抬高，`T210` 还有分子压低共同参与。

- 当前阶段的排序判断：
  1. `B0` 实部结构偏移最可疑，其中最高优先级已经收敛到 strange `term4` 的 `:minus` 分支热权重；
  2. 在 `μ_B=0` 家族内部，决定 `T210/T250` 强弱差异和峰位迁移的主因素是 strange `term4` 的几何/奇点结构，而不是 thermo 权重本身；
  3. 该几何结构的直接放大路径已经更清楚：较小的 `|λ|` 与显著更大的 `regular` 主值积分共同抬高 `|term4|`，而不是 logterm 单独主导；
  4. 在 regular 主值积分内部，这个 family gap 进一步可归结为动量项主导；分布权重项对 `T210` 几乎可忽略，对 `T250` 也只是次级修正；
  5. 动量项内部的直接几何杠杆也已更清楚：主导不是 `E+E0` 分子，而是 `p(E)+p0` 分母更小；
  6. 而这个更小的 `p(E)+p0`，主导也不是 `p(E)` 分支，而是异常小的 `p0`；
  7. 而这个异常小的 `p0`，又可继续收敛到 `k0-(m_s+m_u)` 的阈值超额几乎为零；`m_s-m_u` 只做次级修正；
  8. 而这个高 strange threshold 的内部来源也已有定量归因：light 端约占 `37%`，strange 端约占 `63%`，因此 threshold 抬升更像是 strange sector 仍留在高 constituent-mass 分支；
  9. 当前 tested seed 族里并没有第二个也收敛的低 `m_s` 分支，因此高 `m_s` 不是 selector/seed 误选，而是当前单一稳定 branch 的内生结果；
  10. 当前这条高 `m_s` branch 的自洽来源也更清楚了：`phi_s` 方向主要是 vacuum 压过 thermal，`Φ` 方向则是 Polyakov 势把系统压在较低 loop 值；
  11. 该热权重异常的直接数值表现，是 PNJL `antiquark_distribution` 在 relevant `E0` 处被压得远低于 FD，从而把 `1 - \bar f(E)` 撑到接近 1；
  12. 还需要单独记住一个口径事实：当前 debug 中的 `manual_lit/manual_model` 与 validation literature anchor 不是同一套点位，尤其 `T210 manual_lit` 会把 near-threshold 结构推到更极端；
  13. baseline peak 与 pole/零点的相对位置也已经更清楚：blocking 不是决定峰位的主导量；
  14. 更进一步说，`T210` 的 pre-pole 峰也不是因为 phase/numerator 提前塌掉，而是因为 `|D_s^P|^2` 自身先于 pole 见顶；
  15. 因而峰位问题的下一层更像是传播子层的 `Re/Im denom` 竞争，而不是截面层 kinematic multiplier；
  16. `A_sum` 次之；
  17. `term2` 仍需保留观察，但优先级已低于 `term4`，因为它在相关点位上不表现为同类奇点增强；
  18. `K4567_plus` 映射本身最不可疑；
  19. `term1/term4` 的 hybrid 数值积分器误差基本可排除。

## 6. 建议的下一步

- 继续沿同 `μ_B=0` 家族往下拆，但现在把重点明确锁定为 strange `term4`：
  - “只换权重”与“只换几何”对完整 `σ(usbar→usbar)` 的影响已经拿到，下一步不需要再停留在这一层的因果确认。
  - `p(E)+p0` 的来源拆分也已经拿到，下一步不需要再停留在“`p0` 还是 `p(E)`”这一层。
- 继续追 `antiquark_distribution` 的参数归因：
  - `x` vs `Φ/\bar Φ` 的一阶归因也已经得到；
  - 现在这一步也已经完成；下一步不必再停留在“p(E) vs 权重”这一层。
- 进一步核对文献/旧实现在 `f^-` 口径上的定义和参数取值，确认当前 PNJL `antiquark_distribution` 是否在 relevant `E0` 区域给出了过低的 `\bar f(E)`。
- 对 `A_u`、`A_s`、`G_u`、`K4567_plus` 做来源追踪，判断是 `A` 的口径问题，还是有效耦合映射问题；但这一条当前优先级已经低于 strange `term4` 权重链路本身。
- 若下一轮继续向下做 root-cause 闭环，更值得做的是：
  - `E0-Emin` 的来源已经闭环到 `k0-(m_s+m_u)`；
  - 现在 threshold 本身也已经拆到 `m_u/m_s` 两部分，而 seed audit 又表明高 `m_s` 不是 branch 误选；所以下一步更值得直接追的是：为什么当前 gap 方程会把 `T210` 的 `phi_s`、`Φ/\bar Φ` 自洽地固定在这条高-`m_s` branch 上。
- 在 gap 这条线上，下一步不必再回到 seed/selector；
  - 更值得直接追的是 `vacuum` 与 `thermal` 在 strange condensate 方程里的更细粒度来源，例如：`M_s` 依赖、热积分核、以及它们对 `phi_s` 曲率的分工。
- 继续清理点位口径：
  - 当前 debug 文档里的 `manual_lit/manual_model` 与 validation literature anchor 已确认不是同一套点；
  - 下一轮若继续做跨文献对照，应把 internal debug anchor、model scan peak、validation literature anchor 三套点位拆开命名，避免再混用 `literature_peak`。
- 继续做“模型峰为何贴阈值”的最后一跳：
  - `T210` baseline peak 已确认只比 threshold 高 `~32 MeV`，这说明 near-threshold 不是只由手工 literature 点制造的；
  - baseline peak 相对 `K` pole / `Re denom=0` / blocking 的位置现在也已量化；
  - 而最新链条又表明 phase-space 没有提前塌掉，真正更值得追的是：为什么 `|D_s^P|^2` 会在 `T210` 家族里于 pole 之前先见顶，也就是 `Re/Im denom` 的竞争是如何把 propagator 峰提前的。
- 若后续结论稳定，再把本文档归档到 `docs/dev/archived/`，并补归档元数据。