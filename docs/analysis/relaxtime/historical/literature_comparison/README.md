# Historical Literature Comparison Figure Evidence

## 角色

本目录保存已从 `data/outputs/figures/relaxtime/literature` 迁入的五张历史比较图。它们保留为 `analysis/legacy` evidence，不是当前 `strict` 或 `current_candidate` 正式图，也不构成 external validation gate。

本次迁移只涉及图像文件和本说明；原始结果 CSV、文献 target、绘图脚本和数值语义不变。源文件的 SHA-256 记录在 `docs/analysis/governance/figure_asset_registry_v1/asset_registry.json`、`cleanup_preflight_v1.json` 和 `retirement_execution_v1.json`。

## 图像

- `k_mass_literature_julia_fortranmott_comparison.png`
- `meson_mass_julia_vs_literature_comparison.png`
- `sigma_literature_error_by_process.png`
- `sigma_literature_overlay_by_process.png`
- `tau_literature_comparison.png`

## 输入与生成线索

- K 介子质量比较对应 `data/outputs/results/relaxtime/literature/k_mass_literature_julia_fortranmott_comparison.csv`。
- 介子质量比较对应 `data/outputs/results/relaxtime/literature/meson_mass_julia_vs_literature_comparison.csv`。
- 截面比较使用 `tests/validation/data/raw_long/relaxtime/sigma/relaxtime_sigma_literature_digitized_longtable_v1.csv`、`data/outputs/results/relaxtime/cross_section/literature_compare/` 下的模型扫描和 `sigma_literature_compare_summary.csv`；对应脚本为 `scripts/relaxtime/plot_sigma_literature_comparison.py`。
- 弛豫时间比较使用 `tests/validation/data/raw_long/relaxtime/tau/relaxtime_tau_literature_digitized_longtable_v1.csv`、`data/outputs/results/relaxtime/scan/relaxation_times_vs_T_literature_compare.csv` 和其 summary；对应脚本为 `scripts/relaxtime/plot_tau_literature_comparison.py`。
- 当前 tracked text 中没有找到前两张质量图的专用生成器或完整 manifest，因此它们的生成 provenance 保持为 registry 中记录的历史证据，不在本目录中补写未经验证的来源。

## 解释边界

`plot_sigma_literature_comparison.py` 和 `plot_tau_literature_comparison.py` 会在文献 target 横坐标处对模型扫描作局部线性插值，并据此计算误差统计。该逻辑是历史 comparison/analysis 语义；它不能直接迁移到 `strict`。未来 strict 图只画原始模型 support 点，不把插值后的误差线或连接线当作新的模型 support。

文献 digitization、通道、归一化和模型约定尚未形成统一 external validation 合同。因此这些图可用于历史对照和内部讨论，不应表述为“literature-validated”或当前论文定稿图。
