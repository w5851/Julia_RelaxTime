# Freeze-Out K-/pi- Plot-Review Analysis

## 角色

这是计划从 `freezeout_kminus_piminus_mu_pi_100` 正式图目录迁出的历史 plot-review evidence，目标归入 `docs/analysis`。它用于人工检查 workflow 与 digitized target 的趋势和差异，不是 external validation gate、regression truth source 或 `strict` 论文图。

## 输入与上下文

- target: `data/outputs/results/relaxtime/literature/meson_density_targets/blaschke2019col_kminus_piminus_mu_pi_100_fig4_right.csv`
- workflow source: `data/outputs/results/relaxtime/meson_density/freezeout_validation/blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default/workflow_scan.csv`
- comparison source: `data/outputs/results/relaxtime/meson_density/freezeout_validation/blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default/comparison_vs_target.csv`
- generator: `scripts/relaxtime/build_meson_density_plot_review_case.py`
- source result README: `data/outputs/results/relaxtime/meson_density/plot_review/freezeout_kminus_piminus_mu_pi_100/README.md`

该 workflow 使用 literature target 的 `sqrt(s_NN)` 采样点，经仓库的 freeze-out baseline、path 和化学势 profile 得到模型值。原始结果 CSV 留在 `data/outputs/results`，本次迁移不改写它们。

## 计划迁移的图像证据

- `overlay_kminus_piminus_mu_pi_100.png`：literature target 与 workflow model 的叠图。
- `residual_kminus_piminus_mu_pi_100.png`：absolute difference 与 finite relative difference 的两层诊断图。
- `plot_manifest.json`：原 plot-review manifest，记录 comparison CSV 和图像关系。

现有摘要为 48 个点、最大绝对差 `0.147787`、最大有限相对差 `0.873012`。这些数值只描述该 workflow reproduction / staged comparison，不构成通过或失败的 validation 门槛。

## 迁移后的 manifest 注意事项

当前 `plot_manifest.json` 是源 figure 目录生成的历史 manifest；其中的 figure path 仍指向原始 `data/outputs/figures/...` 路径。本 PR 只建立预检和说明，不执行移动，也不修改 manifest。实际移动获批后，必须由人工确认采用以下哪一种可追溯处理：保留原路径作为历史 provenance，或只更新 manifest 的路径字段并重新记录 manifest hash；两者都不得改动 CSV、图像像素或数值语义。
