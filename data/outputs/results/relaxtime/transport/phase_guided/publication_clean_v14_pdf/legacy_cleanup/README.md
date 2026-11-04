# 旧 phase-guided 数据与图像清理盘点

当前状态：仅盘点并从本工作分支的 current 入口退役旧版。没有删除、移动或新建备份归档。新 v14、full-RS 正式数据、回归基线、数值验证和 Mott 诊断均排除在清理范围之外。

- `working_branch`：`D:/w/meson-freezeout-defaults`，本次 current 更新所在工作树。
- `desktop_checkout`：`D:/Desktop/Julia_RelaxTime`，本次只读盘点，未同步其分支、指针或产物。

共盘点 76 个目录、3872 个文件、877,364,656 字节（约 836.7 MiB）。包含两份 checkout 的重复资产，也包含必须保留的历史收敛证据，不能把总量直接当作可删除量。

| 工作树 | 类别 | 目录 | 文件 | MiB |
| --- | --- | ---: | ---: | ---: |
| working_branch | 历史收敛证据，保留 | 2 | 15 | 3.0 |
| working_branch | 旧图像／显示包 | 32 | 1842 | 373.1 |
| working_branch | 旧数值网格 | 2 | 66 | 40.8 |
| desktop_checkout | 历史收敛证据，保留 | 2 | 15 | 3.0 |
| desktop_checkout | 旧图像／显示包 | 36 | 1868 | 375.9 |
| desktop_checkout | 旧数值网格 | 2 | 66 | 40.8 |

完整绝对路径、逐文件字节数和 SHA-256 见 [inventory.json](inventory.json)。两工作树的活动说明、脚本、测试和配置中有 955 条相关文本匹配，见 [live_references.json](live_references.json)；这些匹配包括历史叙述，不等于同样数量的有效运行依赖。

论文 `main.tex` 仍引用 v13 Figure 1/2 和旧顶部图例 FIG4，三个实际附件的哈希已纳入清单。v14 的形式依据和历史 renderer 也需要可恢复的旧 manifest／代码。旧数值存在已知问题，不意味着诊断、回归或复现证据可以删除。

建议执行范围是：将旧图像／显示包和旧网格从常用产物目录迁出，保留一份经完整恢复校验的历史归档；它们不再作为新论文数据源。历史收敛证据保留于归档证据链。该执行需要明确选定以下清单中的绝对目录和归档位置后进行：

1. 为所选目录生成无损归档，保存仓库、分支、HEAD、路径和逐文件哈希。
2. 解压到新的空目录，逐文件复核，不覆盖当前工作树。
3. 处理实际依赖和论文附件映射，保留新 v14 renderer 所需的旧形式参考。
4. 复核最终绝对目标仍位于所选根目录后，再删除所选冗余副本；不重写任何冻结 CSV、图像或 manifest。

下表是候选盘点，不是已批准的删除 allowlist。

| 工作树 | 目录 | 处理建议 |
| --- | --- | --- |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v10_png_review` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_pdf_review` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_png_review` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v12_png_review` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v13_pdf` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v13_png_review` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v2` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v3` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v4` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v5` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_validated_anchored_prod_v1_convergence` | 保留历史收敛证据 |
| working_branch | `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_xi001_validated_anchored_prod_v1_convergence` | 保留历史收敛证据 |
| working_branch | `data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled` | 归档并完成引用核对后再清理 |
| working_branch | `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v1` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v2` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v3` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v4` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v5` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v1` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v10_png_review` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v11_paper_size_review` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v11_pdf_review` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v11_png_review` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v12_png_review` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v13_pdf` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v13_png_review` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v2` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v2_residual_smoothed` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v3` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4_mechanism_review` | 归档并完成引用核对后再清理 |
| working_branch | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v5` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v1` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v10_png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_pdf_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v12_png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v13_pdf` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v13_png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v2` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v3` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v4` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v5` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/transport_fixed_T_panel_legends_20261021__pdf` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/figures/relaxtime/transport/phase_guided/transport_fixed_T_panel_legends_20261021__png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_validated_anchored_prod_v1_convergence` | 保留历史收敛证据 |
| desktop_checkout | `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_xi001_validated_anchored_prod_v1_convergence` | 保留历史收敛证据 |
| desktop_checkout | `data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/results/relaxtime/transport/phase_guided/transport_fixed_T_panel_legends_20261021__pdf` | 归档并完成引用核对后再清理 |
| desktop_checkout | `data/outputs/results/relaxtime/transport/phase_guided/transport_fixed_T_panel_legends_20261021__png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v1` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v2` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v3` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v4` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_figure_layer_v5` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v1` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v10_png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v11_paper_size_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v11_pdf_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v11_png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v12_png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v13_pdf` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v13_png_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v2` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v2_residual_smoothed` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v3` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4_mechanism_review` | 归档并完成引用核对后再清理 |
| desktop_checkout | `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v5` | 归档并完成引用核对后再清理 |
