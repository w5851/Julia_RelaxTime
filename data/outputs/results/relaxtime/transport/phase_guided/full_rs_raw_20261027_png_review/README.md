# Full-RS 原始 PNG 审阅包

已生成 75 张彩色及 75 张灰度 PNG，600 dpi。图包：data/outputs/figures/relaxtime/transport/phase_guided/full_rs_raw_20261027_png_review

原始 CSV 的 1819 点、12 个量全部保留；所有绘图顶点逐一与源表核对。未平滑、插值、删点或替换值，仅在两处已审计的一阶不连续之间断线。

两张主图、四面板低值细节图及 72 张单图的相对路径与 v13 一一对应。前两张低值子图的下界由 0.58 改为 0.57 fm，以显示新的 αT=1.2 支。

图注、单位及精确温度见 [caption_handoff.md](caption_handoff.md)；机器验收见 [quality_review.json](quality_review.json)。源码与绘图合同的原字节保存在 source/，计算提交为 2e5c411c716d5cbc524c3fc6cfb335fe05b39444。

主要图片：[弛豫时间主图](../../../../../figures/relaxtime/transport/phase_guided/full_rs_raw_20261027_png_review/composites/figure1_relaxation_times_comparison.png)、[输运比值主图](../../../../../figures/relaxtime/transport/phase_guided/full_rs_raw_20261027_png_review/composites/figure2_transport_coefficients_comparison.png)、[低值细节图](../../../../../figures/relaxtime/transport/phase_guided/full_rs_raw_20261027_png_review/details/figure1_low_value_details.png)。

[审图记录](agent_review.json)覆盖三张主要彩色图、全部 72 张单图的缩略总览，以及代表性灰度视图。
μB=900、αT=1.0 的 σ/σ/T 在 ξ≈−0.04 仍有小凸起，αT=1.1 在 ξ≈0.20 也有局部小凸起；
它们来自原始数据，物理成因本轮未判断。[对应原始行](visual_observation_source_rows.csv)保留可核对数值。
17 项定向绘图测试全部通过；机器验收与视觉观察不构成“所有非物理尖峰已消除”的结论。

本包停在用户要求的 PNG 视觉审阅阶段。没有调用 Julia，没有生成 PDF、修改论文或替换当前图层。作者接受尚待记录。
