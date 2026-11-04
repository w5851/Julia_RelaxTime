# publication_clean_v15 PNG 审阅

已完成76张彩色PNG及76张同尺寸灰度PNG，包含72张单图、Figure 1/2/4三张主组合图及1张局部细节图。
数值来自[η局部加密正式版本](../../threshold_cache_upgrade_20261021/full_rs_eta_augmented_20261029_v1/README.md)：
1867点（1810个保留原值、9个原坐标重算、48个新增ξ点），绘图读取原始列，不平滑、插值或值混合。

| 图像 | 入口 |
| --- | --- |
| 全部图与manifest | [PNG目录](../../../../../figures/relaxtime/transport/phase_guided/publication_clean_v15_png_review/) |
| Figure 1：弛豫时间 | [PNG](../../../../../figures/relaxtime/transport/phase_guided/publication_clean_v15_png_review/composites/figure1_relaxation_times_comparison.png) |
| Figure 2：输运比值 | [PNG](../../../../../figures/relaxtime/transport/phase_guided/publication_clean_v15_png_review/composites/figure2_transport_coefficients_comparison.png) |
| Figure 4：固定T | [PNG](../../../../../figures/relaxtime/transport/phase_guided/publication_clean_v15_png_review/composites/figure4_transport_fixed_T.png) |
| Figure 1低值细节 | [PNG](../../../../../figures/relaxtime/transport/phase_guided/publication_clean_v15_png_review/details/figure1_low_value_details.png) |
| 灰度审阅 | [目录](../../../../../figures/relaxtime/transport/phase_guided/publication_clean_v15_png_review/grayscale_review/) |

[形式核验](form_review.json)逐一比较76图的尺寸、面板及轴类型、坐标范围与刻度、标签和单位、
图例位置及锚点、颜色与线型，以及一阶断线端点，全部与v14一致。
复用的v14 renderer SHA-256为 `ae916543fe5ee709361ec2347584320b565cff6b2c2242ff7a7ba77ae4438a17`。
37张彩色PNG与v14逐字节一致，其余反映本次数值选源或补点变化。

[全套质量检查](quality_review.json)通过：33034个绘图顶点与所读源表逐一核对，
全部PNG的字形、线宽、四侧刻度、裁切、遮挡、灰度与来源hash检查通过；最小实测字高2.09621 mm。
[目视记录](visual_review.json)覆盖7张彩色图和3张灰度图，包括三张主图、局部细节图以及三条加密曲线的电导率单图。
已查看的图中标签和图例可读、线型可区分，未见新增接缝；目视检查不替代数值收敛判断。
合并/图像输入合同的[14项测试](implementation_checks.json)通过。

[完整图注、单位及参数映射](caption_handoff.md)与冻结来源随本包保存。
数值回放与派生验证在quote-remote完成，本机仅进行文件核验和Python绘图。
本图包状态为 `png_review`、`manuscript_eligible=false`、`current_publication_layer=false`。
依绘图SOP，新的PNG需要作者视觉接受后才进入PDF交付；当前指针继续指向已接受的v14。
