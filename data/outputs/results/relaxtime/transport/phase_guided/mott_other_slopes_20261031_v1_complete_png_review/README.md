# 其他 Mott 完整参考诊断图（2026-10-31）

13 个窗口各展示 τū、τs̄、η/s、ζ/s、σ/T，共 65 张彩色和 65 张同尺寸灰度 PNG。[图像 manifest](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261031_v1_complete_png_review/plot_manifest.json)绑定冻结输入、绘图源码、样式、全部输出及 SHA-256；[图注](caption.md)说明显示语义。图像保持诊断 PNG 审阅阶段，需按 6.75 英寸宽度阅读。

数据来自 quote-remote 的完整主网格及加严参考：[数值审阅](../../threshold_cache_upgrade_20261021/mott_other_slopes_20261029_v1/scientific_review_v1.json)与[分析说明](../../../../../../../docs/analysis/relaxtime/threshold_cache_upgrade_20261021_analysis/mott_local_bumps_20261028.md)记录结论。K 窗口也存在局部斜率调制；π 最窄斜率的连续极限仍未验收。图面通过不等于积分收敛或介子排他性归因。

| 重点图 | 内容 |
| --- | --- |
| [other01：τs̄](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261031_v1_complete_png_review/other01_tau_sbar.png) | 无 π 粗括区重叠的 K 窗口，奇异味响应 |
| [other03：η/s](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261031_v1_complete_png_review/other03_eta_over_s.png) | K 窗口的末端系数斜率调制 |
| [other07：σ/T](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261031_v1_complete_png_review/other07_sigma_over_T.png) | 混合 σ 窗口的背景趋势与弱细节 |
| [other05：τū](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261031_v1_complete_png_review/other05_tau_ubar.png) · [other08：τū](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261031_v1_complete_png_review/other08_tau_ubar.png) | 两处 π 对照的 ξ 尺度敏感结构 |
| [other13：ζ/s](../../../../../figures/relaxtime/transport/phase_guided/mott_other_slopes_20261031_v1_complete_png_review/other13_zeta_over_s.png) | 有 π 粗括区重叠的 K 窗口，保留半区间标记 |

上面板为直接计算值；下面板为远端表中的 Δln(y)/Δξ，放在各区间中点。连线仅帮助阅读离散值，没有拟合、平滑或新增数值插值。

| 图例 | 含义 |
| --- | --- |
| 蓝线／空心圆 | 主 ξ 细网格，s 表 rtol=10⁻⁵ |
| 橙色空心方块 | 原 ξ 坐标的重算点，仅在上面板显示 |
| 黑叉 | rtol=10⁻⁶ 参考；下面板只比较相同端点 |
| 绿色加号 | 加严参考中新增 ξ 点两侧的半区间割线 |
| 灰色点线 | 该窗口定位到的最低阈值逆传播子零点 |

同端点结果一致时标记会重合，未做横向错位。稀疏参考之间跨度不同的割线不混入精度比较。两个 π 窗口的新增半区间斜率仍显著变化，不能把极窄尖峰幅度作为已收敛结果使用。

[机器验证](validation.json)和再次执行的公共图像校验器通过全部 130 个输出。[图面审阅记录](agent_review.json)披露辅助模型和主代理各自的覆盖范围：DeepSeek 检查 11 张彩色联系表及选定细节，23 项中 21 项结构合格；两项失败由主代理直接补查。主代理通过三张详细／组合预览检查五幅原图，并否决了“绿色加号被裁切”等误报。三张灰度代表获得辅助筛查，其中两张由主代理复核；其余灰度图只做机器检查。原始回复、失败和逐项处置均保存在 [DeepSeek 审阅目录](visual_qa/ds_screen_v1/assessment.json)。

作者尚未审阅这套诊断图；联系表不能替代全部原图在最终尺寸下的作者审阅。本图包没有生成 PDF、替换正式数据或更新 current，v15 PDF 继续按作者要求暂缓。
