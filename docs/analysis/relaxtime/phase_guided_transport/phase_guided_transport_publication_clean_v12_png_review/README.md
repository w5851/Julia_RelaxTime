# publication_clean_v12 PNG 审查层

v12 从相同冻结点表生成 72 张单图、2 张主复合图和 1 张局部视图；
每张都交付彩色及灰度 PNG，共 150 张。所有图按 171.45 mm 宽、600 dpi
导出，保持 v11 两张主图的原生画布尺寸，不使用紧裁切。

所有图改用图外公共图例。复合图的图例和带数学下标的列标题为 13 pt，刻度为
11 pt；最小大写/数字字形必须达到 2 mm，不再使用紧凑字号例外。
主图坐标语义和全部曲线顶点沿用 v11。新增局部图仅给出明确的线性纵轴
窗口，帮助辨认 Figure 1 的低值细节；它不替代完整主图，也不新增数据。

查看 `placement_report.json` 的逐图插入区间和 `caption_handoff.md` 的图注。
推荐按 171.45 mm 宽插入；单图也不能直接缩成单栏。超过推荐宽度时，PNG
有效 DPI 会降低；低于最小宽度时，小字或标记不合格。PDF 后续交付时需要
重新计算其尺寸约束，不用 PNG 分辨率判定矢量曲线质量。

`plot_manifest.json` 位于对应 figure case，集中记录所有输入、生成器、
逐图布局、尺寸限制及彩色/灰度输出 hash。灰度图片已经生成；作者仍需
在最终尺寸逐项检查线型追踪、端点、局部细节和文字。

这是独立绘图任务：manuscript_eligible=false，current=v5，PDF pending。
v5/v10/v11、原始计算和分支断线均保留；不调用求解器或数值收敛门禁。
本包没有验证当前论文整页排版，也不等同于投稿资格或数值生产晋升。

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v12.py --png-review
```
