# publication_clean_v13 PNG 审查层

v13 从相同冻结点表生成 72 张单图、2 张主复合图和 1 张局部视图；
每张都交付彩色及灰度 PNG，共 150 张。所有图按 171.45 mm 宽、600 dpi
导出，保持 v11 两张主图的原生画布尺寸，不使用紧裁切。

两张主图把参数 key 放在 (a)、一阶端点 key 放在 (c)，低值局部图在 (a)
共享参数 key。参数共同标题为 α_T，条目为 1.0 / 1.1 / 1.2；端点共同标题
为 First-order，条目为 restored / broken。标题 α_T 和列标题为 13 pt，
图例条目、一阶标题、复合图刻度为 11 pt；实测字形必须达到 2 mm。
每张单图依次检查声明的四角组合；全部图例通过曲线／端点／互相遮挡检查后
才能使用，无法放下时允许简短顶部备选。实际位置以逐图 manifest 为准。
mode-A 的参数为 α_T；固定 T 的 mode-B 单图标题为 μ_B (MeV)，条目为
0 / 450 / 900。两种模式的参数、颜色／线型映射及端点适用曲线分别保留。
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
v5/v10/v11/v12、原始计算和分支断线均保留；不调用求解器或数值收敛门禁。
本包没有验证当前论文整页排版，也不等同于投稿资格或数值生产晋升。

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v13.py --png-review
```
