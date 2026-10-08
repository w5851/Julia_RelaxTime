# publication_clean_v13 矢量 PDF

v13 PNG 已获作者版本级接受；本包交付 72 张单图、2 张主复合图、1 张局部图，
共 75 份单页彩色矢量 PDF。[图件索引](pdf_index.md)列出全部文件；
[图注交接](caption_handoff.md)保留科学说明与冻结温度映射。

所有图保持 171.45 mm（6.75 in）原生宽度。每张图均通过同一冻结 renderer
重绘与已接受 PNG 的 RGBA 像素逐一比对；坐标、曲线顶点、分支 gap 和图例不变。
PDF 逐文件检查物理尺寸、单页、嵌入字体、非 Type 3 字体和无嵌入栅格图像；
完整逐图测量及文件 hash 存于对应图像目录的唯一 `plot_manifest.json`。

[尺寸报告](placement_report.json)分别给出 PDF 与 PNG 的插入限制。
推荐维持原生宽度；最小实测字形超过 2 mm，未使用紧凑字形例外。
矢量 PDF 不受 PNG 分辨率上限约束；双栏图仍不可直接缩成单栏图。

原 PNG/灰度图及其 manifest 保持冻结，新 manifest 引用原文件。
接受记录仅表明版本级人工接受，不声称每份灰度图均被逐项人工审查。
本包完成 `vector_delivery`，保持 `manuscript_eligible=false`、current=v5；
没有调用求解器、改变数值数据、新增平滑或取得数值生产／论文装配资格。

```powershell
python scripts/analysis/relaxtime/export_phase_guided_publication_clean_v13_pdf.py --check
```

首次生成不带 `--check`；已存在的输出目录拒绝覆盖。
后续独立复核记录见 `verification/`（生成器本身不写入人工审阅结论）。
