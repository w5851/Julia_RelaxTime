# v13 PDF 交付复核

本记录对应本次 PDF 交付与 SOP 精简；v13 PNG 阶段的原始验证包保持冻结。

| 项目 | 结果 |
| --- | --- |
| 交付数量 | 72 单图 + 2 主复合图 + 1 局部图，共 75 PDF |
| 与作者已接受 PNG 的一致性 | 75/75 使用冻结 renderer 在 600 dpi 重绘，RGBA 像素完全一致；panel specs 和所有 legend 声明不变 |
| 实际 PDF | 75/75 单页、尺寸一致、字体嵌入、无 Type 3 字体、无嵌入栅格图像 |
| 尺寸 | 原生宽度 171.45 mm；全包共同 PDF 宽度区间约 163.581–177.8 mm；最小字形约 2.096 mm |
| 原始资产保护 | 802/802 既有文件 SHA-256 不变，包含 v13 PNG 与其验证包、较早图包和 current=v5 指针 |
| 测试 | 11 passed；覆盖接受门槛、源码漂移、像素变化、图例变化、禁止覆盖、全部 PDF 合同和原 PNG/灰度证据 |
| 文档 | docs-consistency、sop-governance 通过；脚本入口检查通过 |
| 实际 PDF 视觉复核 | 两张主图、局部图、mode-B 一阶端点单图、mode-A 顶部 key 单图；未见字体缺失、裁切或图例遮挡 |

SOP 从 379 行精简为 132 行（正文字符减少约 70.9%）。pilot 结果、版本/PR 操作叙述及重复规则移出执行规范，PNJL 三维相图规则放入其案例合同，原文已存入源码快照。
SOP 治理检查保留 16 条通用章节标题建议；它们属于非强制的数值 SOP 模板。本绘图 SOP 使用六个执行章节，不增加收敛性计算、运行记录等无关章节，也未修改检查器或放宽质量阈值。

证据：[汇总](review_checks.json)、[旧文件保护](protected_files_check.json)、[SOP 审阅](sop_review.json)、[实际 PDF 渲染复核](pdf_render_review.json)、[预览缩略图](pdf_previews/contact_sheet.jpg)、[测试日志](tests.log)、[导出日志](pdf_export.log)。
原始 160 dpi PDF 渲染和精简预览保存在 `pdf_previews/`。PNG/PDF 后端的抗锯齿差异只作为诊断记录，不称逐像素 PDF 相同。

本次人工授权为 v13 PNG 版本级接受和对应 PDF 导出；不声称作者逐份审查了 75 份 PDF。`manuscript_eligible=false`、current=v5；本次不运行求解器、不改变数值或分支断线、不验证论文整页排版或物理打印。
