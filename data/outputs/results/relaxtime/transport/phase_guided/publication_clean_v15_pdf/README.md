# v15 正式输运图层

2026-11-05 按作者授权完成 PDF 检查、正式晋升和 current 替换。
[正式化记录](formalization.json)绑定数值资格、授权和视觉复核；
[技术核验](delivery_validation.json)覆盖全部76张PDF。

输入为[η加密正式网格](../../threshold_cache_upgrade_20261021/full_rs_eta_augmented_20261029_v1/README.md)：
1867点，1810个保留值、9个原坐标重算值和48个新增点；不插值、不平滑。
76张彩色和76张灰度PNG直接引用已接受的v15审阅包，全部重绘字节相同，不重复存储。
PDF均为单页、嵌入非Type-3字体、保留矢量曲线和已接受的物理尺寸。
三张论文组合图已从本批PDF重新渲染，经辅助审阅、主代理视觉及几何核对通过。

图件位于 `data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v15_pdf/`。
论文只接入三张组合PDF、两份主表及最小来源记录，独立相图不变；
论文当前科学与文字审阅状态见其 `writing_package/00_status.md`。
旧图层的冻结manifest不回写，旧证据恢复链见工程交付；v15渲染所依赖的源码及PNG继续保留。

复现入口：`scripts/analysis/relaxtime/export_phase_guided_publication_clean_v15_pdf.py`。
该入口拒绝覆盖case，要求绑定PNG的作者授权；若系统字体升级，使用
`font_sources.json`中经过SHA-256核对的同字节保留字体。字体文件不随Git分发。

旧资产退役见[retirement.json](retirement.json)：557张已被替代的PNG退出常用目录，
原字节可从已合并提交`bad5f98ce99ffd81a9ff0f60336a6c3c32a3d44a`恢复，未新增全量归档。
74张仍被v13冻结合同引用的v5图片保留；旧数值表、收敛依据和v10–v15输入继续支持现存复现入口。
论文另清理279份旧附件，来源为paper提交`f8d1bae7377b7ec8976c998a5cba40550a70fdaf`。
当前论文已通过15页编译预览、41项引用和逐页视觉核验；新正文仍为agent_review。
