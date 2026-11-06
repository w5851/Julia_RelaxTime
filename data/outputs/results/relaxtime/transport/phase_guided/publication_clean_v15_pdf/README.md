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
v15渲染所依赖的源码、PNG、正式数值与验收依据继续保留。

复现入口：`scripts/analysis/relaxtime/export_phase_guided_publication_clean_v15_pdf.py`。
该入口拒绝覆盖case，要求绑定PNG的作者授权；若系统字体升级，使用
`font_sources.json`中经过SHA-256核对的同字节保留字体。字体文件不随Git分发。

旧图层、旧raw/显示表、孤立manifest及对应历史入口已整组退役，
恢复来源和最小保留范围见[专题入口](../../../../../../../docs/analysis/relaxtime/phase_guided_transport/README.md)。
v13对v5的历史引用不再要求当前分支保留旧图片；当前数值资格和Mott分析依据继续保留。
退役前版本为`6a3d08663c845587ae457fbe1c90173e0e5dcc25`，不另建全量备份或逐文件退役记录。
论文另清理279份旧附件，来源为paper提交`f8d1bae7377b7ec8976c998a5cba40550a70fdaf`。
当前论文已通过15页编译预览、41项引用和逐页视觉核验；新正文仍为agent_review。
