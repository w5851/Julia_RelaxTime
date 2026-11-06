# Phase-guided transport 绘图合同

当前正式图层是 [publication_clean_v15](publication_clean_current.json)。
公共规则见[绘图 SOP](../../../guides/sop/workflows/figure_production.md)；逐图面板、
输入行、图例、尺寸、字体和输出字节以 v15 manifest 及冻结源码为准。
旧合同的原始字节已在 v15 的 `source/` 中保存，历史资格不由本页重写。

## 数值与曲线

来源为 `full_rs_eta_augmented_20261029_v1/returned/grid` 的正式 1867 点，
mode A/B 分别为 942/925 点；1810 个前驱原值、9 个原坐标重算值和 48 个新增 ξ 点。
读取器核对 package、晋升、完整网格审查及逐点来源，直接读取 12 个同名原始列，
包括已计算的比值列。不平滑、插值、删点、重算比值或掩盖有限质量提示。

只在已有背景审计确认的两处一阶分支间断线：

- mode A：`muB900.0/alpha1.0`，ξ = −0.003/+0.003。
- mode B：`T120.0/muB900.0`，ξ = −0.13/−0.12。

圆圈和方块表示实际采样的手征恢复/破缺侧点，不是外推的共存端点。
新增点的 unknown 相态和缺失迭代数不生成新断线；Mott 或局部尖峰也不自动生成断线。
Matplotlib 路径简化关闭，实际曲线顶点逐图与输入行核对。

## 布局

v15 复用已接受的 v14 renderer：72 张单图、Figure 1/2/4 三张组合图和一张低值细节图。
每幅有彩色、同尺寸灰度 PNG 和单页矢量 PDF。参数身份跨面板保持一致，
蓝/红/绿色为固定参数编码，线型依次为实线、虚线、点划线。

- 原生宽度 6.75 in（171.45 mm）。Figure 1/2 高度分别为 7.1/5.9 in。
  Figure 1 使用 log-y，Figure 2 使用 linear-y。
- Figure 4 为 3×3 固定 T 图：列为 T = 120/160/200 MeV，行为 η/s、ζ/s、σ/T。
  前两行 log-y，第三行 linear-y；参数 key 在 (g)，一阶端点 key 在 (a)，尺寸 6.75×7.4 in。
- 低值细节图对应 Figure 1 的 (a)/(g)/(h)/(i)，线性窗口为
  0.57–0.66、0.57–0.66、0.375–0.43、0.27–0.41 fm。窗口外数据仍保留，视窗不代表物理端点。
- 每份图例分别声明位置和适用曲线，并检查文字、曲线、端点遮挡、图内溢出与裁切。
  使用 strict APS v2 的字号、线宽、marker 和四侧向内刻度要求；实际插入宽度遵守各图 `placement_limits`。

## 复现与正式交付

PNG 入口为 `scripts/analysis/relaxtime/build_phase_guided_publication_clean_v15.py --png-review --case <new-case>`。
PDF 入口为 `scripts/analysis/relaxtime/export_phase_guided_publication_clean_v15_pdf.py`，
提供 `--authorization <record.json>` 和新的 `--case`；已有 case 拒绝覆盖。
作者授权必须绑定实际 PNG 包；PDF 重绘时核对渲染源码、全部原始顶点与 PNG 字节，
检查嵌入的非 Type-3 字体和真实矢量曲线。PDF 直接引用已有 PNG，不再复制。
`--font-sources <map.json>` 仅用于同字节字体换址。

2026-11-05 已完成 v15 PDF、正式晋升、current 和论文图替换，见
[正式交付](../../../../data/outputs/results/relaxtime/transport/phase_guided/publication_clean_v15_pdf/README.md)。
原 PNG manifest 的 review 阶段状态是历史事实，当前资格由正式化记录给出。

v1–v13 生成器与旧 raw/显示包已退役。v14/raw 模块保留为 v15 的渲染依赖，
历史 v14/raw 独立出图需在退役前 checkout 恢复对应参考包。
v14 版式 manifest 是现存 v15 PNG 入口的实际依赖，不要求保留它记录的整套旧图。
恢复版本与保留边界见[专题入口](README.md)。
