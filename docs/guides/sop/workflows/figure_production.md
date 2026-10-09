# 论文级绘图 SOP

状态：`active` · 合同：`figure_production_v2`

## 1. 适用范围与入口

从已冻结的 CSV/JSON 生成可追溯图件，按 **PNG 审阅 → 作者接受 → 矢量 PDF**
交付。本流程只做后处理；数值计算、物理资格晋升、论文装配和历史资产清理另行处理。
运行记录、测试结果和图版演变放在 case 证据包，本页只维护执行规则。

| 职责 | 权威入口 |
| --- | --- |
| 公共样式 | `config/plotting/candidate_aps_v2.toml`、`strict_aps_v2.toml` |
| 样式、导出、测量 | `scripts/plotting/plot_style.py`、`plot_quality.py` |
| 可访问性、两阶段接受链 | `scripts/plotting/plot_accessibility.py`、`plot_delivery.py` |
| 单图／多图合同、来源核验 | `plot_manifest.py`、`plot_bundle.py`、`plot_provenance.py`（均在 `scripts/plotting/`） |
| 产物验证 | `scripts/plotting/validate_plot_artifact.py` |
| 图族数据选择与标签 | 对应生成器及案例合同；入口见 [脚本索引](../../scripts/README.md) |

新图只使用 v2 profile。输入资格先于样式：

| `figure_mode` | 允许的内容与边界 |
| --- | --- |
| `audit` | 可展示失败、unresolved、mask、bracket、raw support 或已披露的旧显示处理；不授予正文资格。 |
| `estimated_midpoint` | 仅 supplement／内部审阅；输入须提供 bracket，声明 midpoint 规则，不称 confirmed CEP。 |
| `strict` | 仅已确认、有限、无重复键且 support 合格的原始物理 series；不插值、外推或添加 connector。 |
| `legacy` | 仅复现历史合同；不作为新图入口，不自动迁移或重绘。 |

`estimated_density` 是 `audit` 中独立的 series 状态，适用于作者明确授权的经验密度
估计，不与有已证明上下界的 `estimated_midpoint` 混用。当前支持
`nearest_subcritical_coexistence_mean`：记录估计坐标、最近亚临界原生行的温度、
两支密度及均值，声明 `interval_semantics=coexistence_pair_not_cep_bracket` 和
`creates_numerical_source_row=false`。源表 CEP 温度／化学势与估计密度分开记录；
闭合段单列为 connector。两者都不能进入 strict，也不构成新数值或误差区间。

## 2. 冻结输入与预检

1. 选择新的 sibling case 目录，图像位于 `data/outputs/figures/<domain>/<figure_family>/`。
   figures case 只保留 PNG/SVG/PDF 和唯一的 `plot_manifest.json`；caption、验收记录、
   源码快照、归档及校验结果放在 `data/outputs/results/<domain>/<figure_family>/` 下
   同名 case，由图像 manifest 引用。迁移已接受包时先保留完整原始字节及 hash，
   再记录路径映射；不得重新渲染图像或把存储迁移解释为新的作者接受。
   已有目录不覆盖；失败重试也使用新目录。
2. 核对字段、有限值、重复键、排序、support、失败状态及来源资格。冻结输入路径、
   bytes、SHA-256、生成器与依赖源码／profile、Git 上下文、实际运行环境和解析字体。
3. 明确源单位、显示单位、转换公式、筛选、mask、分段及已有显示处理。保留全部应绘
   顶点和 gap；不重算、不新增平滑／插值，不修补失败点或改变 marker 物理含义。
4. 选择最终插入宽度、轴类型／范围、共享坐标范围、图例位置和输出阶段。
   所有 case 覆盖均写入 manifest；只在语义和范围一致时共享坐标。

物理状态、绘图质量、作者接受和论文资格分别记录。继承显示数据不等于重新证明
数值收敛；任何派生图都不得仅因排版合格而升级来源证据。

## 3. 排版与最终尺寸门槛

按最终插入尺寸直接绘制，使用固定画布；不以 `bbox_inches="tight"` 改变宽度，
不把缩小后的 PNG 拼接或包进 PDF 冒充矢量曲线。修改论文中的缩放后重新测量。

| 检查项 | 执行规则 |
| --- | --- |
| 画布 | 项目基线：单栏 `3.375 × 2.5 in`，双栏 `6.75 × 4.6 in`；宽度上限 `7 in`。按内容调整高度时记录 `size_override_reason`，实际版心以论文模板为准。 |
| 字体 | 基准 13 pt；测量实际大写／数字字形，包括数学上下标，最终高度至少 **2 mm**。名义字号不能代替字形高度。 |
| 曲线与标记 | 最终线宽至少 **0.5 pt**，landmark 的实际形状宽高至少 **1 mm**。参数须有颜色之外的区分方式，可用形状、直接标签或线型；相态已占用线型时，参数采用独立形状／文字。普通 support marker 按图族需求选择，默认 line-first／landmark-only。 |
| 刻度 | 四侧内向主／小刻度；线性轴默认每对主刻度间一个小刻度，log 轴用相应 locator。更密刻度须说明；数字必须代表真实刻度值，显示精度不代表数值误差。 |
| 标签 | 数学变量、正体单位、圆括号单位；保留反粒子横杠、上下标和数值／单位间距。`mu_q` 与 `mu_B` 不混用，转换须声明。 |
| 相与分支 | 区分 first-order、crossover、spinodal、confirmed CEP、bracket、unresolved 和 nonconverged；不只靠颜色表达。一阶标签同时交代转变和分支含义。 |
| 多面板 | 共用信息只出现一次：列参数、行量纲、底行 x 标签；独立范围在图注说明。网格、log 轴和统一 y 范围按科学目标选择。 |
| 缩放区间 | 用 `placement_limits` 计算字形／线宽／marker 的缩小下限及 PNG 像素／profile 的放大上限。无可用区间时改布局；双栏单图不得直接视为合格单栏图。 |

2 mm 字形、0.5 pt 曲线、1 mm 数据标记和圆括号单位依据
[APS Style Basics](https://journals.aps.org/authors/style-basics)；13 pt、四侧刻度、画布、
7 in 上限、600 dpi PNG 和 PDF 主交付是项目选择。目标期刊要求另行核对。

### 颜色与参数辨识

[PRD 作者指南](https://journals.aps.org/prd/authors) 要求文字或符号提供含义并引用
WCAG。本项目按 [WCAG 2.2 1.4.11](https://www.w3.org/TR/WCAG22/#non-text-contrast)
检查理解内容所需的曲线相对邻近背景至少 3∶1；这是对引用标准的执行口径，
不是 PRD 页面另列的数值条款。颜色不能是唯一参数编码。

新审查记录声明 `rendering.review_contract=prd_review_v1`，并用
`plot_accessibility.py` 写入 `rendering.accessibility`：参数颜色、独立标识、实际
邻近背景及对比度。校验器重算比值，拒绝缺失或不合格记录。透明填色须记录合成后
背景，不可只检查白底；默认调色板合格不意味着任意叠色合格。无额外物理含义的
多层填色优先简化。记录中的参数标识与实际图例／曲线对应仍须视觉复核。

仅 `audit/png_review` 密集复合图可显式登记紧凑字形例外，下限 1.5 mm；该例外
不适用于 `strict`、矢量交付或论文资格，也不能豁免裁切和遮挡检查。

### 图例

优先在图内空白区放一次共享 key，或把参数与端点 key 分置不同 panel。
共同变量／语义放标题，条目去掉重复文字；若图内空白不足，可选简短顶部、轴旁
或独立图例区域，记录理由，不靠不可读的小字解决遮挡。

- 新 case 使用 `legend_policy=declared_geometry_checked`。
- 在 `legend_placements` 逐份声明顺序、`in_axes/outside_axes`、真实 host、
  位置、标题、条目及科学适用范围；局部端点含义不得扩展到全部曲线。
- 实测覆盖全部可见 legend，包括同一 axes 中保留的第二份 key。声明与实测的
  数量、顺序、host、标题和条目一致；`legend_outside` 表示是否全部位于图外。
- 图例与 axes 相交本身允许；遮挡曲线／端点、相互覆盖、图内 key 越出 host、
  文字重叠或裁切均失败。几何 count 与明细一致，缺少证据也不能通过。

## 4. PNG 生成与作者审阅

1. 使用公共样式、刻度、导出和校验 helper，从冻结点表生成 **600 dpi PNG**；
   同时提供同尺寸灰度审查图，记录转换方法、源彩色 hash、像素与 DPI。
2. 检查实际尺寸、最终字形／线宽／marker、四侧刻度、全部图例和文字几何；
   逐图核对输出 hash、有效 DPI、单位、数据选择及分段。
3. 在 manifest 声明 `delivery_stage=png_review`、`manuscript_eligible=false`、
   `current_publication_layer=false`、`vector_delivery_pending=true`，本阶段仅含 PNG。
4. 作者结合 manifest 审图：最终宽度下标签可读、线型可追踪、端点可识别，
   关键低值变化没有因压缩而丢失。记录实际审阅范围和结果；机器通过不等于人工接受。

需要局部视图时，标明源图／panel、轴类型、范围及 viewport 截断含义，保留原顶点
和 gap。局部图补充全范围主图，不生成新数据；曲线离开窗口不代表物理终止。

## 5. 接受后补齐 PDF

1. 记录作者接受的 PNG、manifest 和 hash。再次核对冻结输入、绘图源码及参数；
   使用同一绘图逻辑补出 **单页矢量 PDF**，不重新选择数据或改变坐标、线型、图例、
   平滑、断线及插值语义。独立格式导出器与原 renderer 分别留存 hash。
   新 vector_delivery 使用 `plot_png_acceptance_v1` 接受记录，包含实际作者回复、
   `status=author_accepted`、`pdf_export_authorized=true`、`accepted_png` 和
   `accepted_manifest` 的路径／bytes／hash；manifest 的 `author_review` 引用该记录。
   公共 `plot_delivery.py` 检查接受状态、文件内容、PNG 字节一致、输入内容及角色、
   renderer hash 和两阶段坐标／series／估计／连接／图例／可访问性语义。缺失接受、
   pending 或显示漂移均失败，不以导出成功替代接受；本模块不会生成作者批准。
2. 核验实际 PDF 的页数、物理尺寸、字体嵌入、非 Type 3 字体及矢量曲线；
   线图不得含栅格包裹。把 PDF 渲染为图片复核布局、字体、裁切和 PNG/PDF 一致性。
3. 新阶段记录 `delivery_stage=vector_delivery`、`vector_delivery_pending=false`，
   引用已接受的 PNG 和新增 PDF；原 PNG 阶段 manifest 保持冻结。
4. 灰度 PNG 是审查材料，不替代印刷生产文件。strict 交付须选定彩色／印刷路线；
   选 `color_online_grayscale_print` 时按期刊要求补 EPS/PS。SVG 仅按内部编辑需求提供。

PDF 是本项目的矢量主文件选择；整篇投稿及接受后生产格式按
[目标期刊要求](https://journals.aps.org/authors/web-submission-guidelines-physical-review)核对。
PNG 接受与 PDF 完成交付均不自动修改 `publication_clean_current.json` 或授予
`manuscript_eligible=true`；这两项属于单独授权的晋升动作。

交接论文时同时提供完整 caption、参数／符号映射和 `placement_limits`。仅部分
坐标采用估计时，图例与 caption 指明估计量、来源行选择及显示连接；原生共存区间
不称为 CEP 误差范围或 spinodal 区域。按 APS 当前
[专门 AI 政策](https://journals.aps.org/authors/appropriate-use-ai-tools) 如实说明用于
图件生成的 AI 帮助，不代作者声明未完成的科学验证。

## 6. Manifest、验收与失败处理

每个图包只保留一份 `plot_manifest.json`。单图用 `plot_manifest_v1`；多图用
`plot_manifest_bundle_v1`，共同来源放 `shared`，逐图记录放唯一 `figure_id` 的
`figures[]`。展开后须有完整输入／输出 hash、生成器／环境、计算来源、轴单位与
转换、series 状态、support／mask、插值／连接规则、布局、质量和交付状态。
不生成冗余逐图 sidecar，不允许逐图覆盖共享字段。

```powershell
python scripts/plotting/validate_plot_artifact.py path/to/plot_manifest.json
```

- 验收覆盖总 manifest 的全部记录及实际输出；不能仅凭扩展名、声明的 DPI 或
  `vector=true` 判定通过。PDF 检查需可用的 `pdfinfo`、`pdffonts`、`pdfimages`。
- 公共样式／校验逻辑变化时运行对应 `test_plotting_contract.py`、
  `test_plotting_quality.py` 等公共测试；图族逻辑变化时运行该图族测试。
  文档／入口变化时才运行对应治理检查，普通出图不附带求解器回归或全仓库检查。
- 失败时保留输入与失败原因，在新目录重试；不改数值源、不放宽阈值或隐藏失败项。
- 历史代码／合同按固定 commit 或已登记、hash 匹配的源码快照验证；数据和图像
  按现存字节核对。历史核验不重跑当前 renderer，也不授予当前合同资格。
  新接受链和可访问性要求不追溯改写已冻结 v11/v13 等包；其检查继续使用冻结
  exporter／validator 或 `--snapshot`，不能将旧包缺少新字段解释为已通过新合同。
- 已接受包的存储迁移或资产清理须单独授权，保留原图谱／hash 并证明引用与记录可恢复。

图族专用规则：[输运案例](../../../analysis/relaxtime/phase_guided_transport/plotting_case_contract.md)、
[PNJL 三维相图案例](../../../analysis/pnjl/phase_surface_series/plotting_case_contract.md)。
