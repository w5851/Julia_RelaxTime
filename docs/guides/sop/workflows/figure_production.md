# 论文级绘图资产与生产 SOP

状态：`active`

版本：`figure_production_v2`

最后核验：2026-10-02

## 1. 目的与适用范围

本 SOP 规定如何从已经冻结、可追溯的 CSV/JSON 结果生成论文级图像资产。它统一图像的四层语义、视觉 profile、单位显示、状态表达、输出格式、布局审核和 provenance。

正式图像根目录为 `data/outputs/figures/`。每个新 case 必须写入新的 sibling 目录，并生成 `plot_manifest.json`。图形具体逻辑仍由图族脚本负责，公共 plotting 层只负责样式、合同和验证。

新的 APS/PRD candidate/strict case 采用两阶段交付：先生成 **600 dpi PNG 审查层**，作者审核通过后，再从同一输入 hash、同一生成器和同一绘图逻辑补全 **矢量 PDF**。矢量 PDF 是本项目用于 LaTeX 和最终交付的主文件；PNG 是审查和兼容备用文件。PNG 审查层必须在 manifest 中写明 `delivery_stage=png_review`、`manuscript_eligible=false` 和 `vector_delivery_pending=true`，不能被当作完整投稿交付。SVG 可用于内部编辑，不能默认为期刊交付格式；EPS/PS 按接收后生产要求或 color-online/grayscale-print 路线另行交付。旧 v1 profile 保留以复现历史合同，但已降为 deprecated/legacy，新 case 只能使用 v2 profile。

这是一项项目默认选择，不宣称 APS 强制单独图件为 PDF。2026-10-01 核对的 [APS 网页投稿指南](https://journals.aps.org/authors/web-submission-guidelines-physical-review)要求初投/重投整篇 PDF，接受后请求源文件；[APS Style Basics](https://journals.aps.org/authors/style-basics)仍列 PS/EPS/JPG/PNG 为优选图件格式，并为在线彩色、印刷灰度路线规定 PS/EPS 生产文件。

## 2. 非适用范围

- 不修改 solver、Maxwell、C2、reference、transport 或其他数值语义。
- 不重跑 PNJL，不重生正式 CSV/JSON，不批量重绘、覆盖、改名或删除历史图。
- 不把 `docs/analysis` 诊断图自动升格为正式论文图。
- 不把 `estimated_midpoint` 称为 confirmed CEP。
- 不在 Origin 中改数据、重采样、插值、补点、改单位或改变 marker 的物理含义。
- 不把所有历史绘图脚本强行迁移到单一绘图框架。

## 3. 权威入口

公共合同和验证入口如下。v2 是当前唯一的新图 profile；v1 仅用于历史 manifest/图像复现，不得作为新 case 的默认入口：

- `config/plotting/candidate_aps_v2.toml`
- `config/plotting/strict_aps_v2.toml`
- `scripts/plotting/plot_style.py`
- `scripts/plotting/plot_manifest.py`
- `scripts/plotting/plot_quality.py`
- `scripts/plotting/validate_plot_artifact.py`

图族脚本是数据选择和图形语义的局部权威入口。`scripts/plotting/render_plotting_pilot.py` 仅是历史 v1 compatibility pilot，不替代各正式图族生成器，也不调用数值求解器；新图族应直接使用 v2 profile 和自己的生成器。

## 4. 物理口径、单位与参数约束

- 输入字段、源单位、显示单位和转换公式必须写入 manifest。
- `mu_q` 与 `mu_B` 不得混用；若使用 `mu_B = 3 mu_q`，必须在 axes transform 中明确记录。
- MeV、GeV、fm^-1 和 dimensionless 量必须显式标注；禁止仅凭列名猜单位。
- 数学标题中的数值和单位必须保留可见间距，例如
  `$\mu_B = 450\;\mathrm{MeV}$`，不能依赖过窄的隐式 mathtext 间距。
- 论文轴标签使用数学变量和正体单位，单位置于圆括号中，例如 `$\tau_u\;(\mathrm{fm})$`、`$\mu_B\;(\mathrm{MeV})$`；新论文图禁止方括号单位。保留反夸克横杠及其他物理下标。
- `first_order`、`crossover`、`spinodal`、`cep_confirmed`、`cep_bracket`、`estimated_midpoint`、`unresolved` 和 `nonconverged` 是不同语义，不能只靠颜色区分。
- 缺失 support、失败点和 unresolved 区域必须断线、mask 或排除；不得隐式插值跨越 gap。
- strict 只接受输入侧已经确认、有限、无重复键且 support 合格的 series。

### 4.1 PNJL 三维相图的物理筛选

PNJL 的三维相图必须先做物理状态筛选，再做三角化和视觉排版：

- 原始序参量偏导峰（response peak）只是数值候选，不自动等于 crossover；只有位于同一 `xi` 切片的 CEP 化学势侧、且不落入 Maxwell 一阶区的峰，才可绘制为物理 crossover 面。
- 同一 `(xi, mu_q)` 不能同时标记为 crossover 和 Maxwell。若 `mu_q > mu_CEP`，该处即使仍有偏导峰，也只保留在筛选表/诊断数据中，不绘制为 crossover；该区域的物理面只保留 Maxwell。
- CEP 只有在 strict gate 通过后才能作为 confirmed endpoint。只有温度 bracket 时，图上使用 bracket 或 `estimated_midpoint` 语义，不能把中点写成单值 CEP。
- 高于 CEP 的响应峰不得用灰色叉号混入物理 crossover 面。若为 audit 需要展示，必须置于独立 diagnostic overlay；论文候选图直接排除。
- 端点附近若相邻原生采样点跨过 CEP 筛选边界，保留采样 gap 并在 manifest/派生表记录；不得跨 gap 隐式补线或把两侧面误连为一个面。
- 仅用于作者视觉判断的 `visualization-only closed` 模式可以把 finite/converged 但 geometry/interpolation 未闭合的 Maxwell 行统一绘成同一颜色，并在 manifest 明确声明 display-only 三角化上限；这只是显示连接，不生成缺失数据、不放宽门禁，也不能进入 strict 或 phase-reference promotion。
- 若需要判断空洞是否由三角网格造成，应使用 `diagnostic_no_triangulation` v5：Maxwell/crossover 只绘制原生有序 support 的相邻线段，超过采样门限的 gap 保持断开；该模式不填面、不生成合成点，也不把 unresolved 诊断升级为 Maxwell boundary。

## 5. 输入配置及优先级

绘图输入优先级为：

1. 结果侧 CSV/JSON 及其 source/calculation manifest；
2. 图族脚本声明的字段、筛选、排序和 mask 规则；
3. `figure_mode` 对应的物理资格 gate；
4. `style_profile` 中的公共尺寸、字体、线宽、marker、输出和布局默认值；
5. case-specific 的合法布局覆盖，例如 double-column、external legend，或经人工审查的 in-axes legend。

任何覆盖都必须写入 `plot_manifest.json`。绘图脚本不得把隐藏的单位转换、插值、connector 行或失败点修复放在默认分支中。

## 6. 环境与版本冻结

生成器必须在 manifest 中记录 Python、Matplotlib、平台、可执行文件和解析后的字体。当前 line-first pilot 已在 Python 3.13.2、Matplotlib 3.10.1 和 Times New Roman 解析结果下通过验证；其他机器必须重新记录实际解析结果。

新的 APS v2 profile 统一轴标签、主刻度、图例和标题的基准字号为 13 pt。这是包含 mathtext 上下标在内的保守默认，不是期刊强制字号。实际字体、字号和最终缩放必须通过字形测量；可以用明确的 case 覆盖，但不能用缩小字号解决遮挡。PNG 审查阶段仍必须完成尺寸、字形、刻度、裁切和遮挡测量；只有矢量阶段才执行 PDF 字体嵌入、无栅格包裹等 PDF 专属检查。`pdfinfo`、`pdffonts`、`pdfimages` 是 v2 矢量文件验证依赖，缺失时不得宣称 PDF 验证完成。

样式 profile 的 APS-like 基线为：single-column `3.375 x 2.5 in`，double-column `6.75 x 4.6 in`。宽度作为版式基线，高度是项目当前选择，不宣称为所有 APS 期刊的普遍硬性值。

多面板可以按内容调整高度，但必须记录 `size_override_reason`。新图按最终插入宽度直接绘制，PDF/PNG 保持固定画布，不使用 `bbox_inches="tight"` 静默改变导出宽度。目标宽度原则上不超过项目 APS profile 的 7 in 上限；不同论文模板的实际宽度在交接时再次核验。

## 7. Smoke 预检

在生成图像前执行：

1. 确认目标 sibling 目录不存在；存在时停止，不覆盖。
2. 确认所有 CSV/JSON 输入存在，记录 bytes 和 SHA-256。
3. 检查字段、单位、有限值、重复 key、排序键和失败状态。
4. 确认 profile 可加载，输出格式和尺寸合法。
5. 确认生成器是后处理脚本，不会调用 solver 或产生新的数值输入。
6. 为 strict 预登记 layout policy、legend 位置和是否使用外置 legend。

## 8. 收敛性验证

绘图层不重新证明数值收敛。它必须消费结果侧已经记录的 convergence/support 证据，并把选择规则写入 manifest。

- `audit` 可以展示失败、unresolved、support gap、residual 和 mask。
- `estimated_midpoint` 只能在输入提供明确 bracket/上下界时生成，并记录 midpoint 计算规则。
- `strict` 拒绝 unresolved、nonconverged、未确认 CEP、bracket-only endpoint、隐式 interpolation、外推和 connector。
- `visualization-only closed` 只允许用于诊断全局拓扑：finite/converged 的 Maxwell 行可统一着色，但原始 unresolved 状态必须保留在表格和 manifest 中，且不得被解释为证书通过。
- literature comparison 的模型插值只允许留在 `audit`/`legacy`；strict 只画原始模型 support 点。

## 9. 正式计算命令

本 SOP 不提供数值计算命令。正式图像只从冻结结果生成，不启动 PNJL、Maxwell、C2 或 transport 计算。

代表性 pilot 可使用：

```powershell
python scripts/plotting/render_plotting_pilot.py --only all --suffix __new_review
```

该命令只读取既有 CSV/JSON，并在新 sibling 目录写图。正式图族应使用自己的生成器，但必须复用同一 profile/manifest/validator 合同。

## 10. 输出目录与产物合同

目录命名约定为：

```text
data/outputs/figures/<domain>/<figure_family>/<case_slug>__plotv1__audit/
data/outputs/figures/<domain>/<figure_family>/<case_slug>__plotv1__estimated_midpoint/
data/outputs/figures/<domain>/<figure_family>/<case_slug>__plotv1__strict/
```

四层语义如下：

| mode | 资格和视觉规则 |
| --- | --- |
| `audit` | 内部审计；可以显示 raw support、失败点、unresolved、mask、bracket 和旧 connector，但必须标明 audit。 |
| `estimated_midpoint` | supplement/内部 review；可以显示 bracket 和 midpoint，但不能称为 confirmed CEP。 |
| `strict` | 正文定稿候选；只含已确认、有限、support 合格的物理 series；新的 APS v2 case 默认矢量 PDF + 600 dpi PNG。 |
| `legacy` | 历史兼容；保持原图、原脚本、原单位、原 connector 和原输出语义，不自动迁移。 |

`visualization-only closed` 不是第五种物理状态，而是 `audit` 的一个显示子模式；它只统一 Maxwell 的视觉颜色，不能改变四层语义或晋升资格。

所有新图必须有 `plot_manifest.json`，完整单图记录（多图包展开后）至少包含：`figure_mode`、`style_profile`、输入 hash、generator/hash、Git commit、calculation/postprocess/source provenance、axes 单位和 transform、series state、support/mask 规则、interpolation/connector policy、输出 hash、DPI/vector 和 layout 记录。

每个新图包只交付一份 `plot_manifest.json`。单图可使用 `plot_manifest_v1`；
多图使用 `plot_manifest_bundle_v1`，在 `shared` 中只记录一次共同输入、生成器、
运行环境和 provenance，在 `figures[]` 中以唯一 `figure_id` 保存每张单图及复合图
的轴、series、筛选/断线规则、布局、质量检查和各格式输出 hash。每条记录必须能与
共享字段合并为完整的单图合同；不能只保留文件名/hash 索引，也不能用逐图覆盖共享字段。
独立 sidecar 不是必要的复现信息，默认不生成 `*.plot_manifest.json` 或额外格式 index。
公共 validator 逐记录检查，图族专用 review/preflight 继续保留其资格限制。
新论文图族必须调用公共 `configure_matplotlib`、`configure_axis_ticks`、`export_figure`
和 validator；历史脚本只用于原合同复现，不因命名为 publication 而自动合规。

已接受图包的存储迁移须经作者授权：先保留原 manifest 图谱的逐字节归档和 hash，
再证明总 manifest 可无损还原全部记录，图像、数值表、原验收记录、资格和 current
指针均不改变。旧验收引用只可解析到已登记且 hash 一致的原 manifest 归档；
图像和数值数据仍检查现存字节，不允许使用历史归档掩盖漂移。

### 两阶段交付合同

1. `png_review` 阶段：从冻结输入生成 PNG；使用 `candidate_aps_v2` 的尺寸、字体、单位、刻度和质量合同，但输出 manifest 只声明 PNG，且必须记录 `delivery_stage=png_review`、`manuscript_eligible=false`、`current_publication_layer=false` 和 `vector_delivery_pending=true`。作者审查的对象是这一阶段的 PNG 与其 manifest/hash。
2. `vector_delivery` 阶段：作者明确审核通过后，使用同一输入记录、同一 generator/hash、同一参数和同一绘图逻辑补全 PDF（必要时再补 EPS/PS）。不得在补矢量时重新选择数据、改变平滑/断线/插值规则或重新调用 solver；PDF 必须通过字体嵌入、无栅格包裹、单页和物理尺寸检查。
3. 两阶段都不得自动更新 `publication_clean_current.json` 或把 display derivative 晋升为 manuscript eligible。晋升属于单独的作者接受/正式化动作。

## 11. Regression / Validation 验收

每个新 case 验证本次输入、总 manifest 内全部逐图记录和实际输出，并完成适用的视觉审核：

```powershell
python scripts/plotting/validate_plot_artifact.py <path-to-plot_manifest.json>
```

公共绘图框架、profile、validator 或质量规则变化时，运行受影响的公共测试：

```powershell
python -m pytest -q tests/unit/python/test_plotting_contract.py
python -m pytest -q tests/unit/python/test_plotting_quality.py
```

图族选择、轴类型、单位或布局逻辑变化时选择该图族行为测试。已有生成器的新参数/case
无需重复公共测试。文档链接、SOP registry 或脚本入口变化时才运行对应治理检查；
单纯生成新图不附带全部 Julia 文档与入口检查。

历史快照核验使用 `validate_plot_artifact.py --snapshot --code-ref <full-sha>`，
或已登记图族的专用检查入口。它检查冻结 manifest graph 与数据/输出完整性，
不自动重跑历史 renderer，也不宣称满足新的视觉或数值资格。

新的 APS v2 case 必须检查 PNG 元数据和最终插入宽度下的有效 DPI 均不低于 600，PDF 为单页矢量曲线图、字体嵌入且非 Type 3，导出物理尺寸与测量一致。不能只凭 `.pdf` 扩展名、`dpi=600` 或 manifest 的 `vector=true` 判定矢量输出。输入 hash、forbidden state、跨 gap 连接等 strict 资格要求继续保留。candidate/review 同样执行样式检查，但不因此获得数值资格。

## 12. 失败点、断点续算与重跑

绘图失败时保留输入和失败原因，不修改源 CSV/JSON。若目标目录已经产生部分输出，下一次运行必须使用新的 sibling suffix，或由明确的开发者操作清理未完成的临时目录；不得覆盖已通过验证的 case。

输入和输出 hash、profile 或 source run 不一致时，停止 strict 生成并退回 audit/作者审核。历史
manifest 的 generator/plotting-contract 记录使用 manifest 中的 `git_commit`、登记的
`code_ref` 或保留的 `source_snapshot` 核对；当前脚本或 SOP 后续修改不要求重写已经冻结的
历史 manifest。数据、图像和 manifest 文件本身仍按记录的字节/hash 校验。重跑绘图不等于重跑
数值计算，必须在新 manifest 中记录新的 generator/output hash。

## 13. Diagnostic 与 Formal Production 的边界

`docs/analysis` 是诊断证据区域，允许 C1 unresolved、bracket 和 estimated midpoint。它不能因为图形变得整洁就自动进入 `data/outputs/figures` 的 strict 正式目录。

`strict` 是图像合同通过，不是数值结论升级。只有 source result 的生产、收敛和物理 gate 也已通过，strict 图才可以申请论文定稿。历史 `legacy` 图不因新 profile 存在而失效，也不因新 SOP 自动重绘。

## 14. 后处理与作图

公共视觉规则采用 `line-first / landmark-only`：

- audit 显示 support 点；candidate/strict 默认隐藏普通 support marker；
- confirmed CEP 使用稀疏、醒目的实心圆点；
- `crossover` 使用短虚线，first-order 使用实线，spinodal 使用点线（低 alpha 仅保留在历史 v1 profile）；
- estimated/bracket 优先使用辅助线、开放 landmark 或三维加粗 envelope line；
- support 点可以不画，但其数量、筛选和 mask 规则必须留在 manifest；
- Origin 只能进行最终 panel、字体和尺寸排版，输入 hash 和输出 hash 必须可回溯。

### Strict layout gate

strict 不把 single-column 尺寸和 legend 位置视为不可变硬编码。数据密集时必须在目标尺寸下检查：

1. legend 不遮挡主曲线、CEP、关键边界或误差区域；
2. legend 不占据不合理比例的绘图区，文字不发生裁切或替换；
3. 普通 support 不通过 marker 增加视觉噪声；
4. 若 single-column 不足，优先采用更紧凑的 label、double-column、外置 legend 或独立 legend panel；不得只把字号压到不可读；
5. 最终选择写入 manifest 的 `rendering.legend_policy`、`legend_location`、`legend_outside` 和 column 字段，并由人工视觉审核确认。

`strict_origin_like_v1` 的默认布局策略是 `dense_aware_best_then_external`，允许外置 legend；本次 meson pilot 的两条曲线在 single-column 内部 legend 下通过，但这不替代其他密集图族的 case-level 审核。

默认仍优先使用图外公共图例。只有在密集多面板图族已经通过 case-level
布局审查、且 manifest 明确记录 `legend_policy=shared_in_first_panel_reviewed`、
`shared_in_top_row_panel_reviewed`、`best_in_axes_reviewed` 或
`shared_in_panel_reviewed_geometry_checked` 时，才允许把一次公共图例放入某个 axes。此例外
必须同时记录图例位置、端点 marker 颜色策略、曲线可见性审核和
`legend_outside=false`；validator 不会因此把它误判为默认外置图例。

图例 host、端点标签、字号和 caption 参数映射由图族案例合同规定。
phase-guided transport 的已审查细节见
[案例合同](../../../analysis/relaxtime/phase_guided_transport/plotting_case_contract.md)，
不作为其他图族的默认布局。

### APS v2 最终尺寸检查

1. 测量可见大写字母和数字的字体轮廓高度（包括数学上下标），在预定插入尺寸原则上不低于 2 mm。`1 mm = 72/25.4 = 2.835 pt` 仅是长度换算，字号框/em 高度不等于大写字高。[APS Style Basics](https://journals.aps.org/authors/style-basics) 确实写明最小大写字母和数字的 2 mm 要求，不能把同门的经验替代为期刊规则。只有 `audit`/`png_review` 的密集复合图可以登记 `typography_exception=dense_composite_review_compact_typography`，将紧凑排版的实测下限暂放宽到 1.5 mm；同时记录最小主文本和数学上下标字高。该项目审查例外必须通过无裁切/无文字重叠/无曲线和端点相交检查，不代表 APS 合规，不能进入 `strict`、`vector_delivery` 或 `manuscript_eligible=true` 交付。旧名 `dense_composite_review_compact_legend` 仅兼容此前审查记录。
2. 曲线最终线宽不低于 0.5 pt，landmark 直径不低于 1 mm。参数曲线使用可访问调色板加不同线型；phase 语义不变。v2 spinodal 默认不降低透明度，灰度区分通过线型实现。
3. 四侧 inward 主/小刻度；线性轴默认在相邻主刻度之间只放一个小刻度，log 轴使用对应 log locator。网格、统一 y 范围、log 变换均按科学目标选择，不能作为所有图的强制美化项。若目标图族需要更密的小刻度，必须在 case manifest 中说明理由。
4. 多面板优先使用图外公共图例；对于已按论文/组内博士论文版式审查的密集输运复合图，可以把一次公共图例放在左上角的空白 panel 内。图例 host 必须依据实际曲线占用选择，并记录曲线可见性检查；禁止把图例放在蓝色上升段或其他关键曲线区域。每列只保留一次参数标题、每行一次量纲轴标签、底行显示 x 标签。只在语义和范围一致时共享坐标。独立 y 范围必须在 caption 说明。
5. 默认公共图例不与 axes 相交；若采用上述已审查的 in-axes 例外，必须由 manifest 明确声明并由人工核对曲线、颜色、线型、marker、裁切和文字重叠。validator 只接受已登记的例外，不允许脚本静默绕过检查。
6. 复合图从冻结点表在最终物理尺寸绘制，保持矢量曲线和文本；不得把大单图 PNG 拼接、缩小后再包成 PDF 作为矢量交付。
7. 交接 paper 时记录插入宽度、最小可接受缩放和完整 caption 参数映射；改了 LaTeX 缩放后重新测量。`undecided_review` 表示作者尚未选择彩印/灰印路线，strict 生产交付必须选定。

图族的布局、端点措辞、刻度精度和 log-y 选择放在对应案例文档；
phase-guided transport v11 见
[案例合同](../../../analysis/relaxtime/phase_guided_transport/plotting_case_contract.md)。

## 15. 关联公式、API 和测试

- 样式与 manifest：`scripts/plotting/plot_style.py`、`scripts/plotting/plot_manifest.py`。
- 图像合同验证：`scripts/plotting/validate_plot_artifact.py`。
- 代表性合同测试：`tests/unit/python/test_plotting_contract.py`。
- 科学计算生命周期：`docs/guides/sop/common_scientific_run.md`。
- 图族物理公式、模型 API 和 numerical validation 仍由各自 `docs/reference/`、`docs/api/` 和专题 SOP 权威管理。

## 16. 最后验证记录

2026-08-15 line-first pilot 已由作者完成视觉审核并通过。验证证据：

- Python contract tests：`5 passed`；
- audit、strict、estimated_midpoint 三份 pilot manifest 均通过 `validate_plot_artifact.py`；
- `check_sop_governance.jl`、`check_docs_consistency.jl`、`check_active_docs_governance.jl` 和 `check_script_entrypoints.jl` 通过；
- 历史 strict pilot 输出为 SVG + 600 dpi PNG（v1 合同）；
- 未修改 solver、Maxwell、C2、reference、transport、正式 CSV/JSON 或历史 PNG/PDF/SVG；
- strict 的后续密集图族仍必须执行本节 14 的 layout gate，不能仅因 profile 默认值而跳过人工审核。

## 17. 历史图资产退役与清理

历史 PNG/PDF/SVG 的清理采用 `asset inventory + dry-run + allowlist cleanup`，不是按扩展名、文件名或时间批量删除。

PR A 只运行 `scripts/plotting/inventory_figure_assets.py`，默认扫描 Git 已跟踪的 `data/outputs/figures` 资产，生成：

- `docs/analysis/governance/figure_asset_registry_v1/asset_registry.json`；
- `docs/analysis/governance/figure_asset_registry_v1/cleanup_candidates.csv`。

未跟踪的 C1/C2/pilot 文件默认排除且不修改。`docs/analysis` 诊断证据与正式图像根目录分开治理，不因图形格式相同而自动合并。

registry 只提出 `owner_review_only` 或 `keep_contract_case`，不包含删除、移动和覆盖操作。人工审核必须确认仓库外部引用、canonical case/variant 和历史证据保留策略；未确认项默认保留。

实际退役属于后续 PR B，必须以作者批准的 `path + sha256 + action` allowlist 执行，并再次检查引用、manifest 输出和文档链接。strict 新默认 PDF + 600 dpi PNG 不追溯改变历史 PDF/SVG/PNG 的保留资格。

## 18. 历史验证与新图族迁移

2026-08-15 的 `__plotv1__` pilot 是历史作者审核记录。其未跟踪 sibling 文件不在当前 Git checkout 中时，不能声称仓库内可重新核验；脚本可用新 suffix 重建，本次不会补写旧审核结果。

输运 publication-clean v1–v6 使用独立 renderer，未调用公共 profile/validator，归为保留的历史显示合同；600 dpi PNG 并不等于遵循本 SOP。v7 使用 `candidate_aps_v2` 迁移到 v2 合同，v8 在同一 v5 冻结点表上采用 in-axes legend、5 个 x 轴主刻度、线性轴 `MaxNLocator(nbins=5)`、曲线边距和正值高动态范围 panel 的局部 log-y。v9 在同一冻结点表上收紧为三个 x 轴数字刻度、线性轴一个小刻度、显式单位间距，并将 Figure 1 的 `mu_B=900 MeV` 整列统一为 log-y；仍残留图例遮挡，不是定稿。v10 将参数 key 与带 `1st-order` 标题的端点 key 分置 (a)/(c)，统一收紧复合图字体和留白、采用精确三位小数刻度，并加入 legend–curve / landmark 几何相交检查。v7–v10 都不因 PNG 审查通过而自动获得矢量或投稿资格；v10 的紧凑字形例外仅用于 PNG 审查，strict 交付仍需单独解决最终字高。既有 v7 PDF+PNG 目录保留为迁移前 review artifact，不覆盖、不晋升。v8–v10 保持 `manuscript_eligible=false`、current=v5，不重算、不新增数值、不改 source/raw provenance。

v11 延续 v10 的冻结数值、布局和字形审查边界，只统一弛豫时间图的 log-y、
log 数字标签、线性行小数位，以及 `First-order` 图例标题和 caption 适用范围。
它使用独立 PNG-review sibling，不覆盖 v10，不更新 current=v5，不补矢量文件。

Figure 4 等旧资产的超宽/低 DPI 不追溯修改；若将来进入正文，必须另建 case 重新出图并执行最终尺寸合同。新的 v2 默认不改变旧资产的保留、资格或 hash。
