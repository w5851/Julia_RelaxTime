---
name: plotting-sop
description: 将仓库当前的 figure_production_v2 契约应用于新科研图像，包括冻结输入的来源追踪、PNG 审阅和后续矢量格式交付。不用于求解器或数值生产变更。
---

# 项目绘图 SOP

以已冻结的 CSV/JSON 结果为输入，按[图像生产工作流程](../../../docs/guides/sop/workflows/figure_production.md)创建或审阅论文用图。

## 当前契约

- 新 case 使用 `config/plotting/candidate_aps_v2.toml` 供作者审阅；只有输入资质和最终交付范围足以支持严格状态时，才使用 `config/plotting/strict_aps_v2.toml`。
- v1 配置（`audit_v1`、`candidate_origin_like_v1` 和 `strict_origin_like_v1`）已弃用，仅保留兼容用途。可加载它们来复现历史产物，但新 case 不得选择这些配置。
- 复用 `scripts/plotting/plot_style.py`、`plot_manifest.py`、`plot_quality.py` 和 `validate_plot_artifact.py`。图族脚本负责数据选择和物理标签，共享层负责渲染与产物契约。
- 绘图生产过程中不得调用求解器、修改数值 CSV/JSON 的值、引入隐藏的插值或平滑，或修改论文项目。
- 每个 case 交付一个 `plot_manifest.json`。多图 case 使用 `plot_manifest_bundle_v1`：共享来源信息只记录一次，各图使用唯一标识，并在 `figures[]` 中保留完整的逐图记录。无需另建独立附属记录文件。使用 `scripts/plotting/plot_bundle.py` 提取共享字段并展开记录。

## 两阶段交付

1. 创建新的同级输出目录。拒绝覆盖已有 case。记录每个输入路径、字节数、SHA-256、生成器 hash、Git 提交、运行环境、单位、变换、数据系列状态，以及掩码/插值规则。对保留的历史 case，使用所记录的完整提交或源码快照验证代码与绘图契约的 hash；按记录的字节数/hash 验证当前数据、图像和 manifest。
2. 先生成 `png_review` 阶段。图像格式必须仅为 PNG，并在 manifest 中记录 `delivery_stage=png_review`、`manuscript_eligible=false`、`current_publication_layer=false` 和 `vector_delivery_pending=true`。运行全部非 PDF 质量检查：物理尺寸、有效分辨率、实测字形高度、线宽、四边向内的主/次刻度、裁切、重叠、单位和输出 hash。
3. 停下来等待作者进行视觉审阅。将 PNG 与其 manifest 一同审阅；不得晋升产物或更新 `publication_clean_current.json`。
4. 获得明确接受后，使用相同的冻结输入和绘图逻辑生成 PDF。用 `plot_delivery.py` 校验 `plot_png_acceptance_v1` 接受记录、PNG/manifest hash、冻结输入和两阶段显示语义；公共模块不生成或推断作者批准。运行 PDF 专属检查：单页、嵌入的非 Type-3 字体、曲线未被封装成栅格图、固定物理尺寸和 hash 匹配。只有所选期刊交付路径要求时，才添加 EPS/PS。

## 视觉要求

- 使用数学变量标签，将有量纲的单位放在圆括号内，例如按图族的渲染标签约定使用 `$\tau_u\;(\mathrm{fm})$`；新的论文图不得引入方括号单位。
- 上、下、左、右四边均使用向内刻度。线性轴默认在相邻主刻度之间放置一个次刻度；对数轴使用配置的对数次刻度定位器。若线性轴采用更密的次刻度策略，必须在 manifest 中给出 case 级理由。
- 按预期插入文稿的物理宽度绘制。测量实际字形轮廓高度，而不是名义 em 字号；v2 检查要求实测大写字母/数字字形至少为 2 mm。仅供 PNG 审阅的密集组合图可以使用已明确记录的审阅排版例外，其实测上下标至少为 1.5 mm。这不代表符合 APS 投稿要求，也不能用于晋升产物或授权矢量交付。
- 优先尝试在坐标区未使用的空间放置共享或分组紧凑图例，使用公共标题减少重复标签。几何测量需要时可使用图上方、侧面或独立图例区域。新 v2 case 采用 `declared_geometry_checked`，声明每个图例的承载位置、内容与适用范围；测量包括同一坐标区保留的多个图例在内的所有图例，拒绝遮挡曲线/特征点、超出承载范围及文字/图例重叠。图例与坐标区相交本身不构成失败，几何检查通过不等于作者接受。需要在灰度下区分曲线时，使用独立形状、直接标签或线型。
- 按具体科学问题选择坐标轴尺度、范围、图例承载位置和精度。刻度标签必须标识实际数值；显示格式不能作为数值不确定度的依据。
- 参数用颜色与独立形状／文字／线型区分；相态已使用线型时保留其语义。新 case 声明 `rendering.review_contract=prd_review_v1`，用 `plot_accessibility.py` 记录实际邻近背景和关键曲线的 3∶1 对比度检查，并提供同尺寸灰度图；透明填色按合成背景检查。
- 有明确上下界的中点与经验密度估计分开分类。`estimated_density` 记录来源行、公式及非 bracket 的含义，闭合段单独记录；它们保持 audit 资格。交接时提供完整图注、参数／符号映射和最小插入宽度。
- 图族布局细节保留在对应 case 文档与渲染器中。各版本 phase-guided 布局记录于 `docs/analysis/relaxtime/phase_guided_transport/plotting_case_contract.md`。
- 科学相态标签与分支缺口应遵循输入语义；绘图风格变化不得重新标记强子/夸克结果，也不得重新连接已经审查确认的不连续处。

历史 manifest 是证据记录，不是要求使用当前渲染器重新运行的指令。对保留的快照契约，使用 `validate_plot_artifact.py --snapshot --code-ref <full-sha>`。
后续对现行 SOP 或生成器的修改，可能需要新建 case 或明确执行迁移，但不会使旧图像、输入或 manifest 的字节记录失效。

## 验证

验证 case manifest 中的每条记录及其实际输出。共享绘图代码、配置或校验器变化时，运行公共契约/质量测试；图族逻辑变化时，运行受影响的图族测试。文档与脚本入口检查随对应契约变化而执行。使用已有渲染器生成 case，不要求运行全部共享测试或 Julia 治理检查。

即使所有视觉检查均通过，所得图像产物包仍仅供审阅。正式化、稿件使用资格和当前图层替换，是由作者分别决定的操作。
