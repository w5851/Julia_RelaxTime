# 论文级绘图资产与生产 SOP

状态：`active`

版本：`figure_production_v1`

最后核验：2026-08-15

## 1. 目的与适用范围

本 SOP 规定如何从已经冻结、可追溯的 CSV/JSON 结果生成论文级图像资产。它统一图像的四层语义、视觉 profile、单位显示、状态表达、输出格式、布局审核和 provenance。

正式图像根目录为 `data/outputs/figures/`。每个新 case 必须写入新的 sibling 目录，并生成 `plot_manifest.json`。图形具体逻辑仍由图族脚本负责，公共 plotting 层只负责样式、合同和验证。

本 SOP 的默认交付是 `600 dpi PNG + SVG`。PNG 用于开发者和用户快速查看，SVG 用于论文矢量排版。PDF 只在投稿系统或外部版式流程明确需要时按 case 启用。

## 2. 非适用范围

- 不修改 solver、Maxwell、C2、reference、transport 或其他数值语义。
- 不重跑 PNJL，不重生正式 CSV/JSON，不批量重绘、覆盖、改名或删除历史图。
- 不把 `docs/analysis` 诊断图自动升格为正式论文图。
- 不把 `estimated_midpoint` 称为 confirmed CEP。
- 不在 Origin 中改数据、重采样、插值、补点、改单位或改变 marker 的物理含义。
- 不把所有历史绘图脚本强行迁移到单一绘图框架。

## 3. 权威入口

公共合同和验证入口如下：

- `config/plotting/audit_v1.toml`
- `config/plotting/candidate_origin_like_v1.toml`
- `config/plotting/strict_origin_like_v1.toml`
- `scripts/plotting/plot_style.py`
- `scripts/plotting/plot_manifest.py`
- `scripts/plotting/validate_plot_artifact.py`

图族脚本是数据选择和图形语义的局部权威入口。`scripts/plotting/render_plotting_pilot.py` 是代表性 pilot 生成器，不替代各正式图族生成器，也不调用数值求解器。

## 4. 物理口径、单位与参数约束

- 输入字段、源单位、显示单位和转换公式必须写入 manifest。
- `mu_q` 与 `mu_B` 不得混用；若使用 `mu_B = 3 mu_q`，必须在 axes transform 中明确记录。
- MeV、GeV、fm^-1 和 dimensionless 量必须显式标注；禁止仅凭列名猜单位。
- `first_order`、`crossover`、`spinodal`、`cep_confirmed`、`cep_bracket`、`estimated_midpoint`、`unresolved` 和 `nonconverged` 是不同语义，不能只靠颜色区分。
- 缺失 support、失败点和 unresolved 区域必须断线、mask 或排除；不得隐式插值跨越 gap。
- strict 只接受输入侧已经确认、有限、无重复键且 support 合格的 series。

## 5. 输入配置及优先级

绘图输入优先级为：

1. 结果侧 CSV/JSON 及其 source/calculation manifest；
2. 图族脚本声明的字段、筛选、排序和 mask 规则；
3. `figure_mode` 对应的物理资格 gate；
4. `style_profile` 中的公共尺寸、字体、线宽、marker、输出和布局默认值；
5. case-specific 的合法布局覆盖，例如 double-column 或 external legend。

任何覆盖都必须写入 `plot_manifest.json`。绘图脚本不得把隐藏的单位转换、插值、connector 行或失败点修复放在默认分支中。

## 6. 环境与版本冻结

生成器必须在 manifest 中记录 Python、Matplotlib、平台、可执行文件和解析后的字体。当前 line-first pilot 已在 Python 3.13.2、Matplotlib 3.10.1 和 Times New Roman 解析结果下通过验证；其他机器必须重新记录实际解析结果。

样式 profile 的 APS-like 基线为：single-column `3.375 x 2.5 in`，double-column `6.75 x 4.6 in`。宽度作为版式基线，高度是项目当前选择，不宣称为所有 APS 期刊的普遍硬性值。

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
| `strict` | 正文定稿候选；只含已确认、有限、support 合格的物理 series，默认 SVG + 600 dpi PNG。 |
| `legacy` | 历史兼容；保持原图、原脚本、原单位、原 connector 和原输出语义，不自动迁移。 |

所有新图必须有 `plot_manifest.json`，至少包含：`figure_mode`、`style_profile`、输入 hash、generator/hash、Git commit、calculation/postprocess/source provenance、axes 单位和 transform、series state、support/mask 规则、interpolation/connector policy、输出 hash、DPI/vector 和 layout 记录。

## 11. Regression / Validation 验收

每个新 case 至少运行：

```powershell
python scripts/plotting/validate_plot_artifact.py <path-to-plot_manifest.json>
python -m pytest -q tests/unit/python/test_plotting_contract.py
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. scripts/dev/check_sop_governance.jl
julia --project=. scripts/dev/check_script_entrypoints.jl
```

strict 还必须检查：PNG 实际 DPI 不低于 600、SVG 为 vector、输入 hash 未变、无 forbidden state、无跨 gap 连接、尺寸与 profile 一致。

## 12. 失败点、断点续算与重跑

绘图失败时保留输入和失败原因，不修改源 CSV/JSON。若目标目录已经产生部分输出，下一次运行必须使用新的 sibling suffix，或由明确的开发者操作清理未完成的临时目录；不得覆盖已通过验证的 case。

输入 hash、脚本 hash、profile、Git commit 或 source run 不一致时，停止 strict 生成并退回 audit/作者审核。重跑绘图不等于重跑数值计算，必须在 manifest 中记录新的 generator/output hash。

## 13. Diagnostic 与 Formal Production 的边界

`docs/analysis` 是诊断证据区域，允许 C1 unresolved、bracket 和 estimated midpoint。它不能因为图形变得整洁就自动进入 `data/outputs/figures` 的 strict 正式目录。

`strict` 是图像合同通过，不是数值结论升级。只有 source result 的生产、收敛和物理 gate 也已通过，strict 图才可以申请论文定稿。历史 `legacy` 图不因新 profile 存在而失效，也不因新 SOP 自动重绘。

## 14. 后处理与作图

公共视觉规则采用 `line-first / landmark-only`：

- audit 显示 support 点；candidate/strict 默认隐藏普通 support marker；
- confirmed CEP 使用稀疏、醒目的实心圆点；
- `crossover` 使用短虚线，first-order 使用实线，spinodal 使用低 alpha 点线；
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
- strict pilot 输出为 SVG + 600 dpi PNG；
- 未修改 solver、Maxwell、C2、reference、transport、正式 CSV/JSON 或历史 PNG/PDF/SVG；
- strict 的后续密集图族仍必须执行本节 14 的 layout gate，不能仅因 profile 默认值而跳过人工审核。

## 17. 历史图资产退役与清理

历史 PNG/PDF/SVG 的清理采用 `asset inventory + dry-run + allowlist cleanup`，不是按扩展名、文件名或时间批量删除。

PR A 只运行 `scripts/plotting/inventory_figure_assets.py`，默认扫描 Git 已跟踪的 `data/outputs/figures` 资产，生成：

- `docs/analysis/figure_asset_registry_v1/asset_registry.json`；
- `docs/analysis/figure_asset_registry_v1/cleanup_candidates.csv`。

未跟踪的 C1/C2/pilot 文件默认排除且不修改。`docs/analysis` 诊断证据与正式图像根目录分开治理，不因图形格式相同而自动合并。

registry 只提出 `owner_review_only` 或 `keep_contract_case`，不包含删除、移动和覆盖操作。人工审核必须确认仓库外部引用、canonical case/variant 和历史证据保留策略；未确认项默认保留。

实际退役属于后续 PR B，必须以作者批准的 `path + sha256 + action` allowlist 执行，并再次检查引用、manifest 输出和文档链接。strict 新默认 SVG + 600 dpi PNG 不追溯改变历史 PDF/SVG/PNG 的保留资格。
