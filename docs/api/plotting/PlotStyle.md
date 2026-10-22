# 公共曲线样式与局部渲染上下文

入口为 `scripts/plotting/plot_series.py`、`plot_style.py` 和 `plot_quality.py`。
这些 Python helper 只处理身份、样式和导出；字段单位、数据选择和物理相态由图族定义。
尺寸与交付要求沿用[绘图 SOP](../../guides/sop/workflows/figure_production.md)。

## 曲线身份

先注册完整 `SeriesStyleRegistry(registry_id, styles)`，再用 `select(series_ids)`
抽取或重排。每个 `SeriesStyle(series_id, color, linestyle="-", marker="")` 固定到
物理条件；`plot_kwargs()` 返回颜色、线型、标记及用于实测核对的 Matplotlib `gid`。
注册表不按当前子集重新分配槽位，也不循环颜色。

`SeriesStyle.from_profile(profile, series_id, style_index, linestyle=None, marker="")`
采用从零开始的显式颜色槽位；默认线型采用相同槽位。颜色槽位越界、缺少对应线型、
空白／重复／未知身份及完全相同的编码均抛出 `ValueError`。超出默认样式集合时，
显式声明新的 `SeriesStyle` 或线型；灰度可辨性仍由可访问性合同独立检查。
相态已经使用线型的图族必须保留相态语义，通过 marker／直接标签区分参数。

## 创建、延迟保存与恢复

```python
import matplotlib.pyplot as plt
from pathlib import Path
from scripts.plotting.plot_style import load_profile, style_context, configure_axis_ticks
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
from scripts.plotting.plot_quality import export_figure

profile = load_profile("candidate_aps_v2")
registry = SeriesStyleRegistry("example_conditions_v1", (
    SeriesStyle("xi=-0.2", "#0072B2", "-"),
    SeriesStyle("xi=0.2", "#D55E00", "--"),
))
with style_context(profile) as context:
    figure, ax = plt.subplots(figsize=(6.75, 4.6))
    context.bind(figure)
    for style in registry.select(["xi=0.2", "xi=-0.2"]):
        ax.plot([100, 150, 200], [1, 2, 1], **style.plot_kwargs())
    ax.set_xlabel(r"$T\;(\mathrm{MeV})$")
    ax.set_ylabel("Ratio")
    configure_axis_ticks(ax, profile)

try:
    outputs, quality = export_figure(figure, Path("NEW_CASE/example"), profile, formats=["png"])
    series_styles = registry.record(["xi=0.2", "xi=-0.2"])
finally:
    plt.close(figure)
```

`style_context` 退出及异常退出均恢复调用者的相关 `rcParams`，支持嵌套顺序使用。
`context.bind(figure)` 绑定创建时已解析的字体与设置；`measure_figure` 和
`export_figure` 在延迟绘制／保存时重用它，避免其它主题改变数学字体、PDF 字体类型
或背景。向绑定图传入不同 profile 会报错。图族负责关闭 figure 和固定自身覆盖；
上下文不会重新修改已有 artist。Matplotlib 的全局配置不支持线程并发隔离，
同进程主题按顺序使用。直接调用 `figure.savefig` 时须自行进入
`context.activate()`；公共导出入口自动处理。

`configure_matplotlib(profile)` 保留为独立 CLI 的全局兼容入口，成功后修改当前
进程默认值。字体解析失败时不会先修改配置。v1 保留原全局样式及历史字体解析行为。

## 字体与 manifest

| `font_policy` | 缺失字体行为 |
| --- | --- |
| `error` | 只使用 `font_candidates[0]`；缺失抛出 `ValueError`。当前 strict v2 采用此策略。 |
| `candidate_fallback` | 按声明顺序尝试，替代时发出 `UserWarning`；全部缺失则报错。当前 candidate v2 依次尝试 Times New Roman、Times、STIXGeneral。 |
| `legacy_default` | 仅用于弃用的 v1；保留历史隐式替代行为，不能用于新 v2。 |

`resolve_font`／`context.font` 返回 `plot_font_resolution_v1`：请求、候选、策略、
选择、实际家族、绝对文件路径、bytes、SHA-256、是否替代、缺失候选及数学字体集。
`quality.actual_fonts` 逐份保留实际普通／数学字形的源文件、hash 与角色；
字形高度继续从字体轮廓测量，字体解析成功不能代替 2 mm 检查或 PDF 嵌入检查。

新图在 `rendering.series_styles` 写入 `registry.record(selected_ids)`，在
`rendering.quality` 写入公共测量结果。`build_manifest` 自动将绑定字体放入
`rendering.font`，拒绝与实测字体矛盾的声明；validator 核对实际字体文件及曲线
artist 编码。新字体／身份字段在 PNG 与后续 PDF 阶段必须一致。
未声明这些字段的冻结历史 manifest 保留历史核验方式，不追溯要求迁移。
缺少 `font_policy` 的历史 v2 配置副本可以读取其尺寸等合同字段；若要渲染这类旧图，
使用对应冻结代码。现行 `style_context` 仍拒绝无显式策略的新渲染。
