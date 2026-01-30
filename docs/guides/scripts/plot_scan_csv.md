# plot_scan_csv.py 使用说明（概要）

位置：`scripts/plot_scan_csv.py`

目的：从符合 `scan_csv_v1` 规范的 CSV 文件生成期刊级别的折线图或热力图（支持多格式导出：PDF/EPS/PNG）。

快速示例：

```bash
python scripts/plot_scan_csv.py \
  --csv data/examples/results/tmp_test_scan.csv \
  --mode lines \
  --x X --ys Y1,Y2 \
  --group group \
  --out-dir data/examples/figures/plot_scan_csv \
  --formats pdf,png \
  --dpi 600 \
  --line-style - \
  --check
```

说明与约定：

- 文档位置：该文件为脚本级使用说明，放在 `docs/guides/scripts/` 下以便集中管理所有脚本的使用示例。
- 示例数据与输出：请把用于展示的 CSV 放入 `data/examples/results/`，把生成的输出放入 `data/examples/figures/`，以便提供给 CI、审阅者与新开发者。
- 图像导出要点（已在脚本中实现）：
  - 默认使用 `Times New Roman` 字体、字号 10pt、线宽 1.5pt，刻度方向向内，默认不显示网格，默认不绘制数据点（可通过 `--marker` 指定）。
  - 推荐输出格式为 `pdf`（矢量），并同时生成高分辨率 `png` 作为备份。
  - 使用 `--check` 会对 PNG 文件读取 DPI 元数据并发出警告（若 DPI 低于 `--dpi`）。对于 PDF/EPS，将提醒用户手动验证是否为矢量且字体已嵌入。
  - 导出时使用 `bbox_inches='tight'` 和 `pad_inches=0.05` 以减少多余边距。

更多参数说明与高级用法请参考脚本顶部的说明或直接运行：

```bash
python scripts/plot_scan_csv.py --help
```
