# PNJL 相图查看（步骤 7）

当前绘图入口是 `scripts/pnjl/plot_phase_diagram.py`。它读取已经生成的相结构表，绘制 T–μ 和 T–ρ 图；原来的 Julia 绘图入口与 `plot_curves.py` 已退役。

## 1. 查看 accepted 参考层

在仓库根目录运行以下 PowerShell 命令：

```powershell
python scripts/pnjl/plot_phase_diagram.py --phase-reference-layer accepted --xi 0.0 --output-dir data/outputs/figures/pnjl/phase_preview --format png --no-show
```

默认参考根目录为 `data/reference/pnjl/issue130_phase_reference_v2`。脚本通过 phase-reference adapter 读取所选层；输出包括一阶相变线、可用的 CEP、spinodal 和 crossover。缺失的数据不会补成已验证结果。

`accepted` 是参考层名称，不表示其中每一行均已获得数值认证。查看终端输出的 `uncertified_rows` 诊断；此绘图调用保持 `runtime_consumption=false`。

| 参数 | 作用 |
| --- | --- |
| `--phase-reference-root PATH` | 指定已生成的 reference 根目录 |
| `--phase-reference-layer accepted` | 选择 accepted 参考层；其他层的适用性见脚本 `--help` 与 reference adapter 契约 |
| `--xi VALUE` | 选择 ξ，可重复指定；省略时读取全部可用 ξ |
| `--output-dir PATH` | 图像输出目录 |
| `--format png` | 输出格式，可选 png/pdf/svg |
| `--no-show` | 仅保存文件，不打开交互窗口 |
| `--no-spinodal` / `--no-crossover` | 隐藏对应曲线 |

横轴化学势沿用脚本输入的夸克化学势口径；不能将 μq 数字直接标成 μB。

## 2. 从新扫描生成相结构

原始 TrhoScan CSV 需要先经过 Maxwell、CEP 等相结构分析。使用 [Phase Pipeline](../../api/models/workflows/PhasePipeline.md) 中的 `Models.run_phase_pipeline` 或 `scripts/pnjl/calculate_phase_structure.jl`，核对计算状态与产物路径后再选择对应的绘图输入。

`plot_phase_diagram.py` 不直接接受 `--trho-csv`、`--dump-boundary` 或 `--processed-dir`。显式 CSV 输入使用 `--boundary`、`--cep` 等选项，其列定义必须符合该脚本的输入 schema，不能把任意新扫描表直接替换进去。

## 3. 交付边界

上述命令用于诊断查看，不代表通过正式图像验收。论文图应按[图像生产流程](../sop/workflows/figure_production.md)记录来源、单位、质量检查和审阅状态；保存 PDF/SVG 本身不改变交付状态。
