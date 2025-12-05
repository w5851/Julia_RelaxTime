# PNJL 相图绘制（步骤 7）

本指南介绍如何基于前面步骤生成的 `TrhoScan` 数据，使用 `scripts/pnjl/plot_phase_diagram.jl` 绘制完整的 `T-μ` 相图（包含共存线、自旋odal 与 CEP 标记）。

## 1. 准备输入

1. 已完成 `run_trho_scan`，并在 `data/outputs/results/pnjl/trho_scan.csv` 中存有扫描结果；
2. `TrhoScan` 需要覆盖目标温度区间，并提供足够密的 `ρ` 采样，以便步骤 6（Maxwell）收敛；
3. （可选）提前执行 `scripts/pnjl/plot_curves.py` 以快速检查回环质量。

## 2. 运行脚本

在仓库根目录下执行：

```powershell
julia --project scripts/pnjl/plot_phase_diagram.jl \
    --trho-csv data/outputs/results/pnjl/trho_scan.csv \
    --xi 0.0 \
    --output data/outputs/figures/pnjl/pnjl_phase_diagram.png \
    --dump-boundary data/outputs/results/pnjl/phase_boundary.csv \
    --processed-dir data/processed/results/pnjl
```

关键参数说明：

| 选项 | 默认值 | 说明 |
|------|--------|------|
| `--trho-csv` | `data/outputs/results/pnjl/trho_scan.csv` | `TrhoScan` 输出路径 |
| `--xi` / `--xi-tol` | `0.0 / 1e-6` | 只保留匹配的各向异性参数 |
| `--output` | `data/outputs/figures/pnjl/pnjl_phase_diagram.png` | 图像输出位置（自动创建目录） |
| `--min-samples` | `12` | 每条曲线最少样本数，低于该值的温度被跳过 |
| `--area-tol` | `1e-4` | Maxwell 等面积残差容差 |
| `--candidate-steps` | `64` | μ 初始扫描步数，用于寻找正负残差区间 |
| `--max-iter` | `60` | Maxwell 二分最大迭代次数 |
| `--tmin/--tmax` | `none` | 可选温度窗口（单位 MeV） |
| `--dump-boundary` | `none` | 若提供路径，将同步导出 `T, μ_coex, ρ_g, ρ_l` 表 |
| `--processed-dir` | `data/processed/results/pnjl` | 保存中间结果（曲线重采样、Maxwell 调试、自旋odal、CEP） |
| `--processed-figures-dir` | `data/processed/figures/pnjl` | 输出各类调试图（曲线总览、自旋odal、Maxwell 成功/失败） |

执行成功后，终端会打印 `Saved phase diagram to ...` 并在目标路径下生成 PNG；若指定 `--dump-boundary`，还会得到用于后续分析/绘图的 CSV。同时脚本会在 `--processed-dir` 下写入：

- `curves.csv`：`build_curves` 生成的各温度 `(μ, ρ)` 轨迹；
- `maxwell_results.csv`：每个温度对应的 Maxwell 收敛状态与原因；
- `spinodals.csv`：`detect_s_shape` 得到的上下自旋odal；
- `cep.csv`：CEP 搜索结果（若未找到则只记录 `false`）。
- `curves.png / spinodals.png / maxwell.png`：位于 `--processed-figures-dir`，分别展示所有曲线、上下自旋odal 以及 Maxwell 成功/失败分布，便于快速目视检查。

## 3. 图像内容

- **红色实线**：Maxwell 等面积得到的一阶相变线 `μ_coex(T)`；
- **灰色虚线**：由 `detect_s_shape` 导出的上下自旋odal；
- **蓝色星形**：`CEPFinder.find_cep` 返回的临界点（若存在）；
- 若某些温度缺乏 S 形或样本不足，该温度将被跳过并在日志里提示。

## 4. 常见问题

- **图像为空**：多数情况是 `Trho` 数据未覆盖所选 `ξ` 或温度窗口未命中；通过 `--xi`、`--tmin/--tmax` 检查；
- **Maxwell 不收敛**：提高 `--candidate-steps` 或放宽 `--area-tol`，并确认 ρ 采样足够密；
- **CEP 未标出**：意味着 `find_cep` 没有找到从“有 S 形”到“无 S 形”的温度对，需扩展温度范围或改进扫描；
- **需要导出数据供其他工具使用**：使用 `--dump-boundary` 获取 CSV，再交由自定义可视化脚本或 `plot_curves.py` 叠加显示。

通过该脚本，可以在同一流程内串联步骤 4-7 的结果，快速验证 PNJL 模型的一阶相变结构。