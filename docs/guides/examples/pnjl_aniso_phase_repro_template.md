# 各向异性 PNJL 相图可复现实验模板

本模板用于一键复现“各向异性 PNJL 扫描 + 相结构产出 + 相图绘制”的最小闭环。

## 1. 入口脚本

- 脚本：`scripts/pnjl/run_aniso_phase_template.jl`
- 能力：
  - 对多个 `xi` 依次执行 `run_tmu_scan.jl`
  - 对多个 `xi` 依次执行 `calculate_phase_structure.jl`
  - 可选调用 `plot_phase_diagram.py` 绘图
  - 产出 `manifest.txt` 固化本次实验参数与路径

## 2. 最小可复现命令

在仓库根目录执行：

```powershell
julia --project=. scripts/pnjl/run_aniso_phase_template.jl --profile=smoke --xi-values=0.0,0.2
```

可选：

- 跳过绘图（仅产出数值数据）

```powershell
julia --project=. scripts/pnjl/run_aniso_phase_template.jl --profile=smoke --xi-values=0.0,0.2 --skip-plot
```

- 全量网格（计算耗时明显更长）

```powershell
julia --project=. scripts/pnjl/run_aniso_phase_template.jl --profile=full --xi-values=0.0,0.2,0.4
```

## 3. 输出目录约定

模板输出位于：

- `data/outputs/results/pnjl/experiment_templates/aniso_phase/{tag}/scan/`
  - `tmu_scan_xi{...}.csv`
- `data/outputs/results/pnjl/experiment_templates/aniso_phase/{tag}/phase/`
  - `trho_scan_xi{...}.csv`
  - `boundary.csv`
  - `cep.csv`
  - `spinodals.csv`
  - `crossover.csv`
- `data/outputs/figures/pnjl/experiment_templates/aniso_phase/{tag}/`
  - 相图图片（若未 `--skip-plot`）
- `data/outputs/results/pnjl/experiment_templates/aniso_phase/{tag}/manifest.txt`

其中 `{tag}` 默认为时间戳，可通过 `--tag=...` 固定。

## 4. 验收标准（DoD）

一次模板运行视为通过，需要满足：

1. `manifest.txt` 存在，且记录 `profile/xi_values`。
2. 每个 `xi` 至少有一份 `scan/tmu_scan_xi*.csv`。
3. `phase/` 下存在 `boundary.csv` 与 `cep.csv`（可无 CEP 行，但文件必须存在）。
4. 若未 `--skip-plot`，图目录下存在至少一张相图文件。

## 5. 说明

- `smoke` 用于流程验证与 CI 级快速复现；`full` 用于研究级更密集采样。
- 如本机 Python 命令不是 `python`，可传 `--python=py` 或其它解释器命令。
- 模板不改变既有扫描/相结构实现，仅做稳定编排，便于重复实验与交付追踪。
