# scripts/relaxtime/

本目录包含 relaxtime 主流程脚本与扫描编排脚本。

## 需求固化与防漂移说明

为避免“需求口头更新、脚本行为漂移”，本目录下的需求-能力映射以以下规格文档为准：

- `docs/superpowers/specs/2026-03-27-relaxtime-demand-oriented-scan-config-platform-design.md`
- `docs/superpowers/specs/2026-03-28-mott-transition-temperature-xi-scan-design.md`

当需求与现有脚本能力存在偏差时，优先更新 spec，再实施脚本改造；不要只改命令示例而不更新契约。

## 当前与 Mott 扫描相关的脚本能力（现状声明）

### 1) 计算侧

- `run_gap_meson_mass_scan.jl`
  - 已支持 `(T, muB, xi)` 网格扫描。
  - 可输出 `M_pi`, `M_K`, `Gamma_pi`, `Gamma_K`, `m_u`, `m_d`, `m_s`。
  - 可用于 `muB=0` 固定下的温度-xi 扫描。

### 2) 绘图侧

- `scripts/plot_scan_csv.py`
  - 支持 `--group`：同一图中按分组列画多条曲线（适合“不同 xi 同图”）。
  - 支持 `--split`：按分组列拆分输出目录（适合“不同 xi 分图”）。
  - 支持 `--multi-y`：同一图中画多个因变量（适合“不同因变量同图”）。

## 已冻结的 Mott v1 扫描契约（摘要）

- 固定 `mu_B = 0 MeV`
- `xi_list = [-0.3, -0.15, 0.0, 0.15, 0.3]`
- 温度网格采用 `T_min_MeV`, `T_max_MeV`, `T_step_MeV`
  - 推荐默认值：`120.0`, `260.0`, `2.0`
- 目标因变量：
  - `M_pi`, `M_K`, `Gamma_pi`, `Gamma_K`, `M_u_plus_M_d`, `M_u_plus_M_s`
- 数值连续性策略：
  - 扫描按每个 `xi` 的温度升序执行。
  - 默认启用 meson 连续性初值传递（`meson_seed_state` + `mixed_seed_tracking_state`）以抑制根切换导致的伪跳变。
- 绘图组织：
  - Mode A：不同因变量同图，不同 `xi` 分图
  - Mode B：不同 `xi` 同图，不同因变量分图

完整验收标准与命名规则请以 `2026-03-28` 的 spec 为准。

## 命名约定（生产脚本）

- 本目录面向用户/流程入口的生产脚本统一使用 `run_*.jl` 前缀。
- `scan_*.jl` 保留给底层扫描实现或实验性扫描逻辑，不作为需求固化后的首选入口。
- 对于 Mott v1：
  - `run_mott_phase_scan.jl`
  - `run_mott_phase_derived_csv.jl`
  - `run_mott_phase_plot_modes.jl`

## Mott v1 快速使用

1) 扫描（固定 `mu_B=0`，默认 `xi_list` 与温度网格可通过 CLI 覆盖）

```bash
julia --project=. scripts/relaxtime/run_mott_phase_scan.jl --outdir data/outputs/results/relaxtime/mott_phase/demo --overwrite
```

2) 生成派生因变量（`M_u_plus_M_d`, `M_u_plus_M_s`）

```bash
julia --project=. scripts/relaxtime/run_mott_phase_derived_csv.jl \
  --in data/outputs/results/relaxtime/mott_phase/demo/mott_phase_scan.csv \
  --out data/outputs/results/relaxtime/mott_phase/demo/mott_phase_derived.csv
```

3) 生成 Mode A / Mode B 图

```bash
julia --project=. scripts/relaxtime/run_mott_phase_plot_modes.jl \
  --in data/outputs/results/relaxtime/mott_phase/demo/mott_phase_derived.csv \
  --out-dir data/outputs/figures/relaxtime/mott_phase/demo
```
