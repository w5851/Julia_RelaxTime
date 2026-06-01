# scripts/relaxtime/

本目录包含 relaxtime 主流程脚本与扫描编排脚本。

## 需求固化与防漂移说明

为避免“需求口头更新、脚本行为漂移”，本目录下的需求-能力映射以以下规格文档为准：

- `docs/superpowers/specs/2026-03-27-relaxtime-demand-oriented-scan-config-platform-design.md`
- `docs/superpowers/specs/2026-03-28-mott-transition-temperature-xi-scan-design.md`

当需求与现有脚本能力存在偏差时，优先更新 spec，再实施脚本改造；不要只改命令示例而不更新契约。

## 当前与 Mott 扫描相关的脚本能力（现状声明）

## 当前与 transport 相图邻域扫描相关的脚本能力（2026-05-07 之后）

### 1) 基础规则网格

- `run_gap_transport_scan.jl`
  - 保留为 `(T, mu_B, xi)` 规则网格扫描入口。
  - 继续承担逐点求平衡态 + transport 后处理 + CSV/provenance 落盘。

### 2) phase-guided 上层计划

- `run_phase_guided_transport_scan.jl`
  - 新增 phase-guided transport 上层入口。
  - mode a：固定 `mu_B`，沿 `T/T_phase` 倍率带连续扫描 `xi`。
  - mode b：固定 `T`、离散 `mu_B`、连续扫描 `xi`。
  - 默认 canonical 输出根目录：
    - `data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/`
    - `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/`
  - 当前支持 `--dry-run`，先固化 `sampling_plan.csv`、`README.md`、`effective_config.json`、`run_manifest.json`。

- `run_phase_guided_transport_plots.jl`
  - 面向 canonical case 的 post-processing wrapper。
  - 复用 `scripts/plot_scan_csv.py`，输出目录固定到 `data/outputs/figures/relaxtime/transport/phase_guided/<mode>/<case_name>/`。
  - `fixed-muB-phase-scaled`：按固定 `mu_B` 分图、同图按 `alpha_T` 多线。
  - `fixed-T-sparse-muB`：按固定 `T` 分图、同图按 `mu_B` 多线。
  - 会同步写出 `plot_manifest.json`，并把图层清单追加回 case `README.md`。

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
  - 默认平衡态使用 `stable` 分支策略，即 FixedMu 多初值压力选优；需要亚稳分支追踪时必须显式传 `--equilibrium-branch-mode continuation`。
  - 默认仍启用 meson 连续性初值传递（`meson_seed_state` + `mixed_seed_tracking_state`）以抑制介子根切换导致的伪跳变。
  - 介子 pole 输出采用 `M >= 0` 符号约定；若求解落到负号等价根，CSV 中 `root_sign_flipped_*` 会标记为 `true`。
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
- 对于 phase-guided transport：
  - `run_phase_guided_transport_scan.jl`

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

## Paper P1 图件资产后处理

`build_paper_p1_figure_assets.py` 面向论文 P1 图件，只消费已有主计算产物，不运行求解器：

- `run_gap_meson_mass_scan.jl` 的 `(T, muB, xi)` Mott gap 网格；
- `calculate_phase_structure.jl` 的 phase artifacts；
- `run_isentropic_meson_mass_scan.jl` 的 fixed-`s/nB` 路径 CSV。

输出固定为：

- `mott_lines.csv`
- `isentropic_trajectories.csv`
- `isentropic_mott_crossings.csv`
- `phase_overlay.csv`
- `figure_manifest.json`
- 可选 `figures/p1_mott_phase_diagram.*` 与 `figures/p1_isentropic_mott_paths.*`
- 可选 per-`xi` 子图 `figures/p1_mott_phase_diagram_xi_*.*`

`mott_lines.csv` 会保留 bracket 诊断字段；其中 `bracket_kind=branch_jump` 表示 gap 过零 bracket 同时跨过明显质量分支跳变，不应直接当作普通连续 Mott root 解释。P1 生产入口默认把 Mott 平衡态设为 `stable` 模式，即每个点用 FixedMu 多初值压力选优选择稳定分支；如需复现旧的温度连续亚稳分支扫描，可显式传 `--mott-equilibrium-branch-mode continuation`。

示例：

```bash
python scripts/relaxtime/build_paper_p1_figure_assets.py \
  --mott-grid-csv data/outputs/results/relaxtime/paper_p1_mott_phase_isentropic/mott_grid.csv \
  --isentropic-csv data/outputs/results/relaxtime/paper_p1_mott_phase_isentropic/isentropic_sigma30.csv \
  --phase-dir data/outputs/results/relaxtime/paper_p1_mott_phase_isentropic/phase_xi0 \
  --out-dir data/outputs/figures/relaxtime/paper_p1_mott_phase_isentropic
```

如果 phase 已汇总为正式 reference 产物，可直接消费 `boundary_<tag>.csv`、`spinodals_<tag>.csv`、`crossover_<tag>.csv`、`cep_<tag>.csv`：

```bash
python scripts/relaxtime/build_paper_p1_figure_assets.py \
  --mott-grid-csv data/outputs/results/relaxtime/paper_p1_mott_phase_isentropic/mott_grid.csv \
  --phase-reference-root data/reference/pnjl/paper_p1_mott_phase_isentropic_xi3 \
  --phase-reference-tag paper_p1_mott_phase_isentropic_xi3 \
  --out-dir data/outputs/figures/relaxtime/paper_p1_mott_phase_isentropic
```

若需要把本轮 P1 生产链路统一编排，可使用 `run_paper_p1_pipeline.jl`。各 stage 相互独立，可只跑数据生产或只跑后处理：

```bash
julia --project=. scripts/relaxtime/run_paper_p1_pipeline.jl \
  --stage=assets \
  --mott-grid-csv data/outputs/results/relaxtime/paper_p1_mott_phase_isentropic/mott_grid.csv \
  --phase-reference-root data/reference/pnjl/paper_p1_mott_phase_isentropic_xi3 \
  --phase-tag paper_p1_mott_phase_isentropic_xi3 \
  --isentropic-csv data/outputs/results/relaxtime/paper_p1_mott_phase_isentropic/isentropic_sigma30.csv \
  --figure-dir data/outputs/figures/relaxtime/paper_p1_mott_phase_isentropic
```

当 Mott 主网格在不同 `muB` 区间需要不同温度覆盖时，使用 slice plan CSV，而不是手工拆多批命令。CSV 必须包含 `muB_MeV,T_min_MeV,T_max_MeV,T_step_MeV`：

```csv
muB_MeV,T_min_MeV,T_max_MeV,T_step_MeV
0,100,300,5
200,100,300,5
800,100,300,5
900,30,300,5
1100,30,300,5
```

```bash
julia --project=. scripts/relaxtime/run_paper_p1_pipeline.jl \
  --stage=mott \
  --mott-slice-plan data/outputs/results/relaxtime/paper_p1_mott_phase_isentropic/mott_slice_plan.csv \
  --result-dir data/outputs/results/relaxtime/paper_p1_mott_phase_isentropic/production
```

`stage=mott` 会在合并网格旁写出 `mott_grid_combined_manifest.json`，显式记录 `equilibrium_branch_mode`、`equilibrium_selector_policy`、`equilibrium_selector_tiebreak` 和每个 `muB` slice 的温度范围。`stage=phase` 默认从 `T=30 MeV` 开始，以便低温高密 Mott bracket 能有 phase reference 覆盖；可用 `--phase-tmin` 覆盖。phase 生产默认使用 `--phase-p-num 24 --phase-t-num 8`，与既有正式 PNJL 相图 reference 的 `T-rho` 曲线口径对齐；降低节点数可能在低温 quark-side spinodal 附近产生伪回折。

对应任务单：`docs/dev/active/2026-05-30_P1论文图后端能力任务单.md`。
