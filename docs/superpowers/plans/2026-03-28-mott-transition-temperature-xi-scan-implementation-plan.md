# Mott Transition Temperature-xi Scan Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 固化 `mu_B=0` 的 Mott 温度-xi 扫描需求为可复现脚本流程，并产出两种固定绘图组织（Mode A/Mode B）。

**Architecture:** 采用“扫描脚本 + 派生字段脚本 + 绘图包装脚本 + 配置/元数据契约 + smoke 验证”的最小闭环。计算复用 `run_gap_meson_mass_scan.jl` 的现有能力，新增薄封装以冻结默认参数（`xi_list` 与 `T_min/T_max/T_step`）并统一输出结构。绘图继续复用 `scripts/plot_scan_csv.py`，由专用入口固定参数语义，避免命令漂移。

**Tech Stack:** Julia 1.10+（include-driven scripts）、Python matplotlib 脚本、现有 relaxtime workflow、`tests/unit` + `tests/integration` smoke。

---

## File Structure (Planned)

- Create: `scripts/relaxtime/run_mott_phase_scan.jl`
- Create: `scripts/relaxtime/run_mott_phase_derived_csv.jl`
- Create: `scripts/relaxtime/run_mott_phase_plot_modes.jl`
- Create: `config/workflows/relaxtime/profiles/mott_phase_muB0_xi_scan.toml`
- Create: `tests/unit/relaxtime/test_mott_phase_derived_csv.jl`
- Create: `tests/integration/relaxtime/test_mott_phase_scan_smoke.jl`
- Create: `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`
- Modify: `scripts/relaxtime/README.md`
- Modify: `docs/superpowers/specs/2026-03-28-mott-transition-temperature-xi-scan-design.md` (仅在实现偏差时同步更新)

## Chunk 1: 扫描与配置冻结

### Task 1: 先写 smoke 测试锁定需求契约

**Files:**
- Test: `tests/integration/relaxtime/test_mott_phase_scan_smoke.jl`

- [x] **Step 1: 写失败测试（契约先行）**
  - 断言默认配置下：`muB_MeV == 0.0`、`xi` 为 `[-0.3,-0.15,0.0,0.15,0.3]`、温度网格由 `T_min/T_max/T_step` 生成。
  - 断言输出包含 `M_pi,M_K,Gamma_pi,Gamma_K,m_u,m_d,m_s`。

- [x] **Step 2: 运行单测确认失败**
  - Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_mott_phase_scan_smoke.jl")'`
  - Expected: FAIL（入口脚本尚未实现或输出列不齐）。

- [x] **Step 3: 提交测试基线**
  - `git add tests/integration/relaxtime/test_mott_phase_scan_smoke.jl`
  - `git commit -m "test: add smoke contract for mott phase muB0 xi scan"`

### Task 2: 实现扫描入口与默认配置冻结

**Files:**
- Create: `scripts/relaxtime/run_mott_phase_scan.jl`
- Create: `config/workflows/relaxtime/profiles/mott_phase_muB0_xi_scan.toml`
- Modify: `scripts/relaxtime/README.md`
- Test: `tests/integration/relaxtime/test_mott_phase_scan_smoke.jl`

- [x] **Step 1: 实现最小扫描入口**
  - `run_mott_phase_scan.jl` 解析 `--config --outdir --resume --overwrite`。
  - 默认配置冻结：
    - `muB_MeV=0.0`
    - `xi_list=[-0.3,-0.15,0.0,0.15,0.3]`
    - `T_min_MeV=120.0, T_max_MeV=260.0, T_step_MeV=2.0`
  - 调用现有 meson/gap workflow，限制 meson 为 `pi,K`。

- [x] **Step 2: 输出元数据与列契约**
  - 写 `run_manifest.json`、`effective_config.json`。
  - CSV 至少包含：`T_MeV,muB_MeV,xi,M_pi,M_K,Gamma_pi,Gamma_K,m_u,m_d,m_s,status,error_code,error_message,timestamp_utc`。

- [x] **Step 3: 运行 smoke 测试验证通过**
  - Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_mott_phase_scan_smoke.jl")'`
  - Expected: PASS。

- [x] **Step 4: 提交实现**
  - `git add scripts/relaxtime/run_mott_phase_scan.jl config/workflows/relaxtime/profiles/mott_phase_muB0_xi_scan.toml scripts/relaxtime/README.md tests/integration/relaxtime/test_mott_phase_scan_smoke.jl`
  - `git commit -m "feat: add frozen muB0 mott phase scan entrypoint"`

## Chunk 2: 派生字段与绘图双模式

### Task 3: 先写派生字段单测

**Files:**
- Test: `tests/unit/relaxtime/test_mott_phase_derived_csv.jl`

- [x] **Step 1: 写失败测试**
  - 给定最小输入 CSV，断言新增 `M_u_plus_M_d` 与 `M_u_plus_M_s` 且数值正确。
  - 断言非法数值行（NaN/缺列）保留 `status` 并写错误字段，不中断整批。

- [x] **Step 2: 运行失败测试**
  - Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_mott_phase_derived_csv.jl")'`
  - Expected: FAIL（脚本尚不存在）。

- [x] **Step 3: 提交测试**
  - `git add tests/unit/relaxtime/test_mott_phase_derived_csv.jl`
  - `git commit -m "test: add derived observable csv unit contracts"`

### Task 4: 实现派生字段构建脚本

**Files:**
- Create: `scripts/relaxtime/run_mott_phase_derived_csv.jl`
- Test: `tests/unit/relaxtime/test_mott_phase_derived_csv.jl`

- [x] **Step 1: 实现最小派生逻辑**
  - 读取扫描 CSV。
  - 新增字段：
    - `M_u_plus_M_d = m_u + m_d`
    - `M_u_plus_M_s = m_u + m_s`
  - 保持输入顺序、保留 metadata、补充单位 metadata（`unit.<field>`）。

- [x] **Step 2: 运行单测验证通过**
  - Run: `julia --project=. -e 'include("tests/unit/relaxtime/test_mott_phase_derived_csv.jl")'`
  - Expected: PASS。

- [x] **Step 3: 提交实现**
  - `git add scripts/relaxtime/run_mott_phase_derived_csv.jl tests/unit/relaxtime/test_mott_phase_derived_csv.jl`
  - `git commit -m "feat: add mott derived observables csv builder"`

### Task 5: 先写绘图模式 smoke 测试

**Files:**
- Test: `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`

- [x] **Step 1: 写失败测试（输出布局契约）**
  - Mode A 断言：每个 `xi` 一张图，且每图包含 6 个变量曲线。
  - Mode B 断言：每个变量一张图，且每图包含全部 `xi` 曲线。
  - 断言目录结构与文件命名：
    - `mode_a/mott_mode_a__xi{xi_tag}.png`
    - `mode_b/mott_mode_b__{observable}.png`

- [x] **Step 2: 运行 smoke 预期失败**
  - Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl")'`
  - Expected: FAIL（绘图封装脚本尚未实现）。

- [x] **Step 3: 提交测试**
  - `git add tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`
  - `git commit -m "test: add mott plotting mode contract smoke tests"`

### Task 6: 实现双模式绘图封装

**Files:**
- Create: `scripts/relaxtime/run_mott_phase_plot_modes.jl`
- Test: `tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`

- [x] **Step 1: 封装 Mode A/Mode B 命令**
  - 通过 Julia 调用 `python scripts/plot_scan_csv.py`。
  - Mode A 固定：`--split xi --multi-y`。
  - Mode B 固定：`--group xi`。

- [x] **Step 2: 执行后重命名为契约文件名**
  - 统一输出到 `mode_a/`、`mode_b/`。
  - 规范 `xi_tag`（小数点转 `p`）。

- [x] **Step 3: 运行 smoke 测试通过**
  - Run: `julia --project=. -e 'include("tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl")'`
  - Expected: PASS。

- [x] **Step 4: 提交实现**
  - `git add scripts/relaxtime/run_mott_phase_plot_modes.jl tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl`
  - `git commit -m "feat: add fixed mode A/B plotting wrapper for mott scan"`

## Chunk 3: 文档对齐与验收收口

### Task 7: 同步文档并执行最小验收命令集

**Files:**
- Modify: `scripts/relaxtime/README.md`
- Modify: `docs/superpowers/specs/2026-03-28-mott-transition-temperature-xi-scan-design.md` (若实现细节有偏差)

- [x] **Step 1: 更新 README 的命令示例与输出布局**
  - 补充三段命令：scan / derived / plot。
  - 明确“需求与能力映射”以及权威 spec 路径。

- [x] **Step 2: 运行最小验收命令集**
  - `julia --project=. -e 'include("tests/unit/relaxtime/test_mott_phase_derived_csv.jl")'`
  - `julia --project=. -e 'include("tests/integration/relaxtime/test_mott_phase_scan_smoke.jl")'`
  - `julia --project=. -e 'include("tests/integration/relaxtime/test_mott_phase_plot_modes_smoke.jl")'`
  - `julia --project=. scripts/dev/check_docs_consistency.jl`

- [x] **Step 3: 验证产物完整性**
  - 检查 run 目录中是否存在：`effective_config.json`、`run_manifest.json`、扫描 CSV、派生 CSV、Mode A/B 图目录。

- [x] **Step 4: 提交收口**
  - `git add scripts/relaxtime/README.md docs/superpowers/specs/2026-03-28-mott-transition-temperature-xi-scan-design.md`
  - `git commit -m "docs: align mott scan contracts and usage with implemented workflow"`

## Done Definition (DoD)

- 单命令可完成 `muB=0` + 冻结 `xi_list` + 冻结温度网格的 Mott 扫描。
- 派生字段 `M_u_plus_M_d`、`M_u_plus_M_s` 可重复生成且通过单测。
- Mode A / Mode B 绘图组织满足命名与数量契约并通过 smoke。
- 元数据追溯文件齐全（`effective_config.json`, `run_manifest.json`）。
- `scripts/relaxtime/README.md` 与 `docs/superpowers/specs/2026-03-28-mott-transition-temperature-xi-scan-design.md` 一致。

---

## Execution Notes (2026-03-28)

- 主线实现已完成并提交：`3170002`、`a7cc658`。
- 在扫描入口中新增连续性求解传种（`meson_seed_state` + `mixed_seed_tracking_state`）与诊断列（`residual_*`、`root_quality_*`、`selected_method_*`），用于抑制根切换伪跳变并增强可追溯性。
- 新增回归：`tests/regression/relaxtime/test_mott_phase_continuity_regression.jl`。
- 产物已按用户要求清理为仅保留一份最新结果：
  - `data/outputs/results/relaxtime/mott_phase/reference_100_300_fine`
  - `data/outputs/figures/relaxtime/mott_phase/reference_100_300_fine`
- 当前计划状态：已完成（DoD 满足）。后续若继续推进，建议进入“文献对照报告/阈值对齐”新计划而非继续扩展本计划。
