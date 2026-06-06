# trho_asymmetric 介子数密度正式产物任务单

更新时间：2026-06-06

## 背景

PR #116 接入 `FixedAsymmetricRho` 作为介子数密度后处理的 upstream equilibrium source，并新增 `run_combined_meson_density_scan.jl --path trho_asymmetric`。该路径在 PR #116 中只承诺 smoke / diagnostic 能力，不生产正式高精度产物。

本任务单用于在 PR #116 合并后，从最新 `main` 新开分支执行收敛性测试与正式产物生产。执行时必须遵守 `$formal-production-artifact` skill 的 hard gates：先完成 convergence gate，再按收敛证据生产正式数据和图像。

## Scope Lock

- 物理口径：`FixedAsymmetricRho` density-constrained equilibrium source。
- 扫描入口：`scripts/relaxtime/run_combined_meson_density_scan.jl --path trho_asymmetric`。
- charged profile：先生产 `asymmetric_kplus_over_piplus_signed`，即 `pi+` / `K+`。
- 约束目标：`rho_u/rho_d = 0.876`，`rho_s = 0 fm^-3`。
- 观测量：`n_pi`、`n_K`、`kpi_ratio`、status counts、unsafe Bose diagnostics、constraint diagnostics。
- 结果目录：`data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`。
- 图像目录：`data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`。
- 非目标：不新增求解器；不改变 `FixedAsymmetricRho` 约束语义；不更新 numerical regression baseline；不把 convergence 失败或 diagnostic-only 结果标成正式产物。

## Phase 0: 前置条件

- [ ] PR #116 已合并到 `main`。
- [ ] 从最新 `main` 新开生产分支，例如 `codex/trho-asymmetric-meson-density-production`。
- [ ] 确认 `run_combined_meson_density_scan.jl --path trho_asymmetric`、`render_combined_meson_density_fig3_like.py` 和输出路径治理在 `main` 可用。
- [ ] 明确最终网格范围：`T` 轴、`rho_target` 轴、density regimes。

## Phase 1: Convergence Gate

在 result 目录下创建：

`data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/`

至少运行三档参数：

| tag | 用途 | 建议参数 |
| --- | --- | --- |
| `low` | 粗判形态和失败点 | `stable_q_nodes=96`, `q_nodes=12`, `omega_nodes=12` |
| `mid` | 主判定候选 | `stable_q_nodes=192`, `q_nodes=24`, `omega_nodes=24` |
| `high` | 参考档 | `stable_q_nodes=256`, `q_nodes=36` 或 `48`, `omega_nodes=36` 或 `48` |

每档必须保存：

- [ ] 完整命令。
- [ ] stdout / stderr log。
- [ ] `combined_meson_density_scan.csv`。
- [ ] `README.md`。
- [ ] 参数 manifest。

比较矩阵至少覆盖：

- [ ] `n_pi`、`n_K`、`kpi_ratio` 的 pointwise absolute / relative difference。
- [ ] status counts。
- [ ] `unsafe_bose_count`、`min_E_minus_mu`、`bose_x_min`。
- [ ] `constraint_residual_norm`。
- [ ] `rho_norm`、`rho_u_over_rho_d`、`rho_s_fm3`。
- [ ] NaN / Inf / 负值计数。
- [ ] 失败点、边界点和异常点列表。

收敛判定必须输出三态 verdict：

- `production-grade`
- `diagnostic-only`
- `blocked`

只有 `production-grade` 可以进入 Phase 2。

## Phase 2: 正式生产

使用通过 convergence gate 的最高参数或有证据支持的足够收敛参数，重跑正式 production。

正式 result-side 至少包含：

- [ ] `combined_meson_density_scan.csv`
- [ ] `README.md`
- [ ] `PRODUCTION_AUDIT.md`
- [ ] `manifest.json`
- [ ] `run.stdout.log`
- [ ] `run.stderr.log`
- [ ] `convergence/` 原始证据和摘要

正式 figure-side 至少包含：

- [ ] SVG 或 PNG 主图。
- [ ] `plot_manifest.json`。

result README 必须反向链接 figure 目录和图像文件。

## Phase 3: Audit

`PRODUCTION_AUDIT.md` 必须至少记录：

- [ ] production case。
- [ ] physics scope。
- [ ] non-goals。
- [ ] command log。
- [ ] convergence matrix。
- [ ] selected production parameters。
- [ ] data outputs。
- [ ] figure outputs。
- [ ] validation commands and results。
- [ ] known limitations and residual risks。

已知风险必须显式写入：

- BW 口径仍只作为对照，不能单独作为正式结论依据。
- `trho_asymmetric` 当前是新 path strategy；正式结论依赖 convergence gate 证明。
- 若 phase-shift regime 出现 unsafe Bose domain 或非 `ok` status，不得硬标 production-grade，必须降级为 diagnostic-only 或 blocked。

## Phase 4: 验证与提交前清理

至少运行：

- [ ] `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_combined_meson_density_scan_smoke.jl"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [ ] `julia --project=. scripts/dev/check_relaxtime_script_governance.jl`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_data_output_path_guard.jl`
- [ ] `git diff --check`

若本任务实际修改数值路径、solver、density kernel 或扫描语义，需重新评估 regression / validation；否则只提交正式数据、正式图像、审计文档和必要脚本辅助改动。

