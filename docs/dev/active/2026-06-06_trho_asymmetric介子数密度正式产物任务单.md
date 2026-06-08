# trho_asymmetric 介子数密度正式产物任务单

更新时间：2026-06-08

## 背景

PR #116 接入 `FixedAsymmetricRho` 作为介子数密度后处理的 upstream equilibrium source，并新增 `run_combined_meson_density_scan.jl --path trho_asymmetric`。该路径在 PR #116 中只承诺 smoke / diagnostic 能力，不生产正式高精度产物。

本任务单用于在 PR #116 合并后，从最新 `main` 新开分支执行收敛性测试与正式产物生产。执行时必须遵守 `$formal-production-artifact` skill 的 hard gates：先完成 convergence gate，再按收敛证据生产正式数据和图像。

本机资源策略：高精度 convergence / production 默认不在本机运行，改用手动 GitHub Actions 远程生成 artifact。本机只做 smoke、任务编排、artifact 下载后的审计、正式路径整理和提交。

## 执行状态（2026-06-06）

Verdict: `blocked`

- PR #116 已合并，当前生产分支为 `codex/trho-asymmetric-meson-density-production`。
- 已触发 GitHub Actions remote run `27064914082`，完成 `convergence_low`。
- `convergence_low` 已复制到正式 convergence evidence 路径：
  - `data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low/`
  - `data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low/`
- `convergence_low` 统计显示 880 行中 `ok=245`、`failed=221`、`unsafe_bose_domain=414`；只有 13 个 `(T, rho_target)` 点四种 regime 全部 `ok`。
- blocker 不是积分节点收敛不足，而是当前 `K+ / pi+` 非对称 charged profile 下正 `mu_K` 触发 Bose-domain guard：部分 stable 点已有 `mass <= mu_K`，strict BW 与 phase-shift 在默认 `strict_normal_domain` 下大量遇到 `omega <= mu_K`。
- 因此未运行 `convergence_mid`、`convergence_high` 和 production；不得把本轮结果标成正式 production-grade。
- 已在 result-side 写入 blocked audit：
  - `README.md`
  - `PRODUCTION_AUDIT.md`
  - `manifest.json`
  - `convergence/convergence_summary.json`

后续如需继续正式生产，必须先做物理口径决策：切换到 `K- / pi-`、缩小到严格有效定义域，或显式改成 diagnostic-only density policy。

## 续跑决策（2026-06-08）

Verdict: `blocked pending upstream branch-stability convergence rerun`

- 物理口径决策：复用项目已有 `density_policy=:x_min_cut` / `:bose_x_min_cut`，但只作用于 BU/GBU phase-shift regimes。
- `x_min_cut` 定义为 `(omega - mu_M)/T >= bose_x_min`，等价于 phase-shift energy domain 下界 `omega_lower=max(omega_min, mu_M + bose_x_min*T)`。
- `stable` 与 `strict_bw_stage*` 不采用 x-min 延拓；若触发 `mass <= mu_M` 或 `omega_min <= mu_M` 等 Bose-domain guard，则输出 `unsafe_bose_domain` / `NaN` 行并继续扫描。
- 该策略是显式 policy choice，不是文献唯一处方；正式审计必须记录 `density_policy`、`bose_x_min` 和 policy 作用域，避免后续误认为 stable/BW 也被延拓。
- 由于 combined scan CLI 需要新增 `--bose-x-min` 才能由 GitHub Actions 传递该参数，续跑必须先提交并推送生产分支，再以该分支触发 workflow。
- 远端 `bose_x_min=1e-6` 的 `convergence_low/mid/high` 已在 commit `4c6c1542` 上完成：BU/GBU phase-shift rows 可运行，但 `mid -> high` 的 `n_K` / `kpi_ratio` 最大相对差约 50%，不得进入 production。
- 远端 `bose_x_min=1e-2` 诊断明显降低 cutoff 敏感性，但 `mid -> high` 的 K 相关最大相对差仍约 11%--13%，且 `high36` 与 `custom64` 在相同 `(T, rho_target)` 上出现不同 `FixedAsymmetricRho` equilibrium branch。
- 关键异常例：`T=140 MeV, rho_target=0.150`，`high36` 与 `custom64` 均满足约束残差打印为 0，但给出不同 `muq`、`m_u/m_s`、`m_pi/m_K`；这说明当前 blocker 是上游 density-constrained equilibrium branch instability，不是单纯 density quadrature 节点不足。
- 本分支后续先修 combined scan 的上游遍历策略：`trho_asymmetric` 默认改为温度分组、rho 连续扫描，并采用与 `Models.run_trho_scan` 一致的 `trho_reverse_rho=true` 习惯；修复后仍必须重新跑 convergence gate，不得直接触发 production。

## Scope Lock

- 物理口径：`FixedAsymmetricRho` density-constrained equilibrium source。
- 扫描入口：`scripts/relaxtime/run_combined_meson_density_scan.jl --path trho_asymmetric`。
- 上游分支策略：`trho_asymmetric` 默认使用 `--trho-reverse-rho=true`，每个温度内从高 `rho_target` 到低 `rho_target` 连续传递 equilibrium seed；`--no-trho-reverse-rho` 仅作为顺序敏感诊断。
- 远程入口：`.github/workflows/relaxtime-meson-density-production.yml`，只通过 `workflow_dispatch` 手动触发。
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
- [ ] 确认 `.github/workflows/relaxtime-meson-density-production.yml` 已在默认分支可手动触发。
- [ ] 明确最终网格范围：`T` 轴、`rho_target` 轴、density regimes。
- [ ] 明确是否需要把大网格按 `rho_values` 分片运行；若单次 Action 接近 6 小时上限，必须分片而不是放宽收敛 gate。

## Phase 1: Convergence Gate

优先用 GitHub Actions 远程运行，不在本机跑高精度网格。每个 convergence 档位触发一次 workflow：

```sh
gh workflow run relaxtime-meson-density-production.yml --ref codex/trho-asymmetric-meson-density-production \
  -f case_slug=trho_asymmetric_kplus_piplus_scan_v1 \
  -f run_stage=convergence_low \
  -f resolution_profile=low \
  -f phase_args="--real-axis-mode pv_b0_eta0 --phase-convention arg_inverse_propagator --phase-display fold_0_pi --density-policy x_min_cut --bose-x-min 1e-6 --noanom-policy none --eta 1e-6 --strict-bw-stage stage1"

gh workflow run relaxtime-meson-density-production.yml --ref codex/trho-asymmetric-meson-density-production \
  -f case_slug=trho_asymmetric_kplus_piplus_scan_v1 \
  -f run_stage=convergence_mid \
  -f resolution_profile=mid \
  -f phase_args="--real-axis-mode pv_b0_eta0 --phase-convention arg_inverse_propagator --phase-display fold_0_pi --density-policy x_min_cut --bose-x-min 1e-6 --noanom-policy none --eta 1e-6 --strict-bw-stage stage1"

gh workflow run relaxtime-meson-density-production.yml --ref codex/trho-asymmetric-meson-density-production \
  -f case_slug=trho_asymmetric_kplus_piplus_scan_v1 \
  -f run_stage=convergence_high \
  -f resolution_profile=high \
  -f phase_args="--real-axis-mode pv_b0_eta0 --phase-convention arg_inverse_propagator --phase-display fold_0_pi --density-policy x_min_cut --bose-x-min 1e-6 --noanom-policy none --eta 1e-6 --strict-bw-stage stage1"
```

workflow 为避免 GitHub Actions `workflow_dispatch` 输入数量上限，使用组合式参数：

- `grid_args`：扫描网格，例如 `--tmin 120 --tmax 220 --tstep 10 --rhomin 0.05 --rhomax 1.00 --rhostep 0.05`。
- `physics_args`：约束与 profile，例如 `--asym-ud-ratio-target 0.876 --asym-s-target 0.0 --flavor-profile default --meson-profile asymmetric_kplus_over_piplus_signed`。
- `resolution_profile`：`low` / `mid` / `high` / `custom`。若使用 `custom`，必须填写 `custom_resolution_args`。
- `phase_args`：相移和 BW 相关 policy，例如 `--real-axis-mode pv_b0_eta0 --phase-display fold_0_pi --density-policy x_min_cut --bose-x-min 1e-6 --noanom-policy none`。`density_policy` / `bose_x_min` 在 combined scan 中只作用于 phase-shift current/GBU，stable/BW 继续记录 strict-domain unsafe 行。
- `extra_args`：追加在命令末尾的少量补充参数；参数值不能包含空格。

远程 result-side artifact 内部路径必须保留为：

`data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/<run_stage>/`

远程 figure-side artifact 内部路径必须保留为：

`data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/<run_stage>/`

下载时先落到临时审计目录，确认内容后再复制到仓库正式路径：

```sh
gh run download <run-id> --dir data/outputs/remote_artifacts/trho_asymmetric_kplus_piplus_scan_v1
```

下载后的 artifact 若包含完整 `data/outputs/...` 相对路径，可在审计通过后复制回仓库根目录；不得把 `remote_artifacts/` 本身作为正式产物提交。

至少运行三档参数：

| tag | 用途 | 建议参数 |
| --- | --- | --- |
| `low` | 粗判形态和失败点 | `stable_q_nodes=96`, `q_nodes=12`, `omega_nodes=12` |
| `mid` | 主判定候选 | `stable_q_nodes=192`, `q_nodes=24`, `omega_nodes=24` |
| `high` | 参考档 | `stable_q_nodes=256`, `q_nodes=36` 或 `48`, `omega_nodes=36` 或 `48` |

每档远程 artifact 必须保存：

- [ ] 完整命令。
- [ ] stdout / stderr log。
- [ ] `combined_meson_density_scan.csv`。
- [ ] `README.md`。
- [ ] `remote_run_manifest.json`。
- [ ] 参数 manifest 或等价机器可读记录。

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

使用通过 convergence gate 的最高参数或有证据支持的足够收敛参数，远程触发正式 production。示例：

```sh
gh workflow run relaxtime-meson-density-production.yml --ref codex/trho-asymmetric-meson-density-production \
  -f case_slug=trho_asymmetric_kplus_piplus_scan_v1 \
  -f run_stage=production \
  -f resolution_profile=<selected-profile-or-custom> \
  -f phase_args="--real-axis-mode pv_b0_eta0 --phase-convention arg_inverse_propagator --phase-display fold_0_pi --density-policy x_min_cut --bose-x-min 1e-6 --noanom-policy none --eta 1e-6 --strict-bw-stage stage1" \
  -f custom_resolution_args="<only-when-custom>"
```

production artifact 下载后必须先审计，再整理到仓库正式路径：

- result-side：`data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`
- figure-side：`data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/`

GitHub Actions 只负责生成 artifact，不自动 commit 正式数据。正式入库由生产分支在本机下载、审计、整理后提交。

正式 result-side 至少包含：

- [ ] `combined_meson_density_scan.csv`
- [ ] `README.md`
- [ ] `PRODUCTION_AUDIT.md`
- [ ] `manifest.json`
- [ ] `remote_run_manifest.json`
- [ ] `run.stdout.log`
- [ ] `run.stderr.log`
- [ ] `convergence/` 原始证据和摘要

正式 figure-side 至少包含：

- [x] SVG 或 PNG 主图。`run_combined_meson_density_scan.jl` 的 `trho_asymmetric` 内置 SVG 已使用 `rho_target` 作为横轴；`render_combined_meson_density_fig3_like.py` 可通过 `--x-field rho_target --x-label "rho/rho0" --x-unit ""` 复用正式 CSV 渲染 `T-rho` PNG/SVG；主图已改用 `--color-scale log` 处理跨数量级 `K/pi` 比值。
- [x] `plot_manifest.json`。

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
- 若 phase-shift regime 在 `density_policy=x_min_cut` 下仍出现 unsafe/no-support/非 `ok` status，不得硬标 production-grade，必须降级为 diagnostic-only 或 blocked。
- stable/BW 的 `unsafe_bose_domain` 行是 strict-domain 对照口径的定义域诊断，不单独阻塞 BU/GBU convergence gate；但必须在 audit 中单独统计，不能作为可用 stable/BW 数值解释。

## Phase 4: 验证与提交前清理

至少运行：

- [ ] `julia --project=. -e 'ENV["INTEGRATION_FILES"]="relaxtime/test_combined_meson_density_scan_smoke.jl"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_script_entrypoints.jl`
- [ ] `julia --project=. scripts/dev/check_relaxtime_script_governance.jl`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_data_output_path_guard.jl`
- [ ] `git diff --check`

若本任务实际修改数值路径、solver、density kernel 或扫描语义，需重新评估 regression / validation；否则只提交正式数据、正式图像、审计文档和必要脚本辅助改动。
