# 介子数密度计算 SOP

## 1. 目的与适用范围

本 SOP 用于通过 `Models` 统一 workflow 计算介子数密度与 `K/pi` 比值，并通过组合入口执行 `scan path × density regime` 的可审计扫描。

覆盖：

- stable-particle density；
- strict BW Stage1 reduced 与 Stage2 q-pole；
- phase-shift current 与 GBU reference；
- `tmu` 固定 `mu_q` temperature paths；
- `trho_asymmetric` 的 `FixedAsymmetricRho` density-constrained paths；
- 本地 diagnostic、远程 convergence 和 formal production artifact。

## 2. 非适用范围

- 不把任一 density/phase/no-anomalous policy 表述为文献唯一处方。
- phase-shift density 当前不支持 `xi != 0` 的正式使用。
- 不把 stable/BW 的 unsafe Bose-domain 行替换为有限值。
- 不把 targeted branch repair 等价为新 policy 下的全网格重跑。
- 不负责介子 pressure/EOS；该链由介子热力学 SOP 管理。
- `scripts/analysis/relaxtime/scan_full_kmt_bu_freezeout.jl` 仅是完整 charged
  KMT 耦合在 quark-only BQS 背景上的 A/B 诊断脚本，不属于稳定 workflow，
  不得用来替代下方的正式 production 入口。

## 3. 权威入口

统一 API：

- `Models.solve_meson_density_from_meson_point`
- `Models.solve_gap_and_meson_density_point`
- `Models.solve_strict_bw_meson_density_from_meson_point`
- `Models.solve_gap_and_strict_bw_meson_density_point`
- `Models.solve_phase_shift_meson_density_from_meson_point`
- `Models.solve_gap_and_phase_shift_meson_density_point`

专题组合 CLI：

- `scripts/relaxtime/run_combined_meson_density_scan.jl`
- 轻量参数/输出合同：`scripts/relaxtime/combined_meson_density_scan_contract.jl`

远程长任务入口：

- `.github/workflows/relaxtime-meson-density-production.yml`

组合 CLI 当前是 `domain-candidate`，不是 Phase A 核心稳定白名单。正式长任务优先经 GitHub workflow 生成分层 convergence/production artifacts。

## 4. 物理口径、单位与参数约束

- `T`、`mu_q`、`mu_B`、介子质量和宽度的外部 CSV 字段：MeV。
- `mu_B = 3 mu_q` 仅在对应 equal-flavor 路径口径成立时使用。
- `rho_target/rho_norm`：`rho/rho0`。
- `rho_u/rho_d/rho_s`：`fm^-3`。
- `qmax/omega_min/omega_max`：`fm^-1`。
- `n_pi/n_K`：`fm^-3`。
- `kpi_ratio`、`xi`：无量纲。

Density regimes：

- `stable`；
- `strict_bw_stage1`；
- `strict_bw_stage2`；
- `phase_shift_current`；
- `phase_shift_gbu_reference`。

完整 KMT A/B 诊断：

- `MesonDensity.PhaseShiftInteractionSpec` 保持可审计的耦合与分母因子；
- 传入 `FullKMTInteraction` 时，`pi^±` 使用 `K12`，`K^±` 使用 `K45`；
- 默认仍采用 `2K/(1-4KPi)` 的旧 scalar BU 分母，只用于隔离耦合变化；
- 该路径不把介子密度写回 `Omega`，也不修改 BQS 平衡解。

关键治理边界：

- `current` 是兼容/对照 branch，GBU 是最终比较默认；二者当前都沿用
  `domega/(2pi)` 文献 ratio adapter，尚未完成单电荷 `domega/pi` 迁移，因此均未
  获得 strict 绝对密度 production 授权；
- `finite_eta` 与 `pv_b0_eta0` 必须在 metadata 中可区分；
- `strict_normal_domain` 遇到 Bose unsafe domain 返回 `NaN/status`，不静默 clamp；
- `x_min_cut` 与 `excitation_only_E_gt_mu` 是显式项目延拓策略；
- `low_energy_branch_subtraction` 是 reconstructed diagnostic policy；
- `phase_display=fold_0_pi` 是显示/权重口径选择，不是相移方程本身的唯一事实。

`trho_asymmetric` 默认 `asym_ud_ratio_target=0.876` 是项目历史近似；严格代数参考为 `0.875`。改变 target 需要新 case 和新 convergence gate。

## 5. 输入配置及优先级

组合 CLI 的有效输入来自：

1. `combined_meson_density_scan_contract.jl` 默认值；
2. `config/physics/flavor_chemical/<profile>.toml`；
3. `config/physics/meson_chemical/<profile>.toml`；
4. 显式 CLI 参数覆盖。

远程 workflow 进一步按顺序组合：

1. `grid_args`；
2. `physics_args`；
3. `resolution_profile` 或 custom resolution；
4. `phase_args`；
5. `extra_args` 最后覆盖。

正式 case 必须保存最终 command 和 remote manifest，不能只记录 GitHub 表单默认值。

## 6. 环境与版本冻结

本地预检记录：

```powershell
julia --version
git rev-parse HEAD
git status --short
```

长任务要求：

- CI/production workflow 使用 Julia `1.12.5` 与根项目环境；
- artifact 保持仓库相对目录结构；
- remote manifest 记录 GitHub SHA/ref/run ID、输入参数和退出码；
- 下载后先进入独立 review workspace，完成 convergence/production audit 后再决定是否入库。

## 7. Smoke 预检

轻量 contract smoke：

```powershell
julia --project=. tests/integration/relaxtime/test_combined_meson_density_scan_smoke.jl
```

默认只运行 CLI contract、路径治理、profile 和 branch-policy 静态门禁；慢数值 CLI 由环境变量显式开启：

```powershell
$env:RUN_COMBINED_MESON_DENSITY_CLI_SMOKE = "1"
julia --project=. tests/integration/relaxtime/test_combined_meson_density_scan_smoke.jl
```

独立 diagnostic 命令：

```powershell
julia --project=. scripts/relaxtime/run_combined_meson_density_scan.jl --path tmu --output-dir data/outputs/results/sop_smoke/meson_density --figure-dir data/outputs/figures/sop_smoke/meson_density --overwrite --tmin 150 --tmax 150 --tstep 10 --muq 0 --p-num 4 --t-num 2 --max-iter 6 --stable-q-nodes 6 --qmax 2 --q-nodes 3 --omega-max 2 --omega-nodes 3 --phase-display fold_0_pi
```

低节点 smoke 只验证多 regime workflow、CSV/README、SVG 和 plot manifest，不用于物理结论。

## 8. 收敛性验证

正式 case 至少比较：

- equilibrium `p_num/t_num/max_iter`；
- stable `stable_q_nodes`；
- BW/phase-shift `qmax/q_nodes`；
- `omega_min/omega_max/omega_nodes`；
- real-axis、phase convention/display 和 density policy；
- path T/`mu_q`/rho 步长；
- `FixedAsymmetricRho` branch selection、pressure 和 constraint residual；
- unsafe/failed/status counts。

核心比较量：`n_pi`、`n_K`、`kpi_ratio`、介子质量/宽度、chemical potentials、constraint diagnostics 和 branch stability。

要求使用相邻加密配置，不只比较 coarse 与最高精度。现有 production-grade `trho_asymmetric` case 的 `128 -> 192` 检查在可比较 ok rows 上最大相对差约 1%，但该阈值只属于对应 grid/profile/policy，不是全仓库统一容差。

现有 FIG3-like tmu case 的 `40 -> 48` 局部敏感性仍可达约 6.5%，因此只能按其 audit 中的 diagnostic production 边界解读，不能自动升格为严格 regression baseline。

## 9. 正式计算命令

正式长任务优先使用 GitHub workflow：

```text
.github/workflows/relaxtime-meson-density-production.yml
```

推荐阶段顺序：

1. `convergence_low`；
2. `convergence_mid`；
3. `convergence_high/custom`；
4. 审计相邻精度、branch 和 status；
5. 使用证据支持的 custom resolution 运行 `production`。

本地命令模板：

```powershell
julia --project=. scripts/relaxtime/run_combined_meson_density_scan.jl --path <tmu|trho_asymmetric> --output-dir data/outputs/results/relaxtime/meson_density/<case_id> --figure-dir data/outputs/figures/relaxtime/meson_density/<case_id> --regimes stable,strict_bw_stage1,phase_shift_current,phase_shift_gbu_reference <path/grid/profile arguments> <converged resolution arguments> <explicit phase/density policies> --overwrite
```

`production` 标签本身不构成升格证据；必须同时存在 convergence summary、command/manifest、status audit 和人工 production audit。

## 10. 输出目录与产物合同

组合 CLI 直接输出：

- `combined_meson_density_scan.csv`；
- `README.md`；
- 可选 `combined_meson_density_scan.svg`；
- 可选 `plot_manifest.json`。

远程 workflow 追加：

- `run.command.txt`；
- `run.stdout.log`；
- `run.stderr.log`；
- `run.exitcode`；
- `remote_run_manifest.json`。

Formal production 还应包含：

- `PRODUCTION_AUDIT.md`；
- convergence summary/matrix；
- 必要的合并/修复 manifest；
- figure-side PNG/SVG 与 plot manifest。

CSV 必须保留 path、constraint、chemical profiles、regime、policy、integration、density、unsafe/status/message 和 branch diagnostics。`unsafe_bose_domain` 行必须原样保留，不能当作有限零值。

## 11. Regression / Validation 验收

最小层级：

- Unit kernel：`tests/unit/relaxtime/test_meson_density.jl`；
- Unit workflow：`tests/unit/models/test_meson_density_workflow.jl`；
- Combined contract/integration：`tests/integration/relaxtime/test_combined_meson_density_scan_smoke.jl`；
- Regime regression：`tests/regression/relaxtime/test_meson_density_regimes_regression.jl`；
- Freezeout path regression：`tests/regression/relaxtime/test_meson_density_freezeout_phase_shift_gbu_path_regression.jl`；
- Plot-review regression：`tests/regression/relaxtime/test_meson_density_plot_review_case_regression.jl`；
- Literature target smoke：`tests/validation/relaxtime/test_meson_density_literature_targets_smoke.jl`。

Literature target 或 freezeout regression 只验证其登记的 profile/path/target，不能替代任意 combined production case 的 convergence 和 branch audit。

## 12. 失败点、断点续算与重跑

- Combined CLI 当前不做行级 resume；若输出存在且未传 `--overwrite` 会拒绝运行。
- 新参数、新 profile 或新 policy 必须使用新 case/stage，避免覆盖已审阅产物。
- 远程 workflow 即使失败也应上传已有 result/figure artifacts，并保留 exit code、stdout/stderr 和 remote manifest。
- `status/message`、unsafe count 和 constraint diagnostics 是结果合同，不应通过删行“修复”。
- Branch 异常应先保存候选、pressure 和 residual 证据；targeted repair 必须记录被替换行和来源。
- Targeted repair 不等同于全网格已在新 branch policy 下验证；下一版正式 case 原则上应完整重跑。

## 13. Diagnostic 与 Formal Production 的边界

Diagnostic-only：

- 低节点 CLI smoke；
- 无相邻精度收敛；
- `x_min_cut`、no-anomalous 或 phase display 仍在比较且未写清适用范围；
- branch/order sensitivity 未审计；
- 只有 CSV/图而没有 command/manifest/audit；
- 使用 `production` stage 但未通过 convergence gate。

Formal production 至少需要：

1. 明确 path、profile、constraint、regimes 和 policy；
2. 相邻高精度 convergence evidence；
3. branch/constraint/status audit；
4. command、remote manifest、日志和退出码；
5. regression/validation 中与声明直接相关的覆盖；
6. `PRODUCTION_AUDIT.md` 明确 verdict 与残余风险。

现有 `trho_asymmetric_kplus_piplus_scan_v1` 的 production-grade verdict 只适用于其固定 grid、charged profile、`bose_x_min`、branch repair 和 resolution；不能外推到其他路径或 policy。

## 14. 后处理与作图

后处理只读取冻结 CSV：

- T-only plot：`scripts/analysis/relaxtime/render_combined_meson_density_temperature_scan.py`；
- heatmap：`scripts/analysis/relaxtime/render_combined_meson_density_fig3_like.py`。

坐标约束：

- `tmu` heatmap 默认 x 轴使用 `muq_MeV`；
- `trho_asymmetric` 必须使用 `rho_target`，不能用 equilibrium diagnostic `muq_MeV` 代替扫描轴；
- log color scale 只改变显示，不改变数据。

图像写入 `data/outputs/figures/`，并通过 `plot_manifest.json` 追溯 source CSV、字段、坐标和显示 policy。

## 15. 关联公式、API 和测试

- [Models MesonDensityWorkflow](../../../api/models/workflows/MesonDensityWorkflow.md)
- [MesonDensity API](../../../api/relaxtime/meson_density/MesonDensity.md)
- [领域 workflow 细节](../../../api/relaxtime/workflow/MesonDensityWorkflow.md)
- [公共科学计算生命周期](../common_scientific_run.md)
- `scripts/relaxtime/combined_meson_density_scan_contract.jl`
- `.github/workflows/relaxtime-meson-density-production.yml`

## 16. 最后验证记录

- 验证日期：2026-07-10
- 验证范围：combined contract、CLI/production artifact contract、现有 convergence/production audit
- 执行命令：见第 7 节及 authority map
- 状态：通过；显式开启慢 CLI 后 contract 31/31、numeric CLI 26/26、FixedAsymmetricRho contract 28/28
- 备注：FixedAsymmetricRho 慢数值 branch guard 本批未执行；既有 production audit 仍是长任务权威证据。专题入口保持 domain-candidate，formal verdict 必须逐 case 给出
