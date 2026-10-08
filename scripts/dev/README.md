生成项目依赖图（简要说明）

目的：自动化生成 `src/` 内的模块/文件依赖图，并输出为 Mermaid Markdown，放在 `docs/architecture/dependencies.md`。

快速使用：

```powershell
# 在项目根执行（Windows）
julia --project=. scripts/dev/gen_deps.jl
```

生成依赖审计报告：

```powershell
julia --project=. scripts/dev/analyze_deps.jl
```

文档一致性检查：

```powershell
julia --project=. scripts/dev/check_docs_consistency.jl
```

科学计算 SOP 权威映射与结构检查：

```powershell
julia --project=. scripts/dev/check_sop_governance.jl
```

该门禁以 `config/governance/docs_authority_map.toml` 为机读单一来源，检查 active SOP 的路径、唯一权威范围、稳定入口白名单、必需章节、旧入口模式和复核周期。首期只约束 `docs/guides/sop/`，不一次性阻断全部历史文档。

跨主线 task ledger 检查：

```powershell
julia --project=. scripts/dev/check_task_ledger.jl
julia --project=. scripts/dev/check_task_ledger.jl --preflight
julia --project=. scripts/dev/check_task_ledger.jl --preflight --track rs-transport
julia --project=. scripts/dev/check_task_ledger.jl --preflight --full-paths
```

该检查只读校验 `config/governance/task_tracks.toml` 的状态、依赖、任务文件、evidence
和分支/SHA/run 格式；`--preflight` 报告分支、HEAD、变更路径总数和前 10 项，
`--full-paths` 展开全部路径。有效 ledger 还报告任务文件、阻塞、下一步及最多 5 条
evidence 入口；无效 ledger 仍保留 Git 上下文并失败退出。`--track ID` 选择并校验已有主线。

导出 API 全集索引生成：

```powershell
julia --project=. scripts/dev/generate_api_export_index.jl \
	--module-file src/models/exports_public.jl \
	--module-file src/models/Models.jl \
	--output docs/api/generated/models/ModelsExportIndex.md \
	--title "Models Export API Index"
```

相图主题导出全集生成：

```powershell
julia --project=. scripts/dev/generate_api_export_index.jl \
	--module-file src/models/Models.jl \
	--module-file src/models/entrypoints.jl \
	--include-symbol run_phase_pipeline \
	--include-symbol find_cep \
	--include-symbol build_phase_artifacts \
	--include-symbol resolve_phase_output_target \
	--include-symbol promote_phase_artifacts \
	--include-symbol PhasePipelineResult \
	--include-symbol CEPResult \
	--include-symbol PromotionResult \
	--output docs/api/models/phase/generated/Exports.md \
	--title "Phase Export API Index"
```

工作流主题导出全集生成：

```powershell
julia --project=. scripts/dev/generate_api_export_index.jl \
	--module-file src/models/Models.jl \
	--module-file src/models/entrypoints.jl \
	--include-symbol solve_gap_and_transport \
	--include-symbol solve_transport_from_equilibrium \
	--include-symbol solve_gap_and_meson_point \
	--include-symbol solve_meson_point_from_equilibrium \
	--include-symbol solve_rotation_point \
	--include-symbol solve_gas_liquid_point \
	--include-symbol transport_workflow_module \
	--include-symbol meson_workflow_module \
	--include-symbol workflow_module_for \
	--include-symbol workflow_param_adapters_module \
	--include-symbol pnjl_module \
	--output docs/api/models/workflows/generated/Exports.md \
	--title "Workflow Export API Index"
```

扫描主题导出全集生成：

```powershell
julia --project=. scripts/dev/generate_api_export_index.jl \
	--module-file src/models/entrypoints.jl \
	--include-symbol run_tmu_scan \
	--include-symbol run_trho_scan \
	--include-symbol build_default_rho_grid \
	--output docs/api/models/scans/generated/Exports.md \
	--title "Scan Export API Index"
```

模型变体 magnetic 主题导出全集生成：

```powershell
julia --project=. scripts/dev/generate_api_export_index.jl \
	--module-file src/models/exports_public.jl \
	--module-file src/models/Models.jl \
	--include-symbol PNJLMagneticModel \
	--include-symbol magnetic_gap_residual \
	--include-symbol magnetic_gap_residual_autodiff \
	--include-symbol solve_magnetic_gap \
	--include-symbol alpha_n \
	--include-symbol energy_landau \
	--include-symbol smooth_cutoff \
	--include-symbol resolve_nmax_from_cutoff \
	--include-symbol resolve_magnetic_nmax \
	--include-symbol omega0_flavor_landau \
	--include-symbol omegat_flavor_landau \
	--include-symbol density_flavor_landau \
	--include-symbol MAGNETIC_EB_MIN_MEV2 \
	--include-symbol MAGNETIC_EB_MIN_FM2 \
	--include-symbol validate_magnetic_eB \
	--include-symbol MagneticIMCParams \
	--include-symbol default_imc_params \
	--include-symbol coupling_GB \
	--include-symbol MagneticConfig \
	--include-symbol default_magnetic_config \
	--include-symbol calculate_magnetic_omega_components \
	--include-symbol calculate_magnetic_omega \
	--include-symbol calculate_magnetic_pressure \
	--include-symbol calculate_magnetic_rho \
	--include-symbol calculate_magnetic_number_densities \
	--include-symbol magnetic_nmax_convergence_report \
	--output docs/api/models/variants/magnetic/generated/Exports.md \
	--title "Models Magnetic Variant Export API Index"
```

衍生量 susceptibility 主题导出全集生成：

```powershell
julia --project=. scripts/dev/generate_api_export_index.jl \
	--module-file src/models/exports_public.jl \
	--module-file src/models/Models.jl \
	--include-symbol conserved_charge_susceptibility \
	--include-symbol chi_BQS \
	--include-symbol chi_B \
	--include-symbol chi1_B \
	--include-symbol chi2_B \
	--include-symbol chi3_B \
	--include-symbol chi4_B \
	--include-symbol chi_Q \
	--include-symbol chi1_Q \
	--include-symbol chi2_Q \
	--include-symbol chi3_Q \
	--include-symbol chi4_Q \
	--include-symbol chi_S \
	--include-symbol chi1_S \
	--include-symbol chi2_S \
	--include-symbol chi3_S \
	--include-symbol chi4_S \
	--include-symbol chi11_BQ \
	--include-symbol chi11_BS \
	--include-symbol chi11_QS \
	--include-symbol cumulant_B \
	--include-symbol cumulant_BQS \
	--include-symbol baryon_Ssigma \
	--include-symbol baryon_kappa_sigma2 \
	--output docs/api/models/derived/susceptibility/generated/Exports.md \
	--title "Models Susceptibility Export API Index"
```

衍生量 derivatives 主题导出全集生成：

```powershell
julia --project=. scripts/dev/generate_api_export_index.jl \
	--module-file src/models/exports_public.jl \
	--module-file src/models/Models.jl \
	--include-symbol mass_derivatives \
	--include-symbol thermo_derivatives \
	--include-symbol bulk_derivative_coeffs \
	--include-symbol dP_dT \
	--include-symbol dP_dmu \
	--include-symbol bulk_viscosity_coefficients \
	--include-symbol compute_B_bracket \
	--output docs/api/models/derived/derivatives/generated/Exports.md \
	--title "Models Derivatives Export API Index"
```

relaxtime transport 主题导出全集生成：

```powershell
julia --project=. scripts/dev/generate_api_export_index.jl \
	--module-file src/models/Models.jl \
	--module-file src/relaxtime/TransportCoefficients.jl \
	--module-file src/relaxtime/RelaxationTime.jl \
	--module-file src/relaxtime/AverageScatteringRate.jl \
	--include-symbol transport_provider \
	--include-symbol TransportProvider \
	--include-symbol prepare_transport_provider \
	--include-symbol default_transport_provider \
	--include-symbol TransportIntegrationConfig \
	--include-symbol TransportRequest \
	--include-symbol shear_viscosity \
	--include-symbol electric_conductivity \
	--include-symbol bulk_viscosity \
	--include-symbol bulk_viscosity_isentropic \
	--include-symbol transport_coefficients \
	--include-symbol relaxation_times \
	--include-symbol compute_average_rates \
	--include-symbol REQUIRED_PROCESSES \
	--include-symbol average_scattering_rate \
	--include-symbol CrossSectionCache \
	--include-symbol build_w0cdf_pchip_cache \
	--output docs/api/relaxtime/transport/generated/Exports.md \
	--title "relaxtime Transport Export API Index"
```

核心求解与约束主题导出全集生成：

```powershell
julia --project=. scripts/dev/generate_api_export_index.jl \
	--module-file src/models/exports_public.jl \
	--module-file src/models/Models.jl \
	--include-symbol create_model \
	--include-symbol MeanFieldState \
	--include-symbol meanfield_state \
	--include-symbol state_vector \
	--include-symbol normalize_mu_vec \
	--include-symbol solve_gap \
	--include-symbol ImplicitProblem \
	--include-symbol build_pnjl_fixedmu_problem \
	--include-symbol build_pnjl_flavor_mu_problem \
	--include-symbol build_njl_problem \
	--include-symbol solve \
	--include-symbol solve_multi \
	--include-symbol ConstraintModes \
	--include-symbol FixedMuBConservedCharges \
	--include-symbol flavor_mu_from_bqs \
	--include-symbol conserved_mu_from_flavor \
	--include-symbol conserved_densities_from_flavor \
	--include-symbol GapParams \
	--include-symbol build_conditions \
	--include-symbol build_residual! \
	--include-symbol gap_state_dim \
	--include-symbol gap_residual \
	--include-symbol SeedStrategy \
	--include-symbol DefaultSeed \
	--include-symbol MultiSeed \
	--include-symbol HybridContinuitySeed \
	--include-symbol get_seed \
	--include-symbol get_all_seeds \
	--include-symbol update! \
	--include-symbol solve_with_derivatives \
	--include-symbol solve_pnjl_with_derivatives \
	--include-symbol solve_pnjl_with_flavor_mu_derivatives \
	--include-symbol solve_constraint \
	--output docs/api/models/solver/generated/Exports.md \
	--title "Solver Export API Index"
```

说明：
- 该脚本用于生成“导出 API 全集”这一层视图；
- 现已支持多次传入 `--module-file` 与 `--include-symbol`，可为主题页生成聚合后的公开导出子集；
- 会扫描目标文件中的 `export` 列表，并统计这些符号是否已在 `docs/api/` 的非 generated 页面中被提及；
- 适合与人工维护的“面向用户入口”“算法核心”页面配套使用。

归档元信息补齐（历史文档一次性修复）：

```powershell
julia --project=. scripts/dev/backfill_archived_frontmatter.jl
```

Benchmark 阈值检查：

```powershell
julia --project=. scripts/dev/check_benchmark_thresholds.jl
```

散射截面热点基准（legacy vs fast_path）：

```powershell
julia --project=. scripts/dev/benchmark_total_cross_section_hotpath.jl
```

总截面固定点 baseline 导出：

```powershell
julia --project=. scripts/dev/export_total_cross_section_fixedpoint_baseline.jl
```

PNJL 扫描固定点 baseline 导出：

```powershell
julia --project=. scripts/dev/export_pnjl_scan_fixedpoint_baseline.jl
```

unit skip/deprecated 门禁检查：

```powershell
julia --project=. scripts/dev/check_unit_skip_policy.jl
```

active 文档治理检查（命名 + 归档触发）：

```powershell
julia --project=. scripts/dev/check_active_docs_governance.jl
```

根依赖与 agent 指令治理检查：

```powershell
julia --project=. scripts/dev/check_dependency_policy.jl
julia --project=. scripts/dev/check_agent_instructions.jl
```

- `check_dependency_policy.jl` 强制可选数值 oracle 与根 runtime/test 环境隔离。
- `check_agent_instructions.jl` 强制根 `AGENTS.md` 保持精简，并验证详细命令指南存在且结构完整。

数据与图像输出路径门禁：

```powershell
julia --project=. scripts/dev/check_data_output_path_guard.jl
julia --project=. scripts/dev/check_data_output_path_guard.jl --base origin/main --head HEAD
```

该门禁禁止新增默认输出回流到根目录 `outputs/results`，禁止新增 tracked PNG/SVG/PDF 正式图像资产写入 `data/outputs/results/`，并禁止新增非图像/非 `plot_manifest.json` 资产写入 `data/outputs/figures/`；正式图像应写入对应的 `data/outputs/figures/`，CSV/JSON/README/audit/logs 应留在 `data/outputs/results/`。唯一例外是结果治理中登记的 immutable charged GBU `freezeout_20260907_v2` 哈希绑定历史包；它不得复制为新 case，新的 smoke 运行仍使用 figure-side 输出。

PNJL 迁移门禁检查（限制 src/pnjl 新增核心实现）：

```powershell
julia --project=. scripts/dev/check_pnjl_migration_guard.jl
```

PNJL 裁剪波次门禁（删除前快照 + 删除后 smoke + 失败回滚提示）：

```powershell
julia --project=. scripts/dev/run_prune_wave_gate.jl
```

PNJL 下线前检查表生成（`src/models/scans/*` + `src/pnjl/PNJL.jl`）：

```powershell
julia --project=. scripts/dev/generate_pnjl_decommission_checklist.jl --base HEAD --head HEAD
```

产物路径：
- `outputs/results/pnjl_decommission_checklist_<timestamp>.md`

开发执行记录仅在任务明确要求时按指定文件的现有格式补充，不需要专用追加脚本，也不自动创建台账。记录与历史保留规则见[任务跟踪治理](../../docs/dev/task_tracking_governance.md)。

检查表新增内容：
- `PNJL` 顶层 `scans/*` 默认 include 审计（排除白名单）
- `PNJL` scan loader 可见性清单（`load_*scan*!`）
- `PNJL.solve/PNJL.solve_multi` 外部调用方审计（排除 `src/pnjl/**`）

说明：`run_prune_wave_gate.jl` 会在执行前输出 `models_invokelatest_allowlist.toml` 的增删摘要（若当前 diff 有变更），并将摘要写入快照文件，便于统一入口审阅。

同时会额外生成独立审计产物：
- `outputs/results/models_invokelatest_allowlist_delta_<timestamp>.txt`
- `outputs/results/pnjl_scan_default_include_audit_<timestamp>.txt`

PR 范围裁剪门禁示例：

```powershell
julia --project=. scripts/dev/run_prune_wave_gate.jl --base origin/main --head HEAD
```

失败自动回滚（仅 working-tree 模式）：

```powershell
julia --project=. scripts/dev/run_prune_wave_gate.jl --auto-rollback
```

PR 范围检查示例：

```powershell
julia --project=. scripts/dev/check_pnjl_migration_guard.jl --base origin/main --head HEAD
```

检查项包含：
- `docs/guides/**/*.md` 中历史路径（如 `test_unit/`、`julia server.jl`、`.\\start.bat`、`doc/domain-knowledge/`）
- 当 `README.md` 标注“修复中”时，guides 中是否出现“系统完全就绪/已知问题: 无”等绝对化状态词
- 当 `README.md` 标注“已验证可用”时，guides 中是否残留“修复中”表述
- `src/models/**/*.jl` 新增 `Base.invokelatest(...)` 时必须命中迁移门禁白名单（防止新增分散 world-age 边界）
- 输出 `models-invokelatest-audit` 漂移提示（`observed` vs `allowlist_baseline`）以辅助审计白名单漂移
- 白名单单一来源：`config/ci/models_invokelatest_allowlist.toml`（门禁脚本与文档均以该文件为准）
- 当白名单文件在当前 diff 中变更时，门禁输出 `models-invokelatest-allowlist-diff`（added/removed 摘要）辅助 PR 审阅
- 禁止在 `src/**` 新增 `PNJL.run_tmu_scan/run_trho_scan` 运行时调用（要求统一走 `Models.run_tmu_scan/run_trho_scan`）
- 禁止在 `src/tests/scripts` 新增 `PNJL.TmuScan/PNJL.TrhoScan` 子模块直连依赖（要求统一走 `PNJL.run_*` 顶层入口）
- 禁止在 `src/tests/scripts` 新增对 `src/pnjl/PNJL.jl` 的直接 include 路径依赖（要求统一走 `src/models/Models.jl` 入口）
- 输出 `pnjl-scan-runtime-dependency-audit`，审计 `src/**`（排除 `src/pnjl/**`）中残留 `PNJL.run_*` 调用点
- 禁止在 `src/pnjl/PNJL.jl` 新增 `scans/*` 默认 include（避免默认耦合回流；当前仅允许门禁白名单项）
- 输出 `pnjl-scan-default-include-audit`，审计 `src/pnjl/PNJL.jl` 中 `scans/*` 默认 include 现状（排除白名单）；非 0 将触发门禁失败

可复现安装（推荐）：

```powershell
npm install
```

手动渲染 SVG（可选）：

```powershell
npm run deps:render
```

输出文件：`docs/architecture/dependencies.md`（脚本覆盖）
辅助文件：`docs/architecture/dependencies.mmd`、`docs/architecture/dependencies.svg`

注意：当前脚本为 MVP：
- 只解析 `include("...")` 和 `using .Module` / `import .Module`（带点号的内部模块引用）
- 忽略第三方包依赖（通过 `Pkg.dependencies()` 可单独导出）
- 在发现循环依赖时会在输出中列出强连通分量

后续增强建议：
- 将 `web/js` 的 ES module import 解析加入图
- 在 CI 中加入自动检查并阻止新增循环依赖
- 使用 `docs/architecture/dependency_rules.md` 维护目录级依赖矩阵

---

## 归档开发文档

将已完成或明确取消/被替代的 `docs/dev/active` 直接 Markdown 子文件迁移到 `docs/dev/archived`，添加元信息并保留原文字节。状态与引用责任见[开发文档治理](../../docs/dev/README.md)。

### 交互式归档

```powershell
julia --project=. scripts/dev/archive_docs.jl -i
```

### 批量归档指定文件

```powershell
julia --project=. scripts/dev/archive_docs.jl file1.md file2.md
```

### 使用自定义日期归档

```powershell
julia --project=. scripts/dev/archive_docs.jl -d 2026-01-15 file1.md
```

### 验证已归档文件格式

```powershell
julia --project=. scripts/dev/archive_docs.jl --check 2026-01-15_file1.md
```

### 预览归档操作（不实际执行）

```powershell
julia --project=. scripts/dev/archive_docs.jl --dry-run file1.md
```

默认 `--status completed`；取消或被替代须指定 `--status cancelled` / `--status superseded` 和 `--reason "终止原因"`，预览与执行参数保持一致。

- 帮助、预览和检查不创建目录或文件；路径以脚本所在仓库为准。
- 预览和执行都拒绝已有目标、范围外路径及链接重定向；批量先检查所有计划，执行失败不回滚已完成项。
- 写入临时文件并验证后，通过硬链接发布目标，再核对源文件未变后删除源文件；文件系统不支持硬链接时安全失败。
- `--check <文件...>` 检查本批元数据；省略文件名检查全目录。只解析平面标量 frontmatter，核对字段类型、日期和终止原因，错误返回非零退出码。
- 迁移后由调用者更新有效引用和 ledger 路径；脚本不自动修复历史文档。

## 输运固定点基线候选

[export_transport_fixedpoint_baseline.jl](export_transport_fixedpoint_baseline.jl) 为[输运固定点回归](../../tests/regression/relaxtime/test_transport_fixedpoint_regression.jl)生成候选 CSV。已接受的数据位于 [baseline_transport_fixedpoints_v1.csv](../../tests/baselines/relaxtime/baseline_transport_fixedpoints_v1.csv)；候选的比较与接受遵循[基线管理指南](../../docs/guides/BASELINE_VERSION_MANAGEMENT.md)。

### 计算与输出契约

- 计算脚本预设的 13 个 `(T, mu, xi)` 点，固定各味 `tau=1`、`compute_tau=false`，求得平衡态后计算剪切黏度、电导率和体黏度。这组固定弛豫时间结果用于内部回归，不承担完整弛豫时间计算或正式数值生产。
- CSV 列为 `T,mu,xi,eta,sigma,zeta`。T/mu 使用内部自然单位 fm^-1，xi 无量纲；输运量含义见[输运 API](../../docs/api/relaxtime/transport/README.md)。点位、求解器及积分设置以脚本为准。
- 必须传 `--output <new-candidate.csv>`；既有文件、目录或链接均被拒绝，没有覆盖选项。
- `--backend` 仅支持 `models`。缺失参数、已有目标和 `--help` 在加载数值模型前处理。
- 导出器保留显式 `HADRON_SEED_5` 初值；当前固定点回归使用 `seed_state=nothing`。两者的初值策略并不相同，回归通过不能证明候选一定可以导出。求解失败会报告具体点位；调整初值策略属于数值语义变更，需单独验证。
- 先生成同目录临时文件，检查每点收敛、有限性以及完整列/点位顺序，再通过硬链接发布完整新文件。目标文件系统需支持硬链接；不支持时返回失败，不改用可能覆盖目标的退路。
- 计算或校验失败不会留下目标候选 CSV。临时文件会清理；并发出现同名目标时保留对方文件并失败退出。生成成功仍只是候选，不修改基线消费者。
- 这些文件保护属于本导出器；其他历史导出器的参数、覆盖行为和失败处理需分别核对。

### 生成与验证

生成一次新的临时候选：

```powershell
$candidatePath = Join-Path $env:TEMP ('transport_candidate_' + [guid]::NewGuid().ToString('N') + '.csv')
julia --project=. scripts/dev/export_transport_fixedpoint_baseline.jl --output $candidatePath
```

旧基线的针对性回归：

```powershell
julia --project=. -e 'ENV["REGRESSION_FILES"]="relaxtime/test_transport_fixedpoint_regression.jl"; include("tests/regression/runtests.jl")'
```

文件操作与失败处理的轻量检查：

```powershell
julia --project=. tests/unit/config/test_transport_baseline_export.jl
```

夹具检查不加载物理模型；修改计算参数或真实调用路径时还应运行实际候选导出与对应回归，保持既有数值基线不变。
