# 脚本入口清单

本目录记录当前仓库推荐给用户直接运行的稳定脚本入口。

目标不是把所有脚本都列进来，而是明确：

- 哪些脚本是稳定入口
- 哪些脚本只是分析/开发/排障工具

---

## 1. 推荐稳定入口

全量 `run_*.jl` 脚本功能目录见：

- [run_script_catalog.md](run_script_catalog.md)

该目录用于治理与迁移盘点；本 README 仅保留稳定白名单入口。

### 默认运行规则

稳定 CLI 默认优先通过以下 wrapper 启动：

- `scripts/dev/run_with_sysimage.ps1`
- `scripts/dev/run_with_sysimage.sh`

用途：

- 若本机已有可用 sysimage，则自动追加 `--sysimage=...`
- wrapper 会同时校验 Julia 版本、平台信息与 `git_commit`
- 默认 mismatch policy 为 `rebuild`
- 若 sysimage 缺失、元数据缺失或 `git_commit` 与当前 `HEAD` 不一致，则默认自动重建本地 sysimage
- PowerShell wrapper 可配合 `-MismatchPolicy strict|fallback|rebuild`
- POSIX wrapper 可配合 `--mismatch-policy=strict|fallback|rebuild`
- `-BuildIfMissing` / `--build-if-missing` 仍保留，作为 `rebuild` 别名

如需先获取预构建 sysimage：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/bootstrap_sysimage.ps1
```

```bash
sh scripts/dev/bootstrap_sysimage.sh
```

bootstrap 脚本会根据当前平台、架构和 Julia 版本，解析 GitHub Release 资产名并解包到 `build/`。

### 前端脚本长任务目录

Web 前端的“脚本长任务”面板通过 `GET /api/modules/script-tasks` 获取可解释任务目录，并通过异步 job 路由提交：

- `POST /api/modules/script-tasks/jobs`
- `GET /api/modules/script-tasks/jobs/{job_id}`
- `GET /api/modules/script-tasks/jobs/{job_id}/result`
- `POST /api/modules/script-tasks/jobs/{job_id}/cancel`

前端目录的定位是帮助用户先理解每个 `run_*` 的用途、关键参数、输出与运行风险；默认只推荐 `smoke` preset。`canonical` 和 `custom` 必须显式确认，且仍优先经 `run_with_sysimage` wrapper 启动。为避免前端点击时隐式重建 sysimage，脚本任务服务端使用 wrapper 的 `fallback` mismatch policy：兼容时复用 sysimage，不兼容或缺失时回退到普通 `julia --project=.`。

最小示例（Windows / PowerShell）：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/models/run_unified_scan.jl scan tmu --model_kind=PNJL --T_values=150 --mu_values=0,100 --xi_values=0.0 --output_path=data/outputs/results/tmu_smoke.csv --overwrite=true
```

最小示例（Linux / macOS）：

```bash
sh scripts/dev/run_with_sysimage.sh scripts/models/run_unified_scan.jl scan tmu --model_kind=PNJL --T_values=150 --mu_values=0,100 --xi_values=0.0 --output_path=data/outputs/results/tmu_smoke.csv --overwrite=true
```

### Phase A 默认 wrapper 白名单

以下稳定 CLI 默认应经 wrapper 启动：

| 类别 | 入口 | Windows 默认 | Linux/macOS 默认 |
|---|---|---|---|
| PNJL 相图 | `scripts/pnjl/calculate_phase_structure.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| 统一扫描 | `scripts/models/run_unified_scan.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| 守恒荷易感性 | `scripts/pnjl/run_conserved_charge_susceptibilities.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| Relaxtime 编排 | `scripts/relaxtime/run_relaxtime_orchestrator.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| 输运扫描 | `scripts/relaxtime/run_gap_transport_scan.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| 服务器入口 | `scripts/server/server_full.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |

### PNJL

- [scripts/pnjl/run_conserved_charge_susceptibilities.jl](../../../scripts/pnjl/run_conserved_charge_susceptibilities.jl)
  - 守恒荷广义磁化率、累积量、`Ssigma`、`kappa_sigma2` 的统一单点/小范围扫描入口
- [scripts/models/run_unified_scan.jl](../../../scripts/models/run_unified_scan.jl)
  - 统一扫描入口（`scan tmu|trho`），覆盖 T-μ / T-ρ 网格扫描与 `Models` 主链扫描治理
- [scripts/pnjl/calculate_phase_structure.jl](../../../scripts/pnjl/calculate_phase_structure.jl)
  - 相图自动化产线入口（扫描 -> 判据 -> 报告），支持模板配置 + CLI 覆盖

#### 相图最小单命令产线（PNJL）

在仓库根目录执行：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/calculate_phase_structure.jl --model_kind=PNJL --mode=research --T_min=150 --T_max=150 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --solver_backend=models --output_dir=<your_output_dir>
```

说明：

- 默认会自动加载 `config/models/pnjl/phase_pipeline_default.toml`。
- 可以使用 `--config=<path/to/phase_pipeline.toml>` 指定自定义模板。
- 可以使用 `--preset=smoke` 快速切换到轻量可复现实验参数（随后仍可用 CLI 显式参数覆盖）。
- CLI 显式参数优先级高于模板（同名键会覆盖）。
- production 可显式设置 `crossover_T_max_MeV`、rho 粗细网格几何量门限和温度中点自适应门限；解析后的 `p_num/t_num/iterations` 与这些门限写入 manifest。
- production 默认保持 `rho_refinement_policy=uniform_nested`。`rho_support_cascade` 仅用于显式
  shadow/诊断请求：必须启用 `rho_geometry_convergence`、使用均匀嵌套 coarse rho 网格、
  `rho_refine_levels=1`，并通过 `rho_support_fine_step` 与
  `rho_support_targeted_cap` 记录两层 support/补点合同；它不会自动覆盖旧 reference。
- `rho_support_hybrid` 是更高成本的显式 shadow/候选策略，要求 `rho_refine_levels=4`。
  它按 cascade → 完整 memoized dense → 离散极值 guard 内局部验证的顺序升级；
  guard 取两个 μ 极值外侧首个严格 Stage-B 采样点，不插值、不二分、不使用固定
  padding。Stage-A targeted 总上限仍为 `12`，Stage-C 始终保留完整 Stage-B 曲线；
  没有可靠 guard 时保持 ambiguous，不启动全域 oracle。

最小产物结构（输出目录）：

- `trho_scan.csv`
- `first_order_boundary.csv`
- `spinodal.csv`
- `crossover_line.csv`
- `phase_grid_convergence.csv`
- `phase_summary.json`
- `phase_report.md`
- `run_manifest.json`（记录 argv、config_path、config_hash、git_commit、artifact_paths）

`run_manifest.json` 关键字段：

- `preset`：若使用 `--preset=...`，记录最终采用的预设名。
- `effective_config`：记录本次运行的最终有效参数快照（含被 CLI 覆盖后的值）。

### Server

- [scripts/server/server_full.jl](../../../scripts/server/server_full.jl)
  - Web + API 联调用完整入口
- [scripts/server/server.jl](../../../scripts/server/server.jl)
  - 仅 API 入口

### RelaxTime

- [scripts/relaxtime/run_manual_relaxation_scan_workflow.jl](../../../scripts/relaxtime/run_manual_relaxation_scan_workflow.jl)
  - 手动组合产物入口（`cross_section` / `temperature_scan_muB0_xi0` / `fixed_temperature_xi_scan_muB0`；兼容旧别名 `plan_a` / `plan_b`）
  - 支持 `--base-output-dir` 将结果写到隔离目录（默认 `data/outputs`）
- [scripts/relaxtime/run_phase_guided_transport_scan.jl](../../../scripts/relaxtime/run_phase_guided_transport_scan.jl)
  - phase-guided transport 上层扫描入口
  - 推荐模式名：
    - `fixed-muB-phase-scaled`（兼容短别名 `a`）
    - `fixed-T-sparse-muB`（兼容短别名 `b`）
  - `fixed-muB-phase-scaled`：固定 `mu_B`，沿 `T/T_phase` 倍率带扫描 `xi`
  - `fixed-T-sparse-muB`：固定 `T`、离散 `mu_B`、连续扫描 `xi`
  - 当前 canonical 口径：
    - mode a 固定 `mu_B = 0, 450, 900 MeV`
    - mode b 固定 `T = 120, 160, 200 MeV`
  - 默认 canonical 输出根目录：
    - `data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/`
    - `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/`
  - `compute_bulk` 默认开启；如需更快的预览扫描，可显式传 `--no-compute-bulk`
  - mode a 的一阶区域默认使用 `--phase-anchor-policy direct_coexistence`：旧 phase reference 只作为初始 bracket，再按 `--p-num` / `--t-num`（默认 `12/6`）直接求两分支等热力学势温度；`reference_interpolation` 仅用于显式复现旧计划口径。
  - 一阶可能区每点都使用稳定多初值选优，continuation 只作为候选 seed。bulk 复用主扫描已解析的热力学节点和同一点 equilibrium；导数路径若跨支会显式失败。
  - 直接锚定的 `alpha_T=1` 一阶共存切片不计算严格 `xi=0` 的唯一输运量，而是写入经主热力学节点和独立更高节点共同认证的负/正近零点，分别表示夸克侧/强子侧极限；`sampling_plan.csv` 同步记录 anchor bracket、相别和节点收敛证据。`12/6` 主节点使用 `24/8` 认证；`24/8` 主节点自动使用 `32/10` 认证，禁止同节点自比较。
  - 异常区域诊断可显式传 `--propagator-xi-policy isotropic`，仅让 σ(s)/propagator 使用 `ξ=0`；也可传 `--sigma-cache-policy validated_anchored` 诊断 threshold-subtraction 与 σ-cache 插值放大效应。默认 `match_thermo` / `default` 不变，诊断分支完成复算与收敛检查前不作为正式修复口径
  - GitHub Actions 手动入口 `Relaxtime Phase-Guided Transport Production` 只生成可审阅 artifact；新增或修改该 workflow 后需先合入默认分支让 GitHub 注册，随后才能通过 `workflow_dispatch` 触发。该入口默认 verdict 为 `diagnostic-only`，不会自动把 artifact 晋升为仓库正式数据。
  - 高精度长任务应优先通过 action shard 触发：mode a 可用 `muB_list` + `alpha_t_list` 分片，mode b 可用 `t_list` + `muB_list` 分片，二者都可用 `xi_list` 缩小窗口，并用 `shard_label` 区分 artifact。workflow 会把扫描/绘图日志写入 result artifact；失败或取消时仍尽量上传 partial CSV、`failed_points.csv`、`channel_diagnostics.csv` 和日志，供本地合并与 convergence gate 使用。
  - 当前 phase-guided transport production-grade case：
    - 高 xi 分辨率、同 p128 validated-anchored 积分精度：`first_canonical_v1_p128_xi001_validated_anchored_prod_v1`
      - mode a: `data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/`
      - mode b: `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_xi001_validated_anchored_prod_v1/`
      - shared import gate evidence: `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_xi001_validated_anchored_prod_v1_convergence/`
    - 低 xi 分辨率 anchor / p104-vs-p128 convergence 依据：`first_canonical_v1_p128_validated_anchored_prod_v1`
      - mode a: `data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_validated_anchored_prod_v1/`
      - mode b: `data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_validated_anchored_prod_v1/`
      - shared convergence evidence: `data/outputs/results/relaxtime/transport/phase_guided/first_canonical_v1_p128_validated_anchored_prod_v1_convergence/`
  - 旧 `first_canonical_v1` 已从仓库数据树移除。它缺少 `validated_anchored` sigma-cache policy、显式高精度 tau/sigma 积分参数、channel diagnostics 和全网格 convergence gate；历史对照请通过 git history / PR #122 查看，局部 `xi` 结构的 production-grade 解读应优先使用高 xi p128 production-grade case。
- [scripts/relaxtime/run_phase_guided_transport_plots.jl](../../../scripts/relaxtime/run_phase_guided_transport_plots.jl)
  - canonical case 的 post-processing / plot-review wrapper
  - 图层正式落盘到 `data/outputs/figures/relaxtime/transport/phase_guided/<mode>/<case_name>/`
  - `fixed-muB-phase-scaled`：每个固定 `mu_B` 一张图，图内多条 `alpha_T`
  - `fixed-T-sparse-muB`：每个固定 `T` 一张图，图内多条 `mu_B`
  - 输出 `plot_manifest.json` 并回写 case README

`temperature_scan_muB0_xi0` / `fixed_temperature_xi_scan_muB0` 目录最小溯源产物：

- 扫描 CSV（行级 `run_id`）
- `effective_config.json`（最终有效参数快照）
- `run_manifest.json`（`argv`、`git_commit`、`config_hash`、`artifacts`、`summary`）

其中 `fixed_temperature_xi_scan_muB0_merged.csv` 额外包含：

- `source_file`
- `source_T_MeV`

用于将合并行反向定位到温度分片 CSV。

---

## 2. 不应视为稳定用户入口的目录

- `scripts/dev/`
  - 开发期导出、迁移、一次性工具
- `scripts/analysis/`
  - 后处理与分析脚本
- `scripts/debug/`
  - 排障脚本
- `scripts/perf/`
  - 性能探针与分析

这些目录仍然重要，但默认不应在用户指南里作为正式入口推荐。

---

## 3. 命名约定建议

为了降低“入口太多”的感受，建议以后遵守：

- `run_*.jl`
  - 稳定、可文档化、面向用户的入口
- `calculate_*.jl`
  - 计算型脚本，但不一定是稳定 CLI
- `export_*.jl`
  - 基线/导出/开发辅助
- `diagnose_*.jl`, `debug_*.jl`
  - 排障专用

结论：

- 不是脚本数量本身有问题。
- 问题是稳定入口和内部工具没有层级。
- 这份清单的作用就是把层级固定下来。
