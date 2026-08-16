# Julia RelaxTime

Julia PNJL/NJL 计算平台：相图与扫描、Mott 相变/介子质量、各向异性输运/弛豫时间、介子热力学、`Models` 统一入口与可复现工作流。

[![Latest Release](https://img.shields.io/github/v/release/w5851/Julia_RelaxTime?label=release)](https://github.com/w5851/Julia_RelaxTime/releases/latest)

[中文](#中文) | [面向-llm--agent](#面向-llm--agent)

---

<a name="中文"></a>

## 中文

Latest Release: [releases/latest](https://github.com/w5851/Julia_RelaxTime/releases/latest)

## 1) 三步 Quickstart（脚本链路优先）

下面是新用户默认路径：先跑通相图最小产线，再按文档扩展。

> 预期时长：完成环境准备后，主命令通常在 5-10 分钟内得到最小产物（视机器配置而定）。

<details open>
<summary><b>面向人类用户</b></summary>

### Step 1. 环境准备

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

### Step 2. 运行最小可复现实验

Windows / PowerShell：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/calculate_phase_structure.jl --preset=smoke --output_dir=data/outputs/results/phase_smoke
```

Linux / macOS：

```bash
sh scripts/dev/run_with_sysimage.sh scripts/pnjl/calculate_phase_structure.jl --preset=smoke --output_dir=data/outputs/results/phase_smoke
```

### Step 3. 校验产物

```powershell
julia --project=. -e 'println(isfile("data/outputs/results/phase_smoke/phase_summary.json") && isfile("data/outputs/results/phase_smoke/phase_report.md"))'
```

输出为 `true` 表示最小链路成功。

可选次级入口（非主链路）：

```powershell
julia --project=. scripts/server/server_full.jl
```

Quickstart 运行后若不希望保留示例产物，可清理：

```powershell
Remove-Item -Recurse -Force "data/outputs/results/phase_smoke"
```

</details>

<details>
<summary><b>面向 LLM / Agent</b></summary>

请在执行任何修改前先遵循以下规则（与“面向人类”说明分离）：

1. 先读取仓库协作约束：`AGENTS.md`、`.github/copilot-instructions.md`。
2. 统一入口优先：`Models` 与 `src/models/entrypoints.jl`。
3. non-fixedmu 求解模式按“展平联合求解”治理，不要引入新的分层默认路径。
4. PR76 后契约对齐：mixed-meson 治理与 non-fixedmu 联合求解语义默认保持不变；若需变更必须在任务范围内显式声明并补回归证据。
5. 目录治理：
   - 分析脚本放 `scripts/analysis/`
   - 性能探针放 `scripts/perf/`
   - 非测试脚本不要放入 `tests/`
6. 测试执行顺序优先 smoke profile；测试分层保持 `unit/integration/regression/validation`。
7. 可优先使用分层 wrapper 入口：`test/unit.jl`、`test/integration.jl`、`test/regression.jl`、`test/validation.jl`。
8. 运行稳定 CLI 时，优先走 `scripts/dev/run_with_sysimage.ps1`，以便在本机存在 sysimage 时稳定复用冷启动优化。
9. 稳定公共入口变更需同步更新 `docs/api/`；新增核心模块必须补 unit tests。
10. 若工作区有用户已有改动：不要覆盖/回滚无关改动；仅提交本任务相关文件。

建议最小验证命令（agent 默认基线）：

```powershell
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. scripts/dev/check_models_entry_contract.jl
```

</details>

## 2) 稳定脚本入口矩阵（白名单）

稳定用户入口以 `docs/guides/scripts/README.md` 为准；README 仅列核心白名单。

### 这个项目能做哪些物理计算（以及用哪个脚本）

下面按“核心稳定入口 + 专题能力入口”给出可执行能力矩阵。

#### A. 核心稳定入口（默认优先）

| 计算能力 | 典型用途 | 推荐脚本入口 | 最小用法（示例） |
|---|---|---|---|
| PNJL 相结构 / 相图产线 | 生成 boundary/spinodal/crossover/CEP 与报告 | `scripts/pnjl/calculate_phase_structure.jl` | `powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/calculate_phase_structure.jl --preset=smoke --output_dir=data/outputs/results/phase_smoke` |
| PNJL T-μ / T-ρ 扫描 | `Models` 主链统一网格扫描、单点/批量求解 | `scripts/models/run_unified_scan.jl` | `powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/models/run_unified_scan.jl scan tmu --model_kind=PNJL --T_values=150 --mu_values=0,100 --xi_values=0.0 --output_path=data/outputs/results/tmu_smoke.csv --overwrite=true` |
| 守恒荷易感性与累积量 | `chi_BQS` / cumulant / `Ssigma` / `kappa_sigma2` | `scripts/pnjl/run_conserved_charge_susceptibilities.jl` | `powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/run_conserved_charge_susceptibilities.jl --help` |
| 各向异性 PNJL 输运系数扫描（PNJL_aniso） | 平衡求解 + 弛豫时间 + RTA 输运系数批量计算 | `scripts/relaxtime/run_gap_transport_scan.jl` | `powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/relaxtime/run_gap_transport_scan.jl --help` |
| RelxTime 工作流编排 | 统一触发 `transport` / `cross-section` 产线 | `scripts/relaxtime/run_relaxtime_orchestrator.jl` | `powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/relaxtime/run_relaxtime_orchestrator.jl transport --help` |
| 模型服务/API 调用 | 通过 HTTP 服务调用模型求解能力 | `scripts/server/server_full.jl` | `powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/server/server_full.jl` |

#### B. 专题能力入口（研究/后处理常用）

| 计算能力 | 典型用途 | 脚本入口 | 最小用法（示例） |
|---|---|---|---|
| Mott 相变扫描 | 介子质量与 Mott 阈值随 `T/xi` 变化扫描 | `scripts/relaxtime/run_mott_phase_scan.jl` | `julia --project=. scripts/relaxtime/run_mott_phase_scan.jl --help` |
| Phase-guided transport 邻域扫描 | 围绕相变线 / 固定温度稀疏窗口生成 canonical transport case | `scripts/relaxtime/run_phase_guided_transport_scan.jl` | `julia --project=. scripts/relaxtime/run_phase_guided_transport_scan.jl --mode fixed-T-sparse-muB --xi-list -0.5,-0.2,0.0,0.2 --muB-list 0,450,900 --T-list 120,160,200 --dry-run` |
| Gap + 介子质量联合扫描 | 生成 Mott 相关基础数据（mass/width/threshold） | `scripts/relaxtime/run_gap_meson_mass_scan.jl` | `julia --project=. scripts/relaxtime/run_gap_meson_mass_scan.jl --help` |
| Mott 派生 CSV / 可视化模式 | 从主扫描结果生成派生字段与绘图输入 | `scripts/relaxtime/run_mott_phase_derived_csv.jl` / `scripts/relaxtime/run_mott_phase_plot_modes.jl` | `julia --project=. scripts/relaxtime/run_mott_phase_derived_csv.jl --help` |
| 各向异性相图模板实验 | 按 `xi` 批量跑扫描 + 相结构 + 可选绘图 | `scripts/pnjl/run_aniso_phase_template.jl` | `julia --project=. scripts/pnjl/run_aniso_phase_template.jl --profile=smoke --xi-values=0.0,0.2` |
| 磁场 PNJL 单点/扫描 | `eB` 依赖的热力学与密度计算 | `scripts/pnjl/run_magnetic_point.jl` / `scripts/pnjl/run_magnetic_eb_scan.jl` | `julia --project=. scripts/pnjl/run_magnetic_eb_scan.jl` |
| 手动工作流产物编排 | 人工控制 `cross_section/temperature_scan_muB0_xi0/fixed_temperature_xi_scan_muB0` 产物生成（兼容旧别名 `plan_a/plan_b`） | `scripts/relaxtime/run_manual_relaxation_scan_workflow.jl` | `julia --project=. scripts/relaxtime/run_manual_relaxation_scan_workflow.jl --help` |

说明：
- 稳定白名单以 `docs/guides/scripts/README.md` 为准。
- `run_*.jl` 全量能力目录见 `docs/guides/scripts/run_script_catalog.md`。
- 若在 Windows / PowerShell 环境运行稳定 CLI，默认优先使用 `scripts/dev/run_with_sysimage.ps1`。
- 若在 Linux / macOS 环境运行稳定 CLI，默认优先使用 `scripts/dev/run_with_sysimage.sh`。
- 两个 wrapper 都会在本机 sysimage 可用且与当前 `HEAD` 匹配时自动接管冷启动优化；默认若 sysimage 缺失、元数据缺失或 commit 漂移，则自动重建本地 sysimage。
- fresh clone / 新机器如需获取匹配的预构建 sysimage，可先运行 `scripts/dev/bootstrap_sysimage.ps1` 或 `scripts/dev/bootstrap_sysimage.sh`。
- wrapper mismatch policy 为 `fallback | strict | rebuild`，默认 `rebuild`。

### 从稳定入口到深层文档

- 相图主产线：`scripts/pnjl/calculate_phase_structure.jl`
  - 执行与复现 SOP：`docs/guides/sop/workflows/pnjl_phase_structure.md`
  - 用户说明：`docs/guides/scripts/README.md`
  - API 入口：`docs/api/models/phase/README.md`
- T-μ / T-ρ 扫描：`scripts/models/run_unified_scan.jl`
  - 用户说明：`docs/guides/scripts/README.md`
  - API 入口：`docs/api/models/scans/README.md`
- `Models` 工作流与统一编排
  - API 入口：`docs/api/models/workflows/README.md`
- 守恒荷 susceptibility / cumulant
  - API 入口：`docs/api/models/derived/susceptibility/README.md`
  - 当前导数后端口径：单方向 `B/Q/S` 默认走 TaylorDiff fast path，mixed BQS 走内部 multivariate Taylor jet；旧 `:forwarddiff` susceptibility fallback 已下线。
- 各向异性输运 / Relaxtime 工作流
  - 执行与复现 SOP：`docs/guides/sop/workflows/relaxtime_transport.md`
  - 用户说明：`docs/guides/scripts/README.md`
  - API 入口：`docs/api/relaxtime/transport/README.md`
  - phase-guided production-grade asset（mode a）：`data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi005_validated_anchored_prod_v1/`
  - phase-guided production-grade 图层（mode a）：`data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi005_validated_anchored_prod_v1/`
  - phase-guided production-grade asset（mode b）：`data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_xi005_validated_anchored_prod_v1/`
  - phase-guided production-grade 图层（mode b）：`data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_xi005_validated_anchored_prod_v1/`
- 介子热力学 / 介子数密度
  - 介子热力学 SOP：`docs/guides/sop/workflows/meson_thermodynamics.md`
  - 介子数密度 SOP：`docs/guides/sop/workflows/meson_density.md`
  - API 入口：`docs/api/models/workflows/MesonThermoWorkflow.md`、`docs/api/models/workflows/MesonDensityWorkflow.md`
- Web/API 服务入口：`scripts/server/server_full.jl`
  - 当前状态：`docs/guides/STATUS.md`

## 3) 能力边界与非目标

- 稳定入口采用白名单治理；未进入白名单的脚本默认不作为用户入口。
- `scripts/dev/`、`scripts/analysis/`、`scripts/debug/`、`scripts/perf/` 主要用于开发/分析/排障/性能探针。
- 历史路径与旧调用方式仅作兼容参考，不建议作为新流程基线。
- 前端页面当前定位为框架与联调载体，不作为本仓库“主计算能力入口”。

## 4) 推荐验证命令

### 产物级验证（面向运行者）

```powershell
julia --project=. -e 'println(isfile("data/outputs/results/phase_smoke/phase_summary.json") && isfile("data/outputs/results/phase_smoke/phase_report.md"))'
```

### 维护者检查（面向协作/提交流程）

```powershell
julia --project=. scripts/dev/check_docs_consistency.jl
julia --project=. scripts/dev/check_sop_governance.jl
julia --project=. scripts/dev/check_script_entrypoints.jl
julia --project=. scripts/dev/check_models_entry_contract.jl
julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'
```

## 5) 文档导航（按目标）

### 用户/运行

- `docs/guides/QUICKSTART.md`
- `docs/guides/USER_GUIDE.md`
- `docs/guides/STATUS.md`
- `docs/guides/scripts/README.md`
- `docs/guides/sop/README.md`

### 开发/治理

- `docs/dev/README.md`
- `docs/dev/active/`
- `docs/dev/backlog/`
- `docs/dev/archived/`
- `docs/guides/RELEASE_GOVERNANCE_v0.1.1.md`

### API/实现细节

- `docs/api/README.md`
- `docs/architecture/`
- `docs/reference/`

## 6) 仓库结构（极简图）

```text
src/            Julia 核心实现（models/relaxtime/simulation/utils）
scripts/        可执行脚本入口（pnjl/server 为稳定入口核心）
tests/          unit/integration/regression/validation 分层测试
docs/           guides/api/architecture/reference/dev 文档体系
config/         模型与物理配置（TOML）
data/outputs/   默认运行产物目录
web/            前端框架与静态资源
```

架构与原理请看 `docs/architecture/` 与 `docs/reference/`；README 只保留可执行索引。

## 7) 贡献与治理入口

- 贡献指南：`.github/CONTRIBUTING.md`
- 行为准则：`.github/CODE_OF_CONDUCT.md`
- 安全策略：`.github/SECURITY.md`
- 文档治理：`scripts/dev/check_docs_consistency.jl`、`scripts/dev/check_sop_governance.jl` 与 `scripts/dev/check_script_entrypoints.jl`

协作原则：README 保持“入口与边界”，深内容以 `docs/` 为准。
