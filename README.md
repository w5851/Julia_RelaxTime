# Julia RelaxTime

Julia PNJL/NJL 计算平台：相图与扫描、各向异性输运/弛豫时间、`Models` 统一入口与可复现工作流。

Latest Release: [releases/latest](https://github.com/w5851/Julia_RelaxTime/releases/latest)

## 1) 三步 Quickstart（脚本链路优先）

下面是新用户默认路径：先跑通相图最小产线，再按文档扩展。

> 预期时长：完成环境准备后，主命令通常在 5-10 分钟内得到最小产物（视机器配置而定）。

### Step 1. 环境准备

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

### Step 2. 运行最小可复现实验

```powershell
julia --project=. scripts/pnjl/calculate_phase_structure.jl --preset=smoke --output_dir=data/outputs/results/phase_smoke
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

## 2) 稳定脚本入口矩阵（白名单）

稳定用户入口以 `docs/guides/scripts/README.md` 为准；README 仅列核心白名单。

### PNJL

- `scripts/pnjl/run_conserved_charge_susceptibilities.jl`
- `scripts/pnjl/run_tmu_scan.jl`
- `scripts/pnjl/calculate_phase_structure.jl`

### Server/API

- `scripts/server/server_full.jl`
- `scripts/server/server.jl`

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

### 开发/治理

- `docs/dev/active/`
- `docs/dev/archived/`
- `docs/guides/RELEASE_GOVERNANCE_v0.1.1.md`

### API/实现细节

- `docs/api/`
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
- 文档治理：`scripts/dev/check_docs_consistency.jl` 与 `scripts/dev/check_script_entrypoints.jl`

协作原则：README 保持“入口与边界”，深内容以 `docs/` 为准。
