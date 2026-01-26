# Julia RelaxTime

## 【重要单位约定】
```
=============================================================================
本模块统一使用自然单位制 (ℏ = c = 1)，所有物理量单位为 fm⁻¹：
- 温度 T: fm⁻¹ 单位 (转换：T[MeV] = T[fm⁻¹] × ℏc = T[fm⁻¹] × 197.327 MeV·fm)
- 化学势 μ: fm⁻¹ 单位 (μ[MeV] = μ[fm⁻¹] × 197.327)  
- 质量 m: fm⁻¹ 单位 (m[MeV] = m[fm⁻¹] × 197.327)
- 四动量 k₀, k: fm⁻¹ 单位 (p[MeV] = p[fm⁻¹] × 197.327)
- 能量变量 z: fm⁻¹ 单位 (E[MeV] = E[fm⁻¹] × 197.327)
- 偏振函数 Π: fm² 单位 (真空极化振幅，量纲 = 1/k²)
- BPM积分 B: fm² 单位 (泡泡图积分振幅)
- Polyakov环参数 A: 无量纲
- 截断参数 Λf: fm⁻¹ 单位 (由常数模块提供，Λf ≈ 3.05 fm⁻¹)
=============================================================================
```

## 项目说明

本项目用于计算弛豫时间相关的物理量，并包含：散射运动学/矩阵元计算、PNJL 平衡求解与扫描、以及若干用于对比与可视化的脚本与网页前端。

协作与规范文件位于 `.github/`：
- 贡献指南：[`CONTRIBUTING.md`](.github/CONTRIBUTING.md)
- 行为准则：[`CODE_OF_CONDUCT.md`](.github/CODE_OF_CONDUCT.md)
- 安全策略：[`SECURITY.md`](.github/SECURITY.md)

## 快速开始（Quickstart）

完整环境复现说明见 [INSTALL.md](INSTALL.md)。这里提供最短可跑路径：

1. 初始化 Julia 环境：

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

2. 启动服务端（API + 前端）：

```powershell
julia --project=. scripts/server/server_full.jl
```

3. （可选）生成依赖图与 SVG：

```powershell
npm install
julia --project=. scripts/dev/gen_deps.jl
```

> Python 脚本环境与前端开发工具的可复现安装步骤请见 [INSTALL.md](INSTALL.md)。

## 当前功能概览

- **散射运动学**：`src/simulation/MomentumMapping.jl` 提供 2→2 运动学求解、Mandelstam 变量与椭球包络，并在 `scripts/server/server_full.jl` 中通过 `/compute` 端点暴露；`tests/unit/test_momentum_mapping.jl`、`test_frame_transformations.jl` 已覆盖核心校验。
- **散射矩阵元（当前可用）**：`src/relaxtime/ScatteringAmplitude.jl` 及依赖模块（`Polarization*`, `EffectiveCouplings`, `MesonPropagator` 等）仍是稳定入口，可提供 Σ|M|² 结果给外部积分器；相关推导见 `docs/api/relaxtime/scattering/ScatteringAmplitude.md` 与 `docs/reference/formula`。
- **截面/弛豫时间链路（修复中）**：`DifferentialCrossSection.jl`, `TotalCrossSection.jl`, `RelaxationTime*.jl` 等仍包含已知缺陷（阈值处理、归一因子、输运系数整合尚未校对），目前默认不在服务器或前端中暴露，仅供研究性参考。
- **积分与数值工具**：`src/integration/` 提供 Cauchy 主值与 Gauss-Legendre 节点，`src/utils/` 集中常用校验、数值辅助；`QuarkDistribution*.jl` 暴露各向同性/各向异性分布函数。
- **HTTP + 前端**：`scripts/server/server_full.jl` 同时提供 API 与静态资源，`web/index.html` + `web/js` 展示 3D 椭球、输入面板与健康检查指示灯；`web/simple_test.html` 适合最小交互验证。
- **文档与流程**：`docs/guides/QUICKSTART.md`、`docs/guides/USER_GUIDE.md`、`docs/guides/STATUS.md` 说明部署/排错；`docs/process/` 保留 prompt 与计划；`docs/reference/` 存放公式与推导。
- **数据与结果**：`data/outputs/` 用于收集服务器或批处理输出（例如 `data/outputs/results/relaxtime/` 下的扫描 CSV），便于跨语言/跨实现对比。
- **PNJL（求解器 + 扫描）**：`src/pnjl/` 提供 PNJL 平衡求解与扫描能力。
	- 推荐入口：`PNJL.solve(...)` + seed 策略（`MultiSeed/PhaseAwareContinuitySeed` 等），见 `docs/api/pnjl/PNJL.md` 与 `docs/api/pnjl/SeedStrategies.md`。
	- 扫描脚本：`scripts/relaxtime/run_gap_transport_scan.jl` 可批量输出平衡量与输运相关派生量到 CSV。
	- HTTP 端（实验性）：`scripts/server/server_full.jl` 提供 `POST /api/modules/pnjl-gap/run` 单点调用（请求体仍使用 `T_mev`/`mu_mev` 这类 MeV 输入字段；内部会换算到自然单位）。
- **平均散射率（实验性）**：`src/relaxtime/AverageScatteringRate.jl` 基于 Gauss-Legendre (p=32, 角度=4) 计算各向异性平均散射率，散射截面支持预计算+插值缓存。

> ⚠️ **状态说明**：截面/弛豫时间链路仍标记为“修复中”（阈值处理、归一因子与输运系数整合仍需校对）。PNJL 求解与扫描链路已可用，但物理与数值精度仍建议通过 `docs/` 下的对比报告持续验证。

## 计算链路概览（各向异性输运）

1. **能隙求解 → 序参量/有效质量**：在各向异性 PNJL 下先解能隙方程，得到序参量与三味夸克有效质量、粒子数密度（`src/pnjl/`）。
2. **弛豫时间近似 (RTA)**：输运系数采用 RTA，需要弛豫时间 \(\tau\)。
3. **弛豫时间依赖平均散射率**：\(\tau\) 由各散射过程的平均散射率 \(\Gamma\) 决定，\(\Gamma\) 需要在入射动量可行域上对散射截面加权积分。
4. **总截面需要微分截面积分**：`TotalCrossSection.jl` 对 \(\mathrm{d}\sigma/\mathrm{d}t\) 在 \(t\) 或 \(\theta^*\) 空间积分，包含各向异性分布的阻塞因子（`quark_distribution_aniso(...,\cos\theta^*)`）。
5. **微分截面依赖散射矩阵元**：`DifferentialCrossSection.jl` 使用 \(|\mathcal{M}|^2\) 与 Källén 函数的动量因子；\(|\mathcal{M}|^2\) 由 `ScatteringAmplitude.jl` 计算。
6. **矩阵元依赖介子传播子**：各散射道调用 `MesonPropagator`/`TotalPropagator` 获取传播子，按道合成振幅。
7. **传播子依赖极化函数**：`Polarization*` 使用单圈积分 `OneLoopIntegrals*`（含各向异性积分）与夸克分布函数求取极化张量。

> 当前位置：截面/弛豫时间链路仍标记为“修复中”，上述流程为实际调用路径，便于排查或扩展时定位入口。

## 最近更新

### 2026-01-25: Parameter Struct Migration
完成参数结构体迁移，提供类型安全的参数接口：
- ✅ 新增 `QuarkParams` 和 `ThermoParams` 结构体
- ✅ 双接口模式：支持结构体和 NamedTuple 两种格式
- ✅ 零性能开销：内联归一化确保无运行时损耗
- ✅ 向后兼容：现有 NamedTuple 代码无需修改
- ✅ 100%测试通过：所有模块支持双接口
- 📖 完整文档：`docs/guides/PARAMETER_STRUCT_MIGRATION.md`
- 📖 API 文档：`docs/api/PARAMETER_TYPES_API.md`

**快速使用**：
```julia
using Main.ParameterTypes: QuarkParams, ThermoParams

# 使用结构体（推荐）
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t = ThermoParams(0.15, 0.5, 0.5, 0.0)
result = relaxation_times(q, t, K_coeffs; densities=densities)

# 使用 NamedTuple（向后兼容）
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
result = relaxation_times(q_nt, t_nt, K_coeffs; densities=densities)
```

**支持的模块**：
- RelaxationTime.jl
- AverageScatteringRate.jl
- TotalCrossSection.jl
- ScatteringAmplitude.jl
- DifferentialCrossSection.jl
- TotalPropagator.jl
- ParticleSymbols.jl

### 2025-11-17: 极化函数缓存模块
新增 `PolarizationCache` 模块，通过哈希表缓存优化极化函数计算性能：
- ✅ 自动缓存相同参数的极化函数，避免重复计算
- ✅ 典型加速：10-30倍（大规模输运系数计算）
- ✅ 100%测试通过（34/34测试用例）
- 📖 完整文档：`docs/api/relaxtime/polarization/PolarizationCache.md`

**快速使用**：
```julia
using .PolarizationCache

reset_cache!()  # 开始新计算
Π = polarization_aniso_cached(...)  # 自动缓存
stats = get_cache_stats()  # 查看统计
reset_cache!()  # 释放内存
```

## 安装与使用

1. 激活项目环境：
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

2. 使用模块：
```julia
include("src/relaxtime/relaxtime.jl")
```

3. 常用入口：

- 启动本地服务器（API + 前端）：

```powershell
julia --project=. scripts/server/server_full.jl
```

- 运行 PNJL + 输运相关扫描（示例命令见脚本 Usage）：

```powershell
julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --help
```

## 项目结构

更完整的目录职责与文档/源码对齐规则见：[docs/dev/项目结构约定.md](docs/dev/%E9%A1%B9%E7%9B%AE%E7%BB%93%E6%9E%84%E7%BA%A6%E5%AE%9A.md)。

当前顶层目录及作用：

- `src/`：核心 Julia 源码，按 `integration/`、`relaxtime/`、`simulation/`、`utils/` 等子模块拆分。
- `web/`：前端页面与静态资源，`web/js/` 含 `api.js`、`ui.js`、`visualization.js`。
- `scripts/server/`：所有启动脚本与 HTTP 服务端代码（`server.jl`、`server_full.jl`、`test_server.jl`、`test_minimal_server.jl`、`start.bat`）。
- `tests/`：测试资产；`tests/unit/` 存放原 `test_unit`，`tests/analysis/` 存放原 `test_other` 中的调试脚本与诊断报告。
- `docs/`：文档中心。
	- `docs/guides/`：用户/开发手册（`README` 补充材料、`QUICKSTART.md`、`USER_GUIDE.md`、`FRONTEND_DEBUG.md`、`FIXES.md` 等）。
	- `docs/reference/`：公式、Mathematica、domain-knowledge 等原 `doc/` 内容。
	- `docs/process/prompt`, `docs/process/plans`：原 `prompt/` 与 `plans/`，用于流程记录与规范。
	- `docs/guides/examples/`：原 `examples/` 下的示例说明。
- `docs/api/`：面向外部的 API/公式描述。
- `data/outputs/`：运行结果与缓存输出（原 `results/` 位于 `data/outputs/results/`）。
- `scripts/`、`tests/`、`docs/` 之外的根文件：`Project.toml`、`Manifest.toml`、`README.md` 等项目元数据。

## 目录迁移指南

| 旧位置 | 新位置 | 说明 |
| --- | --- | --- |
| `server*.jl`, `start.bat` | `scripts/server/` | 所有后端/启动脚本集中到单一目录，`start.bat` 会自动回到仓库根目录再拉起 `server_full.jl`。 |
| `test_unit/` | `tests/unit/` | 原全部单元测试未改名，只调整路径；引用 `../../src/...` 即可。 |
| `test_other/` | `tests/analysis/` | 各类性能分析、调试脚本、诊断报告集中。 |
| `results/` | `data/outputs/results/` | 将运行产物与原始数据分离，方便清理或忽略。 |
| `doc/`（公式、domain-knowledge 等） | `docs/reference/` | 文档分类更明确。 |
| `prompt/` | `docs/process/prompt/` | 规范类文档整合到流程档案。 |
| `plans/` | `docs/process/plans/` | 规划和想法集中管理。 |
| `examples/` | `docs/guides/examples/` | 所有示例写在 Guides 下便于索引。 |
| `FIXES.md`, `FRONTEND_DEBUG.md`, `QUICKSTART.md`, `STATUS.md`, `USER_GUIDE.md` | `docs/guides/` | 一致收纳在指南目录，README 仅保留入口链接。 |

> **启动方式更新**：
> - Windows：执行 `scripts/server/start.bat`（会先 `Pkg.instantiate`，随后运行 `scripts/server/server_full.jl` 并自动打开浏览器）。
> - CLI：`julia scripts/server/server_full.jl [port]` 或 `julia scripts/server/server.jl`（仅 API）。
