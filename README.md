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

本项目用于计算弛豫时间相关的物理量。

## 当前功能概览

- **散射运动学**：`src/simulation/MomentumMapping.jl` 提供 2→2 运动学求解、Mandelstam 变量与椭球包络，并在 `scripts/server/server_full.jl` 中通过 `/compute` 端点暴露；`tests/unit/test_momentum_mapping.jl`、`test_frame_transformations.jl` 已覆盖核心校验。
- **散射矩阵元（当前可用）**：`src/relaxtime/ScatteringAmplitude.jl` 及依赖模块（`Polarization*`, `EffectiveCouplings`, `MesonPropagator` 等）仍是稳定入口，可提供 Σ|M|² 结果给外部积分器；相关推导见 `api/ScatteringAmplitude.md` 与 `docs/reference/formula`。
- **截面/弛豫时间链路（修复中）**：`DifferentialCrossSection.jl`, `TotalCrossSection.jl`, `RelaxationTime*.jl` 等仍包含已知缺陷（阈值处理、归一因子、输运系数整合尚未校对），目前默认不在服务器或前端中暴露，仅供研究性参考。
- **积分与数值工具**：`src/integration/` 提供 Cauchy 主值与 Gauss-Legendre 节点，`src/utils/` 集中常用校验、数值辅助；`QuarkDistribution*.jl` 暴露各向同性/各向异性分布函数。
- **HTTP + 前端**：`scripts/server/server_full.jl` 同时提供 API 与静态资源，`web/index.html` + `web/js` 展示 3D 椭球、输入面板与健康检查指示灯；`web/simple_test.html` 适合最小交互验证。
- **文档与流程**：`docs/guides/QUICKSTART.md`、`USER_GUIDE.md`、`STATUS.md` 说明部署/排错，`docs/process/*` 保留 prompt 与计划，`docs/reference` 存放公式与 Mathematica 推导。
- **数据与结果**：`data/outputs/results/` 用于收集服务器或批处理输出，便于与 PNJL 结果对比；尚未与 PNJL 求解器联通（见“下一步”）。

> ⚠️ **未决功能**：PNJL 能隙方程与 Excel-like UI 尚未合入；新增 API/前端模式切换在 plan 中但无实现代码，集成前请确认挂钩方案。

## 最近更新

### 2025-11-17: 极化函数缓存模块
新增 `PolarizationCache` 模块，通过哈希表缓存优化极化函数计算性能：
- ✅ 自动缓存相同参数的极化函数，避免重复计算
- ✅ 典型加速：10-30倍（大规模输运系数计算）
- ✅ 100%测试通过（34/34测试用例）
- 📖 完整文档：`api/PolarizationCache.md`

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

## 项目结构

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
- `api/`：面向外部的 API/公式描述（保持不变）。
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
