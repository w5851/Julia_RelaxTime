# 测试组织与入口规范（Julia_RelaxTime）

更新时间：2026-03-04

本文件用于统一项目内测试的放置规则、命名规则与入口策略，目标是：
- 默认测试入口"稳定、确定、快速"（CI/本地都不折腾）。
- 五层分类清晰：Unit / Integration / Regression / Validation / Benchmark。
- 非测试制品（分析脚本、性能探针）与测试严格分离。

---

## 目录分层（必须遵守）

| 层级 | 目录 | 用途 | CI 层 |
|------|------|------|-------|
| Unit | `tests/unit/` | 确定性、快速、无外部依赖的单元测试 | smoke / full |
| Numerics | `tests/unit/numerics/` | 核心数值算法测试（Gauss-Legendre、Cauchy PV 等） | smoke |
| Integration | `tests/integration/` | 跨模块端到端正确性 | smoke / nightly |
| Regression | `tests/regression/` | 内部 baseline 数值回归，对比仓库内 CSV 基线 | smoke / full |
| Validation | `tests/validation/` | 外部参考值验证（Fortran / Mathematica / 文献） | nightly / milestone |
| Benchmark | `benchmark/` | PkgBenchmark 性能基准，独立 `Project.toml` | on-demand |
| 分析脚本 | `scripts/analysis/` | 一次性分析/诊断/扫描脚本（非测试入口） | — |
| 性能脚本 | `scripts/perf/` | Profile/timing 探针脚本（非测试入口） | — |
| 交互脚本 | `scripts/` | 调试/可视化/开发辅助脚本 | — |

### 已废弃目录

- ~~`tests/analysis/`~~ → 已迁移至 `scripts/analysis/`（2026-03-04）
- ~~`tests/perf/`~~ → benchmark 文件迁至 `benchmark/`，profile 脚本迁至 `scripts/perf/`（2026-03-04）
- ~~`tests/unit/integration/`~~ → 已重命名为 `tests/unit/numerics/`（2026-03-04）

## 五层分类标准

### Unit（单元测试）
准入条件——必须全部满足：
1. **确定性**：固定随机种子；结果不依赖运行时机/线程调度/浮点噪声放大。
2. **快速**：默认配置下单文件 < 1–2s，整套 smoke < 1min。
3. **无外部依赖**：不依赖网络、不依赖大数据文件、不依赖用户本地环境（除 `--project=.` 之外）。
4. **不混杂**：不把脚本、性能探针、交互诊断代码放进 `tests/unit`。

### Integration（集成测试）
- 跨模块端到端正确性验证：gap → densities → transport 链路。
- 可以较慢（单文件 < 30s），但仍须确定性。

### Regression（回归测试）
- 使用仓库内固定 baseline CSV，验证代码变更不会引入静默数值漂移。
- 入口统一为 `tests/regression/runtests.jl`。
- smoke 运行小规模固定点；full 运行更大覆盖面或 nightly 基线。
- 导出脚本保留在 `scripts/dev/`，但门禁逻辑必须落到 `@testset`。

### Validation（验证测试）
- Fortran / Mathematica 参考值对照。
- 使用固定参考数据文件，容差在 `isapprox` 层显式标注。
- 若某个 digitized literature target 经多份 legacy 实现交叉验证后被判定为 disputed，则保留原始数据文件不改写，在 acceptance target 中禁用该点，并以单独记录 provenance 的 legacy consensus target 替补。
- 允许更宽松的运行时间。

### Benchmark（性能基准）
- 使用 `BenchmarkTools.jl` / `PkgBenchmark.jl`。
- 独立 `benchmark/Project.toml`，不污染主项目依赖。
- 入口：`benchmark/benchmarks.jl`（PkgBenchmark 标准布局）。

## 子系统分目录

按子系统管理 unit 测试：
- `tests/unit/pnjl/`：PNJL 求解器、热力学、积分缓存等
- `tests/unit/relaxtime/`：RelaxTime/散射率/截面等
- `tests/unit/numerics/`：核心数值算法（Gauss-Legendre、Cauchy PV、动量映射）
- `tests/unit/models/`：模型工厂、适配器、后端桥接
- `tests/unit/config/`：配置加载与校验

根目录 `tests/unit/*.jl` 只保留：
- `tests/unit/runtests.jl`（入口）
- 少量无法归属子系统的通用测试

## 单测命名

- 测试文件统一：`test_*.jl`
- Benchmark 文件：`bench_*.jl`（`benchmark/` 下）或 `benchmark_*.jl`
- 非测试脚本严禁命名为 `test_*.jl`

## 单测入口策略

入口文件：`tests/unit/runtests.jl`

运行档：
- `UNIT_PROFILE=smoke`（默认）：精选、稳定、确定性的测试集合，长期保持绿色。
- `UNIT_PROFILE=full`：更大范围 include，逐步修复；允许短期不全绿。

环境变量（以入口实现为准）：
- `UNIT_PROFILE=smoke|full`
- `UNIT_INCLUDE_PERF=1`：full 下允许 include performance 相关测试（默认关闭）。
- `UNIT_INCLUDE_SLOW=1`：允许 include 标记为 slow 的测试。
- `UNIT_INCLUDE_WIP=1`：允许 include 标记为 WIP 的测试。
- `UNIT_FILES=path1,path2,...`：仅运行指定文件。

## 回归测试入口策略

入口文件：`tests/regression/runtests.jl`

环境变量（以入口实现为准）：
- `REGRESSION_PROFILE=smoke|full`
- `REGRESSION_FILES=path1,path2,...`：仅运行指定回归文件
- `REGRESSION_MAGNETIC_SCOPE=smoke|nightly`：切换 PNJL magnetic baseline 口径

## PNJL 求解器相关约定

- 随机采样测试必须**确定性**（固定 RNG seed），默认样本数保守。
- 参数空间扫描放 `scripts/pnjl/`，不放 unit。

## 后端对比测试约定

> 注：`solver_backend` 参数已弃用（P3.5, 2026-03-04），现统一走 models 路径。

- 基线使用固定点集（固定 `T/μ/xi`、固定积分节点与配置），统一存放于 `tests/baselines/<domain>/`。
- 容差在测试断言层（`isapprox(...; rtol=..., atol=...)`），不通过被测函数 keyword 注入。
- Transport smoke 口径：
	- 固定点：`(T,μ)=(0.75,0.00),(0.90,0.00),(1.05,0.00),(0.90,0.15)`
	- 比较量：`η/σ/ζ`
	- 容差：`rtol=8e-2`, `atol=1e-6`
- 容差调整须在 `docs/dev/active` 记录变更原因与对比。

## Regression 与 Validation 的边界

- `tests/regression/`：当前代码 vs 仓库内部 baseline。
- `tests/validation/`：当前代码 vs 外部参考实现或文献数据。
- 不再将内部 baseline 回归继续堆入 `tests/validation/`。

## smoke → full 迁移触发条件

当测试不再满足 smoke 的"高价值、低成本、稳定"原则时应迁出：

- 移入 `full`：
	- 单文件运行时长稳定超过 ~10s；
	- 依赖较重初始化或更大输入规模。
- 移入 `benchmark/`：
	- 主要目标是性能/缩放趋势比较，非功能正确性。
- 触发后动作：
	- 在 `runtests.jl` 分组注释中说明迁移原因；
	- 在 `docs/dev/active` 记录回归影响；
	- 补一条替代 smoke 维持覆盖。

## Benchmark 基础设施

```
benchmark/
├── Project.toml       # 独立依赖（BenchmarkTools）
├── benchmarks.jl      # PkgBenchmark 入口
├── README.md          # 说明与运行方法
├── relaxtime/         # 输运相关 benchmark
│   ├── bench_*.jl
│   └── benchmark_*.jl
└── pnjl/              # PNJL 相关 benchmark
    └── *.jl
```

运行方式：
```julia
using PkgBenchmark
benchmarkpkg(".")              # 运行全部
judge(benchmarkpkg("."), "main")  # 对比分支
```
